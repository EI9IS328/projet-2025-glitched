//************************************************************************
//   proxy application v.0.0.1
//
//  semproxy.cpp: the main interface of  proxy application
//
//************************************************************************

#include "sem_proxy.h"

#include <cartesian_struct_builder.h>
#include <cartesian_unstruct_builder.h>
#include <math.h>
#include <measure.h>
#include <sem_solver_acoustic.h>
#include <source_and_receiver_utils.h>

#include <algorithm>
#include <cstdint>  // For uint32_t (RLE)
#include <cxxopts.hpp>
#include <filesystem>
#include <filesystem>  // For creating directories
#include <fstream>
#include <iostream>
#include <limits>   // For numeric_limits
#include <numeric>  // For std::accumulate (mean)
#include <ostream>
#include <string>
#include <valarray>  // For FFT

using namespace SourceAndReceiverUtils;

const double PI = 3.141592653589793238460;
typedef std::complex<double> Complex;
typedef std::valarray<Complex> CArray;

// https://stackoverflow.com/a/37729648/5353851
// /!\ FFT implementation was done using gemini-cli
void fft(CArray& x)
{
  const size_t N = x.size();
  if (N <= 1) return;
  // divide
  CArray even = x[std::slice(0, N / 2, 2)];
  CArray odd = x[std::slice(1, N / 2, 2)];
  // conquer
  fft(even);
  fft(odd);
  // combine
  for (size_t k = 0; k < N / 2; ++k)
  {
    Complex t = std::polar(1.0, -2 * PI * k / N) * odd[k];
    x[k] = even[k] + t;
    x[k + N / 2] = even[k] - t;
  }
}

SEMproxy::SEMproxy(const SemProxyOptions& opt)
{
  const int order = opt.order;
  nb_elements_[0] = opt.ex;
  nb_elements_[1] = opt.ey;
  nb_elements_[2] = opt.ez;
  nb_nodes_[0] = opt.ex * order + 1;
  nb_nodes_[1] = opt.ey * order + 1;
  nb_nodes_[2] = opt.ez * order + 1;

  const float spongex = opt.boundaries_size;
  const float spongey = opt.boundaries_size;
  const float spongez = opt.boundaries_size;
  const std::array<float, 3> sponge_size = {spongex, spongey, spongez};
  src_coord_[0] = opt.srcx;
  src_coord_[1] = opt.srcy;
  src_coord_[2] = opt.srcz;

  domain_size_[0] = opt.lx;
  domain_size_[1] = opt.ly;
  domain_size_[2] = opt.lz;

  bool isModelOnNodes = opt.isModelOnNodes;
  isElastic_ = opt.isElastic;
  cout << boolalpha;
  bool isElastic = isElastic_;

  snapshot_folder_ = opt.snapshot_folder_path;
  snapshot_iterations_interval_ = opt.snapshot_interval;
  snapshot_in_situ_ = opt.snapshot_in_situ;
  snapshot_colormap_ = colormapFromString(opt.snapshot_colormap);
  snapshot_slice_axis_ = getSliceAxis(opt.snapshot_slice_axis);
  if (opt.snapshot_folder_path.length() > 0)
  {
    should_snapshot_ = true;
    // Create snapshot directory if it doesn't exist
    std::filesystem::create_directories(snapshot_folder_);
  }

  in_situ_stats_ = opt.in_situ_stats;
  in_situ_folder_ = opt.in_situ_folder;
  if (in_situ_stats_)
  {
    std::filesystem::create_directories(in_situ_folder_);
  }

  const SolverFactory::methodType methodType = getMethod(opt.method);
  const SolverFactory::implemType implemType = getImplem(opt.implem);
  const SolverFactory::meshType meshType = getMesh(opt.mesh);
  const SolverFactory::modelLocationType modelLocation =
      isModelOnNodes ? SolverFactory::modelLocationType::OnNodes
                     : SolverFactory::modelLocationType::OnElements;
  const SolverFactory::physicType physicType =
      SolverFactory::physicType::Acoustic;

  float lx = domain_size_[0];
  float ly = domain_size_[1];
  float lz = domain_size_[2];
  int ex = nb_elements_[0];
  int ey = nb_elements_[1];
  int ez = nb_elements_[2];
  rcvs_size_ = opt.rcvs.size();

  if (meshType == SolverFactory::Struct)
  {
    switch (order)
    {
      case 1: {
        model::CartesianStructBuilder<float, int, 1> builder(
            ex, lx, ey, ly, ez, lz, isModelOnNodes);
        m_mesh = builder.getModel();
        break;
      }
      case 2: {
        model::CartesianStructBuilder<float, int, 2> builder(
            ex, lx, ey, ly, ez, lz, isModelOnNodes);
        m_mesh = builder.getModel();
        break;
      }
      case 3: {
        model::CartesianStructBuilder<float, int, 3> builder(
            ex, lx, ey, ly, ez, lz, isModelOnNodes);
        m_mesh = builder.getModel();
        break;
      }
      default:
        throw std::runtime_error(
            "Order other than 1 2 3 is not supported (semproxy)");
    }
  }
  else if (meshType == SolverFactory::Unstruct)
  {
    model::CartesianParams<float, int> param(order, ex, ey, ez, lx, ly, lz,
                                             isModelOnNodes);
    model::CartesianUnstructBuilder<float, int> builder(param);
    m_mesh = builder.getModel();
  }
  else
  {
    throw std::runtime_error("Incorrect mesh type (SEMproxy ctor.)");
  }

  // time parameters
  if (opt.autodt)
  {
    float cfl_factor = (order == 2) ? 0.5 : 0.7;
    dt_ = find_cfl_dt(cfl_factor);
  }
  else
  {
    dt_ = opt.dt;
  }
  timemax_ = opt.timemax;
  num_sample_ = timemax_ / dt_;

  m_solver = SolverFactory::createSolver(methodType, implemType, meshType,
                                         modelLocation, physicType, order);
  m_solver->computeFEInit(*m_mesh, sponge_size, opt.surface_sponge,
                          opt.taper_delta);

  // watched receivers list
  if (rcvs_size_ > m_mesh->getNumberOfElements())
  {
    throw std::runtime_error(
        "trying to define more receivers than there are elements");
  }
  for (int i = 0; i < rcvs_size_; i++)
  {
    auto rcv = opt.rcvs[i];
    auto x = std::get<0>(rcv);
    auto y = std::get<1>(rcv);
    auto z = std::get<2>(rcv);
    if (x < 0 || y < 0 || z < 0 || x > m_mesh->domainSize(0) ||
        y > m_mesh->domainSize(1) || z > m_mesh->domainSize(2))
    {
      std::ostringstream errorSS;
      errorSS << "trying to allocate coordinates outside of the domain: x=" << x
              << "y=" << y << "z=" << z;
      throw std::runtime_error(errorSS.str());
    }
    rcvs_coord_.push_back({x, y, z});
  }

  // watched reveivers output
  if (!opt.watchedReceiversOutputPath.empty())
  {
    saveWatchedReceiversOutput = true;
    watchedReceiversOutputPath = opt.watchedReceiversOutputPath;
    watchedReceiversOutputFormat =
        opt.watchedReceiversOutputFormat == "bin" ? BIN : PLAIN;
  }

  snapshot_format = opt.snapshot_format == "bin" ? BIN : PLAIN;

  bool set_quantize_level = false;
  if (opt.compression_method_ == "none")
  {
    compression_method_ = None;
  }
  else if (opt.compression_method_ == "rle")
  {
    compression_method_ = RLE;
  }
  else if (opt.compression_method_ == "quant")
  {
    compression_method_ = Quant;
    set_quantize_level = true;
  }
  else if (opt.compression_method_ == "quant_rle")
  {
    compression_method_ = QuantRLE;
    set_quantize_level = true;
  }
  else
  {
    throw std::runtime_error("Unsupported compression method " +
                             opt.compression_method_);
  }

  if (set_quantize_level)
  {
    switch (opt.quant_level_)
    {
      case 1:
        quant_level_ = OneByte;
        break;
      case 2:
        quant_level_ = TwoByte;
        break;
      default:
        throw std::runtime_error("unsupported quantization level: " +
                                 std::to_string(opt.quant_level_));
    }
  }

  if (!opt.saveReport.empty())
  {
    saveReportPath = opt.saveReport;
  }

  initFiniteElem();

  std::cout << "Number of node is " << m_mesh->getNumberOfNodes() << std::endl;
  std::cout << "Number of element is " << m_mesh->getNumberOfElements()
            << std::endl;
  std::cout << "Launching the Method " << opt.method << ", the implementation "
            << opt.implem << " and the mesh is " << opt.mesh << std::endl;
  std::cout << "Model is on " << (isModelOnNodes ? "nodes" : "elements")
            << std::endl;
  std::cout << "Physics type is " << (isElastic ? "elastic" : "acoustic")
            << std::endl;
  std::cout << "Order of approximation will be " << order << std::endl;
  std::cout << "Time step is " << dt_ << "s" << std::endl;
  std::cout << "Simulated time is " << timemax_ << "s" << std::endl;
}

void SEMproxy::run()
{
  Measure metrics;
  SEMsolverDataAcoustic solverData(i1, i2, myRHSTerm, pnGlobal, rhsElement,
                                   rhsWeights);

  metrics.startClock(Global);
  for (int indexTimeSample = 0; indexTimeSample < num_sample_;
       indexTimeSample++)
  {
    metrics.startClock(Kernel);
    m_solver->computeOneStep(dt_, indexTimeSample, solverData);
    metrics.stopClockAndAppend(Kernel);

    if (indexTimeSample % 50 == 0)
    {
      m_solver->outputSolutionValues(indexTimeSample, i1, rhsElement[0],
                                     pnGlobal, "pnGlobal");
    }

    if (should_snapshot_ &&
        indexTimeSample % snapshot_iterations_interval_ == 0)
    {
      // in-situ stats for snapshots
      if (in_situ_stats_)
      {
        this->i2 = i2;
        computeInSituSnapshotStats(solverData.m_pnGlobal, indexTimeSample);
      }

      metrics.startClock(MakeSnapshots);
      // create path string
      std::string extension = (snapshot_format == BIN) ? ".bin" : ".txt";
      std::ostringstream stringStream;
      stringStream << snapshot_folder_;
      stringStream << "/snapshot";
      stringStream << indexTimeSample;
      if (snapshot_in_situ_)
      {
        if (snapshot_colormap_ == COLORMAP_GRAYSCALE)
        {
          stringStream << ".pgm";
        }
        else
        {
          stringStream << ".ppm";
        }
      }
      else
      {
        stringStream << ".bin";
      }
      std::string snapshot_file_path = stringStream.str();

      std::cout << "snapshoting at " << snapshot_file_path << std::endl;

      // open snapshot file
      ofstream snapshot_file;
      snapshot_file.open(snapshot_file_path);
      int dim = m_mesh->getOrder() + 1;
      int ex_ = nb_elements_[0];
      int ey_ = nb_elements_[1];
      int ez_ = nb_elements_[2];
      int order_ = m_mesh->getOrder();

      if (!snapshot_in_situ_)
      {
        if (snapshot_format == BIN)
        {
          struct snapshot_header_t
          {
            int ex;
            int ey;
            int ez;
            int order;
            enum CompressionMethod
                compression_method;  // 0 = raw, 1 = RLE encoded
          };

          struct snapshot_header_t hdr = {
              .ex = ex_,
              .ey = ey_,
              .ez = ez_,
              .order = order_,
              .compression_method = compression_method_,
          };
          snapshot_file.write(reinterpret_cast<char*>(&hdr), sizeof(hdr));

          if (compression_method_ == RLE)
          {
            // RLE encoding: write (count, value) pairs
            const float* data = solverData.m_pnGlobal.data();
            size_t size = solverData.m_pnGlobal.size();
            size_t i = 0;
            while (i < size)
            {
              float value = data[i];
              uint32_t count = 1;
              while (i + count < size && data[i + count] == value &&
                     count < UINT32_MAX)
              {
                count++;
              }
              snapshot_file.write(reinterpret_cast<char*>(&count),
                                  sizeof(uint32_t));
              snapshot_file.write(reinterpret_cast<const char*>(&value),
                                  sizeof(float));
              i += count;
            }
          }
          else if (compression_method_ == Quant)
          {
            // /!\ header additions: 1 byte for quant level (1 or 2) and 4 bytes
            // for the scale
            const float* original_data = solverData.m_pnGlobal.data();
            size_t nb_elems = solverData.m_pnGlobal.size();
            // first we find the absolute max value
            float abs_max = 0.0f;
            for (size_t i = 0; i < nb_elems; i++)
            {
              float cur_v = original_data[i];
              if (std::isfinite(cur_v))
              {
                float val = std::abs(cur_v);
                if (abs_max < val) abs_max = val;
              }
            }
            if (abs_max == 0.0f) abs_max = 1.0f;  // avoid divide by 0
            // get max possible value depending on quant level
            float quant_type_max =
                (quant_level_ == OneByte) ? 127.0f : 32767.0f;
            float scale = quant_type_max / abs_max;
            // write the header additons
            snapshot_file.write(reinterpret_cast<char*>(&quant_level_),
                                sizeof(uint8_t));
            snapshot_file.write(reinterpret_cast<char*>(&scale), sizeof(float));

            // Variable to accumulate squared errors
            double sum_squared_error = 0.0;
            size_t valid_count = 0;

            // we quantize and save the content
            if (quant_level_ ==
                OneByte)  // ik its ugly but I don't know c++ very well
            {
              std::vector<int8_t> buffer(nb_elems);
              for (size_t cur = 0; cur < nb_elems; cur++)
              {
                float cur_v = original_data[cur];
                if (!std::isfinite(cur_v)) cur_v = 0.0f;
                buffer[cur] = static_cast<int8_t>(std::max(
                    -quant_type_max,
                    std::min(quant_type_max, std::round(cur_v * scale))));

                // Compute error for RMSE
                float quantized_value = static_cast<float>(buffer[cur]) / scale;
                float error = cur_v - quantized_value;
                sum_squared_error += error * error;
                valid_count++;
              }
              snapshot_file.write(reinterpret_cast<char*>(buffer.data()),
                                  buffer.size() * sizeof(int8_t));
            }
            else
            {
              std::vector<int16_t> buffer(nb_elems);
              for (size_t cur = 0; cur < nb_elems; cur++)
              {
                float cur_v = original_data[cur];
                if (!std::isfinite(cur_v)) cur_v = 0.0f;
                buffer[cur] = static_cast<int16_t>(std::max(
                    -quant_type_max,
                    std::min(quant_type_max, std::round(cur_v * scale))));

                float quantized_value = static_cast<float>(buffer[cur]) / scale;
                float error = cur_v - quantized_value;
                sum_squared_error += error * error;
                valid_count++;
              }
              snapshot_file.write(reinterpret_cast<char*>(buffer.data()),
                                  buffer.size() * sizeof(int16_t));
            }

            // Compute and output RMSE
            double rmse = std::sqrt(sum_squared_error / valid_count);
            std::cout << "Quantization RMSE: " << rmse << std::endl;
          }
          else if (compression_method_ ==
                   QuantRLE)  // asked Gemini to fuse Quant and RLE, would be
                              // better to extract everything into proper
                              // functions and do it by hand but time is running
                              // out
          {
            const float* original_data = solverData.m_pnGlobal.data();
            size_t nb_elems = solverData.m_pnGlobal.size();

            // --- STEP 1: CALCULATE SCALE (Same as your Quant code) ---
            float abs_max = 0.0f;
            for (size_t i = 0; i < nb_elems; i++)
            {
              float val = std::abs(original_data[i]);
              if (std::isfinite(val) && abs_max < val) abs_max = val;
            }
            if (abs_max == 0.0f) abs_max = 1.0f;

            // Determine max value based on bit depth
            float quant_type_max =
                (quant_level_ == OneByte) ? 127.0f : 32767.0f;
            float scale = quant_type_max / abs_max;

            // Write Header
            snapshot_file.write(reinterpret_cast<char*>(&quant_level_),
                                sizeof(uint8_t));
            snapshot_file.write(reinterpret_cast<char*>(&scale), sizeof(float));

            // --- RMSE TRACKING ---
            double sum_squared_error = 0.0;
            size_t valid_count = 0;

            // --- STEP 2: HELPER LAMBDA ---
            // A small helper to quantize a specific index.
            // This keeps the loop logic below clean.
            auto get_quantized_val = [&](size_t idx) -> int32_t {
              float cur_v = original_data[idx];
              if (!std::isfinite(cur_v)) cur_v = 0.0f;
              // The quantization math:
              int32_t quantized = static_cast<int32_t>(std::max(
                  -quant_type_max,
                  std::min(quant_type_max, std::round(cur_v * scale))));

              // Compute error for RMSE
              float quantized_value = static_cast<float>(quantized) / scale;
              float error = cur_v - quantized_value;
              sum_squared_error += error * error;
              valid_count++;

              return quantized;
            };

            // --- STEP 3: RLE LOOP ON QUANTIZED DATA ---
            if (nb_elems > 0)
            {
              // Initialize with the first element
              int32_t current_run_val = get_quantized_val(0);
              uint32_t current_count = 1;

              for (size_t i = 1; i < nb_elems; i++)
              {
                int32_t next_val = get_quantized_val(i);

                // If it matches and we haven't overflowed the counter
                if (next_val == current_run_val && current_count < UINT32_MAX)
                {
                  current_count++;
                }
                else
                {
                  // Value changed: Write the previous run
                  snapshot_file.write(reinterpret_cast<char*>(&current_count),
                                      sizeof(uint32_t));

                  // Write the value (cast back to int8 or int16)
                  if (quant_level_ == OneByte)
                  {
                    int8_t v = static_cast<int8_t>(current_run_val);
                    snapshot_file.write(reinterpret_cast<char*>(&v),
                                        sizeof(int8_t));
                  }
                  else
                  {
                    int16_t v = static_cast<int16_t>(current_run_val);
                    snapshot_file.write(reinterpret_cast<char*>(&v),
                                        sizeof(int16_t));
                  }

                  // Reset for the new run
                  current_run_val = next_val;
                  current_count = 1;
                }
              }

              // --- STEP 4: FLUSH FINAL RUN ---
              snapshot_file.write(reinterpret_cast<char*>(&current_count),
                                  sizeof(uint32_t));
              if (quant_level_ == OneByte)
              {
                int8_t v = static_cast<int8_t>(current_run_val);
                snapshot_file.write(reinterpret_cast<char*>(&v),
                                    sizeof(int8_t));
              }
              else
              {
                int16_t v = static_cast<int16_t>(current_run_val);
                snapshot_file.write(reinterpret_cast<char*>(&v),
                                    sizeof(int16_t));
              }
            }

            // --- STEP 5: COMPUTE AND OUTPUT RMSE ---
            if (valid_count > 0)
            {
              double rmse = std::sqrt(sum_squared_error / valid_count);
              std::cout << "Quantization+RLE RMSE: " << rmse << std::endl;
            }
          }
          else
          {
            snapshot_file.write(
                reinterpret_cast<char*>(solverData.m_pnGlobal.data()),
                solverData.m_pnGlobal.size() * sizeof(float));
          }
        }
        else
        {
          snapshot_file << ex_ << ',' << ey_ << ',' << ez_ << ',' << order_
                        << '\n';
          for (int elementNumber = 0;
               elementNumber < m_mesh->getNumberOfElements(); elementNumber++)
          {
            for (int i = 0; i < m_mesh->getNumberOfPointsPerElement(); ++i)
            {
              int x = i % dim;
              int z = (i / dim) % dim;
              int y = i / (dim * dim);
              int const globalIdx =
                  m_mesh->globalNodeIndex(elementNumber, x, y, z);
              snapshot_file << solverData.m_pnGlobal(globalIdx, i2);

              if (i != m_mesh->getNumberOfPointsPerElement() -
                           1)  // if not last point of the element
              {
                snapshot_file << ",";
              }
            }
            snapshot_file << std::endl;
          }
        }
      }
      else  // IN SITU
      {
        // slice along the configured axis (0=X, 1=Y, 2=Z)
        int sliced_dim = snapshot_slice_axis_;
        int slice_pos_along_sliced_dim = domain_size_[sliced_dim] / 2;

        // determine the two axes for the image (the ones we're NOT slicing)
        // img_axis1 and img_axis2 are the remaining dimensions
        int img_axis1 = (sliced_dim + 1) % 3;
        int img_axis2 = (sliced_dim + 2) % 3;

        // Write PGM (grayscale) or PPM (color) header
        if (snapshot_colormap_ == COLORMAP_GRAYSCALE)
        {
          snapshot_file.write("P5\n", 3);
          snapshot_file << nb_nodes_[img_axis1] << " " << nb_nodes_[img_axis2]
                        << std::endl;
          snapshot_file.write("255\n", 4);
        }
        else
        {
          snapshot_file.write("P6\n", 3);
          snapshot_file << nb_nodes_[img_axis1] << " " << nb_nodes_[img_axis2]
                        << std::endl;
          snapshot_file.write("255\n", 4);
        }

        float max_pressure;
        bool first_max_pressure = false;
        for (int elementNumber = 0;
             elementNumber < m_mesh->getNumberOfElements(); elementNumber++)
        {
          for (int i = 0; i < m_mesh->getNumberOfPointsPerElement(); ++i)
          {
            int x = i % dim;
            int z = (i / dim) % dim;
            int y = i / (dim * dim);

            const int globalIdx =
                m_mesh->globalNodeIndex(elementNumber, x, y, z);

            float global_coords[3];
            global_coords[0] = m_mesh->nodeCoord(globalIdx, 0);
            global_coords[1] = m_mesh->nodeCoord(globalIdx, 1);
            global_coords[2] = m_mesh->nodeCoord(globalIdx, 2);

            if (global_coords[sliced_dim] != slice_pos_along_sliced_dim)
            {  // only get data from the slice
              continue;
            }
            if (!first_max_pressure)
            {
              max_pressure = solverData.m_pnGlobal(globalIdx, i2);
              first_max_pressure = true;
            }
            if (solverData.m_pnGlobal(globalIdx, i2) > max_pressure)
            {
              max_pressure = solverData.m_pnGlobal(globalIdx, i2);
            }
          }
        }

        printf("max_pressure=%f\n", max_pressure);
        float _step[3] = {
            domain_size_[0] / (float)(nb_nodes_[0] - 1),
            domain_size_[1] / (float)(nb_nodes_[1] - 1),
            domain_size_[2] / (float)(nb_nodes_[2] - 1),
        };
        printf("step_x = %f, step_y = %f, step_z = %f\n", _step[0], _step[1],
               _step[2]);

        // Allocate image buffer (1 byte for grayscale, 3 bytes for color)
        int bytes_per_pixel =
            (snapshot_colormap_ == COLORMAP_GRAYSCALE) ? 1 : 3;
        char* img = (char*)calloc(
            nb_nodes_[img_axis1] * nb_nodes_[img_axis2] * bytes_per_pixel, 1);

        for (int elementNumber = 0;
             elementNumber < m_mesh->getNumberOfElements(); elementNumber++)
        {
          for (int i = 0; i < m_mesh->getNumberOfPointsPerElement(); ++i)
          {
            int x = i % dim;
            int z = (i / dim) % dim;
            int y = i / (dim * dim);

            const int globalIdx =
                m_mesh->globalNodeIndex(elementNumber, x, y, z);

            float global_coords[3];
            global_coords[0] = m_mesh->nodeCoord(globalIdx, 0);
            global_coords[1] = m_mesh->nodeCoord(globalIdx, 1);
            global_coords[2] = m_mesh->nodeCoord(globalIdx, 2);

            if (global_coords[sliced_dim] != slice_pos_along_sliced_dim)
            {  // only get data from the slice
              continue;
            }
            int pixel_index =
                ((int)(global_coords[img_axis1] / _step[img_axis1])) *
                    nb_nodes_[img_axis2] +
                ((int)(global_coords[img_axis2] / _step[img_axis2]));

            uint8_t grayscale_value;
            if (max_pressure > 0.0 &&
                solverData.m_pnGlobal(globalIdx, i2) > 0.0)
            {
              grayscale_value =
                  (solverData.m_pnGlobal(globalIdx, i2) * 255.0) / max_pressure;
            }
            else
            {
              grayscale_value = 0;
            }

            if (snapshot_colormap_ == COLORMAP_GRAYSCALE)
            {
              img[pixel_index] = grayscale_value;
            }
            else
            {
              // Apply colormap and write RGB values
              RGB rgb = applyColormap(grayscale_value, snapshot_colormap_);
              int img_offset = pixel_index * 3;
              img[img_offset] = rgb.r;
              img[img_offset + 1] = rgb.g;
              img[img_offset + 2] = rgb.b;
            }
          }
        }
        snapshot_file.write(
            img, nb_nodes_[img_axis1] * nb_nodes_[img_axis2] * bytes_per_pixel);
        free(img);
      }

      snapshot_file.close();
      metrics.stopClockAndAppend(MakeSnapshots);
      metrics.measureIO(stringStream.str());
    }

    // Save pressure for every receiver
    const int order = m_mesh->getOrder();

    metrics.startClock(MakeSismos);
    for (int rcvIdx = 0; rcvIdx < rcvs_size_; rcvIdx++)
    {
      float varnp1 = 0.0;
      for (int i = 0; i < order + 1; i++)
      {
        for (int j = 0; j < order + 1; j++)
        {
          for (int k = 0; k < order + 1; k++)
          {
            int nodeIdx =
                m_mesh->globalNodeIndex(rhsElementRcv[rcvIdx], i, j, k);
            int globalNodeOnElement =
                i + j * (order + 1) + k * (order + 1) * (order + 1);
            varnp1 += pnGlobal(nodeIdx, i2) *
                      rhsWeightsRcv(rcvIdx, globalNodeOnElement);
          }
        }
      }
      pnAtReceiver(rcvIdx, indexTimeSample) = varnp1;
    }
    metrics.stopClockAndAppend(MakeSismos);

    swap(i1, i2);

    auto tmp = solverData.m_i1;
    solverData.m_i1 = solverData.m_i2;
    solverData.m_i2 = tmp;
  }

  // handling save of watched receiver data:
  if (saveWatchedReceiversOutput)
  {
    metrics.startClock(OutputSismos);
    if (watchedReceiversOutputFormat == BIN)
    {
      save_watched_receivers_output_bin(metrics);
    }
    else
    {
      save_watched_receivers_output_plain(metrics);
    }
    metrics.stopClockAndAppend(OutputSismos);
  }

  if (in_situ_stats_)
  {
    exportInSituStats();
  }

  metrics.stopClockAndAppend(Global);

  cout << metrics;
  if (saveReportPath)
  {
    std::ofstream output(saveReportPath.value(),
                         std::ios::trunc | std::ios::out);
    output << metrics;
    output.close();
  }
}

// Initialize arrays
void SEMproxy::init_arrays()
{
  cout << "Allocate host memory for source and pressure values ..." << endl;

  rhsElement = allocateVector<vectorInt>(myNumberOfRHS, "rhsElement");
  rhsWeights = allocateArray2D<arrayReal>(
      myNumberOfRHS, m_mesh->getNumberOfPointsPerElement(), "RHSWeight");
  myRHSTerm = allocateArray2D<arrayReal>(myNumberOfRHS, num_sample_, "RHSTerm");
  pnGlobal =
      allocateArray2D<arrayReal>(m_mesh->getNumberOfNodes(), 2, "pnGlobal");
  pnAtReceiver =
      allocateArray2D<arrayReal>(rcvs_size_, num_sample_, "pnAtReceiver");
  // Receivers
  rhsElementRcv = allocateVector<vectorInt>(rcvs_size_, "rhsElementRcv");
  rhsWeightsRcv = allocateArray2D<arrayReal>(
      rcvs_size_, m_mesh->getNumberOfPointsPerElement(), "RHSWeightRcv");
}

// Initialize sources
void SEMproxy::init_source()
{
  arrayReal myRHSLocation = allocateArray2D<arrayReal>(1, 3, "RHSLocation");
  // std::cout << "All source are currently are coded on element 50." <<
  // std::endl;
  std::cout << "All source are currently are coded on middle element."
            << std::endl;
  int ex = nb_elements_[0];
  int ey = nb_elements_[1];
  int ez = nb_elements_[2];

  int lx = domain_size_[0];
  int ly = domain_size_[1];
  int lz = domain_size_[2];

  // Get source element index

  int source_index = floor((src_coord_[0] * ex) / lx) +
                     floor((src_coord_[1] * ey) / ly) * ex +
                     floor((src_coord_[2] * ez) / lz) * ey * ex;

  for (int i = 0; i < 1; i++)
  {
    rhsElement[i] = source_index;
  }

  // Get coordinates of the corners of the sourc element
  float cornerCoords[8][3];
  int I = 0;
  int nodes_corner[2] = {0, m_mesh->getOrder()};
  for (int k : nodes_corner)
  {
    for (int j : nodes_corner)
    {
      for (int i : nodes_corner)
      {
        int nodeIdx = m_mesh->globalNodeIndex(rhsElement[0], i, j, k);
        cornerCoords[I][0] = m_mesh->nodeCoord(nodeIdx, 0);
        cornerCoords[I][2] = m_mesh->nodeCoord(nodeIdx, 2);
        cornerCoords[I][1] = m_mesh->nodeCoord(nodeIdx, 1);
        I++;
      }
    }
  }

  // initialize source term
  vector<float> sourceTerm =
      myUtils.computeSourceTerm(num_sample_, dt_, f0, sourceOrder);
  for (int j = 0; j < num_sample_; j++)
  {
    myRHSTerm(0, j) = sourceTerm[j];
    if (j % 100 == 0)
      cout << "Sample " << j << "\t: sourceTerm = " << sourceTerm[j] << endl;
  }

  // get element number of source term
  myElementSource = rhsElement[0];
  cout << "Element number for the source location: " << myElementSource << endl
       << endl;

  int order = m_mesh->getOrder();

  switch (order)
  {
    case 1:
      SourceAndReceiverUtils::ComputeRHSWeights<1>(cornerCoords, src_coord_,
                                                   rhsWeights, 0);
      break;
    case 2:
      SourceAndReceiverUtils::ComputeRHSWeights<2>(cornerCoords, src_coord_,
                                                   rhsWeights, 0);
      break;
    case 3:
      SourceAndReceiverUtils::ComputeRHSWeights<3>(cornerCoords, src_coord_,
                                                   rhsWeights, 0);
      break;
    default:
      throw std::runtime_error("Unsupported order: " + std::to_string(order));
  }

  // preparing every receiver
  for (int i = 0; i < rcvs_size_; i++)
  {
    // Receiver computation
    int receiver_index =
        floor((std::get<0>(rcvs_coord_[i]) * ex) / lx) +
        floor((std::get<1>(rcvs_coord_[i]) * ey) / ly) * ex +
        floor((std::get<2>(rcvs_coord_[i]) * ez) / lz) * ey * ex;

    rhsElementRcv[i] = receiver_index;

    // Get coordinates of the corners of the receiver element
    float cornerCoordsRcv[8][3];
    I = 0u;
    for (int k : nodes_corner)
    {
      for (int j : nodes_corner)
      {
        for (int i : nodes_corner)
        {
          int nodeIdx = m_mesh->globalNodeIndex(rhsElementRcv[i], i, j, k);
          cornerCoordsRcv[I][0] = m_mesh->nodeCoord(nodeIdx, 0);
          cornerCoordsRcv[I][2] = m_mesh->nodeCoord(nodeIdx, 2);
          cornerCoordsRcv[I][1] = m_mesh->nodeCoord(nodeIdx, 1);
          I++;
        }
      }
    }

    switch (order)
    {
      case 1:
        SourceAndReceiverUtils::ComputeRHSWeights<1>(
            cornerCoordsRcv, rcvs_coord_[i], rhsWeightsRcv, i);
        break;
      case 2:
        SourceAndReceiverUtils::ComputeRHSWeights<2>(
            cornerCoordsRcv, rcvs_coord_[i], rhsWeightsRcv, i);
        break;
      case 3:
        SourceAndReceiverUtils::ComputeRHSWeights<3>(
            cornerCoordsRcv, rcvs_coord_[i], rhsWeightsRcv, i);
        break;
      default:
        throw std::runtime_error("Unsupported order: " + std::to_string(order));
    }
  }
}

int SEMproxy::getSliceAxis(string axisArg)
{
  if (axisArg == "X" || axisArg == "x") return 0;
  if (axisArg == "Y" || axisArg == "y") return 1;
  if (axisArg == "Z" || axisArg == "z") return 2;

  throw std::invalid_argument("Slice axis must be X, Y, or Z. Got: " + axisArg);
}

SolverFactory::implemType SEMproxy::getImplem(string implemArg)
{
  if (implemArg == "makutu") return SolverFactory::MAKUTU;
  if (implemArg == "shiva") return SolverFactory::SHIVA;

  throw std::invalid_argument(
      "Implentation type does not follow any valid type.");
}

SolverFactory::meshType SEMproxy::getMesh(string meshArg)
{
  if (meshArg == "cartesian") return SolverFactory::Struct;
  if (meshArg == "ucartesian") return SolverFactory::Unstruct;

  std::cout << "Mesh type found is " << meshArg << std::endl;
  throw std::invalid_argument("Mesh type does not follow any valid type.");
}

SolverFactory::methodType SEMproxy::getMethod(string methodArg)
{
  if (methodArg == "sem") return SolverFactory::SEM;
  if (methodArg == "dg") return SolverFactory::DG;

  throw std::invalid_argument("Method type does not follow any valid type.");
}

float SEMproxy::find_cfl_dt(float cfl_factor)
{
  float sqrtDim3 = 1.73;  // to change for 2d
  float min_spacing = m_mesh->getMinSpacing();
  float v_max = m_mesh->getMaxSpeed();

  float dt = cfl_factor * min_spacing / (sqrtDim3 * v_max);

  return dt;
}

void SEMproxy::save_watched_receivers_output_bin(Measure& metrics)
{
  /*
   * <HEADER>
   * <rcvs_coord_.dump><pnAtReceiver.dump>
   *
   * with <HEADER> being two integers, nb_receivers and
   * nb_samples_per_receiver.
   */
  std::ofstream watchedReceiversOutput(
      watchedReceiversOutputPath,
      std::ios::trunc | std::ios::out | std::ios::binary);
  // first we write the header
  watchedReceiversOutput.write(reinterpret_cast<char*>(&rcvs_size_),
                               sizeof(int));
  watchedReceiversOutput.write(reinterpret_cast<char*>(&num_sample_),
                               sizeof(int));
  // then we dump rcvs_coord_
  watchedReceiversOutput.write(
      reinterpret_cast<char*>(rcvs_coord_.data()),
      sizeof(std::array<float, 3>) * rcvs_coord_.size());
  // and finally the sample array
  watchedReceiversOutput.write(
      reinterpret_cast<char*>(pnAtReceiver.data()),
      sizeof(float) * pnAtReceiver.size());  // pnAtReceiver.size() is the full
                                             // array size, not one single dim
  watchedReceiversOutput.close();
  metrics.measureIO(watchedReceiversOutputPath);
  metrics.getTotalBytes();
}

void SEMproxy::save_watched_receivers_output_plain(Measure& metrics)
{ /*
   * plaintext format will be fairly simply:
   * nb_receivers;nb_samples_per_receiver
   * coords_rcv_1
   * result_rcv_1_1;result_rcv_1_2;result_rcv_1_3(...)result_rcv_1_{nb_samples_per_receiver}
   * coords_rcv_2
   * result_rcv_2_1;result_rcv_2_2;result_rcv_2_3(...)result_rcv_2_{nb_samples_per_receiver}
   * (...)
   * coords_rcv_{nb_receivers}
   * result_rcv_{nb_receivers}_1;result_rcv_{nb_receivers}_2;(...)result_rcv_{nb_receivers}_{nb_samples_per_receiver}
   */
  std::ofstream watchedReceiversOutput(watchedReceiversOutputPath,
                                       std::ios::trunc | std::ios::out);
  watchedReceiversOutput << rcvs_size_ << ";" << num_sample_ << std::endl;

  for (int i = 0; i < rcvs_size_; i++)
  {
    auto rcv_coord = rcvs_coord_[i];
    watchedReceiversOutput << rcv_coord[0] << ";" << rcv_coord[1] << ";"
                           << rcv_coord[2] << std::endl;
    for (int j = 0; j < num_sample_; j++)
    {
      watchedReceiversOutput << pnAtReceiver(i, j);
      if (j + 1 < num_sample_) watchedReceiversOutput << ";";
      // we always add `\n` even if it's the last receiver, as POSIX
      // compliance is the key for an healthy life
    }
  }
  watchedReceiversOutput.close();
  metrics.measureIO(watchedReceiversOutputPath);
  metrics.getTotalBytes();
}

void SEMproxy::computeInSituSnapshotStats(const arrayReal& pressure, int index)
{
  float min_val = std::numeric_limits<float>::max();
  float max_val = std::numeric_limits<float>::lowest();
  double sum = 0.0;

  int num_nodes = m_mesh->getNumberOfNodes();

  for (int i = 0; i < num_nodes; ++i)
  {
    float val = pressure(i, i2);
    if (val < min_val) min_val = val;
    if (val > max_val) max_val = val;
    sum += val;
  }
  float mean_val = static_cast<float>(sum / num_nodes);
  snapshot_stats_.push_back({index, min_val, max_val, mean_val});
}

void SEMproxy::exportInSituStats()
{
  // export sismos
  std::string sismos_path = in_situ_folder_ + "/sismos_data.csv";
  std::ofstream sismos_file(sismos_path);
  sismos_file << "time";
  for (int i = 0; i < rcvs_size_; ++i) sismos_file << ",recv_" << i;
  sismos_file << "\n";

  // pnAtReceiver(rcv, time)
  for (int t = 0; t < num_sample_; ++t)
  {
    float time = t * dt_;
    sismos_file << time;
    for (int i = 0; i < rcvs_size_; ++i)
    {
      sismos_file << "," << pnAtReceiver(i, t);
    }
    sismos_file << "\n";
  }
  sismos_file.close();

  // export FFT
  int n = num_sample_;
  int N = 1;
  while (N < n) N *= 2;

  std::string fft_path = in_situ_folder_ + "/fft_data.csv";
  std::ofstream fft_file(fft_path);
  fft_file << "freq";
  for (int i = 0; i < rcvs_size_; ++i) fft_file << ",recv_" << i;
  fft_file << "\n";

  std::vector<std::vector<double>> fft_magnitudes(rcvs_size_);

  for (int i = 0; i < rcvs_size_; ++i)
  {
    CArray data(N);
    for (int t = 0; t < n; ++t)
    {
      data[t] = Complex(pnAtReceiver(i, t), 0);
    }

    fft(data);

    fft_magnitudes[i].resize(N / 2);
    for (int k = 0; k < N / 2; ++k)
    {
      fft_magnitudes[i][k] = std::abs(data[k]);
    }
  }

  double df = (1.0 / dt_) / N;
  for (int k = 0; k < N / 2; ++k)
  {
    fft_file << k * df;
    for (int i = 0; i < rcvs_size_; ++i)
    {
      fft_file << "," << fft_magnitudes[i][k];
    }
    fft_file << "\n";
  }
  fft_file.close();

  // export
  std::string snap_path = in_situ_folder_ + "/snapshot_stats.csv";
  std::ofstream snap_file(snap_path);
  snap_file << "index,min,max,mean\n";
  for (const auto& s : snapshot_stats_)
  {
    snap_file << s.index << "," << s.min_val << "," << s.max_val << ","
              << s.mean_val << "\n";
  }
  snap_file.close();
}
