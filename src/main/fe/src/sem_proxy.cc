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
#include <cxxopts.hpp>
#include <fstream>
#include <iostream>
#include <ostream>
#include <string>

using namespace SourceAndReceiverUtils;

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
  if (opt.snapshot_folder_path.length() > 0)
  {
    should_snapshot_ = true;
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
      metrics.startClock(MakeSnapshots);
      // create path string
      std::ostringstream stringStream;
      stringStream << snapshot_folder_;
      stringStream << "/snapshot";
      stringStream << indexTimeSample;
      if (snapshot_in_situ_)
      {
        stringStream << ".pgm";
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

      if (!snapshot_in_situ_)
      {
        if (snapshot_format == BIN)
        {
          snapshot_file.write(
              reinterpret_cast<char*>(solverData.m_pnGlobal.data()),
              solverData.m_pnGlobal.size() * sizeof(float));
        }
        else
        {
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
        // lock x for a slice view
        int sliced_dim = 0;  // x
        int slice_pos_along_sliced_dim = domain_size_[sliced_dim] / 2;

        snapshot_file.write("P5\n", 3);
        snapshot_file << nb_nodes_[1] << " " << nb_nodes_[2] << std::endl;
        snapshot_file.write("255\n", 4);

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
        printf("step_x = %f, step_y = %f, step_z = %f\n", _step[0], _step[1], _step[2]);

        char *img = (char*)calloc(nb_nodes_[1] * nb_nodes_[2], 1);

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
            int img_offset = ((int)(global_coords[1] / _step[1])) * nb_nodes_[2] +
                             ((int)(global_coords[2] / _step[2]));
            if (max_pressure > 0.0 &&
                solverData.m_pnGlobal(globalIdx, i2) > 0.0)
            {
              img[img_offset] =
                  (solverData.m_pnGlobal(globalIdx, i2) * 255.0) / max_pressure;
            }
            else
            {
              img[img_offset] = 0;
            }
          }
        }
        snapshot_file.write(img, nb_nodes_[1] * nb_nodes_[2]);
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
