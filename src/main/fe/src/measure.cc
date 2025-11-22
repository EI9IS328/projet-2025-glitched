#include <measure.h>

#include <cstddef>
#include <filesystem>
#include <iostream>
#include <ostream>
#include <sstream>
#include <stdexcept>

Measure::Measure()
{
  m_clock_measures =
      std::map<MeasurePoint,
               std::chrono::time_point<std::chrono::system_clock>>();
  m_io_measures = std::map<std::string, std::size_t>();
  m_running_clocks = std::map<
      MeasurePoint,
      std::optional<std::chrono::time_point<std::chrono::system_clock>>>();

  m_clock_measures.insert(
      {Global, std::chrono::time_point<std::chrono::system_clock>()});
  m_clock_measures.insert(
      {Kernel, std::chrono::time_point<std::chrono::system_clock>()});
  m_clock_measures.insert(
      {MakeSnapshots, std::chrono::time_point<std::chrono::system_clock>()});
  m_clock_measures.insert(
      {MakeSismos, std::chrono::time_point<std::chrono::system_clock>()});
  m_clock_measures.insert(
      {OutputSismos, std::chrono::time_point<std::chrono::system_clock>()});

  m_running_clocks.insert({Global, std::nullopt});
  m_running_clocks.insert({Kernel, std::nullopt});
  m_running_clocks.insert({MakeSnapshots, std::nullopt});
  m_running_clocks.insert({MakeSismos, std::nullopt});
  m_running_clocks.insert({OutputSismos, std::nullopt});
}

void Measure::measureIO(std::string closedFilePath)
{
  if (!std::filesystem::exists(closedFilePath))
  {
    std::stringstream message;
    message << "couldn't find file " << closedFilePath;
    throw std::runtime_error(message.str());
  }
  m_io_measures.insert(
      {closedFilePath, std::filesystem::file_size(closedFilePath)});
}

void Measure::startClock(MeasurePoint measuredPoint)
{
  if (!m_running_clocks[measuredPoint].has_value())
  {
    m_running_clocks[measuredPoint] = std::chrono::system_clock::now();
  }
}

void Measure::stopClockAndAppend(MeasurePoint measuredPoint)
{
  if (m_running_clocks[measuredPoint].has_value())
  {
    m_clock_measures[measuredPoint] += std::chrono::system_clock::now() -
                                       m_running_clocks[measuredPoint].value();
    m_running_clocks[measuredPoint] = std::nullopt;
  }
}

float Measure::getTimeMs(MeasurePoint measuredPoint)
{
  return std::chrono::time_point_cast<std::chrono::microseconds>(
             m_clock_measures[measuredPoint])
      .time_since_epoch()
      .count();
}

const std::map<std::string, std::size_t> Measure::getDetailedBytes()
{
  return std::map<std::string, std::size_t>(m_io_measures);
}

std::size_t Measure::getTotalBytes()
{
  std::size_t ttl = 0;
  for (auto kv : m_io_measures)
  {
    ttl += kv.second;
  }
  return ttl;
}

std::ostream& operator<<(std::ostream& os, Measure& metrics)
{
  float kerneltime_ms = metrics.getTimeMs(Global);
  float simtime_ms = metrics.getTimeMs(Kernel);
  float make_snapshots_ms = metrics.getTimeMs(MakeSnapshots);
  float make_sismos_ms = metrics.getTimeMs(MakeSismos);
  float make_output_sismos = metrics.getTimeMs(OutputSismos);
  float totalBytes = metrics.getTotalBytes();
  auto detailedBytes = metrics.getDetailedBytes();

  os << "------------------------------------------------" << std::endl;
  os << "---- Time Kernel Total : " << kerneltime_ms / 1E6 << " seconds."
     << std::endl;
  os << "------------------------------------------------" << std::endl
     << std::endl;

  os << "------------------------------------------------" << std::endl;
  os << "---- Time Spent Simulating : " << simtime_ms / 1E6 << " seconds."
     << std::endl;
  os << "------------------------------------------------" << std::endl
     << std::endl;

  os << "------------------------------------------------" << std::endl;
  os << "---- Time Making and Saving Snapshots : " << make_snapshots_ms / 1E6
     << " seconds." << std::endl;
  os << "------------------------------------------------" << std::endl
     << std::endl;

  os << "------------------------------------------------" << std::endl;
  os << "---- Time Making Sismos : " << make_sismos_ms / 1E6 << " seconds."
     << std::endl;
  os << "------------------------------------------------" << std::endl
     << std::endl;

  os << "------------------------------------------------" << std::endl;
  os << "---- Time Saving Outputs : " << make_output_sismos / 1E6 << " seconds."
     << std::endl;
  os << "------------------------------------------------" << std::endl
     << std::endl;

  os << "------------------------------------------------" << std::endl;
  os << "---- Total written data: " << totalBytes << " Bytes." << std::endl;
  os << "------------------------------------------------" << std::endl
     << std::endl;

  os << "------------------------------------------------" << std::endl;
  os << "---- Detail:" << std::endl;
  os << "------------------------------------------------" << std::endl;
  for (auto kv : detailedBytes)
  {
    auto key = kv.first;
    auto value = kv.second;

    os << "--- " << key << ": " << value << " Bytes" << std::endl;
  }
  return os;
}
