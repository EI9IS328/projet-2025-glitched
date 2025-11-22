#include <metrics.h>

#include <chrono>
#include <cstddef>
#include <filesystem>
#include <iostream>
#include <optional>
#include <sstream>
#include <stdexcept>

#define TODO(msg) (throw std::runtime_error(msg " not implemented"));

Metrics::Metrics()
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

void Metrics::measureIO(std::string closedFilePath)
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

void Metrics::startClock(MeasurePoint measuredPoint)
{
  if (!m_running_clocks[measuredPoint].has_value())
  {
    m_running_clocks[measuredPoint] = std::chrono::system_clock::now();
  }
}

void Metrics::stopClockAndAppend(MeasurePoint measuredPoint)
{
  if (m_running_clocks[measuredPoint].has_value())
  {
    m_clock_measures[measuredPoint] += std::chrono::system_clock::now() -
                                       m_running_clocks[measuredPoint].value();
    m_running_clocks[measuredPoint] = std::nullopt;
  }
}

float Metrics::getTimeMs(MeasurePoint measuredPoint)
{
  return std::chrono::time_point_cast<std::chrono::microseconds>(
             m_clock_measures[measuredPoint])
      .time_since_epoch()
      .count();
}

const std::map<std::string, std::size_t> Metrics::getDetailedBytes()
{
  return std::map<std::string, std::size_t>(m_io_measures);
}

std::size_t Metrics::getTotalBytes()
{
  std::cout << " wtf " << m_io_measures.size() << std::endl;
  std::size_t ttl = 0;
  for (auto kv : m_io_measures)
  {
    ttl += kv.second;
  }
  return ttl;
}
