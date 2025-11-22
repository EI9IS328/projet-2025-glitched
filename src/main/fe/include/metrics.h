#include <ostream>
#ifndef METRICS_HPP_
#include <chrono>
#include <cstddef>
#include <map>
#include <optional>
#include <string>
#define METRICS_HPP_

enum MeasurePoint
{
  Global,
  Kernel,
  MakeSnapshots,
  MakeSismos,
  OutputSismos,
};

class Metrics
{
 public:
  Metrics();
  void measureIO(std::string closedFilePath);
  void startClock(MeasurePoint measuredPoint);
  void stopClockAndAppend(MeasurePoint measuredPoint);
  float getTimeMs(MeasurePoint measuredPoint);
  std::size_t getTotalBytes();
  const std::map<std::string, std::size_t> getDetailedBytes();
  friend std::ostream& operator<<(std::ostream& os, Metrics& metrics);

 private:
  std::map<MeasurePoint, std::chrono::time_point<std::chrono::system_clock>>
      m_clock_measures;
  std::map<MeasurePoint,
           std::optional<std::chrono::time_point<std::chrono::system_clock>>>
      m_running_clocks;
  std::map<std::string, std::size_t> m_io_measures;
};

#endif
