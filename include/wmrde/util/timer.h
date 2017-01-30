#ifndef _WMRDE_TIMER_H_
#define _WMRDE_TIMER_H_

#include <sys/time.h> //Linux only

namespace wmrde
{
class Timer
{
private:
  timeval start_time_;
  timeval stop_time_;

public:
  void start() { gettimeofday(&start_time_, nullptr); }
  void stop() { gettimeofday(&stop_time_, nullptr); }

  double elapsedTimeSec() const
  {
    return (stop_time_.tv_sec - start_time_.tv_sec) +
        (stop_time_.tv_usec - start_time_.tv_usec)/1000000.0;
  }

  double elapsedTimeMs() const
  {
    return (stop_time_.tv_sec - start_time_.tv_sec)*1000.0 +
        (stop_time_.tv_usec - start_time_.tv_usec)/1000.0;
  }
};

} //namespace

#endif
