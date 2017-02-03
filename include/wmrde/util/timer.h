#ifndef _WMRDE_TIMER_H_
#define _WMRDE_TIMER_H_

//Linux only includes:
#include <sys/time.h> //for gettimeofday
#include <time.h> //for nanosleep

namespace wmrde
{

//TODO, remove this?
class TimedWait
{
public:
  TimedWait(const double sec)
  {
    struct timespec tim, tim2;
    double floor_sec = floor(sec);
    tim.tv_sec = floor_sec;
    tim.tv_nsec = (sec - floor_sec)*1e9;
    if(nanosleep(&tim , &tim2) < 0 )
    {
       printf("nanosleep system call failed!\n");
       return;
    }
    printf("Waited for %f sec.\n", sec);
  }
};

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
