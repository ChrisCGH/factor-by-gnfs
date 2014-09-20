#ifndef TIMINGS_H
#define TIMINGS_H

#ifndef WIN32
#include <sys/times.h>
#include <unistd.h>
#endif
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <fstream>
#include <string>
#include <map>

extern "C"
{
   void timing_file(const char* filename);
   void timing_start(const char* str);
   void timing_stop();
   void timing_summary();
}
class Timing
{
   public:
      Timing(const char* filename, bool summary_only = true);
      ~Timing();
      void start(const char* str);
      void stop();
      void summary();
      void reset();

   private:
      static int clk_tck_;
      static char* time_str();
      std::fstream* timing_file_;
      bool summary_only_;
      struct Timing_
      {
         Timing_();
#ifndef WIN32
         struct tms start_tms_;
         struct tms stop_tms_;
#endif
         clock_t start_time_;
         clock_t stop_time_;
         double elapsed_real_;
         double elapsed_user_;
         double elapsed_system_;
         double total_elapsed_real_;
         double total_elapsed_user_;
         double total_elapsed_system_;
         std::string message_;
      };
      std::map<std::string, Timing_*> timing_map_;
      Timing_* timing_;
      double grand_total_elapsed_real_;
      double grand_total_elapsed_user_;
      double grand_total_elapsed_system_;
};

class Timer
{
public:
	Timer();
	~Timer();
	void start();
	double stop();
	double total_elapsed() const;
    static void reset();
private:
    static int clk_tck_;
#ifndef WIN32
    struct tms start_tms_;
    struct tms stop_tms_;
#endif
    clock_t start_time_;
    clock_t stop_time_;
    double elapsed_real_;
    double elapsed_user_;
    double elapsed_system_;
	static time_t creation_time_;
};
#endif
