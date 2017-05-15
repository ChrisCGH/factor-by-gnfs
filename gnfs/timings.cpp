#include "timings.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <unordered_map>
#ifndef WIN32
#include <sys/times.h>
#include <unistd.h>
#endif
#include <stdio.h>
#include <time.h>
#include <string.h>

#ifndef WIN32
static struct tms Start_tms;
static struct tms Stop_tms;
#endif
static clock_t Start_time = 0;
static clock_t Stop_time = 0;
static int Clk_tck = 0;

static char Message[1024];
static FILE* Timing_file = 0;

static double Elapsed_real = 0;
static double Elapsed_user = 0.0;
static double Elapsed_system = 0.0;

static double Total_elapsed_real = 0;
static double Total_elapsed_user = 0.0;
static double Total_elapsed_system = 0.0;

static char* time_str()
{
   time_t the_time = time(0);
   static char the_time_str[132];
   strftime(the_time_str, 132, "%b %d %H:%M:%S", localtime(&the_time));
   return the_time_str;
}

void timing_file(const char* filename)
{
   Timing_file = fopen(filename, "a");
}

void timing_start(const char* str)
{
#ifndef WIN32
   if (Clk_tck == 0)
   {
      Clk_tck = sysconf(_SC_CLK_TCK);
   }

   strcpy(Message, str);

   Start_time = times(&Start_tms);
#else
   Clk_tck = CLOCKS_PER_SEC;
   Start_time = clock();
#endif
}

void timing_stop()
{
#ifndef WIN32
   Stop_time = times(&Stop_tms);

   Elapsed_real = (double)(Stop_time - Start_time)/(double)Clk_tck;
   Elapsed_user = (double)(Stop_tms.tms_utime - Start_tms.tms_utime)/(double)Clk_tck;
   Elapsed_system = (double)(Stop_tms.tms_stime - Start_tms.tms_stime)/(double)Clk_tck;

   Total_elapsed_real = Total_elapsed_real + Elapsed_real;
   Total_elapsed_user = Total_elapsed_user + Elapsed_user;
   Total_elapsed_system = Total_elapsed_system + Elapsed_system;

   if (Timing_file)
   {
      fprintf(Timing_file, "%s,%s,Real=%4.2f,User=%4.2f,System=%4.2f\n",
              time_str(), Message, Elapsed_real, Elapsed_user, Elapsed_system);
      fflush(Timing_file);
   }
   else
   {
      fprintf(stderr, "%s,%s,Real=%4.2f,User=%4.2f,System=%4.2f\n",
              time_str(), Message, Elapsed_real, Elapsed_user, Elapsed_system);
      fflush(stderr);
   }
#else
   Stop_time = clock();

   Elapsed_real = (double)(Stop_time - Start_time)/(double)Clk_tck;
   Elapsed_user = 0.0;
   Elapsed_system = 0.0;

   Total_elapsed_real = Total_elapsed_real + Elapsed_real;
   Total_elapsed_user = Total_elapsed_user + Elapsed_user;
   Total_elapsed_system = Total_elapsed_system + Elapsed_system;

   if (Timing_file)
   {
      fprintf(Timing_file, "%s,%s,Real=%4.2f,User=%4.2f,System=%4.2f\n",
              time_str(), Message, Elapsed_real, Elapsed_user, Elapsed_system);
      fflush(Timing_file);
   }
   else
   {
      fprintf(stderr, "%s,%s,Real=%4.2f,User=%4.2f,System=%4.2f\n",
              time_str(), Message, Elapsed_real, Elapsed_user, Elapsed_system);
      fflush(stderr);
   }
#endif
}

void timing_summary()
{
   if (Timing_file)
   {
      fprintf(Timing_file, "%s,Total,Real=%4.2f,User=%4.2f,System=%4.2f\n",
              time_str(), Total_elapsed_real, Total_elapsed_user, Total_elapsed_system);
      fflush(Timing_file);
   }
   else
   {
      fprintf(stderr, "%s,Total,Real=%4.2f,User=%4.2f,System=%4.2f\n",
              time_str(), Total_elapsed_real, Total_elapsed_user, Total_elapsed_system);
      fflush(stderr);
   }
}

int Timer::clk_tck_ = 0;
time_t Timer::creation_time_ = time(0);

Timer::Timer()
: elapsed_real_(0.0), elapsed_user_(0.0), elapsed_system_(0.0)
{
#if defined(WIN32) || defined(linux)
   if (clk_tck_ == 0)
   {
      clk_tck_ = CLOCKS_PER_SEC;
   }
#else
   if (clk_tck_ == 0)
   {
      clk_tck_ = sysconf(_SC_CLK_TCK);
   }
#endif
}

Timer::~Timer()
{
}

void Timer::start()
{
#ifndef WIN32
   //struct tms stms;
   //clock_t st = times(&stms);
   clock_t st = clock();
#else
   clock_t st = clock();
#endif
   start_time_ = st;
#ifndef WIN32
   start_time_ = st;
#endif
}

double Timer::stop()
{
#ifndef WIN32
   //stop_time_ = times(&stop_tms_);
   stop_time_ = clock();
   elapsed_real_ = (double)(stop_time_ - start_time_);
   //std::cout << "stop_time_ = " << stop_time_ << ", start_time_ = " << start_time_ << ", elapsed_real_ = " << elapsed_real_ << std::endl;
   elapsed_user_ = (double)(stop_tms_.tms_utime - start_tms_.tms_utime);
   elapsed_system_ = (double)(stop_tms_.tms_stime - start_tms_.tms_stime);
#else
   stop_time_ = clock();
   elapsed_real_ = (double)(stop_time_ - start_time_);
   elapsed_user_ = 0.0;
   elapsed_system_ = 0.0;
#endif
   return elapsed_real_ / (double)clk_tck_;
}

double Timer::total_elapsed() const
{
    time_t now = time(0);
	return ((double)now - (double)creation_time_);
}

void Timer::reset()
{
    creation_time_ = time(0);
}

int Timing::clk_tck_ = 0;

Timing::Timing(const char* filename, bool summary_only) : timing_file_(0), summary_only_(summary_only)
{
   timing_file_ = new std::fstream(filename, std::ios::out|std::ios::app);
   timing_map_.clear();
#if defined(WIN32) || defined(linux)
   if (clk_tck_ == 0)
   {
      clk_tck_ = CLOCKS_PER_SEC;
   }
#else
   if (clk_tck_ == 0)
   {
      clk_tck_ = sysconf(_SC_CLK_TCK);
   }
#endif
   timing_ = 0;
   grand_total_elapsed_real_ = 0.0;
   grand_total_elapsed_user_ = 0.0;
   grand_total_elapsed_system_ = 0.0;
}

Timing::~Timing()
{
   summary();
   reset();
   delete timing_file_;
}

char* Timing::time_str()
{
   time_t the_time = time(0);
   static char the_time_str[132];
   strftime(the_time_str, 132, "%b %d %H:%M:%S", localtime(&the_time));
   return the_time_str;
}

void Timing::start(const char* str)
{
#if defined(WIN32) || defined(linux)
   clock_t st = clock();
#else
   struct tms stms;
   clock_t st = times(&stms);
#endif
   std::string s(str);
   if (timing_map_.find(s) == timing_map_.end())
   {
      timing_map_[s] = new Timing_;
      timing_map_[s]->message_ = s;
   }
   timing_ = timing_map_[s];
   timing_->start_time_ = st;
#if defined(WIN32) || defined(linux)
#else
   timing_->start_tms_ = stms;
#endif
}

void Timing::stop()
{
#if defined(WIN32) || defined(linux)
   timing_->stop_time_ = clock();
   //timing_->elapsed_real_ = (double)(timing_->stop_time_ - timing_->start_time_)/(double)clk_tck_;
   timing_->elapsed_real_ = (double)(timing_->stop_time_ - timing_->start_time_);
   timing_->elapsed_user_ = 0.0;
   timing_->elapsed_system_ = 0.0;

   timing_->total_elapsed_real_ = timing_->total_elapsed_real_ + timing_->elapsed_real_;
   timing_->total_elapsed_user_ = timing_->total_elapsed_user_ + timing_->elapsed_user_;
   timing_->total_elapsed_system_ = timing_->total_elapsed_system_ + timing_->elapsed_system_;
   grand_total_elapsed_real_ = grand_total_elapsed_real_ + timing_->elapsed_real_;
   grand_total_elapsed_user_ = grand_total_elapsed_user_ + timing_->elapsed_user_;
   grand_total_elapsed_system_ = grand_total_elapsed_system_ + timing_->elapsed_system_;
#else
   timing_->stop_time_ = times(&timing_->stop_tms_);
   timing_->elapsed_real_ = (double)(timing_->stop_time_ - timing_->start_time_);
   timing_->elapsed_user_ = (double)(timing_->stop_tms_.tms_utime - timing_->start_tms_.tms_utime);
   timing_->elapsed_system_ = (double)(timing_->stop_tms_.tms_stime - timing_->start_tms_.tms_stime);

   timing_->total_elapsed_real_ = timing_->total_elapsed_real_ + timing_->elapsed_real_;
   timing_->total_elapsed_user_ = timing_->total_elapsed_user_ + timing_->elapsed_user_;
   timing_->total_elapsed_system_ = timing_->total_elapsed_system_ + timing_->elapsed_system_;
   grand_total_elapsed_real_ = grand_total_elapsed_real_ + timing_->elapsed_real_;
   grand_total_elapsed_user_ = grand_total_elapsed_user_ + timing_->elapsed_user_;
   grand_total_elapsed_system_ = grand_total_elapsed_system_ + timing_->elapsed_system_;
#endif
   if (summary_only_) return;
   if (timing_file_)
   {
      *timing_file_ << time_str() << "," << timing_->message_ << ",Real=" << std::setprecision(6) << timing_->elapsed_real_/(double)clk_tck_ << ",User=" << timing_->elapsed_user_/(double)clk_tck_ << ",System=" << timing_->elapsed_system_/(double)clk_tck_ << std::endl;
   }
   else
   {
      std::cerr << time_str() << "," << timing_->message_ << ",Real=" << std::setprecision(6) << timing_->elapsed_real_/(double)clk_tck_ << ",User=" << timing_->elapsed_user_/(double)clk_tck_ << ",System=" << timing_->elapsed_system_/(double)clk_tck_ << std::endl;
   }
}

void Timing::summary()
{
   Timing_* save = timing_;
   std::unordered_multimap<double, std::string> summary_map;
   for (auto& t: timing_map_)
   {
      timing_ = t.second;
      std::stringstream ss;
      double pct = 0.0;
      ss << time_str() << "," << timing_->message_ << ",Total,Real=" << std::setprecision(6)
      << timing_->total_elapsed_real_/(double)clk_tck_;
#if defined(WIN32) || defined(linux)
      pct = (100.0 * timing_->total_elapsed_real_/(double)clk_tck_) / (grand_total_elapsed_real_/(double)clk_tck_);
      ss << " (" << pct << "%),User=" << timing_->total_elapsed_user_/(double)clk_tck_ << ",System=";
#else
      pct = (100.0 * timing_->total_elapsed_user_/(double)clk_tck_) / (grand_total_elapsed_user_/(double)clk_tck_);
      ss << ",User=" << timing_->total_elapsed_user_/(double)clk_tck_ << " (" << pct << "%),System=";
#endif
      ss << timing_->total_elapsed_system_/(double)clk_tck_;
      summary_map.insert(std::unordered_multimap<double, std::string>::value_type(pct, ss.str()));
   }

   if (timing_file_)
   {
      *timing_file_ << time_str() << "================ Summary ==================" << std::endl;
   }
   else
   {
      std::cerr << time_str() << "================ Summary ==================" << std::endl;
   }
   for (auto& s: summary_map)
   {
      if (timing_file_)
      {
         *timing_file_ << s.second << std::endl;
      }
      else
      {
         std::cerr << s.second << std::endl;
      }
   }
   if (timing_file_)
   {
      *timing_file_ << time_str() << ",Grand Total,Real=" << std::setprecision(6) << grand_total_elapsed_real_/(double)clk_tck_ << ",User=" << grand_total_elapsed_user_/(double)clk_tck_ << ",System=" << grand_total_elapsed_system_/(double)clk_tck_ << std::endl;
   }
   else
   {
      std::cerr << time_str() << ",Grand Total,Real=" << std::setprecision(6) << grand_total_elapsed_real_/(double)clk_tck_ << ",User=" << grand_total_elapsed_user_/(double)clk_tck_ << ",System=" << grand_total_elapsed_system_/(double)clk_tck_ << std::endl;
   }
   timing_ = save;
}

void Timing::reset()
{
   for (auto& t: timing_map_)
   {
      timing_ = t.second;
      delete timing_;
   }
   timing_map_.clear();
}

Timing::Timing_::Timing_() :
      start_time_(0), stop_time_(0), elapsed_real_(0), elapsed_user_(0), elapsed_system_(0),
      total_elapsed_real_(0), total_elapsed_user_(0), total_elapsed_system_(0), message_("")
{
#if defined(WIN32) || defined(linux)
#else
   start_tms_.tms_utime = 0;
   start_tms_.tms_stime = 0;
   stop_tms_.tms_utime = 0;
   stop_tms_.tms_stime = 0;
#endif
}
