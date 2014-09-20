#ifndef LOGGER_H
#define LOGGER_H
#include <string>
#include <sstream>
#include <vector>

class Logger
{
    public:
        enum level { debug = 0, info, warning, error };

        Logger(const std::string& logfile = "");
        virtual ~Logger();
        void set_level(level l);
        void log(const std::string& msg, level l);

    private:
        std::string logfile_;
        level level_;
        bool use_stderr_;
        std::ostream* logstream_;
};

class LogManager
{
    public:
        static LogManager& instance();
        typedef int log_id_type;
        log_id_type start_logging(const std::string& logfile, bool switch_logging = true);
        void switch_logging(log_id_type log_id);
        void stop_logging(log_id_type log_id);
        void log(const std::string& msg, Logger::level l);
        void set_level(Logger::level l);

    private:
        LogManager();
        ~LogManager();
        log_id_type current_logger_;
        std::vector<Logger*> logger_list_;
};

#define LOG_ERROR(msg) \
{ \
    std::ostringstream oss; \
    oss << msg << " : " << __FUNCTION__ << " " << __FILE__ << ":" << __LINE__; \
    LogManager::instance().log(oss.str(), Logger::error); \
}  
#define LOG_DEBUG(msg) \
{ \
    std::ostringstream oss; \
    oss << msg << " : " << __FUNCTION__ << " " << __FILE__ << ":" << __LINE__; \
    LogManager::instance().log(oss.str(), Logger::debug); \
}  
#endif
