#include "Logger.h"
#include <string>
#include <iostream>
#include <fstream>
#include <ctime>

namespace
{
std::string timestamp()
{
    char buf[100];
    std::time_t t = std::time(0);
    std::strftime(buf, sizeof(buf), "%d %b %Y %H:%M:%S", std::localtime(&t));
    return std::string(buf);
}

std::string level_to_string(Logger::level l)
{
    switch (l)
    {
        case Logger::error:
            return "ERROR  ";
        case Logger::warning:
            return "WARNING";
        case Logger::info:
            return "INFO   ";
        case Logger::debug:
            return "DEBUG  ";
    }
    return "NOTICE ";
}
}

Logger::Logger(const std::string& logfile)
    : logfile_(logfile), level_(error), use_stderr_(false), logstream_(0)
{
    if (logfile.empty())
    {
        // use stderr
        use_stderr_ = true;
        logstream_ = &std::cerr;
    }
    else
    {
        logstream_ = new std::fstream(logfile.c_str(), std::ios::out | std::ios::app);
    }
}

Logger::~Logger()
{
    if (!use_stderr_)
    {
        std::fstream* f = dynamic_cast<std::fstream*>(logstream_);
        if (f)
        {
            f->close();
        }
        delete logstream_;
    }
}

void Logger::set_level(level l)
{
    level_ = l;
}

void Logger::log(const std::string& msg, level l)
{
    if (l >= level_)
    {
        *logstream_ << timestamp() << " : " << level_to_string(l) << " : " << msg << std::endl;
    }
}

LogManager::LogManager()
    : current_logger_(-1)
{
}

LogManager::~LogManager()
{
    for (auto& i: logger_list_)
    {
        delete i;
    }
    logger_list_.clear();
}

LogManager& LogManager::instance()
{
    static LogManager log_manager;
    return log_manager;
}

LogManager::log_id_type LogManager::start_logging(const std::string& logfile, bool switch_logging)
{
    Logger* logger = new Logger(logfile);
    log_id_type log_id = 0;
    for (auto it = logger_list_.begin();
         it != logger_list_.end();
         ++it, ++log_id)
    {
        if (*it == 0)
        {
            *it = logger;
        }
    }
    if (log_id >= (log_id_type)logger_list_.size())
    {
        logger_list_.push_back(logger);
    }
    current_logger_ = log_id;
    return log_id;
}

void LogManager::switch_logging(LogManager::log_id_type log_id)
{
    if (log_id >= 0 && log_id < (log_id_type)logger_list_.size() && logger_list_[log_id])
    {
        current_logger_ = log_id;
    }
}

void LogManager::stop_logging(LogManager::log_id_type log_id)
{
    delete logger_list_[log_id];
    logger_list_[log_id] = 0;
    if (current_logger_ == log_id)
    {
        log_id_type id = 0;
        for (auto it = logger_list_.begin();
             it != logger_list_.end();
             ++it, ++log_id)
        {
            if (*it)
            {
               current_logger_ = id;
               return;
            }
        }
        current_logger_ = -1;
    }
}

void LogManager::log(const std::string& msg, Logger::level l)
{
    if (current_logger_ >= 0 && current_logger_ < (log_id_type)logger_list_.size())
    {
        Logger* logger = logger_list_[current_logger_]; 
        logger->log(msg, l);
    }
}

void LogManager::set_level(Logger::level l)
{
    if (current_logger_ >= 0 && current_logger_ < (log_id_type)logger_list_.size())
    {
        Logger* logger = logger_list_[current_logger_];
        logger->set_level(l);
    }
}
