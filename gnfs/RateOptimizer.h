#ifndef RATE_OPTIMIZER_H
#define RATE_OPTIMIZER_H

#include <vector>
#include <deque>
#include <cmath>
#include <algorithm>

// Class to track relation generation rates and optimize sieving parameters
class RateOptimizer
{
public:
    struct RateSample
    {
        double rate;           // Relations per second
        double time;           // Time when sampled
        int relations_count;   // Number of relations in this sample
        
        // Parameters at time of sampling
        long int B1;
        long int B2;
        long int sieve_bound_adj1;
        long int sieve_bound_adj2;
        long int initial_cutoff;
        
        RateSample() : rate(0.0), time(0.0), relations_count(0),
                       B1(0), B2(0), sieve_bound_adj1(0), 
                       sieve_bound_adj2(0), initial_cutoff(0) {}
    };
    
    RateOptimizer(bool enable_auto_tuning = false, size_t history_size = 20)
        : enable_auto_tuning_(enable_auto_tuning),
          history_size_(history_size),
          sample_count_(0),
          best_rate_(0.0),
          stable_count_(0),
          adjustment_count_(0)
    {
    }
    
    // Record a new rate sample
    void record_sample(double rate, double time, int relations_count,
                      long int B1, long int B2, 
                      long int sieve_bound_adj1, long int sieve_bound_adj2,
                      long int initial_cutoff)
    {
        RateSample sample;
        sample.rate = rate;
        sample.time = time;
        sample.relations_count = relations_count;
        sample.B1 = B1;
        sample.B2 = B2;
        sample.sieve_bound_adj1 = sieve_bound_adj1;
        sample.sieve_bound_adj2 = sieve_bound_adj2;
        sample.initial_cutoff = initial_cutoff;
        
        rate_history_.push_back(sample);
        if (rate_history_.size() > history_size_)
        {
            rate_history_.pop_front();
        }
        
        sample_count_++;
        
        if (rate > best_rate_)
        {
            best_rate_ = rate;
            best_sample_ = sample;
        }
    }
    
    // Get moving average of recent rates
    double get_moving_average(size_t window = 5) const
    {
        if (rate_history_.empty()) return 0.0;
        
        size_t n = std::min(window, rate_history_.size());
        double sum = 0.0;
        auto it = rate_history_.rbegin();
        for (size_t i = 0; i < n; i++, ++it)
        {
            sum += it->rate;
        }
        return sum / n;
    }
    
    // Get rate trend (positive = improving, negative = degrading)
    double get_rate_trend(size_t window = 10) const
    {
        if (rate_history_.size() < 2) return 0.0;
        
        size_t n = std::min(window, rate_history_.size());
        if (n < 2) return 0.0;
        
        // Simple linear regression on recent samples
        // Use forward iteration so older samples have lower x values
        double sum_x = 0.0, sum_y = 0.0, sum_xy = 0.0, sum_x2 = 0.0;
        size_t start_idx = rate_history_.size() - n;
        for (size_t i = 0; i < n; i++)
        {
            double x = static_cast<double>(i);
            double y = rate_history_[start_idx + i].rate;
            sum_x += x;
            sum_y += y;
            sum_xy += x * y;
            sum_x2 += x * x;
        }
        
        double slope = (n * sum_xy - sum_x * sum_y) / (n * sum_x2 - sum_x * sum_x);
        return slope;
    }
    
    // Check if rate is stable (low variance)
    bool is_rate_stable(double threshold = 0.1) const
    {
        if (rate_history_.size() < 5) return false;
        
        double avg = get_moving_average(5);
        if (avg < 0.001) return false;
        
        double variance = 0.0;
        auto it = rate_history_.rbegin();
        for (size_t i = 0; i < 5; i++, ++it)
        {
            double diff = it->rate - avg;
            variance += diff * diff;
        }
        variance /= 5;
        
        double std_dev = std::sqrt(variance);
        double coeff_var = std_dev / avg;
        
        return coeff_var < threshold;
    }
    
    // Suggest parameter adjustments
    // Returns true if adjustments are recommended
    bool suggest_adjustments(long int& B1_adj, long int& B2_adj,
                            long int& sieve_bound_adj1, long int& sieve_bound_adj2,
                            long int& initial_cutoff_adj)
    {
        if (!enable_auto_tuning_) return false;
        if (rate_history_.size() < 10) return false;  // Need enough data
        
        // Don't adjust too frequently
        if (sample_count_ - adjustment_count_ < 5) return false;
        
        double trend = get_rate_trend(10);
        double current_rate = rate_history_.back().rate;
        double avg_rate = get_moving_average(5);
        
        // If rate is improving and stable, don't change
        if (trend > 0 && is_rate_stable(0.15))
        {
            stable_count_++;
            if (stable_count_ > 3)
            {
                return false;  // Keep current parameters
            }
        }
        else
        {
            stable_count_ = 0;
        }
        
        // If rate is degrading significantly, try adjustments
        if (trend < -0.5 || current_rate < best_rate_ * 0.8)
        {
            // Try to move towards best known parameters
            const RateSample& current = rate_history_.back();
            
            // Gradual adjustment (10% steps)
            double factor = 0.1;
            B1_adj = static_cast<long int>((best_sample_.B1 - current.B1) * factor);
            B2_adj = static_cast<long int>((best_sample_.B2 - current.B2) * factor);
            sieve_bound_adj1 = static_cast<long int>((best_sample_.sieve_bound_adj1 - current.sieve_bound_adj1) * factor);
            sieve_bound_adj2 = static_cast<long int>((best_sample_.sieve_bound_adj2 - current.sieve_bound_adj2) * factor);
            initial_cutoff_adj = static_cast<long int>((best_sample_.initial_cutoff - current.initial_cutoff) * factor);
            
            // If no best sample yet or very close, try exploration
            if (best_rate_ < 0.001 || std::abs(B1_adj) < 100)
            {
                // Small exploratory changes
                B1_adj = (sample_count_ % 2 == 0) ? 1000 : -1000;
                sieve_bound_adj1 = (sample_count_ % 3 == 0) ? 1 : -1;
                initial_cutoff_adj = (sample_count_ % 4 == 0) ? 5 : -5;
            }
            
            adjustment_count_ = sample_count_;
            return true;
        }
        
        return false;
    }
    
    // Get statistics for reporting
    double get_best_rate() const { return best_rate_; }
    size_t get_sample_count() const { return sample_count_; }
    const RateSample& get_best_sample() const { return best_sample_; }
    const std::deque<RateSample>& get_history() const { return rate_history_; }
    
    // Enable/disable auto-tuning
    void set_auto_tuning(bool enable) { enable_auto_tuning_ = enable; }
    bool is_auto_tuning_enabled() const { return enable_auto_tuning_; }
    
private:
    bool enable_auto_tuning_;
    size_t history_size_;
    std::deque<RateSample> rate_history_;
    size_t sample_count_;
    double best_rate_;
    RateSample best_sample_;
    size_t stable_count_;
    size_t adjustment_count_;
};

#endif // RATE_OPTIMIZER_H
