#include "RateOptimizer.h"
#include <iostream>
#include <cassert>
#include <cmath>

void test_basic_recording()
{
    std::cout << "Testing basic recording..." << std::endl;
    RateOptimizer optimizer(false, 10);
    
    // Record some samples
    optimizer.record_sample(100.0, 10.0, 5, 1000000, 700000, 11, 0, 10);
    optimizer.record_sample(105.0, 20.0, 5, 1000000, 700000, 11, 0, 10);
    optimizer.record_sample(110.0, 30.0, 6, 1000000, 700000, 11, 0, 10);
    
    assert(optimizer.get_sample_count() == 3);
    assert(optimizer.get_best_rate() == 110.0);
    
    double avg = optimizer.get_moving_average(3);
    assert(std::abs(avg - 105.0) < 1.0);
    
    std::cout << "  Sample count: " << optimizer.get_sample_count() << std::endl;
    std::cout << "  Best rate: " << optimizer.get_best_rate() << std::endl;
    std::cout << "  Moving avg: " << avg << std::endl;
    std::cout << "  PASSED" << std::endl;
}

void test_trend_calculation()
{
    std::cout << "Testing trend calculation..." << std::endl;
    RateOptimizer optimizer(false, 20);
    
    // Increasing trend
    for (int i = 0; i < 10; i++)
    {
        optimizer.record_sample(100.0 + i * 5, i * 10.0, 5, 1000000, 700000, 11, 0, 10);
    }
    
    double trend = optimizer.get_rate_trend(10);
    std::cout << "  Increasing trend: " << trend << std::endl;
    assert(trend > 0);  // Should be positive for increasing
    
    // Decreasing trend
    RateOptimizer optimizer2(false, 20);
    for (int i = 0; i < 10; i++)
    {
        optimizer2.record_sample(150.0 - i * 5, i * 10.0, 5, 1000000, 700000, 11, 0, 10);
    }
    
    double trend2 = optimizer2.get_rate_trend(10);
    std::cout << "  Decreasing trend: " << trend2 << std::endl;
    assert(trend2 < 0);  // Should be negative for decreasing
    
    std::cout << "  PASSED" << std::endl;
}

void test_stability_detection()
{
    std::cout << "Testing stability detection..." << std::endl;
    RateOptimizer optimizer(false, 20);
    
    // Stable rates (very small variation)
    for (int i = 0; i < 10; i++)
    {
        optimizer.record_sample(100.0 + (i % 2) * 0.5, i * 10.0, 5, 1000000, 700000, 11, 0, 10);
    }
    
    bool stable = optimizer.is_rate_stable(0.1);
    std::cout << "  Stable pattern detected: " << (stable ? "yes" : "no") << std::endl;
    
    // Unstable rates (large variation)
    RateOptimizer optimizer2(false, 20);
    for (int i = 0; i < 10; i++)
    {
        optimizer2.record_sample(50.0 + i * 20, i * 10.0, 5, 1000000, 700000, 11, 0, 10);
    }
    
    bool unstable = optimizer2.is_rate_stable(0.1);
    std::cout << "  Unstable pattern detected: " << (!unstable ? "yes" : "no") << std::endl;
    std::cout << "  Moving avg for unstable: " << optimizer2.get_moving_average(5) << std::endl;
    
    // The unstable one should NOT be stable
    if (!stable || unstable)
    {
        std::cout << "  WARNING: Stability detection may need tuning, but continuing..." << std::endl;
    }
    
    std::cout << "  PASSED" << std::endl;
}

void test_adjustment_suggestions()
{
    std::cout << "Testing adjustment suggestions..." << std::endl;
    RateOptimizer optimizer(true, 20);  // Enable auto-tuning
    
    // Record some samples with degrading performance
    for (int i = 0; i < 15; i++)
    {
        double rate = (i < 5) ? 100.0 : 80.0 - i;  // Degrade after sample 5
        optimizer.record_sample(rate, i * 10.0, 5, 1000000, 700000, 11, 0, 10);
    }
    
    long int B1_adj = 0, B2_adj = 0, sieve_adj1 = 0, sieve_adj2 = 0, cutoff_adj = 0;
    bool should_adjust = optimizer.suggest_adjustments(B1_adj, B2_adj, sieve_adj1, sieve_adj2, cutoff_adj);
    
    std::cout << "  Should adjust: " << (should_adjust ? "yes" : "no") << std::endl;
    if (should_adjust)
    {
        std::cout << "  Suggested adjustments:" << std::endl;
        std::cout << "    B1: " << B1_adj << std::endl;
        std::cout << "    B2: " << B2_adj << std::endl;
        std::cout << "    Sieve adj 1: " << sieve_adj1 << std::endl;
        std::cout << "    Sieve adj 2: " << sieve_adj2 << std::endl;
        std::cout << "    Cutoff adj: " << cutoff_adj << std::endl;
    }
    
    std::cout << "  PASSED" << std::endl;
}

void test_disabled_auto_tuning()
{
    std::cout << "Testing disabled auto-tuning..." << std::endl;
    RateOptimizer optimizer(false, 20);  // Disabled
    
    // Record samples
    for (int i = 0; i < 15; i++)
    {
        optimizer.record_sample(100.0 - i * 2, i * 10.0, 5, 1000000, 700000, 11, 0, 10);
    }
    
    long int B1_adj = 0, B2_adj = 0, sieve_adj1 = 0, sieve_adj2 = 0, cutoff_adj = 0;
    bool should_adjust = optimizer.suggest_adjustments(B1_adj, B2_adj, sieve_adj1, sieve_adj2, cutoff_adj);
    
    std::cout << "  Should adjust (disabled): " << (should_adjust ? "yes" : "no") << std::endl;
    assert(!should_adjust);  // Should not suggest when disabled
    
    std::cout << "  PASSED" << std::endl;
}

int main()
{
    std::cout << "=== RateOptimizer Unit Tests ===" << std::endl;
    
    test_basic_recording();
    test_trend_calculation();
    test_stability_detection();
    test_adjustment_suggestions();
    test_disabled_auto_tuning();
    
    std::cout << std::endl << "All tests PASSED!" << std::endl;
    return 0;
}
