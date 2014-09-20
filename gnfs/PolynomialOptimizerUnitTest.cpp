#include "VeryLong.h"
#include "Polynomial.h"
#include "PolynomialOptimizer.h"
#include "MPFloat.h"
#include "pselect.h"
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/ui/text/TestRunner.h>
#include <iomanip>

class PolynomialOptimizerTest : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE(PolynomialOptimizerTest);
    CPPUNIT_TEST(test1);
    CPPUNIT_TEST_SUITE_END();

    public:
    void setUp() {}
    void tearDown() {}

    void test1()
    {
        VeryLong N("308268017314015502198864479573232259026111539470765498704205704273225473378735381489745954136566971705214233108249933319714229109554576259335341");
        Polynomial<VeryLong> f1 = Polynomial<VeryLong>::read_polynomial("-16612792256731454888860424973312 + 1498121907742379786368701480 X - 71616325394087869671190 X^2 - 436671758051847283 X^3 + 25011687432982 X^4 + 23794992 X^5");
        std::cout << "f1 = " << f1 << std::endl;
        MPFloat::set_precision(100);
        VeryLong a("11774494093137707");
        VeryLong b("1669121445503958099638459425"); 
        VeryLong m("254841105695668414977236991947256294566541981105672266845065891716312859988279394572523515845653926731492992971798027921010177494642224601803274");
        VeryLong best_s(53378L);
        double I_F_S = PolynomialOptimizer::average_log_size(f1, best_s);
        MPFloat I_F_S_mpf(0.0);
        VeryLong new_b;
        VeryLong new_m;

        VeryLong check = f1.evaluate_homogeneous(b, a) % N;
        CPPUNIT_ASSERT(check == 0L);
        double alpha = PolynomialOptimizer::alpha_F(f1, 2000, 200);
        double E_F = I_F_S + alpha;
	    std::cout << "f1 = " << f1 << std::endl;
	    std::cout << "best_s = " << best_s << std::endl;
	    std::cout << "alpha = " << alpha << std::endl;
	    std::cout << "I_F_S = " << std::setprecision(20) << I_F_S << std::endl;
	    std::cout << "E(F) = " << E_F << std::endl;

        double s1 = PolynomialOptimizer::minimize_I_over_s<double>(Polynomial<VeryLong>::convert_to_double<double>(f1), best_s);
        std::cout << "s1 = " << s1 << std::endl;

        VeryLong adjusted_s;
        std::vector<PolynomialOptimizer::Poly_info> poly_list;
        {
            std::fstream config_file("ut_skewed.cfg", std::ios::out);
            config_file << "N = 308268017314015502198864479573232259026111539470765498704205704273225473378735381489745954136566971705214233108249933319714229109554576259335341" << std::endl;
            config_file << "DEGREE = 5" << std::endl;
            config_file << "OUTPUT_FILE = skewed.out" << std::endl;
            config_file << "MIN_AD = 10000" << std::endl;
            config_file << "MAX_AD = 100000000" << std::endl;
            config_file << "C_START = 2310" << std::endl;
            config_file << "C_FACTOR = 0" << std::endl;
            config_file << "C_RESTART = 33989340" << std::endl;
            config_file << "GOOD_M_CUTOFF = 20.0" << std::endl;
            config_file << "MAX_FRACTION = 0.50" << std::endl;
            config_file << "MAX_ALS = 60.0" << std::endl;
            config_file << "MAX_J0 = 50000" << std::endl;
            config_file << "MAX_J1 = 20" << std::endl;
            config_file << "MAX_SMALL_PRIME = 1000" << std::endl;
            config_file << "ALPHA_CUTOFF = -4" << std::endl;
            config_file << "PRINTING_BOUND = 41.0" << std::endl;
            config_file << "REPEAT_CUTOFF = 39.0" << std::endl;
        }
        Skewed_selection_config config("ut_skewed.cfg");
        Polynomial<VeryLong> adjusted_f1 = adjust_root_properties(config, f1, a, b, adjusted_s, I_F_S, poly_list);

        double new_I_F_S(0.0);
        Polynomial<VeryLong> f2 = PolynomialOptimizer::minimize_I<double>(f1, a, b, m, best_s, new_I_F_S, new_b, new_m);
        alpha = PolynomialOptimizer::alpha_F(f2, 2000, 200);
        E_F = I_F_S + alpha;
        CPPUNIT_ASSERT(new_I_F_S < I_F_S);
	    std::cout << "f2 = " << f2 << std::endl;
	    std::cout << "best_s = " << best_s << std::endl;
	    std::cout << "I_F_S = " << new_I_F_S << std::endl;
	    std::cout << "alpha = " << alpha << std::endl;
	    std::cout << "E(F) = " << E_F << std::endl;
	    std::cout << "new_b = " << new_b << std::endl;
	    std::cout << "new_m = " << new_m << std::endl;
        check = f2.evaluate_homogeneous(new_b, a) % N;
        std::cout << "check = " << check << std::endl;
        CPPUNIT_ASSERT(check == 0L);

        VeryLong better_b;    
        VeryLong better_s;    
        Polynomial<VeryLong> translated_f2 = PolynomialOptimizer::translate(f2, a, new_b, best_s, better_b, better_s);
        alpha = PolynomialOptimizer::alpha_F(translated_f2, 2000, 200);
        double translated_I_F_S = PolynomialOptimizer::average_log_size(translated_f2, better_s);
        double translated_E_F = translated_I_F_S + alpha;
        CPPUNIT_ASSERT(translated_E_F < E_F);
	    std::cout << "translated_f2 = " << translated_f2 << std::endl;
	    std::cout << "better_s = " << better_s << std::endl;
	    std::cout << "I_F_S = " << translated_I_F_S << std::endl;
	    std::cout << "alpha = " << alpha << std::endl;
	    std::cout << "E(F) = " << translated_E_F << std::endl;
	    std::cout << "better_b = " << better_b << std::endl;
        check = translated_f2.evaluate_homogeneous(better_b, a) % N;
        std::cout << "check = " << check << std::endl;
        CPPUNIT_ASSERT(check == 0L);

        Polynomial<VeryLong> f3 = PolynomialOptimizer::minimize_I<MPFloat>(f1, a, b, m, best_s, I_F_S_mpf, new_b, new_m);
        alpha = PolynomialOptimizer::alpha_F(f3, 2000, 200);
        E_F = I_F_S_mpf + alpha;
        CPPUNIT_ASSERT((double)I_F_S_mpf < I_F_S);
	    std::cout << "f3 = " << f2 << std::endl;
	    std::cout << "best_s = " << best_s << std::endl;
	    std::cout << "I_F_S_mpf = " << I_F_S_mpf << std::endl;
	    std::cout << "alpha = " << alpha << std::endl;
	    std::cout << "E(F) = " << E_F << std::endl;
	    std::cout << "new_b = " << new_b << std::endl;
	    std::cout << "new_m = " << new_m << std::endl;
        check = f3.evaluate_homogeneous(new_b, a) % N;
        std::cout << "check = " << check << std::endl;
        CPPUNIT_ASSERT(check == 0L);

    }
};

int main()
{
    CppUnit::TextUi::TestRunner runner;
    runner.addTest(PolynomialOptimizerTest::suite());
    runner.run();

    return 0;
}

