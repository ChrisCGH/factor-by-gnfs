#include <fstream>
#include "LatticeSiever.h"
#include "SieveUtils.h"
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/ui/text/TestRunner.h>
#include <iomanip>
#include <cstdlib>

namespace
{
//LatticeSiever siever(std::string("ut_sieve_fixed.cfg"));
LatticeSiever* sieverp = 0;
bool verbose()
{
    if (std::getenv("GNFS_TEST_VERBOSE"))
    {
        return true;
    }
    return false;
}
}

class LatticeSieverTest : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE(LatticeSieverTest);
    CPPUNIT_TEST(test_parallelogram);
    CPPUNIT_TEST(test_fixed_sieve_region);
    CPPUNIT_TEST_SUITE_END();

public:
    void setUp()
    {
        make_config();
        VeryLong::generate_prime_table();
    }

    void tearDown()
    {
        cleanup();
    }

    void cleanup()
    {
        ::unlink("ut_sieve.cfg");
        //unlink("C144_104_45_UT.fb.dat");
        //unlink("C144_104_45_UT.rb.dat");
        ::unlink("sieve.tim");
        ::unlink("qs.tim");
    }

    void make_config()
    {
        {
            std::fstream config_file("ut_sieve.cfg", std::ios::out);
            config_file << "SIEVE_ID = C144_104_45_UT" << std::endl;
            config_file << "N = 308268017314015502198864479573232259026111539470765498704205704273225473378735381489745954136566971705214233108249933319714229109554576259335341" << std::endl;
            config_file << "m = 254841105695668414977236991947256294566541981105672266845065891716312859988279394572523515845653926731492992971798027921010177494642224601803274" << std::endl;
            config_file << "SKEWEDNESS = 53378" << std::endl;
            config_file << "f1 = -16612792256731454888860424973312 + 1498121907742379786368701480 X - 71616325394087869671190 X^2 - 436671758051847283 X^3 + 25011687432982 X^4 + 23794992 X^5" << std::endl;
            config_file << "f2 = -1669121445503958099638459425 + 11774494093137707 X" << std::endl;
            config_file << "B1 = 16000000" << std::endl;
            config_file << "L1 = 300000000" << std::endl;
            config_file << "LP1 = 2" << std::endl;
            config_file << "B2 = 16000000" << std::endl;
            config_file << "L2 = 300000000" << std::endl;
            config_file << "LP2 = 2" << std::endl;
            config_file << "SIEVE_BOUND_ADJUSTMENT1 = 20" << std::endl;
            config_file << "SIEVE_BOUND_ADJUSTMENT2 = 2" << std::endl;
            config_file << "SMALL_PRIME_BOUND1 = 500" << std::endl;
            config_file << "SMALL_PRIME_BOUND2 = 300" << std::endl;
            config_file << "MIN_A = -2600000000" << std::endl;
            config_file << "MAX_A = 2600000000" << std::endl;
            config_file << "MIN_B = 1" << std::endl;
            config_file << "MAX_B = 80000" << std::endl;
            config_file << "INITIAL_CUTOFF = 10" << std::endl;
            config_file << "RELATION_FILE = C144_104_45_UT.relations" << std::endl;
        }

        {
            std::fstream config_file("ut_sieve_fixed.cfg", std::ios::out);
            config_file << "SIEVE_ID = C144_104_45_UT" << std::endl;
            config_file << "N = 308268017314015502198864479573232259026111539470765498704205704273225473378735381489745954136566971705214233108249933319714229109554576259335341" << std::endl;
            config_file << "m = 254841105695668414977236991947256294566541981105672266845065891716312859988279394572523515845653926731492992971798027921010177494642224601803274" << std::endl;
            config_file << "SKEWEDNESS = 53378" << std::endl;
            config_file << "f1 = -16612792256731454888860424973312 + 1498121907742379786368701480 X - 71616325394087869671190 X^2 - 436671758051847283 X^3 + 25011687432982 X^4 + 23794992 X^5" << std::endl;
            config_file << "f2 = -1669121445503958099638459425 + 11774494093137707 X" << std::endl;
            config_file << "B1 = 16000000" << std::endl;
            config_file << "L1 = 300000000" << std::endl;
            config_file << "LP1 = 2" << std::endl;
            config_file << "B2 = 16000000" << std::endl;
            config_file << "L2 = 300000000" << std::endl;
            config_file << "LP2 = 2" << std::endl;
            config_file << "SIEVE_BOUND_ADJUSTMENT1 = 20" << std::endl;
            config_file << "SIEVE_BOUND_ADJUSTMENT2 = 2" << std::endl;
            config_file << "SMALL_PRIME_BOUND1 = 500" << std::endl;
            config_file << "SMALL_PRIME_BOUND2 = 300" << std::endl;
            config_file << "MIN_A = -2600000000" << std::endl;
            config_file << "MAX_A = 2600000000" << std::endl;
            config_file << "MIN_B = 1" << std::endl;
            config_file << "MAX_B = 80000" << std::endl;
            config_file << "INITIAL_CUTOFF = 10" << std::endl;
            config_file << "RELATION_FILE = C144_104_45_fixed_UT.relations" << std::endl;
            //config_file << "FIXED_SIEVE_REGION = true" << std::endl;
        }
    }

    void test_parallelogram()
    {
        const long int min_c(-8192L);
        const long int max_c(8191L);
        const long int min_d(0L);
        const long int max_d(4095L);
        Point<double> p1(min_c, min_d);
        Point<double> p2(max_c, min_d);
        Point<double> p3(min_c, max_d);
        Parallelogram c_region(p1, p2, p3);

        {
            std::pair<long int, long int> e1(11L, 32L);
            std::pair<long int, long int> e2(-14L, 5L);
            Parallelogram E_region(c_region, e1, e2);
            if (verbose())
            {
                c_region.display(std::cerr);
                E_region.display(std::cerr);
            }
            int32_t e_min = static_cast<int32_t>(std::ceil(E_region.min_x()));
            int32_t e_max = static_cast<int32_t>(std::floor(E_region.max_x()));
            if (verbose())
            {
                std::cerr << "e_min = " << e_min << std::endl;
                std::cerr << "e_max = " << e_max << std::endl;
            }
            int32_t f_min = 0L;
            int32_t f_max = 0L;
            int32_t e(e_min);
            while (e <= e_max)
            {
                while (!E_region.y_limits1(e, f_min, f_max) && e <= e_max) ++e;
                if (e <= e_max)
                {
                    long int c1 = e * e1.first + f_min * e2.first;
                    long int d1 = e * e1.second + f_min * e2.second;
                    long int c2 = e * e1.first + f_max * e2.first;
                    long int d2 = e * e1.second + f_max * e2.second;
                    if (verbose())
                    {
                        std::cerr << "e = " << e << ", f_min = " << f_min << ", f_max = " << f_max << ", (c1,d1) = (" << c1 << "," << d1 << "), (c2,d2) = (" << c2 << "," << d2 << ")" << std::endl;
                    }

                    for (long int f = f_min; f <= f_max; ++f)
                    {
                        long int c = e * e1.first + f * e2.first;
                        long int d = e * e1.second + f * e2.second;
                        if (verbose())
                        {
                            if (c < min_c || c > max_c || d < min_d || d > max_d)
                            {
                                std::cerr << "(c,d) = (" << c << "," << d << ")" << std::endl;
                            }
                        }
                        CPPUNIT_ASSERT(c >= min_c && c <= max_c);
                        CPPUNIT_ASSERT(d >= min_d && d <= max_d);
                    }

                    ++e;
                }
            }
        }
        {
            std::pair<long int, long int> e1(16L, 25L);
            std::pair<long int, long int> e2(19L, -9L);
            Parallelogram E_region(c_region, e1, e2);
            if (verbose())
            {
                c_region.display(std::cerr);
                E_region.display(std::cerr);
            }
            int32_t e_min = static_cast<long int>(std::ceil(E_region.min_x()));
            int32_t e_max = static_cast<long int>(std::floor(E_region.max_x()));
            if (verbose())
            {
                std::cerr << "e_min = " << e_min << std::endl;
                std::cerr << "e_max = " << e_max << std::endl;
            }
            int32_t f_min = 0L;
            int32_t f_max = 0L;
            int32_t e(e_min);
            while (e <= e_max)
            {
                while (!E_region.y_limits1(e, f_min, f_max) && e <= e_max) ++e;
                if (e <= e_max)
                {
                    long int c1 = e * e1.first + f_min * e2.first;
                    long int d1 = e * e1.second + f_min * e2.second;
                    long int c2 = e * e1.first + f_max * e2.first;
                    long int d2 = e * e1.second + f_max * e2.second;
                    if (verbose())
                    {
                        std::cerr << "e = " << e << ", f_min = " << f_min << ", f_max = " << f_max << ", (c1,d1) = (" << c1 << "," << d1 << "), (c2,d2) = (" << c2 << "," << d2 << ")" << std::endl;
                    }

                    for (long int f = f_min; f <= f_max; ++f)
                    {
                        long int c = e * e1.first + f * e2.first;
                        long int d = e * e1.second + f * e2.second;
                        if (verbose())
                        {
                            if (c < min_c || c > max_c || d < min_d || d > max_d)
                            {
                                std::cerr << "(e,f) = (" << e << "," << f << "), (c,d) = (" << c << "," << d << ")" << std::endl;
                            }
                        }
                        CPPUNIT_ASSERT(c >= min_c && c <= max_c);
                        CPPUNIT_ASSERT(d >= min_d && d <= max_d);
                    }
                    ++e;
                }
            }
        }
        {
            std::pair<long int, long int> e1(-1943L, -1854L);
            std::pair<long int, long int> e2(1057L, -955L);
            Parallelogram E_region(c_region, e1, e2);
            if (verbose())
            {
                c_region.display(std::cerr);
                E_region.display(std::cerr);
            }
            int32_t e_min = static_cast<long int>(std::ceil(E_region.min_x()));
            int32_t e_max = static_cast<long int>(std::floor(E_region.max_x()));
            if (verbose())
            {
                std::cerr << "e_min = " << e_min << std::endl;
                std::cerr << "e_max = " << e_max << std::endl;
            }
            int32_t f_min = 0L;
            int32_t f_max = 0L;
            int32_t e(e_min);
            while (e <= e_max)
            {
                while (!E_region.y_limits1(e, f_min, f_max) && e <= e_max) ++e;
                if (e <= e_max)
                {
                    long int c1 = e * e1.first + f_min * e2.first;
                    long int d1 = e * e1.second + f_min * e2.second;
                    long int c2 = e * e1.first + f_max * e2.first;
                    long int d2 = e * e1.second + f_max * e2.second;
                    if (verbose())
                    {
                        std::cerr << "e = " << e << ", f_min = " << f_min << ", f_max = " << f_max << ", (c1,d1) = (" << c1 << "," << d1 << "), (c2,d2) = (" << c2 << "," << d2 << ")" << std::endl;
                    }

                    for (long int f = f_min; f <= f_max; ++f)
                    {
                        long int c = e * e1.first + f * e2.first;
                        long int d = e * e1.second + f * e2.second;
                        if (verbose())
                        {
                            if (c < min_c || c > max_c || d < min_d || d > max_d)
                            {
                                std::cerr << "(e,f) = (" << e << "," << f << "), (c,d) = (" << c << "," << d << ")" << std::endl;
                            }
                        }
                        CPPUNIT_ASSERT(c >= min_c && c <= max_c);
                        CPPUNIT_ASSERT(d >= min_d && d <= max_d);
                    }
                    ++e;
                }
            }
        }
    }

    void test_fixed_sieve_region()
    {
        sieverp = new LatticeSiever(std::string("ut_sieve_fixed.cfg"));
        LatticeSiever& siever(*sieverp);
        CPPUNIT_ASSERT(siever.f1_ == Polynomial<VeryLong>::read_polynomial("-16612792256731454888860424973312 + 1498121907742379786368701480 X - 71616325394087869671190 X^2 - 436671758051847283 X^3 + 25011687432982 X^4 + 23794992 X^5"));
        CPPUNIT_ASSERT(siever.f2_ == Polynomial<VeryLong>::read_polynomial("-1669121445503958099638459425 + 11774494093137707 X"));
        CPPUNIT_ASSERT(siever.B1_ == 16000000L);
        CPPUNIT_ASSERT(siever.L1_ == 300000000L);
        CPPUNIT_ASSERT(siever.B2_ == 16000000L);
        CPPUNIT_ASSERT(siever.L2_ == 300000000L);
        CPPUNIT_ASSERT(siever.LP1_ == 2L);
        CPPUNIT_ASSERT(siever.LP2_ == 2L);
        CPPUNIT_ASSERT(siever.MIN_A_ == -2600000000.0);
        CPPUNIT_ASSERT(siever.MAX_A_ == +2600000000.0);
        CPPUNIT_ASSERT(siever.MIN_B_ == 1L);
        CPPUNIT_ASSERT(siever.MAX_B_ == 80000L);
        CPPUNIT_ASSERT(siever.SIEVE_BOUND_ADJUSTMENT1_ == 20L);
        CPPUNIT_ASSERT(siever.SIEVE_BOUND_ADJUSTMENT2_ == 2L);
        CPPUNIT_ASSERT(siever.SMALL_PRIME_BOUND1_ == 500L);
        CPPUNIT_ASSERT(siever.SMALL_PRIME_BOUND2_ == 300L);
        CPPUNIT_ASSERT(siever.INITIAL_CUTOFF_ == 10L);
        CPPUNIT_ASSERT(siever.SKEWEDNESS_ == 53378.0);

        long int q = 10000019L;
        std::vector<LongModular> q_roots;
        find_roots_mod_p<VeryLong, long int, LongModular>(siever.f1_, q, q_roots);
        CPPUNIT_ASSERT(!q_roots.empty());
        CPPUNIT_ASSERT(q_roots.size() == 2);

        long int s = q_roots[0].get_long();
        CPPUNIT_ASSERT(s == 6621906L);

        siever.generate_lattice(q, s);
        CPPUNIT_ASSERT(siever.c1_.first == -1900593L);
        CPPUNIT_ASSERT(siever.c1_.second == -32L);
        CPPUNIT_ASSERT(siever.c2_.first == 134320L);
        CPPUNIT_ASSERT(siever.c2_.second == -3L);

        s = q_roots[1].get_long();
        CPPUNIT_ASSERT(s == 5918925L);

        siever.generate_lattice(q, s);
        CPPUNIT_ASSERT(siever.c1_.first == -621535L);
        CPPUNIT_ASSERT(siever.c1_.second == -17L);
        CPPUNIT_ASSERT(siever.c2_.first == 405432L);
        CPPUNIT_ASSERT(siever.c2_.second == -5L);

        CPPUNIT_ASSERT(siever.allocate_c_d_region());

        CPPUNIT_ASSERT_DOUBLES_EQUAL(siever.c_region_.p1_.x, siever.min_c, 1.0e-2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(siever.c_region_.p1_.y, siever.min_d, 1.0e-2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(siever.c_region_.p2_.x, siever.min_c, 1.0e-2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(siever.c_region_.p2_.y, siever.max_d, 1.0e-2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(siever.c_region_.p3_.x, siever.max_c, 1.0e-2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(siever.c_region_.p3_.y, siever.min_d, 1.0e-2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(siever.c_region_.p4_.x, siever.max_c, 1.0e-2);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(siever.c_region_.p4_.y, siever.max_d, 1.0e-2);

#ifdef USE_FIXED_SIEVE_REGION
        {
            std::pair<long int, long int> cd;
            cd.first = -4000L;
            cd.second = 2345L;
            size_t offset = siever.c_d_to_offset(cd);
            std::pair<long int, long int> cd1 = siever.offset_to_c_d(offset);
            CPPUNIT_ASSERT(cd == cd1);

        }
#endif

        //std::cout << "sieve_array_size_ = " << siever.sieve_array_size_ << std::endl;
        //std::cout << "max_d_span_ = " << siever.max_d_span_ << std::endl;
        //CPPUNIT_ASSERT(siever.sieve_array_size_ == 40022001L);
        //CPPUNIT_ASSERT(siever.max_d_span_ == 26667L);

        SieveConfig config("ut_sieve_fixed.cfg");
        std::string relfile = config.RELATION_FILE();
        ::unlink(relfile.c_str());


        CPPUNIT_ASSERT(siever.sieve(q));

        CPPUNIT_ASSERT(SieveUtils::checkRelations(relfile.c_str(), config));

    }
};


int main()
{
    CppUnit::TextUi::TestRunner runner;
    runner.addTest(LatticeSieverTest::suite());
    runner.run();

    return 0;
}

