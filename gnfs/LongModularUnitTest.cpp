#include <iostream>
#include <sstream>
#include "LongModular.h"
#include "gcd.h"
#include "crt.h"
#ifdef USING_CPPUNIT
#include <limits>
#include <complex>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/ui/text/TestRunner.h>
#else
#include "UnitTest.h"
#endif

#ifdef USING_CPPUNIT
class LongModularTest : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE(LongModularTest);
    CPPUNIT_TEST(test);
    CPPUNIT_TEST_SUITE_END();

    private:

    public:
    void setUp()
    {
    }

    void tearDown()
    {
    }

    void test()
    {
        const LongModular zero;
        CPPUNIT_ASSERT(zero.get_long() == 0UL);
        
        LongModular::set_default_modulus(11L);
        const LongModular lm1(13L);
        CPPUNIT_ASSERT(lm1.get_long() == 2UL);

        long long int lli = 11LL * 0x100000000LL + 9LL;
        const LongModular lm2(lli);
        CPPUNIT_ASSERT(lm2.get_long() == 9UL);

        const LongModular lm3(lm2);
        CPPUNIT_ASSERT(lm3.get_long() == 9UL);

        const LongModular lm4(29L, 57L);
        CPPUNIT_ASSERT(lm4.get_long() == 28UL);

        const LongModular lm200("1");
        CPPUNIT_ASSERT(lm200.get_long() == 1UL);

        CPPUNIT_ASSERT(lm2 == lm3);
        CPPUNIT_ASSERT(lm2 != lm4);

        CPPUNIT_ASSERT(!(lm2 < lm2));
        CPPUNIT_ASSERT(!(lm2 > lm2));
        CPPUNIT_ASSERT(!(lm2 < 0L));
   
        LongModular::set_default_modulus(11L);
   
        const LongModular lm201("21");
        CPPUNIT_ASSERT(lm201.get_long() == 10UL);

        LongModular lm5;
        lm5 = lm2;
        CPPUNIT_ASSERT(lm2 == lm5);
        lm5 += lm1;
        CPPUNIT_ASSERT(lm5.get_long() == 0);

        lm5 -= lm2;
        CPPUNIT_ASSERT(lm5.get_long() == 2);

        lm5 *= lm2;
        CPPUNIT_ASSERT(lm5.get_long() == 7);

        lm5 /= lm2;
        CPPUNIT_ASSERT(lm5.get_long() == 2);

        LongModular lm6 = -lm5;
        CPPUNIT_ASSERT(lm6.get_long() == 9);

        lm6.exp(5L);
        CPPUNIT_ASSERT(lm6.get_long() == 1);

        LongModular lm7 = exp(lm6, lm1);
        CPPUNIT_ASSERT(lm7.get_long() == 1);

        lm7 = exp(lm2, lm1);
        CPPUNIT_ASSERT(lm7.get_long() == 4);

        LongModular lm8 = lm1 + lm2;
        CPPUNIT_ASSERT(lm8.get_long() == 0);

        lm8 = lm1 - lm2;
        CPPUNIT_ASSERT(lm8.get_long() == 4);

        lm8 = lm1 * lm2;
        CPPUNIT_ASSERT(lm8.get_long() == 7);

        lm8 = lm1 / lm2;
        CPPUNIT_ASSERT(lm8.get_long() == 10);

        CPPUNIT_ASSERT(zero.is_zero());
        CPPUNIT_ASSERT(!lm8.is_zero());

        CPPUNIT_ASSERT_THROW_MESSAGE("", LongModular lm9 = lm1.square_root(), std::string);

        LongModular lm10(5L);
        LongModular lm9;
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", lm9 = lm10.square_root());

        CPPUNIT_ASSERT(lm9.get_long() == 4);

        LongModular::set_default_modulus(15999251L);
        const LongModular lm100(2230183161LL);
        CPPUNIT_ASSERT(lm100.get_long() == 6287272L);

        // tests for crt
        std::vector<long int> p(4);
        std::vector<long long> b(4);
        p[0] = 1201;
        b[0] = 1030;
        p[1] = 1171;
        b[1] = 593;
        p[2] = 83;
        b[2] = 30;
        p[3] = 1409;
        b[3] = 523;
        long long x = crt<long long, long int>(b, p);
        CPPUNIT_ASSERT(x == 58131719393LL);

        p[0] = 1201;
        b[0] = 67;
        p[1] = 1171;
        b[1] = 593;
        p[2] = 83;
        b[2] = 30;
        p[3] = 1409;
        b[3] = 523;
        x = crt<long long, long int>(b, p);
        CPPUNIT_ASSERT(x == -19379114949LL);

        // Tests for add_product
        LongModular::set_default_modulus(15999251L);
        LongModular lm301(12345678L);
        LongModular lm302(2345678L);
        LongModular lm303(345678L);

        LongModular lm304(lm301);
        lm304 *= lm302;
        lm304 *= lm303;
        LongModular lm400(15999250L);
        lm400 += lm304;
        LongModular lm305(15999250L);
        lm305.add_product(lm301, lm302, lm303);
        CPPUNIT_ASSERT(lm400 == lm305);

    }

};
#endif

int main()
{
#ifdef USING_CPPUNIT
   CppUnit::TextUi::TestRunner runner;
   runner.addTest(LongModularTest::suite());
   runner.run();
#else
   UnitTest t;
   std::cout << "Testing LongModular ..." << std::endl;

   const LongModular zero;
   t.check(zero.get_long() == 0UL, "zero should be 0");

   LongModular::set_default_modulus(11L);
   const LongModular lm1(13L);
   t.check(lm1.get_long() == 2UL, "lm1 should be 2");

   long long int lli = 11LL * 0x100000000LL + 9LL;
   const LongModular lm2(lli);
   t.check(lm2.get_long() == 9UL, "lm2 should be 9");

   const LongModular lm3(lm2);
   t.check(lm3.get_long() == 9UL, "lm3 should be 9");

   const LongModular lm4(29L, 57L);
   t.check(lm4.get_long() == 28UL, "lm4 should be 28");

   t.check(lm2 == lm3, "lm2 should equal lm3");
   t.check(lm2 != lm4, "lm2 should not equal lm4");

   t.check(!(lm2 < lm2), "operator< should always return false");
   t.check(!(lm2 > lm2), "operator> should always return false");
   t.check(!(lm2 < 0L), "operator< should always return false");

   LongModular::set_default_modulus(11L);
   LongModular lm5;
   lm5 = lm2;
   t.check(lm2 == lm5, "lm2 should equal lm5");
   lm5 += lm1;
   t.check(lm5.get_long() == 0, "lm5 should equal 0");

   lm5 -= lm2;
   t.check(lm5.get_long() == 2, "lm5 should equal 2");

   lm5 *= lm2;
   t.check(lm5.get_long() == 7, "lm5 should equal 7");

   lm5 /= lm2;
   t.check(lm5.get_long() == 2, "lm5 should equal 2");

   LongModular lm6 = -lm5;
   t.check(lm6.get_long() == 9, "lm6 should equal 9");

#if 0
   lm6.inverse();
   t.check(lm6.get_long() == 5, "lm6 should equal 5");
#endif

   lm6.exp(5L);
   t.check(lm6.get_long() == 1, "lm6 should equal 1");

   LongModular lm7 = exp(lm6, lm1);
   t.check(lm7.get_long() == 1, "lm7 should equal 1");

   lm7 = exp(lm2, lm1);
   t.check(lm7.get_long() == 4, "lm7 should equal 4");

#if 0
   lm7 = exp(lm2, 3L);
   t.check(lm7.get_long() == 3, "lm7 should equal 3");
#endif

   LongModular lm8 = lm1 + lm2;
   t.check(lm8.get_long() == 0, "lm8 should equal 0");

   lm8 = lm1 - lm2;
   t.check(lm8.get_long() == 4, "lm8 should equal 4");

   lm8 = lm1 * lm2;
   t.check(lm8.get_long() == 7, "lm8 should equal 7");

   lm8 = lm1 / lm2;
   t.check(lm8.get_long() == 10, "lm8 should equal 10");

   t.check(zero.is_zero(), "zero.is_zero() should be true");
   t.check(!lm8.is_zero(), "lm8.is_zero() should be false");

   bool thrown = false;
   try
   {
   	LongModular lm9 = lm1.square_root();
   }
   catch (std::string& s)
   {
	thrown = true;
   }
   t.check(thrown, "square_root() should have thrown");

   const LongModular lm10(5L);
   LongModular lm9 = lm10.square_root();
   t.check(lm9.get_long() == 4, "4 should be square root of 5 mod 11");

   LongModular::set_default_modulus(15999251L);
   const LongModular lm100(2230183161LL);
   t.check(lm100.get_long() == 6287272L, "2230183161 mod 15999251 should be 6287272L");

   // tests for crt
   std::vector<long int> p(4);
   std::vector<long long> b(4);
   p[0] = 1201;
   b[0] = 1030;
   p[1] = 1171;
   b[1] = 593;
   p[2] = 83;
   b[2] = 30;
   p[3] = 1409;
   b[3] = 523;
   long long x = crt<long long, long int>(b, p);
   t.check(x == 58131719393LL, "result of crt should be 58131719393");

   p[0] = 1201;
   b[0] = 67;
   p[1] = 1171;
   b[1] = 593;
   p[2] = 83;
   b[2] = 30;
   p[3] = 1409;
   b[3] = 523;
   x = crt<long long, long int>(b, p);
   t.check(x == -19379114949LL, "result of crt should be -19379114949");

   t.test_summary();
#endif
   return 0;
}
