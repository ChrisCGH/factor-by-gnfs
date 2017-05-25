#include <iostream>
#include "VeryLongModular.h"
#ifdef USING_CPPUNIT
#include <limits>
#include <complex>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/ui/text/TestRunner.h>
#else
#include "UnitTest.h"
#endif

#ifdef USING_CPPUNIT
class VeryLongModularTest : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE(VeryLongModularTest);
    CPPUNIT_TEST(test);
    CPPUNIT_TEST_SUITE_END();

public:
    void setUp()
    {
    }

    void tearDown()
    {
    }

    void test()
    {
        const VeryLongModular zero;
        CPPUNIT_ASSERT(zero.get_very_long() == VeryLong(0UL));

        VeryLongModular::set_default_modulus(11L);
        const VeryLongModular lm1(13L);
        CPPUNIT_ASSERT(lm1.get_very_long() == VeryLong(2UL));

        long long int lli = 11LL * 0x100000000LL + 9LL;
        const VeryLongModular lm2(lli);
        CPPUNIT_ASSERT(lm2.get_very_long() == VeryLong(9UL));

        const VeryLongModular lm3(lm2);
        CPPUNIT_ASSERT(lm3.get_very_long() == VeryLong(9UL));

        const VeryLongModular lm4(29L, 57L);
        CPPUNIT_ASSERT(lm4.get_very_long() == VeryLong(28UL));

        CPPUNIT_ASSERT(lm2 == lm3);
        CPPUNIT_ASSERT(lm2 != lm4);

        CPPUNIT_ASSERT(!(lm2 < lm2));
        CPPUNIT_ASSERT(!(lm2 > lm2));
        CPPUNIT_ASSERT(!(lm2 < 0L));

        VeryLongModular::set_default_modulus(11L);
        VeryLongModular lm5;
        lm5 = lm2;
        CPPUNIT_ASSERT(lm2 == lm5);
        lm5 += lm1;
        CPPUNIT_ASSERT(lm5.get_very_long() == VeryLong(0L));

        lm5 -= lm2;
        CPPUNIT_ASSERT(lm5.get_very_long() == VeryLong(2L));

        lm5 *= lm2;
        CPPUNIT_ASSERT(lm5.get_very_long() == VeryLong(7L));

        lm5 /= lm2;
        CPPUNIT_ASSERT(lm5.get_very_long() == VeryLong(2L));

        VeryLongModular lm6 = -lm5;
        CPPUNIT_ASSERT(lm6.get_very_long() == VeryLong(9L));

        lm6.exp(5L);
        CPPUNIT_ASSERT(lm6.get_very_long() == VeryLong(1L));

        VeryLongModular lm7 = exp(lm6, lm1);
        CPPUNIT_ASSERT(lm7.get_very_long() == VeryLong(1L));

        lm7 = exp(lm2, lm1);
        CPPUNIT_ASSERT(lm7.get_very_long() == VeryLong(4L));

        VeryLongModular lm8 = lm1 + lm2;
        CPPUNIT_ASSERT(lm8.get_very_long() == VeryLong(0L));

        lm8 = lm1 - lm2;
        CPPUNIT_ASSERT(lm8.get_very_long() == VeryLong(4L));

        lm8 = lm1 * lm2;
        CPPUNIT_ASSERT(lm8.get_very_long() == VeryLong(7L));

        lm8 = lm1 / lm2;
        CPPUNIT_ASSERT(lm8.get_very_long() == VeryLong(10L));

        CPPUNIT_ASSERT(zero.is_zero());
        CPPUNIT_ASSERT(!lm8.is_zero());

        VeryLongModular::set_default_modulus(15999251L);
        const VeryLongModular lm100(2230183161LL);
        CPPUNIT_ASSERT(lm100.get_very_long() == VeryLong(6287272L));

        // Tests for add_product
        VeryLongModular::set_default_modulus(15999251L);
        VeryLongModular lm301(12345678L);
        VeryLongModular lm302(2345678L);
        VeryLongModular lm303(345678L);

        VeryLongModular lm304(lm301);
        lm304 *= lm302;
        lm304 *= lm303;
        VeryLongModular lm400(45678L);
        lm400 += lm304;
        VeryLongModular lm305(45678L);
        lm305.add_product(lm301, lm302, lm303);
        CPPUNIT_ASSERT(lm400 == lm305);

        VeryLongModular::set_default_modulus(5L);
        VeryLongModular lm500(1L);
        VeryLongModular lm501(4L);
        CPPUNIT_ASSERT(lm500 < lm501);
        CPPUNIT_ASSERT((!(lm501 < lm500)));
    }
};
#endif

int main()
{
#ifdef USING_CPPUNIT
    CppUnit::TextUi::TestRunner runner;
    runner.addTest(VeryLongModularTest::suite());
    runner.run();
#else
    UnitTest t;
    std::cout << "Testing VeryLongModular ..." << std::endl;

    const VeryLongModular zero;
    t.check(zero.get_very_long() == 0UL, "zero should be 0");

    VeryLongModular::set_default_modulus(11L);
    const VeryLongModular lm1(13L);
    t.check(lm1.get_very_long() == 2UL, "lm1 should be 2");

    long long int lli = 11LL * 0x100000000LL + 9LL;
    const VeryLongModular lm2(lli);
    t.check(lm2.get_very_long() == 9UL, "lm2 should be 9");

    const VeryLongModular lm3(lm2);
    t.check(lm3.get_very_long() == 9UL, "lm3 should be 9");

    const VeryLongModular lm4(29L, 57L);
    t.check(lm4.get_very_long() == 28UL, "lm4 should be 28");

    t.check(lm2 == lm3, "lm2 should equal lm3");
    t.check(lm2 != lm4, "lm2 should not equal lm4");

    t.check(!(lm2 < lm2), "operator< should always return false");
    t.check(!(lm2 > lm2), "operator> should always return false");
    t.check(!(lm2 < 0L), "operator< should always return false");

    VeryLongModular::set_default_modulus(11L);
    VeryLongModular lm5;
    lm5 = lm2;
    t.check(lm2 == lm5, "lm2 should equal lm5");
    lm5 += lm1;
    t.check(lm5.get_very_long() == 0L, "lm5 should equal 0");

    lm5 -= lm2;
    t.check(lm5.get_very_long() == 2L, "lm5 should equal 2");

    lm5 *= lm2;
    t.check(lm5.get_very_long() == 7L, "lm5 should equal 7");

    lm5 /= lm2;
    t.check(lm5.get_very_long() == 2L, "lm5 should equal 2");

    VeryLongModular lm6 = -lm5;
    t.check(lm6.get_very_long() == 9L, "lm6 should equal 9");

#if 0
    lm6.inverse();
    t.check(lm6.get_very_long() == 5L, "lm6 should equal 5");
#endif

    lm6.exp(5L);
    t.check(lm6.get_very_long() == 1L, "lm6 should equal 1");

    VeryLongModular lm7 = exp(lm6, lm1);
    t.check(lm7.get_very_long() == 1L, "lm7 should equal 1");

    lm7 = exp(lm2, lm1);
    t.check(lm7.get_very_long() == 4L, "lm7 should equal 4");

#if 0
    lm7 = exp(lm2, 3L);
    t.check(lm7.get_very_long() == 3, "lm7 should equal 3");
#endif

    VeryLongModular lm8 = lm1 + lm2;
    t.check(lm8.get_very_long() == 0L, "lm8 should equal 0");

    lm8 = lm1 - lm2;
    t.check(lm8.get_very_long() == 4L, "lm8 should equal 4");

    lm8 = lm1 * lm2;
    t.check(lm8.get_very_long() == 7L, "lm8 should equal 7");

    lm8 = lm1 / lm2;
    t.check(lm8.get_very_long() == 10L, "lm8 should equal 10");

    t.check(zero.is_zero(), "zero.is_zero() should be true");
    t.check(!lm8.is_zero(), "lm8.is_zero() should be false");

    VeryLongModular::set_default_modulus(15999251L);
    const VeryLongModular lm100(2230183161LL);
    t.check(lm100.get_very_long() == 6287272L, "2230183161 mod 15999251 should be 6287272L");

    t.test_summary();
#endif
    return 0;
}
