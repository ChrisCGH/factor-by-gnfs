#include "Quotient.h"
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/ui/text/TestRunner.h>
#include "VeryLong.h"

class QuotientTest : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE(QuotientTest);
    CPPUNIT_TEST(testQuotientLongInt);
    CPPUNIT_TEST(testQuotientVeryLong);
    CPPUNIT_TEST_SUITE_END();

public:
    void setUp()
    {
    }
    void tearDown()
    {
    }

    void testQuotientLongInt()
    {
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", Quotient<long int>());
        Quotient<long int> q1;
        CPPUNIT_ASSERT(q1.is_zero());
        CPPUNIT_ASSERT(q1.numerator() == 0L);
        CPPUNIT_ASSERT(q1.denominator() == 1L);

        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", Quotient<long int>(10L));
        Quotient<long int> q2(10L);
        CPPUNIT_ASSERT(!q2.is_zero());
        CPPUNIT_ASSERT(q2.numerator() == 10L);
        CPPUNIT_ASSERT(q2.denominator() == 1L);

        CPPUNIT_ASSERT_THROW_MESSAGE("", Quotient<long int>(44L*1000234L, 0L), std::string);
        CPPUNIT_ASSERT_THROW_MESSAGE("", Quotient<long int>(0L, 0L), std::string);

        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", Quotient<long int>(44L*1000234L, 14L*1000234L));
        Quotient<long int> q3(44L*1000234L, 14L*1000234L);
        CPPUNIT_ASSERT(!q3.is_zero());
        CPPUNIT_ASSERT(q3.numerator() == 22L);
        CPPUNIT_ASSERT(q3.denominator() == 7L);

        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", Quotient<long int>(q3));
        Quotient<long int> q4(q3);
        CPPUNIT_ASSERT(!q4.is_zero());
        CPPUNIT_ASSERT(q3 == q4);
        CPPUNIT_ASSERT(q4.numerator() == 22L);
        CPPUNIT_ASSERT(q4.denominator() == 7L);

        Quotient<long int> q5(5L, 5L);
        CPPUNIT_ASSERT(q5.is_one());
        CPPUNIT_ASSERT(q5 != q4);

        Quotient<long int> q6(10L, 3L);
        q3 = q6;
        Quotient<long int> q7(q6);
        CPPUNIT_ASSERT(q3 == q7); 

        Quotient<long int> q8(5L, 8L);
        Quotient<long int> q9(2L, 3L);
        CPPUNIT_ASSERT(q8 < q9);
        CPPUNIT_ASSERT(q9 > q8);
        CPPUNIT_ASSERT(q8 <= q9);
        CPPUNIT_ASSERT(q9 >= q8);

        Quotient<long int> q10 = -q9;
        CPPUNIT_ASSERT(q9 + q10 == Quotient<long int>(0L));
        Quotient<long int> q11(q10);
        CPPUNIT_ASSERT(q10 - q11 == Quotient<long int>(0L));

        q10 -= q11;
        CPPUNIT_ASSERT(q10 == Quotient<long int>(0L));

        q10 += q11;
        CPPUNIT_ASSERT(q10 == q11);

        Quotient<long int> q12(13L, 23L);
        Quotient<long int> q13(23L, 13L);

        q12 *= q13;
        CPPUNIT_ASSERT(q12.is_one());

        q13 *= 26L;
        CPPUNIT_ASSERT(q13 == Quotient<long int>(46L));

        Quotient<long int> q14(117L, 13L);
        Quotient<long int> q15(13L, 9L);
        Quotient<long int> q16 = q14 * q15;
        CPPUNIT_ASSERT(q16 == Quotient<long int>(13L));

        CPPUNIT_ASSERT(q14 * Quotient<long int>(0L) == Quotient<long int>(0L));
        CPPUNIT_ASSERT(Quotient<long int>(0L) * q15 == Quotient<long int>(0L));
        CPPUNIT_ASSERT(Quotient<long int>(10L, 10L) * q15 == q15);
        CPPUNIT_ASSERT(q15 == Quotient<long int>(10L, 10L) * q15);

        CPPUNIT_ASSERT_THROW_MESSAGE("", (Quotient<long int>(10L, 3L) / Quotient<long int>(0L)), std::string);
        CPPUNIT_ASSERT(q14 / Quotient<long int>(5L, 5L) == q14);

        Quotient<long int> q17(111L, 5L);
        Quotient<long int> q18(222L, 5L);
        CPPUNIT_ASSERT(q17 / q18 == Quotient<long int>(1L, 2L));

/*
        std::istringstream iss("102325/2302938");
        Quotient<long int> q19;
        iss >> q19;
        CPPUNIT_ASSERT(q19 == Quotient<long int>(102325L, 2302938L));
        std::ostringstream oss;
        oss << q19;
        CPPUNIT_ASSERT(oss.str() == "102325/2302938");
*/

        Quotient<long int> q20(5L, 2L);
        Quotient<long int> q21 = exp(q20, 3L);
        CPPUNIT_ASSERT(q21 == Quotient<long int>(125L, 8L));

        CPPUNIT_ASSERT_THROW_MESSAGE("", (q21 /= Quotient<long int>(0L)), std::string);
        Quotient<long int> q22(0L);
        CPPUNIT_ASSERT_THROW_MESSAGE("", (q22 /= Quotient<long int>(0L)), std::string);
        q21 /= q21;
        CPPUNIT_ASSERT(q21 == 1L);
        Quotient<long int> q23(100L, 3L);
        q23 /= Quotient<long int>(10L, 3L);
        CPPUNIT_ASSERT(q23 == 10L);
    }

    void testQuotientVeryLong()
    {
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", Quotient<VeryLong>());
        Quotient<VeryLong> q1;
        CPPUNIT_ASSERT(q1.is_zero());
        CPPUNIT_ASSERT(q1.numerator() == VeryLong(0L));
        CPPUNIT_ASSERT(q1.denominator() == VeryLong(1L));

        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", Quotient<VeryLong>(VeryLong(10L)));
        Quotient<VeryLong> q2(VeryLong(10L));
        CPPUNIT_ASSERT(!q2.is_zero());
        CPPUNIT_ASSERT(q2.numerator() == VeryLong(10L));
        CPPUNIT_ASSERT(q2.denominator() == VeryLong(1L));

        CPPUNIT_ASSERT_THROW_MESSAGE("", Quotient<VeryLong>(VeryLong(44L*1000234L), VeryLong(0L)), std::string);
        CPPUNIT_ASSERT_THROW_MESSAGE("", Quotient<VeryLong>(VeryLong(0L), VeryLong(0L)), std::string);
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", Quotient<VeryLong>(VeryLong(44L*1000234L), VeryLong(14L*1000234L)));
        Quotient<VeryLong> q3(VeryLong(44L*1000234L), VeryLong(14L*1000234L));
        CPPUNIT_ASSERT(!q3.is_zero());
        CPPUNIT_ASSERT(q3.numerator() == VeryLong(22L));
        CPPUNIT_ASSERT(q3.denominator() == VeryLong(7L));

        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", Quotient<VeryLong>(q3));
        Quotient<VeryLong> q4(q3);
        CPPUNIT_ASSERT(!q4.is_zero());
        CPPUNIT_ASSERT(q3 == q4);
        CPPUNIT_ASSERT(q4.numerator() == VeryLong(22L));
        CPPUNIT_ASSERT(q4.denominator() == VeryLong(7L));

        Quotient<VeryLong> q5(VeryLong(5L), VeryLong(5L));
        CPPUNIT_ASSERT(q5.is_one());
        CPPUNIT_ASSERT(q5 != q4);

        Quotient<VeryLong> q6(VeryLong(10L), VeryLong(3L));
        q3 = q6;
        Quotient<VeryLong> q7(q6);
        CPPUNIT_ASSERT(q3 == q7);

        Quotient<VeryLong> q8(VeryLong(5L), VeryLong(8L));
        Quotient<VeryLong> q9(VeryLong(2L), VeryLong(3L));
        CPPUNIT_ASSERT(q8 < q9);
        CPPUNIT_ASSERT(q9 > q8);
        CPPUNIT_ASSERT(q8 <= q9);
        CPPUNIT_ASSERT(q9 >= q8);

        Quotient<VeryLong> q10 = -q9;
        CPPUNIT_ASSERT(q9 + q10 == Quotient<VeryLong>(VeryLong(0L)));
        Quotient<VeryLong> q11(q10);
        CPPUNIT_ASSERT(q10 - q11 == Quotient<VeryLong>(VeryLong(0L)));

        q10 -= q11;
        CPPUNIT_ASSERT(q10 == Quotient<VeryLong>(VeryLong(0L)));

        q10 += q11;
        CPPUNIT_ASSERT(q10 == q11);

        Quotient<VeryLong> q12(VeryLong(13L), VeryLong(23L));
        Quotient<VeryLong> q13(VeryLong(23L), VeryLong(13L));

        q12 *= q13;
        CPPUNIT_ASSERT(q12.is_one());

        q13 *= 26L;
        CPPUNIT_ASSERT(q13 == Quotient<VeryLong>(VeryLong(46L)));

        Quotient<VeryLong> q14(VeryLong(117L), VeryLong(13L));
        Quotient<VeryLong> q15(VeryLong(13L), VeryLong(9L));
        Quotient<VeryLong> q16 = q14 * q15;
        CPPUNIT_ASSERT(q16 == Quotient<VeryLong>(VeryLong(13L)));

        CPPUNIT_ASSERT(q14 * Quotient<VeryLong>(VeryLong(0L)) == Quotient<VeryLong>(VeryLong(0L)));
        CPPUNIT_ASSERT(Quotient<VeryLong>(VeryLong(0L)) * q15 == Quotient<VeryLong>(VeryLong(0L)));
        CPPUNIT_ASSERT(Quotient<VeryLong>(VeryLong(10L), VeryLong(10L)) * q15 == q15);
        CPPUNIT_ASSERT(q15 == Quotient<VeryLong>(VeryLong(10L), VeryLong(10L)) * q15);

        CPPUNIT_ASSERT_THROW_MESSAGE("", (Quotient<VeryLong>(VeryLong(10L), VeryLong(3L)) / Quotient<VeryLong>(VeryLong(0L))), std::string);
        CPPUNIT_ASSERT(q14 / Quotient<VeryLong>(VeryLong(5L), VeryLong(5L)) == q14);

        Quotient<VeryLong> q17(VeryLong(111L), VeryLong(5L));
        Quotient<VeryLong> q18(VeryLong(222L), VeryLong(5L));
        CPPUNIT_ASSERT(q17 / q18 == Quotient<VeryLong>(VeryLong(1L), VeryLong(2L)));

/*
        std::istringstream iss("102325/2302938");
        Quotient<VeryLong> q19;
        iss >> q19;
        CPPUNIT_ASSERT(q19 == Quotient<VeryLong>(VeryLong(102325L), VeryLong(2302938L)));
        std::ostringstream oss;
        oss << q19;
        CPPUNIT_ASSERT(oss.str() == "102325/2302938");
*/

        Quotient<VeryLong> q20(VeryLong(5L), VeryLong(2L));
        Quotient<VeryLong> q21 = exp(q20, 3L);
        CPPUNIT_ASSERT(q21 == Quotient<VeryLong>(VeryLong(125L), VeryLong(8L)));

        CPPUNIT_ASSERT_THROW_MESSAGE("", (q21 /= Quotient<VeryLong>(VeryLong(0L))), std::string);
        Quotient<VeryLong> q22(VeryLong(0L));
        CPPUNIT_ASSERT_THROW_MESSAGE("", (q22 /= Quotient<VeryLong>(VeryLong(0L))), std::string);
        q21 /= q21;
        CPPUNIT_ASSERT(q21 == VeryLong(1L));
        Quotient<VeryLong> q23(VeryLong(100L), VeryLong(3L));
        q23 /= Quotient<VeryLong>(VeryLong(10L), VeryLong(3L));
        CPPUNIT_ASSERT(q23 == VeryLong(10L));
    }
};
    
int main()
{
    CppUnit::TextUi::TestRunner runner;
    runner.addTest(QuotientTest::suite());
    runner.run();

    return 0;
}
