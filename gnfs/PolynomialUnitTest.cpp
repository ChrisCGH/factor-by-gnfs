#include <iostream>
#include <sstream>
#include "Polynomial.h"
#include "VeryLong.h"
#include "LongModular.h"
#include "VeryLongModular.h"
#include "Polynomial.inl"
#include "MPFloat.h"
#include <limits>
#include <complex>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/ui/text/TestRunner.h>

class PolynomialTest : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE(PolynomialTest);
    CPPUNIT_TEST(testConstructors);
    CPPUNIT_TEST(testArithmetic);
    CPPUNIT_TEST(testFactorisation);
    CPPUNIT_TEST(testRootFinding);
    CPPUNIT_TEST(testEvaluate);
    CPPUNIT_TEST_SUITE_END();

    public:
    void setUp()
    {
    }

    void tearDown()
    {
    }

    void testConstructors()
    {
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", Polynomial<long int>());
        Polynomial<long int> p;
        CPPUNIT_ASSERT(p.deg() == -1);
        CPPUNIT_ASSERT(p.content() == 0L);
        CPPUNIT_ASSERT(p == 0L);
        CPPUNIT_ASSERT(!(p == 1L));
        CPPUNIT_ASSERT(p.is_zero());
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", p.coefficient(0));
        CPPUNIT_ASSERT(p.coefficient(0) == 0L);

        long int li(10L);
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", Polynomial<long int>(li));
        Polynomial<long int> p1(li);
        CPPUNIT_ASSERT(p1.deg() == 0);
        CPPUNIT_ASSERT(p1.content() == 10L);
        CPPUNIT_ASSERT(p1.coefficient(0) == 10L);
        CPPUNIT_ASSERT(p1 == 10L);
        CPPUNIT_ASSERT(!p1.is_zero());
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", p1.coefficient(1));
        CPPUNIT_ASSERT(p1.coefficient(1) == 0L);

        static long int c2[8] = { 1, 2, 3, 4, 5, 6, 7, 8 };
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", Polynomial<long int>(8, c2));
        Polynomial<long int> p2(8, c2);
        CPPUNIT_ASSERT(p2.deg() == 7);
        CPPUNIT_ASSERT(p2.content() == 1L);
        CPPUNIT_ASSERT(p2.coefficient(0) == 1L);
        CPPUNIT_ASSERT(p2 != p1);
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", p2.coefficient(8));
        CPPUNIT_ASSERT(p1.coefficient(8) == 0L);

        std::vector<long int> c3;
        c3.push_back(2);
        c3.push_back(4);
        c3.push_back(6);
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", Polynomial<long int>(c3));
        Polynomial<long int> p3(c3);
        CPPUNIT_ASSERT(p3.deg() == 2);
        CPPUNIT_ASSERT(p3.content() == 2L);
        CPPUNIT_ASSERT(p3.coefficient(0) == 2L);
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", p3.coefficient(3));
        CPPUNIT_ASSERT(p3.coefficient(3) == 0L);

        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", Polynomial<long int>(p3));
        Polynomial<long int> p4(p3);
        CPPUNIT_ASSERT(p4.deg() == 2);
        CPPUNIT_ASSERT(p4.content() == 2L);
        CPPUNIT_ASSERT(p4.coefficient(0) == 2L);
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", p4.coefficient(3));
        CPPUNIT_ASSERT(p3 == p4);
     
        Polynomial<long int> p5; 
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", p5 = p3);
        CPPUNIT_ASSERT(p5.deg() == 2);
        CPPUNIT_ASSERT(p5.content() == 2L);
        CPPUNIT_ASSERT(p5.coefficient(0) == 2L);
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", p5.coefficient(3));
        CPPUNIT_ASSERT(p3 == p5);

        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", -p5);
        Polynomial<long int> p6 = -p5;
        CPPUNIT_ASSERT(p6.deg() == 2);
        CPPUNIT_ASSERT(p6.content() == 2L);
        CPPUNIT_ASSERT(p6.coefficient(0) == -2L);
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", p6.coefficient(3));
       
        CPPUNIT_ASSERT(!(p5 < 1));
        // p2 = 1 + 2 X + 3 X^2 + 4 X^3 + 5 X^4 + 6 X^5 + 7 X^6 + 8 X^7
        // p5 = 2 + 4 X + 6 X^2
        // p6 = -2 - 4 X - 6 X^2

        CPPUNIT_ASSERT(p5 < p2);
        CPPUNIT_ASSERT(p6 < p5);

        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", Polynomial<VeryLong> p7 = Polynomial<VeryLong>::read_polynomial("1 - 2 X + 3 X^2"));
        Polynomial<VeryLong> p7 = Polynomial<VeryLong>::read_polynomial("-3 + 2 X + 1 X^2");
        Polynomial<VeryLong> p8 = Polynomial<VeryLong>::read_polynomial("3 + 1 X") * Polynomial<VeryLong>::read_polynomial("-1 + 1 X");
        CPPUNIT_ASSERT(p7 == p8);

        std::istringstream iss("-3 + 2 X + 1 X^2");
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", Polynomial<VeryLong> p9 = Polynomial<VeryLong>::read_polynomial(iss));
        std::istringstream iss1("-3 + 2 X + 1 X^2");
        Polynomial<VeryLong> p9 = Polynomial<VeryLong>::read_polynomial(iss1);
        CPPUNIT_ASSERT(p9 == p8);

        std::ostringstream oss;
        oss << "2395827 + 295827 X - 2934517951 X^2 + 2501953 X^3 + 63609 X^4 - 20024092 X^5";
        std::istringstream iss2(oss.str());
        Polynomial<VeryLong> p10 = Polynomial<VeryLong>::read_polynomial(iss2);
        std::ostringstream oss1;
        oss1 << p10;
        CPPUNIT_ASSERT(oss.str() == oss1.str());

        Polynomial<VeryLong> p11 = Polynomial<VeryLong>::read_polynomial("");
        CPPUNIT_ASSERT(p11 == Polynomial<VeryLong>());

        Polynomial<VeryLong> p12 = Polynomial<VeryLong>::read_polynomial("   ");
        CPPUNIT_ASSERT(p12 == Polynomial<VeryLong>());

        CPPUNIT_ASSERT_THROW_MESSAGE("", Polynomial<VeryLong>::read_polynomial("+"), std::string);
        CPPUNIT_ASSERT_THROW_MESSAGE("", Polynomial<VeryLong>::read_polynomial("-"), std::string);

        Polynomial<VeryLong> p13 = Polynomial<VeryLong>::read_polynomial("+ 1234   ");
        CPPUNIT_ASSERT(p13 == Polynomial<VeryLong>(1234L));

        Polynomial<VeryLong> p14 = Polynomial<VeryLong>::read_polynomial("+ X   ");
        std::vector<VeryLong> c15;
        c15.push_back(0L);
        c15.push_back(1L);
        CPPUNIT_ASSERT(p14 == Polynomial<VeryLong>(c15));

        Polynomial<VeryLong> p16 = Polynomial<VeryLong>::read_polynomial(" - X   ");
        CPPUNIT_ASSERT(p16 == -Polynomial<VeryLong>(c15));

        Polynomial<VeryLong> p17 = Polynomial<VeryLong>::read_polynomial(" + 3456 X ");
        std::vector<VeryLong> c18;
        c18.push_back(0L);
        c18.push_back(3456L);
        CPPUNIT_ASSERT(p17 == Polynomial<VeryLong>(c18));
        
        Polynomial<VeryLong> p19 = Polynomial<VeryLong>::read_polynomial(" 3456 X ");
        CPPUNIT_ASSERT(p19 == Polynomial<VeryLong>(c18));

        CPPUNIT_ASSERT_THROW_MESSAGE("", Polynomial<VeryLong>::read_polynomial("-1234X^"), std::string);
        CPPUNIT_ASSERT_THROW_MESSAGE("", Polynomial<VeryLong>::read_polynomial("-1234X^ "), std::string);
        CPPUNIT_ASSERT_THROW_MESSAGE("", Polynomial<VeryLong>::read_polynomial("-1234X1"), std::string);
        CPPUNIT_ASSERT_THROW_MESSAGE("", Polynomial<VeryLong>::read_polynomial("-1234X 1"), std::string);
        CPPUNIT_ASSERT_THROW_MESSAGE("", Polynomial<VeryLong>::read_polynomial("-1234/X 1"), std::string);
        CPPUNIT_ASSERT_THROW_MESSAGE("", Polynomial<VeryLong>::read_polynomial("-12/34/X 1"), std::string);

        Polynomial<VeryLong> p20 = Polynomial<VeryLong>::read_polynomial(" 3456 X^5 ");
        std::vector<VeryLong> c21;
        c21.push_back(0L);
        c21.push_back(0L);
        c21.push_back(0L);
        c21.push_back(0L);
        c21.push_back(0L);
        c21.push_back(3456L);
        CPPUNIT_ASSERT(p20 == Polynomial<VeryLong>(c21));

        Polynomial<VeryLong> p22 = Polynomial<VeryLong>::read_polynomial(" 3456 X^5 + 123 - 456 X^3 - X");
        std::vector<VeryLong> c23;
        c23.push_back(123L);
        c23.push_back(-1L);
        c23.push_back(0L);
        c23.push_back(-456L);
        c23.push_back(0L);
        c23.push_back(3456L);
        CPPUNIT_ASSERT(p22 == Polynomial<VeryLong>(c23));

        Polynomial<VeryLong> p24 = Polynomial<VeryLong>::read_polynomial(" 3456 X^5 + 123 - 456 X^3 - X +456X^3");
        std::vector<VeryLong> c25;
        c25.push_back(123L);
        c25.push_back(-1L);
        c25.push_back(0L);
        c25.push_back(0L);
        c25.push_back(0L);
        c25.push_back(3456L);
        CPPUNIT_ASSERT(p24 == Polynomial<VeryLong>(c25));

        CPPUNIT_ASSERT_THROW_MESSAGE("", Polynomial<VeryLong>::read_polynomial("23f2oowie"), std::string);
        Polynomial<VeryLong> p26 = Polynomial<VeryLong>::read_polynomial("a^5 - a^4");
        std::vector<VeryLong> c27;
        c27.push_back(0L);
        c27.push_back(0L);
        c27.push_back(0L);
        c27.push_back(0L);
        c27.push_back(-1L);
        c27.push_back(1L);
        CPPUNIT_ASSERT(p26 == Polynomial<VeryLong>(c27));
        CPPUNIT_ASSERT_THROW_MESSAGE("", Polynomial<VeryLong>::read_polynomial("a^5 - b^4"), std::string);
        CPPUNIT_ASSERT_THROW_MESSAGE("", Polynomial<VeryLong>::read_polynomial("1 + a^5 - b^4"), std::string);

        Polynomial<long int> p30 = Polynomial<long int>::read_polynomial("1 - 2 X + X^2");
        std::vector<long int> c31;
        c31.push_back(1L);
        c31.push_back(-2L);
        c31.push_back(1L);
        CPPUNIT_ASSERT(p30 == Polynomial<long int>(c31));

        Polynomial<double> p32 = Polynomial<double>::read_polynomial("1 - 2 X + X^2");
        std::vector<double> c33;
        c33.push_back(1);
        c33.push_back(-2);
        c33.push_back(1L);
        CPPUNIT_ASSERT(p32 == Polynomial<double>(c33));

        Polynomial<double> p34 = Polynomial<double>::read_polynomial("1.0 - 2234.345 X + 0.123e+10 X^2");
        Polynomial<double> p34a = Polynomial<double>::read_polynomial("1.0 - 2234.345 X + 0.123E+10 X^2");
        std::vector<double> c35;
        c35.push_back(1.0);
        c35.push_back(-2234.345);
        c35.push_back(0.123e10);
        CPPUNIT_ASSERT(p34 == Polynomial<double>(c35));
        CPPUNIT_ASSERT(p34a == Polynomial<double>(c35));

        Polynomial<double> p36 = Polynomial<double>::read_polynomial("1.0 - 2234.345 e + 0.123e^2");
        std::vector<double> c37;
        c37.push_back(1.0);
        c37.push_back(-2234.345);
        c37.push_back(0.123);
        CPPUNIT_ASSERT(p36 == Polynomial<double>(c37));

        LongModular::set_default_modulus(12347L);
        Polynomial<LongModular> p38 = Polynomial<LongModular>::read_polynomial("1 + X + X^2");
        std::vector<LongModular> c39;
        c39.push_back(1L);
        c39.push_back(1L);
        c39.push_back(1L);
        CPPUNIT_ASSERT(p38 == Polynomial<LongModular>(c39));

        Polynomial<LongModular> p40 = Polynomial<LongModular>::read_polynomial("12342326 + 23423X + 23526262  X^2");
        std::vector<LongModular> c41;
        c41.push_back(12342326L);
        c41.push_back(23423L);
        c41.push_back(23526262L);
        CPPUNIT_ASSERT(p40 == Polynomial<LongModular>(c41));

        VeryLongModular::set_default_modulus(VeryLong("12356985701"));
        Polynomial<VeryLongModular> p42 = Polynomial<VeryLongModular>::read_polynomial("1 + X + X^2");
        std::vector<VeryLongModular> c43;
        c43.push_back(1L);
        c43.push_back(1L);
        c43.push_back(1L);
        CPPUNIT_ASSERT(p42 == Polynomial<VeryLongModular>(c43));

        Polynomial<VeryLongModular> p44 = Polynomial<VeryLongModular>::read_polynomial("-1 - X - 234 X^2");
        std::vector<VeryLongModular> c45;
        c45.push_back(-1L);
        c45.push_back(-1L);
        c45.push_back(-234L);
        CPPUNIT_ASSERT(p44 == Polynomial<VeryLongModular>(c45));
        
        Polynomial<Quotient<long int> > p46 = Polynomial<Quotient<long int> >::read_polynomial("1 + X + X^2");
        std::vector<Quotient<long int> > c47;
        c47.push_back(1L);
        c47.push_back(1L);
        c47.push_back(1L);
        CPPUNIT_ASSERT(p46 == Polynomial<Quotient<long int> >(c47));

        Polynomial<Quotient<long int> > p48 = Polynomial<Quotient<long int> >::read_polynomial("1/2 + X + 2/3 X^2");
        std::vector<Quotient<long int> > c49;
        c49.push_back(Quotient<long int>(1L, 2L));
        c49.push_back(1L);
        c49.push_back(Quotient<long int>(2L, 3L));
        CPPUNIT_ASSERT(p48 == Polynomial<Quotient<long int> >(c49));

        CPPUNIT_ASSERT_THROW_MESSAGE("", Polynomial<VeryLong>::read_polynomial("1/2 + Y - 2/ Y"), std::string);
        CPPUNIT_ASSERT_THROW_MESSAGE("", Polynomial<VeryLong>::read_polynomial("1/2 + Y - 2/Y"), std::string);

        Polynomial<Quotient<VeryLong> > p50 = Polynomial<Quotient<VeryLong> >::read_polynomial("12247572935827555/235036038027 + X + 25080803805845555/320580285025028058555 X^2");
        std::vector<Quotient<VeryLong> > c51;
        c51.push_back(Quotient<VeryLong>(VeryLong("12247572935827555"), VeryLong("235036038027"))); 
        c51.push_back(Quotient<VeryLong>(1L));
        c51.push_back(Quotient<VeryLong>(VeryLong("5016160761169111"), VeryLong("64116057005005611711")));
        CPPUNIT_ASSERT(p50 == Polynomial<Quotient<VeryLong> >(c51));

        Polynomial<Quotient<VeryLong> > p52 = Polynomial<Quotient<VeryLong> >::read_polynomial("1/7 + 32730265/42 a + 275165447/105 a^2 + 7409138/35 a^3 + 12818 a^4");
        std::vector<Quotient<VeryLong> > c53;
        c53.push_back(Quotient<VeryLong>(VeryLong("1"), VeryLong("7"))); 
        c53.push_back(Quotient<VeryLong>(VeryLong("32730265"), VeryLong("42"))); 
        c53.push_back(Quotient<VeryLong>(VeryLong("275165447"), VeryLong("105"))); 
        c53.push_back(Quotient<VeryLong>(VeryLong("7409138"), VeryLong("35"))); 
        c53.push_back(Quotient<VeryLong>(VeryLong("12818"))); 
        CPPUNIT_ASSERT(p52 == Polynomial<Quotient<VeryLong> >(c53));

        Polynomial<Quotient<VeryLong> > p54 = Polynomial<Quotient<VeryLong> >::read_polynomial("40458 alpha + 2691780 alpha^2");
        std::vector<Quotient<VeryLong> > c55;
        c55.push_back(Quotient<VeryLong>(VeryLong(0L))); 
        c55.push_back(Quotient<VeryLong>(VeryLong("40458"))); 
        c55.push_back(Quotient<VeryLong>(VeryLong("2691780"))); 
        CPPUNIT_ASSERT(p54 == Polynomial<Quotient<VeryLong> >(c55));

        Polynomial<VeryLong> p59 = Polynomial<VeryLong>::read_polynomial("4 + X");
        Polynomial<VeryLong> p60 = Polynomial<VeryLong>::read_polynomial("1 + 4 X");
        CPPUNIT_ASSERT(p60 < p59 || p59 < p60 || p59 == p60);
        CPPUNIT_ASSERT((!(p60 < p59)));

        VeryLongModular::set_default_modulus(VeryLong("5"));
        Polynomial<VeryLongModular> p56 = Polynomial<VeryLongModular>::read_polynomial("4 + X");
        Polynomial<VeryLongModular> p57 = Polynomial<VeryLongModular>::read_polynomial("1 + 4 X");
        Polynomial<VeryLongModular> p58 = gcd(p56, p57);
        std::cout << "p58 = " << p58 << std::endl;
        CPPUNIT_ASSERT(p58 == p57);
    }

    void testArithmetic()
    {
        std::vector<long int> c1;
        c1.push_back(1L);
        c1.push_back(2L);
        c1.push_back(3L);
        c1.push_back(4L);
        Polynomial<long int> p1(c1);

        Polynomial<long int> p1a(p1);

        std::vector<long int> c2;
        c2.push_back(4L);
        c2.push_back(3L);
        c2.push_back(2L);
        c2.push_back(1L);
        Polynomial<long int> p2(c2);

        std::vector<long int> c3;
        c3.push_back(5L);
        c3.push_back(5L);
        c3.push_back(5L);
        c3.push_back(5L);
        Polynomial<long int> p3(c3);

        p1 += p2;
        CPPUNIT_ASSERT(p1 == p3);

        p1 -= p2;
        CPPUNIT_ASSERT(p1 == p1a);
        
        std::vector<long int> c4;
        c4.push_back(0L);
        c4.push_back(2L);
        c4.push_back(3L);
        c4.push_back(4L);
        Polynomial<long int> p4(c4);

        p1 -= p4;
        CPPUNIT_ASSERT(p1.deg() == 0);
        CPPUNIT_ASSERT(p1 == 1L);

        p1 = p1a;

        Polynomial<long int> p5 = p1 * p2;
        CPPUNIT_ASSERT(p5.deg() == 6);
        
        std::vector<long int> c7;
        c7.push_back(4L);
        c7.push_back(11L);
        c7.push_back(20L);
        c7.push_back(30L);
        c7.push_back(20L);
        c7.push_back(11L);
        c7.push_back(4L);
        Polynomial<long int> p7(c7);
        
        CPPUNIT_ASSERT(p5 == p7);

        Polynomial<long int> p8(p1);
        p8 *= p2;
        CPPUNIT_ASSERT(p8 == p7);

        Polynomial<long int> p9;
        p9 = p1 + p2;
        CPPUNIT_ASSERT(p9 == p3);

        Polynomial<long int> p10;
        p10 = p9 - p2;
        CPPUNIT_ASSERT(p10 == p1);

        Polynomial<long int> p11(p1);
        p11 -= 1L;
        std::vector<long int> c12;
        c12.push_back(0L);
        c12.push_back(2L);
        c12.push_back(3L);
        c12.push_back(4L);
        Polynomial<long int> p12(c12);
        CPPUNIT_ASSERT(p11 == p12);
    
        p11 += 1L;    
        CPPUNIT_ASSERT(p11 == p1);

        Polynomial<long int> p13;
        p13 -= 1L;
        CPPUNIT_ASSERT(p13.deg() == 0);
        p13 += 1L;
        CPPUNIT_ASSERT(p13 == 0L);

        Polynomial<long int> p14;
        p14 = p1 * 2L;
        std::vector<long int> c15;
        c15.push_back(2L);
        c15.push_back(4L);
        c15.push_back(6L);
        c15.push_back(8L);
        Polynomial<long int> p15(c15);
        CPPUNIT_ASSERT(p15 == p14);

        p14 = 2L * p1;
        CPPUNIT_ASSERT(p15 == p14);
        
        Polynomial<long int> p16 = p14 / 2L;
        CPPUNIT_ASSERT(p16 == p1);

        std::vector<long int> c17;
        c17.push_back(1L);
        c17.push_back(2L);
        c17.push_back(1L);
        Polynomial<long int> p17(c17);

        std::vector<long int> c18;
        c18.push_back(1L);
        c18.push_back(1L);
        Polynomial<long int> p18(c18);

        Polynomial<long int> p19 = p17 / p18;
        CPPUNIT_ASSERT(p19 == p18);

        Polynomial<long int> p20(p17);
        p20 /= p18;
        CPPUNIT_ASSERT(p19 == p18);

        std::vector<long int> c21;
        c21.push_back(1L);
        c21.push_back(2L);
        c21.push_back(2L);
        Polynomial<long int> p21(c21);
        Polynomial<long int> p22 = p21 / p18;
        std::vector<long int> c23;
        c23.push_back(0L);
        c23.push_back(2L);
        Polynomial<long int> p23(c23);
        CPPUNIT_ASSERT(p22 == p23);

        p20 = p21;
        p20 /= p18;
        CPPUNIT_ASSERT(p20 == p23);

        Polynomial<long int> p24 = p18 / p21;
        CPPUNIT_ASSERT(p24 == 0L);

        std::vector<VeryLong> c25;
        c25.push_back(1L);
        c25.push_back(2L);
        c25.push_back(3L);
        c25.push_back(4L);
        c25.push_back(5L);
        c25.push_back(6L);
        Polynomial<VeryLong> p25(c25);

        std::vector<VeryLong> c26;
        c26.push_back(1L);
        c26.push_back(2L);
        c26.push_back(3L);
        Polynomial<VeryLong> p26(c26);

        Polynomial<VeryLong> q27;
        Polynomial<VeryLong> r27;
        pseudo_divide(p25, p26, q27, r27);
        CPPUNIT_ASSERT(q27 == Polynomial<VeryLong>::read_polynomial("48 + 36 X + 27 X^2 + 162 X^3"));
        CPPUNIT_ASSERT(r27 == Polynomial<VeryLong>::read_polynomial("33 + 30 X"));

        Polynomial<long int> p30;
        p30.divide_by_X();
        CPPUNIT_ASSERT(p30 == 0L);
        Polynomial<VeryLong> p31(1234L);
        p31.divide_by_X();
        CPPUNIT_ASSERT(p31 == 0L);

        Polynomial<VeryLong> p32 = Polynomial<VeryLong>::read_polynomial("2 + 3 X - 4 X^2");
        p32.divide_by_X();
        CPPUNIT_ASSERT(p32 == Polynomial<VeryLong>::read_polynomial("3 - 4 X"));

        Polynomial<VeryLong> p33 = Polynomial<VeryLong>::read_polynomial("2395827 + 295827 X - 2934517951 X^2 + 2501953 X^3 + 63609 X^4 - 20024092 X^5");
        Polynomial<VeryLong> p34 = Polynomial<VeryLong>::read_polynomial("123 + 227 X - 934951 X^2");

        CPPUNIT_ASSERT_THROW_MESSAGE("", p33 % p34, std::string);
        Polynomial<VeryLong> p35 = remainder(p33, p34);
        CPPUNIT_ASSERT(p35 == Polynomial<VeryLong>::read_polynomial("1535682626909597333575400626731 - 318118391007329613026988361751 X"));

        Polynomial<VeryLong> p36 = Polynomial<VeryLong>::read_polynomial("123 + 227 X - 455093 X^2");
        CPPUNIT_ASSERT_THROW_MESSAGE("", p33 % p36, std::string);

        Polynomial<Quotient<VeryLong> > p37 = convert_to_quotient_field(p33);
        Polynomial<Quotient<VeryLong> > p38 = convert_to_quotient_field(p34);
        Polynomial<Quotient<VeryLong> > p39 = p37 % p38;
        std::vector<Quotient<VeryLong> > c40;
        c40.push_back(Quotient<VeryLong>("1535682626909597333575400626731","764109152745145348504801"));
        c40.push_back(Quotient<VeryLong>("-318118391007329613026988361751","764109152745145348504801"));
        CPPUNIT_ASSERT(p39 == Polynomial<Quotient<VeryLong> >(c40));

        VeryLongModular::set_default_modulus(1234577L);
        const VeryLong p41(1234577L);
        Polynomial<VeryLongModular> p42 = convert_to_F_p<VeryLong, VeryLong, VeryLongModular>(p33, p41);

        std::vector<VeryLongModular> c43;
        c43.push_back(VeryLongModular("1161250"));
        c43.push_back(VeryLongModular("295827"));
        c43.push_back(VeryLongModular("71578"));
        c43.push_back(VeryLongModular("32799"));
        c43.push_back(VeryLongModular("63609"));
        c43.push_back(VeryLongModular("963717"));
        Polynomial<VeryLongModular> p43(c43);
        CPPUNIT_ASSERT(p42 == p43);

        VeryLong power("12345");

        Polynomial<VeryLongModular> p44 = convert_to_F_p<VeryLong, VeryLong, VeryLongModular>(p34, p41);

        Polynomial<VeryLongModular> p45 = powmodf_vl(power, p42, p44);
        Polynomial<VeryLong> p46 = Polynomial<VeryLong>::read_polynomial("1078629 + 1023579 X");
        Polynomial<VeryLong> p47 = lift<VeryLong, VeryLongModular>(p45);
        CPPUNIT_ASSERT(p46 == p47);

        Polynomial<VeryLongModular> p48 = powpowmodf(power, 127L, p42, p44);
        Polynomial<VeryLong> p49 = Polynomial<VeryLong>::read_polynomial("284392 + 10173 X");
        CPPUNIT_ASSERT(p49 == (lift<VeryLong, VeryLongModular>(p48)));
    }

    void testFactorisation()
    {
        std::vector<VeryLong> c1;
        c1.push_back(1L);
        c1.push_back(2L);
        c1.push_back(1L);
        Polynomial<VeryLong> p1(c1);
        std::vector<Polynomial<VeryLong> > factors1;
        VeryLong cont1;
        Polynomial<VeryLong>::factor(p1, factors1, cont1);
        CPPUNIT_ASSERT(factors1.size() == 2);
        CPPUNIT_ASSERT(cont1 == VeryLong(1L));

        std::vector<VeryLong> c2a;
        c2a.push_back(1L);
        c2a.push_back(1L); 
        Polynomial<VeryLong> p2a(c2a); // X + 1

        CPPUNIT_ASSERT(factors1[0] == p2a);
        CPPUNIT_ASSERT(factors1[1] == p2a);

        std::vector<VeryLong> c2;
        c2.push_back(-1L);
        c2.push_back(1L); 
        Polynomial<VeryLong> p2(c2); // X - 1

        std::vector<VeryLong> c3;
        c3.push_back(-2L);
        c3.push_back(1L); 
        Polynomial<VeryLong> p3(c3); // X - 2

        std::vector<VeryLong> c4;
        c4.push_back(-23482374L);
        c4.push_back(583469L); 
        Polynomial<VeryLong> p4(c4); // 583469 X - 23482374 

        std::vector<VeryLong> c5;
        c5.push_back(1L);
        c5.push_back(1L);
        c5.push_back(1L);
        Polynomial<VeryLong> p5(c5); // X^2 + X + 1 

        Polynomial<VeryLong> p6 = p2 * p3 * p4 * p5;
        // -46964748 + 24649312 X - 583469 X^2 + 46964748 X^3 - 24649312 X^4 + 583469 X^5

        std::vector<VeryLong> c7;
        c7.push_back("-46964748");
        c7.push_back("24649312");
        c7.push_back("-583469");
        c7.push_back("46964748");
        c7.push_back("-24649312");
        c7.push_back("583469");
        Polynomial<VeryLong> p7(c7);
        
        CPPUNIT_ASSERT(p6 == p7);

        std::vector<Polynomial<VeryLong> > factors6;
        VeryLong cont6;
        Polynomial<VeryLong>::factor(p6, factors6, cont6);

        CPPUNIT_ASSERT(factors6.size() == 4);
  
        CPPUNIT_ASSERT(cont6 == VeryLong(1L));
        CPPUNIT_ASSERT(factors6[0] == p3);
        CPPUNIT_ASSERT(factors6[1] == p4);
        CPPUNIT_ASSERT(factors6[2] == p2);
        CPPUNIT_ASSERT(factors6[3] == p5);
       
        Polynomial<VeryLong> p33 = Polynomial<VeryLong>::read_polynomial("2395827 + 295827 X - 2934517951 X^2 + 2501953 X^3 + 63609 X^4 - 20024092 X^5");

        std::vector<std::pair<int, Polynomial<VeryLongModular> > > Ai;
        squarefree_factorisation<VeryLong, VeryLongModular>(convert_to_F_p<VeryLong, VeryLong, VeryLongModular>(p33, VeryLong(7L)), VeryLong(7L), Ai);
        CPPUNIT_ASSERT(Ai.size() == 2);
        CPPUNIT_ASSERT(Ai[0].first == 1);
        CPPUNIT_ASSERT(lift<VeryLong>(Ai[0].second) == Polynomial<VeryLong>::read_polynomial("4 + 1 X"));
        CPPUNIT_ASSERT(Ai[1].first == 2);
        // Need to fix read_polynomial()
        CPPUNIT_ASSERT(lift<VeryLong>(Ai[1].second) == Polynomial<VeryLong>::read_polynomial("5 X + 1 X^2"));
    }

    void testRootFinding()
    {
        Polynomial<VeryLong> p1 = Polynomial<VeryLong>::read_polynomial("23487236 + 92487 X - 982375 X^2 + 29592 X^3");
        long int p = 127L;
        std::vector<LongModular> roots;
        find_roots_mod_p<VeryLong, long int, LongModular>(p1, p, roots); 
        CPPUNIT_ASSERT(roots.size() == 1);
        CPPUNIT_ASSERT(roots[0].get_long() == 31L);

        std::vector<Polynomial<VeryLongModular> > factors;
        factor_over_F_p<VeryLong, VeryLong, VeryLongModular>(p1, p, factors);
        CPPUNIT_ASSERT(factors.size() == 2);
        std::vector<VeryLongModular> c2;
        VeryLongModular::set_default_modulus(127L);
        c2.push_back(96L);
        c2.push_back(1L);
        CPPUNIT_ASSERT(factors[0] == Polynomial<VeryLongModular>(c2));

        std::vector<VeryLongModular> c3;
        c3.push_back(62L);
        c3.push_back(1L);
        c3.push_back(1L);
        CPPUNIT_ASSERT(factors[1] == Polynomial<VeryLongModular>(c3));

        p = 12347L;
        roots.clear();
        find_roots_mod_p<VeryLong, long int, LongModular>(p1, p, roots); 
        CPPUNIT_ASSERT(roots.size() == 3);
        CPPUNIT_ASSERT(roots[0].get_long() == 7032L);
        CPPUNIT_ASSERT(roots[1].get_long() == 2511L);
        CPPUNIT_ASSERT(roots[2].get_long() == 4207L);

        factors.clear();
        factor_over_F_p<VeryLong, VeryLong, VeryLongModular>(p1, p, factors);
        CPPUNIT_ASSERT(factors.size() == 3);
        std::vector<VeryLongModular> c4;
        VeryLongModular::set_default_modulus(12347L);
        c4.push_back(8758L);
        c4.push_back(10719L);
        CPPUNIT_ASSERT(factors[0] == Polynomial<VeryLongModular>(c4));

        std::vector<VeryLongModular> c5;
        c5.push_back(10017L);
        c5.push_back(7298L);

        std::vector<VeryLongModular> c6;
        c6.push_back(11500L);
        c6.push_back(792L);

        CPPUNIT_ASSERT(factors[1] == Polynomial<VeryLongModular>(c6));
        CPPUNIT_ASSERT(factors[2] == Polynomial<VeryLongModular>(c5));

        Polynomial<VeryLong> p2000 = Polynomial<VeryLong>::read_polynomial("9401401242022932575419604204700 + 623872110646578368801362410 X - 8591401659640532521423 X^2 - 222405543007291322 X^3 + 1249913278668 X^4 + 2691780 X^5");
        static long int ppp[] = { 
2924023, 2924177, 2924191, 2924321, 2924371, 2924431, 2924437, 2924497, 2924609,
2924827, 2924927, 2924953, 2925001, 2925077, 2925127, 2925137, 2925193, 2925359,
2925367, 2925389, 2925437, 2925463, 2925511, 2925521, 2925523, 2925551, 2925851,
2926039, 2926111, 2926349, 2926421, 2926661, 2926687, 2926699, 2926799, 2926801,
2926823, 2926897, 2926919, 2927131, 2927219, 2927293, 2927339, 2927531, 2927663,
2927713, 2927777, 2927801, 2927867, 2928017, 2928113, 2928137, 2928227, 2928307,
2928517, 2928713, 2928817, 2928823, 2928829, 2928929, 2929373, 2929523, 2929571,
2929603, 2929657, 2929669, 2929939, 2929943, 2929961, 2929963, 2929973, 2930119,
2930153, 2930173, 2930449, 2930507, 2930519, 2930563, 2930699, 2930723, 2930737,
2931083, 2931127, 2931427, 2931443, 2931647, 2931827, 2931829, 2931883, 2931919,
2931941, 2932003, 2932309, 2932373, 2932429, 2932493, 2932511, 2932591, 2932691,
2932757, 2932777, 2932871, 2932933, 2932981, 2933083, 2933093, 2933171, 2933279,
2933347, 2933453, 2933759, 2933803, 2933849, 2934073, 2934233, 2934301, 2934473,
2934629, 2934641, 2934671, 2934703, 2934881, 2934917, 2935013, 2935027, 2935057,
2935189, 2935379, 2935417, 2935531, 2935727, 2935769, 2935993, 2936023, 2936159,
2936177, 2936273, 2936369, 2936407, 2936441, 2936491, 2936551, 2936711, 2936719,
2936807, 2936821, 2936977, 2936993, 2937173, 2937217, 2937283, 2937463, 2937521,
2937523, 2937581, 2937679, 2937731, 2937751, 2937791, 2937817, 2937967, 2938009,
2938037, 2938279, 2938493, 2938517, 2938543, 2938601, 2938609, 2938627, 2938667,
2938801, 2938853, 2938987, 2939081, 2939171, 2939219, 2939473, 2939549, 2939597,
2939633, 2939693, 2939707, 2939749, 2939753, 2939941, 2940061, 2940073, 2940083,
2940263, 2940271, 2940521, 2940601, 2940689, 2940863, 2940911, 2940989, 2941007,
2941031, 2941039, 2941207, 2941319, 2941339, 2941387, 2941573, 2941733, 2941843,
2941849, 2941859, 2941877, 2941943, 2942011, 2942113, 2942141, 2942257, 2942399,
2942441, 2942477, 2942519, 2942767, 2942809, 2943121, 2943191, 2943289, 2943301,
2943361, 2943403, 2943467, 2943503, 2943653, 2943671, 2943767, 2943781, 2943821,
2943859, 2943869, 2943887, 2943911, 2943929, 2944133, 2944219, 2944243, 2944289,
2944301, 2944321, 2944429, 2944433, 2944499, 2944547, 2944589, 2944673, 2944807,
2944999, 2945021, 2945167, 2945177, 2945191, 2945611, 2945707, 2945749, 2945869,
2945959, 2945983, 2945993, 2946211, 2946241, 2946257, 2946259, 2946331, 2946679,
2946683, 2946883, 2947097, 2947111, 2947117, 2947123, 2947159, 2947169, 2947339,
2947471, 2947499, 2947537, 2947709, 2947723, 2947729, 2947823, 2947853, 2947939,
2948009, 2948063, 2948333, 2948369, 2948389, 2948411, 2948527, 2948683, 2948711,
2948779, 2948887, 2948951, 2948977, 2948987, 2949013, 2949059, 2949169, 2949223,
2949343, 2949361, 2949539, 2949577, 2949643, 2949679, 2949733, 2949803, 2949839,
2949841, 2949901, 2949953, 2950037, 2950393, 2950447, 2950511, 2950579, 2950589,
2950609, 2950667, 2950697, 2950817, 2950853, 2950873, 2950939, 2950949, 2951087,
2951161, 2951233, 2951279, 2951287, 2951297, 2951413, 2951483, 2951497, 2951537,
2951579, 2951717, 2951719, 2951777, 2951779, 2951801, 2951813, 2951849, 2951983,
2952017, 2952029, 2952133, 2952143, 2952179, 2952241, 2952311, 2952317, 2952319,
2952371, 2952377, 2952421, 2952427, 2952601, 2952617, 2952629, 2952749, 2952769,
2952863, 2953117, 2953129, 2953231, 2953411, 2953459, 2953501, 2953589, 2953609,
2953619, 2953697, 2953931, 2953961, 2953981, 2954069, 2954093, 2954113, 2954129,
2954333, 2954353, 2954389, 2954437, 2954599, 2954681, 2954737, 2954863, 2954947,
2955107, 2955131, 2955307, 2955319, 2955397, 2955497, 2955541, 2955629, 2955683,
2955761, 2955781, 2955817, 2955859, 2955907, 2955929, 2956021, 2956061, 2956111,
2956237, 2956363, 2956553, 2956651, 2956763, 2956769, 2956813, 2956819, 2956931,
2956969, 2957189, 2957191, 2957231, 2957237, 2957261, 2957321, 2957413, 2958023,
2958173, 2958331, 2958359, 2958413, 2958521, 2958539, 2958563, 2958887, 2958941,
2958971, 2959031, 2959081, 2959219, 2959393, 2959549, 2959597, 2959609, 2959643,
2959753, 2959771, 2959841, 2959861, 2959951, 2960017, 2960101, 2960231, 2960273,
2960327, 2960341, 2960393, 2960437, 2960603, 2960609, 2960809, 2960821, 2960869,
2960891, 2960957, 2961053, 2961113, 2961253, 2961271, 2961533, 2961583, 2961589,
2961619, 2961677, 2961737, 2961773, 2961787, 2961821, 2961877, 2961923, 2961979,
2961983, 2962033, 2962039, 2962093, 2962159, 2962163, 2962207, 2962277, 2962327,
2962381, 2962411, 2962441, 2962451, 2962549, 2962871, 2962891, 2962931, 2963203,
2963269, 2963299, 2963333, 2963377, 2963399, 2963453, 2963599, 2963669, 2963671,
2963809, 2963833, 2963911, 2963977, 2964041, 2964097, 2964289, 2964331, 2964347,
2964349, 2964433, 2964461, 2964719, 2964827, 2964869, 2964971, 2965037, 2965093,
2965141, 2965223, 2965247, 2965307, 2965351, 2965379, 2965421, 2965537, 2965561,
2965681, 2965717, 2965733, 2965763, 2965777, 2965793, 2965901, 2965951, 2965961,
2966017, 2966219, 2966329, 2966407, 2966507, 2966531, 2966633, 2966681, 2966773,
2966911, 2966987, 2967011, 2967049, 2967187, 2967199, 2967247, 2967269, 2967277,
2967373, 2967383, 2967389, 2967409, 2967709, 2967749, 2967953, 2967977, 2968193,
2968261, 2968337, 2968391, 2968517, 2968519, 2968711, 2968781, 2968859, 2968891,
2969023, 2969257, 2969333, 2969699, 2969797, 2969809, 2969971, 2970047, 2970157,
2970193, 2970277, 2970467, 2970599, 2970661, 2970701, 2970769, 2970893, 2971081,
2971093, 2971121, 2971123, 2971183, 2971193, 2971223, 2971259, 2971469, 2971489,
2971559, 2971597, 2971603, 2971663, 2971699, 2971879, 2971987, 2972027, 2972069,
2972119, 2972219, 2972287, 2972429, 2972491, 2972503, 2972663, 2972687, 2972791,
2972803, 2972881, 2972899, 2973017, 2973029, 2973031, 2973059, 2973073, 2973119,
2973121, 2973151, 2973331, 2973359, 2973433, 2973437, 2973449, 2973497, 2973563,
2973617, 2973667, 2973703, 2973787, 2973833, 2973857, 2973889, 2974033, 2974121,
2974129, 2974159, 2974187, 2974219, 2974289, 2974351, 2974369, 2974561, 2974649,
2974669, 2974817, 2974859, 2975057, 2975129, 2975149, 2975279, 2975339, 2975383,
2975417, 2975513, 2975677, 2975801, 2975813, 2976089, 2976097, 2976199, 2976209,
2976257, 2976401, 2976409, 2976433, 2976529, 2976541, 2976551, 2976643, 2976719,
2976971, 2977067, 2977171, 2977283, 2977367, 2977391, 2977417, 2977549, 2977697,
2977729, 2977817, 2977831, 2977841, 2977937, 2977963, 2978119, 2978219, 2978233,
2978251, 2978273, 2978323, 2978401, 2978419, 2978681, 2978797, 2978923, 2978959,
2979049, 2979107, 2979149, 2979239, 2979299, 2979311, 2979343, 2979359, 2979407,
2979491, 2979497, 2979611, 2979623, 2979679, 2979701, 2979737, 2979877, 2979919,
2980007, 2980027, 2980051, 2980139, 2980319, 2980333, 2980391, 2980441, 2980697,
2980721, 2980759, 2980799, 2980819, 2980883, 2980921, 2980993, 2981023, 2981057,
2981059, 2981087, 2981123, 2981179, 2981263, 2981333, 2981339, 2981371, 2981387,
2981471, 2981527, 2981773, 2981837, 2981897, 2981981, 2982017, 2982139, 2982149,
2982223, 2982247, 2982379, 2982431, 2982457, 2982461, 2982491, 2982527, 2982559,
2982583, 2982643, 2982647, 2982649, 2982799, 2982841, 2982869, 2982899, 2982971,
2983021, 2983081, 2983289, 2983301, 2983367, 2983411, 2983457, 2983459, 2983769,
2983777, 2983831, 2983837, 2983943, 2984209, 2984221, 2984239, 2984291, 2984297,
2984327, 2984581, 2984609, 2984689, 2984713, 2984743, 2984749, 2984771, 2984831,
2984899, 2985007, 2985061, 2985133, 2985179, 2985187, 2985253, 2985259, 2985293,
2985329, 2985391, 2985491, 2985613, 2985677, 2985809, 2985817, 2986003, 2986129,
2986171, 2986253, 2986259, 2986349, 2986409, 2986411, 2986429, 2986453, 2986561,
2986567, 2986661, 2986741, 2986799, 2986811, 2987057, 2987147, 2987287, 2987437,
2987497, 2987591, 2987827, 2988101, 2988523, 2988647, 2988841, 2988911, 2989003,
2989069, 2989087, 2989169, 2989201, 2989319, 2989367, 2989417, 2989451, 2989487,
2989537, 2989757, 2989823, 2989849, 2989927, 2989963, 2990059, 2990123, 2990137,
2990191, 2990227, 2990231, 2990443, 2990501, 2990587, 2990629, 2990791, 2990797,
2990917, 2990983, 2991199, 2991229, 2991293, 2991463, 2991631, 2991673, 2991701,
2991773, 2991803, 2991853, 2991917, 2992013, 2992123, 2992201, 2992229, 2992261,
2992279, 2992333, 2992361, 2992373, 2992567, 2992607, 2992651, 2992667, 2992679,
2992723, 2992853, 2992859, 2993047, 2993297, 2993303, 2993311, 2993399, 2993447,
2993663, 2993671, 2993687, 2993723, 2994083, 2994127, 2994311, 2994317, 2994403,
2994643, 2994701, 2994857, 2994923, 2995379, 2995463, 2995469, 2995667, 2995763,
2995799, 2995807, 2995823, 2995999, 2996041, 2996047, 2996087, 2996113, 2996183,
2996207, 2996237, 2996333, 2996387, 2996479, 2996611, 2996717, 2997067, 2997101,
2997263, 2997329, 2997509, 2997641, 2997667, 2997767, 2997773, 2997893, 2997919,
2998013, 2998097, 2998153, 2998169, 2998291, 2998319, 2998357, 2998393, 2998423,
2998573, 2998649, 2998883, 2998921, 2998991, 2999047, 2999063, 2999069, 2999099,
2999107, 2999263, 2999299, 2999459, 2999509, 2999539, 2999611, 2999629, 2999831,
2999833, 2999957
};
        const size_t ppp_count = sizeof(ppp) / sizeof(long int);
        for (size_t i = ppp_count; i < ppp_count; ++i)
        { 
            long int pp = ppp[i];
            factors.clear();
            factor_over_F_p<VeryLong, VeryLong, VeryLongModular>(p2000, VeryLong(pp), factors);
            std::vector<Polynomial<LongModular> > factors2;
            factor_over_F_p<VeryLong, long int, LongModular>(p2000, pp, factors2);
            CPPUNIT_ASSERT(factors.size() == 1);
            CPPUNIT_ASSERT(factors.size() == factors2.size());
            CPPUNIT_ASSERT(factors[0].deg() == factors2[0].deg());
            for (int j = 0; j <= factors[0].deg(); ++j)
            {
               CPPUNIT_ASSERT((long)factors[0].coefficient(j).get_very_long().get_long() == (long)factors2[0].coefficient(j).get_long());
            }
        }

        MPFloat::set_precision(100);
        Polynomial<VeryLong> p7 = Polynomial<VeryLong>::read_polynomial("-46964748 + 24649312 X - 583469 X^2 + 46964748 X^3 - 24649312 X^4 + 583469 X^5");
        std::vector<complex<MPFloat > > co(p7.deg() + 1, complex<MPFloat>(0.0, 0.0));
        for (int i = 0; i <= p7.deg(); i++)
        {
           co[i] = complex<MPFloat>((MPFloat)p7.coefficient(i).get_double(), MPFloat(0.0));
        }
        Polynomial<complex<MPFloat > > cp7(co);

        std::vector<complex<MPFloat> > roots7(0, complex<MPFloat>(0.0, 0.0));
        CPPUNIT_ASSERT(find_roots_over_C_q1<MPFloat>(cp7, roots7) == 5);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(roots7[0].real(), 1.0, std::numeric_limits<float>::epsilon());
        CPPUNIT_ASSERT_DOUBLES_EQUAL(roots7[0].imag(), 0.0, std::numeric_limits<float>::epsilon());
        CPPUNIT_ASSERT_DOUBLES_EQUAL(roots7[1].real(), -0.5, std::numeric_limits<float>::epsilon());
        CPPUNIT_ASSERT_DOUBLES_EQUAL(roots7[1].imag(), -0.8660254037, std::numeric_limits<float>::epsilon());
        CPPUNIT_ASSERT_DOUBLES_EQUAL(roots7[2].real(), -0.5, std::numeric_limits<float>::epsilon());
        CPPUNIT_ASSERT_DOUBLES_EQUAL(roots7[2].imag(), 0.8660254037, std::numeric_limits<float>::epsilon());
        CPPUNIT_ASSERT_DOUBLES_EQUAL(roots7[3].real(), 2.0, std::numeric_limits<float>::epsilon());
        CPPUNIT_ASSERT_DOUBLES_EQUAL(roots7[3].imag(), 0.0, std::numeric_limits<float>::epsilon());
        CPPUNIT_ASSERT_DOUBLES_EQUAL(roots7[4].real(), 40.24613818, std::numeric_limits<float>::epsilon());
        CPPUNIT_ASSERT_DOUBLES_EQUAL(roots7[4].imag(), 0.0, std::numeric_limits<float>::epsilon());

        std::vector<complex<MPFloat> > roots7a(0, complex<MPFloat>(0.0, 0.0));
        CPPUNIT_ASSERT(find_roots_over_C_q<MPFloat>(cp7, roots7a) == 5);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(roots7a[0].real(), 1.0, std::numeric_limits<float>::epsilon());
        CPPUNIT_ASSERT_DOUBLES_EQUAL(roots7a[0].imag(), 0.0, std::numeric_limits<float>::epsilon());
        CPPUNIT_ASSERT_DOUBLES_EQUAL(roots7a[3].real(), -0.5, std::numeric_limits<float>::epsilon());
        CPPUNIT_ASSERT_DOUBLES_EQUAL(roots7a[3].imag(), 0.8660254037, std::numeric_limits<float>::epsilon());
        CPPUNIT_ASSERT_DOUBLES_EQUAL(roots7a[2].real(), -0.5, std::numeric_limits<float>::epsilon());
        CPPUNIT_ASSERT_DOUBLES_EQUAL(roots7a[2].imag(), -0.8660254037, std::numeric_limits<float>::epsilon());
        CPPUNIT_ASSERT_DOUBLES_EQUAL(roots7a[1].real(), 2.0, std::numeric_limits<float>::epsilon());
        CPPUNIT_ASSERT_DOUBLES_EQUAL(roots7a[1].imag(), 0.0, std::numeric_limits<float>::epsilon());
        CPPUNIT_ASSERT_DOUBLES_EQUAL(roots7a[4].real(), 40.24613818, std::numeric_limits<float>::epsilon());
        CPPUNIT_ASSERT_DOUBLES_EQUAL(roots7a[4].imag(), 0.0, std::numeric_limits<float>::epsilon());

    }

    void testEvaluate()
    {
        VeryLong N("308268017314015502198864479573232259026111539470765498704205704273225473378735381489745954136566971705214233108249933319714229109554576259335341");
        Polynomial<VeryLong> f1 = Polynomial<VeryLong>::read_polynomial("-16612792256731454888860424973312 + 1498121907742379786368701480 X - 71616325394087869671190 X^2 - 436671758051847283 X^3 + 25011687432982 X^4 + 23794992 X^5");
        VeryLong a("11774494093137707");
        VeryLong b("1669121445503958099638459425"); 
        VeryLong m("254841105695668414977236991947256294566541981105672266845065891716312859988279394572523515845653926731492992971798027921010177494642224601803274");
        VeryLong result = f1.evaluate(m);
        CPPUNIT_ASSERT(result == VeryLong("25576020225388922064853273145846255804849998901928725992047982109147879231589564510262602638681713677104160277689381673293778048742285490401843226866638559137918685465209312134776963902072912189733245953513534217052245699388176744848535835322610130182780545587319034110560331020018964681187978547542075241694007992046904143150228032302996988010974417913117845614176968678454827965057250461237216969863626699205576595366239232326711891741989091443782009008279751055516044123694729224164767894647401843285109077445727384062253206087700664270610333488664282669686405978476533854335897013017962978948887925767623561306960933446068556410818685397429809083372970419320796107755350042252135204232913546920330533664000081254886627616"));
        CPPUNIT_ASSERT(result % N == VeryLong(0L));
        result = f1.evaluate_homogeneous(b, a);
        CPPUNIT_ASSERT(result == N);
        CPPUNIT_ASSERT(result % N == VeryLong(0L));
    }
};

int main()
{
    CppUnit::TextUi::TestRunner runner;
    runner.addTest(PolynomialTest::suite());
    runner.run();
    
    return 0;

}

