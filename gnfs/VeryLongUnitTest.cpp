#include <iostream>
#include <iomanip>
#include "VeryLong.h"
#include "squfof.h"
#ifdef USING_CPPUNIT
#include <limits>
#include <complex>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/ui/text/TestRunner.h>
#else
#include "UnitTest.h"
#endif

#ifdef USING_CPPUNIT
class VeryLongTest : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE( VeryLongTest );
    CPPUNIT_TEST( testConstructors );
    CPPUNIT_TEST( testFastVeryLong );
    CPPUNIT_TEST( testArithmetic );
    CPPUNIT_TEST( testMaths );
    CPPUNIT_TEST( testFactorisation );
    CPPUNIT_TEST_SUITE_END();
private:
    VeryLong v23;
    VeryLong v24;
    VeryLong v27;
    VeryLong v28;
    VeryLong v32;
public:
    void setUp()
    {
        VeryLong::clearDebug();
        v23 = VeryLong("1234567890123456789012345678901234567890");
        v24 = VeryLong("9876543210987654321098765432109876543210");
        v27 = VeryLong("1232355834509287492857913847");
        v28 = VeryLong("1521426942491462939710131061863585936186206450924386713067892572830");
        v32 = VeryLong("123123");
    }

    void tearDown()
    {
    }

    void testConstructors()
    {
        const VeryLong vl1;
        CPPUNIT_ASSERT(vl1.get_long() == 0L);
        CPPUNIT_ASSERT(vl1.is_zero());

        const VeryLong vl2("1234567890123456789");
        CPPUNIT_ASSERT(vl2.get_long_long() == 1234567890123456789LL);

        const VeryLong vl3(vl2);
        CPPUNIT_ASSERT(vl3.get_long_long() == 1234567890123456789LL);

        CPPUNIT_ASSERT(vl2 == vl3);

        std::string s("123451234512345");
        const VeryLong vl4(s);
        CPPUNIT_ASSERT(vl4.get_long_long() == 123451234512345LL);

        const VeryLong v15(-234567L);
        CPPUNIT_ASSERT(v15.get_long() == -234567L);

        const VeryLong v16(234567UL);
        CPPUNIT_ASSERT(v16.get_long() == 234567UL);

        const VeryLong v17(234567123456ULL);
        CPPUNIT_ASSERT(v17.get_long_long() == 234567123456ULL);

        const VeryLong v17a = 2230183161LL;
        CPPUNIT_ASSERT(v17a.get_long_long() == 2230183161LL);

        const VeryLong v18(-234567123456LL);
        CPPUNIT_ASSERT(v18.get_long_long() == -234567123456LL);

        const VeryLong v19(-1234567.23456);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(v19.get_double(), -1234567, std::numeric_limits<float>::epsilon());

        const VeryLong v20(-1234567798798719287192371928392.0);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(v20.get_long_double(), -1234567798798719287192371928392.0, std::numeric_limits<float>::epsilon());

        const VeryLong v21(1L);
        CPPUNIT_ASSERT(v21.is_one());

        CPPUNIT_ASSERT(v20 != v21);
        CPPUNIT_ASSERT(v20 < v21);
        CPPUNIT_ASSERT(v20 <= v21);
        CPPUNIT_ASSERT(v20 < VeryLong(-123456L));
        CPPUNIT_ASSERT(v21 > v20);
        CPPUNIT_ASSERT(v21 >= v20);

        VeryLong v22;
        v22 = v20;
        CPPUNIT_ASSERT_DOUBLES_EQUAL(v22.get_long_double(), -1234567798798719287192371928392.0, std::numeric_limits<float>::epsilon());

        long long int lli96(0x7FFFFFFFFFFFFFFFLL);
        const VeryLong v97(lli96);
        CPPUNIT_ASSERT(lli96 == v97.get_long_long());

        long long int lli96m(0x8000000000000000LL);
        const VeryLong v97m(lli96m);
        CPPUNIT_ASSERT(lli96m == v97m.get_long_long());

        long int li98(0x7FFFFFFFL);
        const VeryLong v99(li98);
        CPPUNIT_ASSERT(li98 == v99.get_long());

        long int li98m(0x80000000L);
        const VeryLong v99m(li98m);
        CPPUNIT_ASSERT(li98m == v99m.get_long());
    }

    void testFastVeryLong()
    {
        FastVeryLong fvl1;
        CPPUNIT_ASSERT(fvl1.get_long() == 0L);

        const VeryLong v106(0x7FFFFFFFL);
        const FastVeryLong fvl2(v106);
        CPPUNIT_ASSERT(fvl2.is_faster());
        CPPUNIT_ASSERT(fvl2.is_fast());

        const FastVeryLong fvl3(fvl2);
        CPPUNIT_ASSERT(fvl2 == fvl3);

        CPPUNIT_ASSERT(!(fvl2 != v106));
    }

    void testArithmetic()
    {
        VeryLong vv1("1");
        VeryLong vv2("4");
        CPPUNIT_ASSERT(vv1 < vv2);
        CPPUNIT_ASSERT((!(vv2 < vv1)));

        v23 = VeryLong("1234567890123456789012345678901234567890");
        v24 = VeryLong("9876543210987654321098765432109876543210");
        VeryLong v25(v23);
        v25 += v24;
        const VeryLong v26("11111111101111111110111111111011111111100");
        CPPUNIT_ASSERT(v25 == v26);

        v25 -= v24;
        CPPUNIT_ASSERT(v25 == v23);

        v27 = VeryLong("1232355834509287492857913847");

        v25 *= v27;
        v28 = VeryLong("1521426942491462939710131061863585936186206450924386713067892572830");
        CPPUNIT_ASSERT(v25 == v28);

        v25 /= v27;
        CPPUNIT_ASSERT(v25 == v23);

        v25 /= 12345L;
        const VeryLong v29("100005499402467135602458135188435363");
        CPPUNIT_ASSERT(v25 == v29);

        VeryLong v30 = -v25;
        const VeryLong v31("-100005499402467135602458135188435363");
        CPPUNIT_ASSERT(v30 == v31);

        v30.negate();
        CPPUNIT_ASSERT(v30 == v29);

        v32 = VeryLong("123123");
        const VeryLong v33("23");
        VeryLong v34 = exp(v32, v33);
        const VeryLong v35("1196193243742734313072657062937753550797883517849819344399014727405544587087316341405542835997613669890217861709298267");
        CPPUNIT_ASSERT(v34 == v35);

        double d36 = ln(v35);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(d36, 269.581600097816, std::numeric_limits<float>::epsilon());

        double d37 = log10(v35);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(d37, 117.07780134513, std::numeric_limits<float>::epsilon());

        VeryLong v38 = v23 + v24;
        CPPUNIT_ASSERT(v38 == v26);

        VeryLong v39 = v38 - v24;
        CPPUNIT_ASSERT(v39 == v23);

        VeryLong v40 = v39 * v27;
        CPPUNIT_ASSERT(v40 == v28);

        VeryLong v41 = v40 / v27;
        CPPUNIT_ASSERT(v41 == v23);

        VeryLong v42 = v41 / 12345L;
        CPPUNIT_ASSERT(v42 == v29);

        v42 %= v32;
        const VeryLong v43("53367");
        CPPUNIT_ASSERT(v42 == v43);

        VeryLong v44 = v29 % v32;
        CPPUNIT_ASSERT(v44 == v43);

        long int l45 = v29 % 123123L;
        CPPUNIT_ASSERT(l45 == 53367);

        std::complex<long double> alpha(-13501.183274494388);
        //std::complex<long double> alpha_ = std::conj(alpha);
        const VeryLong a("-434022535");
        const VeryLong b("32147");
        long double delta = 1.0;

        long double aa = a.get_long_double();
        long double bb = b.get_long_double();
        std::complex<long double> x = std::complex<long double>(aa) - alpha * bb;
        std::complex<long double> x_ = std::conj(x);
        long double mod2 = (x * x_).real();
        long double answer = (long double)log((double)mod2) / 2.0 + (long double)log(delta);

        CPPUNIT_ASSERT_DOUBLES_EQUAL(answer, 1.3151127783, std::numeric_limits<float>::epsilon());
    }

    void testMaths()
    {
        int i46 = jacobi_symbol(v23, v32);
        CPPUNIT_ASSERT(i46 == 0);

        int i47 = jacobi_symbol(VeryLong("5"), VeryLong("11"));
        CPPUNIT_ASSERT(i47 == 1);

        int i48 = jacobi_symbol(VeryLong("2"), VeryLong("11"));
        CPPUNIT_ASSERT(i48 == -1);

        const VeryLong v49("12345");
        VeryLong v50(v49);
        v50 *= v50;
        v50 *= v50;
        v50 *= v50;
        v50 *= v49;
        VeryLong v51 = v50.nth_root(9);
        CPPUNIT_ASSERT(v51 == v49);

        v50 += v49;
        VeryLong v52 = v50.nth_root(9);
        CPPUNIT_ASSERT(v52 == v49);

        VeryLong v53 = gcd(v27, v28);
        CPPUNIT_ASSERT(v53 == v27);

        VeryLong v54;
        VeryLong v55;
        VeryLong v56 = extended_gcd(v23, v24, v54, v55);
        CPPUNIT_ASSERT(v54 == VeryLong(-8L));
        CPPUNIT_ASSERT(v55 == VeryLong(1L));
        CPPUNIT_ASSERT(v56 == VeryLong("90000000009000000000900000000090"));

        CPPUNIT_ASSERT((v23 * v54 + v24 * v55) == v56);

        bool b57 = v23.is_probable_prime();
        CPPUNIT_ASSERT(!b57);

        bool b58 = v27.is_probable_prime();
        CPPUNIT_ASSERT(!b58);

        const VeryLong v59("1232355834509287492857913897");
        bool b60 = v59.is_probable_prime();
        CPPUNIT_ASSERT(b60);

        VeryLong v61(v59);
        v61 *= v61;
        v61 *= v61;
        v61 *= v61;
        v61 *= v59;
        bool b62 = v61.is_prime_power();
        CPPUNIT_ASSERT(b62);

        v61 -= 1L;
        bool b63 = v61.is_prime_power();
        CPPUNIT_ASSERT(!b63);

        v61 += 1L;
        bool b64 = v61.is_square();
        CPPUNIT_ASSERT(!b64);

        v61 *= v59;
        bool b65 = v61.is_square();
        CPPUNIT_ASSERT(b65);

        const VeryLong v100("123346620910319027308203983058475023948203498720394");
        const VeryLong v101("504930582034509586048345983759385739857394520570157");

        VeryLong v102 = VeryLong::random(v100, v101);
        CPPUNIT_ASSERT(v102 >= v100 && v102 <= v101);

        const VeryLong v103(-v102);
        VeryLong v104 = abs(v103);
        CPPUNIT_ASSERT(v104 == v102);

        VeryLong v105 = v103.abs();
        CPPUNIT_ASSERT(v105 == v102);

        VeryLong::clear_prime_table();
        VeryLong::generate_prime_table();
        long int p = VeryLong::firstPrime();
        CPPUNIT_ASSERT(p == 2L);

        CPPUNIT_ASSERT(VeryLong::nextPrime() == 3L);
        CPPUNIT_ASSERT(VeryLong::nextPrime() == 5L);
        CPPUNIT_ASSERT(VeryLong::nextPrime() == 7L);
        CPPUNIT_ASSERT(VeryLong::nextPrime() == 11L);
        CPPUNIT_ASSERT(VeryLong::nextPrime() == 13L);
        CPPUNIT_ASSERT(VeryLong::nextPrime() == 17L);

        VeryLong::set_max_prime(97L);
        CPPUNIT_ASSERT(VeryLong::get_max_prime() == 97L);
    }

    void testFactorisation()
    {
        VeryLong::clear_prime_table();
        VeryLong::generate_prime_table();

        VeryLong v62(v27);
        VeryLong v63;
        VeryLong v64;
        bool b66 = v62.factorise_pollardrho(&v63, &v64);
        CPPUNIT_ASSERT(b66);
        CPPUNIT_ASSERT(v63 * v64 == v62);

        VeryLong v67;
        VeryLong v68;
        VeryLong v69(0x7ffffffffffffffdLL);
        bool b70 = v69.factorise_qs(&v67, &v68);
        CPPUNIT_ASSERT(b70);
        CPPUNIT_ASSERT(v67 * v68 == v69);

        VeryLong v71;
        VeryLong v72;
        VeryLong v73(12347LL * 34583LL);
        bool b74 = v73.factorise_squfof(&v71, &v72);
        CPPUNIT_ASSERT(b74);
        CPPUNIT_ASSERT(v71 * v72 == v73);

        VeryLong v75;
        VeryLong v76;
        bool b77 = v73.factorise_p_minus_1(&v75, &v76);
        CPPUNIT_ASSERT(b77);
        CPPUNIT_ASSERT(v75 * v76 == v73);

        VeryLong v77;
        std::vector<VeryLong> factors;
        CPPUNIT_ASSERT(v23.factorise_trial_division(&factors, &v77));
        CPPUNIT_ASSERT(factors.size() == 12);
        VeryLong v78(v77);
        for (auto& f: factors)
        {
            v78 *= f;
        }
        CPPUNIT_ASSERT(v78 == v23);

        VeryLong v79;
        std::vector<long int> ifactors;
        CPPUNIT_ASSERT(v23.factorise_trial_division(&ifactors, &v79));
        CPPUNIT_ASSERT(ifactors.size() == 12);

        VeryLong v80(v79);
        for (auto& f: ifactors)
        {
            v80 *= f;
        }
        CPPUNIT_ASSERT(v80 == v23);

        VeryLong v81;
        std::vector<long int> iifactors;
        CPPUNIT_ASSERT(v23.factorise_trial_division(&iifactors, &v81, 10000));
        CPPUNIT_ASSERT(iifactors.size() == 10);

        VeryLong v82(v81);
        for (auto& iif: iifactors)
        {
            v82 *= iif;
        }
        CPPUNIT_ASSERT(v82 == v23);

        VeryLong v83;
        VeryLong v84;
        CPPUNIT_ASSERT(v23.factorise_trial_division(&v83, &v84));
        CPPUNIT_ASSERT(v83 * v84 == v23);

        VeryLong v85;
        VeryLong v86;
        CPPUNIT_ASSERT(v73.factorise_ecm(&v85, &v86));
        CPPUNIT_ASSERT(v85 * v86 == v73);

        VeryLong v87;
        VeryLong v88;
        CPPUNIT_ASSERT(v73.factorise_ecm(&v87, &v88));
        CPPUNIT_ASSERT(v87 * v88 == v73);

        std::vector<VeryLong> factors89;
        CPPUNIT_ASSERT(v23.factorise(&factors89));
        CPPUNIT_ASSERT(factors89.size() == 13);
        VeryLong v90(1L);
        for (auto& f: factors89)
        {
            v90 *= f;
        }
        CPPUNIT_ASSERT(v90 == v23);

        std::vector<VeryLong> factors91;
        CPPUNIT_ASSERT(v73.factorise_no_trial(&factors91));
        CPPUNIT_ASSERT(factors91.size() == 2);
        if (factors91.size() == 2)
        {
            CPPUNIT_ASSERT(factors91[0] * factors91[1] == v73);
        }

        long long int lli92 = v73.get_long_long();
        std::vector<long int> ifactors93;
        CPPUNIT_ASSERT(VeryLong::factorise_no_trial_ll(lli92, &ifactors93));
        CPPUNIT_ASSERT(ifactors93.size() == 2);
        if (ifactors93.size() == 2)
        {
            CPPUNIT_ASSERT(ifactors93[0] * (long long int)ifactors93[1] == lli92);
        }

        long int li94 = 1237L * 2351L;
        std::vector<long int> ifactors95;
        CPPUNIT_ASSERT(VeryLong::factorise_no_trial_l(li94, &ifactors95));
        CPPUNIT_ASSERT(ifactors95.size() == 2);
        if (ifactors95.size() == 2)
        {
            CPPUNIT_ASSERT(ifactors95[0] * ifactors95[1] == li94);
        }

        const VeryLong v2000("24304511608315849");
        std::vector<VeryLong> v2000_factors;
        VeryLong v2000_unfactored;
        CPPUNIT_ASSERT(v2000.factorise_no_trial(&v2000_factors, v2000_unfactored));

        VeryLong v2001("68454066430929935312919949146284654064784560463655945766546376196659857165127911441901256275347366429914574181910564173262954907243781843526369433234400129532781555616626107827873263556651889521214284201244347354065119176141583589197364828686588669194047973922985132474269362219900353120931603959821770405789500417740915047315436720600435558247495565875952517451686427355414486482150319474507487");
        VeryLong v2001_factor;
        VeryLong v2001_unfactored;
        CPPUNIT_ASSERT(v2001.factorise_ecm(&v2001_factor, &v2001_unfactored));

        const VeryLong v200("5544779380905324760346515880849056979247549397556131607090256471929448430375360826794001758303136680823080508734755698034299347486746329325635924091986410492155306004946714734057734348088803051218357020300792135679274653267468270724986551123613682204717885887761795730415818339811928602795459920745563402868949533837014118832550374368635280218047140835952153913586600615788573405054175877435106447");
        std::vector<VeryLong> factors200;

        VeryLong v200_unfactored;
        CPPUNIT_ASSERT(v200.factorise(&factors200, v200_unfactored));
        VeryLong v200_product(1L);

        for (size_t i = 0; i < factors200.size(); ++i)
        {
            v200_product *= factors200[i];
            //std::cout << factors200[i] << std::endl;
        }
        //std::cout << "Unfactored : " << v200_unfactored << std::endl;
        CPPUNIT_ASSERT(v200_unfactored == VeryLong(1L));
        CPPUNIT_ASSERT(v200 == v200_product);
    }
};
#endif
int main()
{
#ifdef USING_CPPUNIT
    CppUnit::TextUi::TestRunner runner;
    runner.addTest(VeryLongTest::suite());
    runner.run();
#else
    UnitTest t;

    VeryLong::clearDebug();


    const VeryLong vl1;
    t.check((vl1.get_long() == 0L), "vl1 should be 0");

    t.check(vl1.is_zero(), "vl1 should be 0");

    const VeryLong vl2("1234567890123456789");
    t.check((vl2.get_long_long() == 1234567890123456789LL), "vl2 should be 1234567890123456789");

    const VeryLong vl3(vl2);
    t.check((vl3.get_long_long() == 1234567890123456789LL), "vl3 should be 1234567890123456789");

    t.check(vl2 == vl3, "vl2 should be equal to vl3");

    std::string s("123451234512345");
    const VeryLong vl4(s);
    t.check((vl4.get_long_long() == 123451234512345LL), "vl4 should be 123451234512345");

    const VeryLong v15(-234567L);
    t.check((v15.get_long() == -234567L), "vl5 should be -234567L");

    const VeryLong v16(234567UL);
    t.check((v16.get_long() == 234567UL), "vl6 should be 234567UL");

    const VeryLong v17(234567123456ULL);
    t.check((v17.get_long_long() == 234567123456ULL), "vl7 should be 234567123456ULL");

    const VeryLong v17a = 2230183161LL;
    t.check((v17a.get_long_long() == 2230183161LL), "v17a should be 2230183161LL");

    const VeryLong v18(-234567123456LL);
    t.check((v18.get_long_long() == -234567123456LL), "vl8 should be -234567123456LL");

    const VeryLong v19(-1234567.23456);
    t.check(UnitTest::compare_double(v19.get_double(), -1234567), "v19 should be -1234567");

    const VeryLong v20(-1234567798798719287192371928392.0);
    t.check(UnitTest::compare_double(v20.get_long_double(), -1234567798798719287192371928392.0), "v20 should be -1234567798798719287192371928392");

    const VeryLong v21(1L);
    t.check(v21.is_one(), "v21 should be 1L");

    t.check(v20 != v21, "v20 should not be equal to v21");
    t.check(v20 < v21, "v20 should be less than v21");
    t.check(v20 <= v21, "v20 should be less than or equal to v21");
    t.check(v20 < -123456L, "v20 should be less than or equal to -123456L");
    t.check(v21 > v20, "v21 should be greater than v20");
    t.check(v21 >= v20, "v21 should be greater than or equal to v20");

    VeryLong v22;
    v22 = v20;
    t.check(UnitTest::compare_double(v22.get_long_double(), -1234567798798719287192371928392.0), "v22 should be -1234567798798719287192371928392");

    const VeryLong v23("1234567890123456789012345678901234567890");
    const VeryLong v24("9876543210987654321098765432109876543210");
    VeryLong v25(v23);
    v25 += v24;
    const VeryLong v26("11111111101111111110111111111011111111100");
    t.check(v25 == v26, "v25 should equal v26");

    v25 -= v24;
    t.check(v25 == v23, "v25 should equal v23");

    const VeryLong v27("1232355834509287492857913847");

    v25 *= v27;
    const VeryLong v28("1521426942491462939710131061863585936186206450924386713067892572830");
    t.check(v25 == v28, "v25 should equal v28");

    v25 /= v27;
    t.check(v25 == v23, "v25 should equal v23");

    v25 /= 12345L;
    const VeryLong v29("100005499402467135602458135188435363");
    t.check(v25 == v29, "v25 should equal v29");

    VeryLong v30 = -v25;
    const VeryLong v31("-100005499402467135602458135188435363");
    t.check(v30 == v31, "v30 should equal v31");

    v30.negate();
    t.check(v30 == v29, "v30 should equal v29");

    const VeryLong v32("123123");
    const VeryLong v33("23");
    VeryLong v34 = exp(v32, v33);
    const VeryLong v35("1196193243742734313072657062937753550797883517849819344399014727405544587087316341405542835997613669890217861709298267");
    t.check(v34 == v35, "v34 should equal v35");

    double d36 = ln(v35);
    t.check(UnitTest::compare_double(d36, 269.581600097816), "ln(v35) should equal 269.581600097816");

    double d37 = log10(v35);
    t.check(UnitTest::compare_double(d37, 117.07780134513), "log10(v35) should equal 117.07780134513");

    VeryLong v38 = v23 + v24;
    t.check(v38 == v26, "v38 should equal v26");

    VeryLong v39 = v38 - v24;
    t.check(v39 == v23, "v39 should equal v23");

    VeryLong v40 = v39 * v27;
    t.check(v40 == v28, "v40 should equal v28");

    VeryLong v41 = v40 / v27;
    t.check(v41 == v23, "v41 should equal v23");

    VeryLong v42 = v41 / 12345L;
    t.check(v42 == v29, "v42 should equal v29");

    v42 %= v32;
    const VeryLong v43("53367");
    t.check(v42 == v43, "v42 should equal v43");

    VeryLong v44 = v29 % v32;
    t.check(v44 == v43, "v44 should equal v43");

    long int l45 = v29 % 123123L;
    t.check(l45 == 53367, "l45 should equal 53367");

    int i46 = jacobi_symbol(v23, v32);
    t.check(i46 == 0, "jacobi_symbol(v23, v32) should equal 0");

    int i47 = jacobi_symbol(VeryLong("5"), VeryLong("11"));
    t.check(i47 == 1, "jacobi_symbol(5, 11) should equal 1");

    int i48 = jacobi_symbol(VeryLong("2"), VeryLong("11"));
    t.check(i48 == -1, "jacobi_symbol(2, 11) should equal -1");

    const VeryLong v49("12345");
    VeryLong v50(v49);
    v50 *= v50;
    v50 *= v50;
    v50 *= v50;
    v50 *= v49;
    VeryLong v51 = v50.nth_root(9);
    t.check(v51 == v49, "v51 should equal v49");

    v50 += v49;
    VeryLong v52 = v50.nth_root(9);
    t.check(v52 == v49, "v52 should equal v49");

    VeryLong v53 = gcd(v27, v28);
    t.check(v53 == v27, "gcd(v27, v28) should equal v28");

    VeryLong v54;
    VeryLong v55;
    VeryLong v56 = extended_gcd(v23, v24, v54, v55);
    t.check(v54 == VeryLong(-8L), "v54 should be -8");
    t.check(v55 == VeryLong(1L), "v55 should be 1");
    t.check(v56 == VeryLong("90000000009000000000900000000090"), "v56 should equal 90000000009000000000900000000090");

    t.check(v23 * v54 + v24 * v55 == v56, "Should have v23 * v54 + v24 * v55 == v56");

    bool b57 = v23.is_probable_prime();
    t.check(!b57, "v23 is not a prime");

    bool b58 = v27.is_probable_prime();
    t.check(!b58, "v27 is not a prime");

    const VeryLong v59("1232355834509287492857913897");
    bool b60 = v59.is_probable_prime();
    t.check(b60, "v59 is a prime");

    VeryLong v61(v59);
    v61 *= v61;
    v61 *= v61;
    v61 *= v61;
    v61 *= v59;
    bool b62 = v61.is_prime_power();
    t.check(b62, "v61 is a prime power");

    v61 -= 1L;
    bool b63 = v61.is_prime_power();
    t.check(!b63, "v61 is not a prime power");

    v61 += 1L;
    bool b64 = v61.is_square();
    t.check(!b64, "v61 is not a square");

    v61 *= v59;
    bool b65 = v61.is_square();
    t.check(b65, "v61 is a square");

    VeryLong v62(v27);
    VeryLong v63;
    VeryLong v64;
    bool b66 = v62.factorise_pollardrho(&v63, &v64);
    t.check(b66, "factorise_pollardrho should succeed");
    t.check(v63 * v64 == v62, "should have v63 * v64 == v62");

    VeryLong v67;
    VeryLong v68;
    VeryLong v69(0x7ffffffffffffffdLL);
    bool b70 = v69.factorise_qs(&v67, &v68);
    t.check(b70, "factorise_qs should succeed");
    t.check(v67 * v68 == v69, "should have v67 * v68 == v69");

    VeryLong v71;
    VeryLong v72;
    VeryLong v73(12347LL * 34583LL);
    bool b74 = v73.factorise_squfof(&v71, &v72);
    t.check(b74, "factorise_squfof should succeed");
    t.check(v71 * v72 == v73, "should have v71 * v72 == v73");

    VeryLong v75;
    VeryLong v76;
    bool b77 = v73.factorise_p_minus_1(&v75, &v76);
    t.check(b77, "factorise_p_minus_1 should succeed");
    t.check(v75 * v76 == v73, "should have v75 * v76 == v73");

    VeryLong v77;
    std::vector<VeryLong> factors;
    t.check(v23.factorise_trial_division(&factors, &v77), "factorise_trial_division should succeed");
    t.check(factors.size() == 12, "factors.size() should be 12");
    VeryLong v78(v77);
    for (auto& f: factors)
    {
        v78 *= f;
    }
    t.check(v78 == v23, "v78 should equal v23");

    VeryLong v79;
    std::vector<long int> ifactors;
    t.check(v23.factorise_trial_division(&ifactors, &v79), "factorise_trial_division should succeed");
    t.check(ifactors.size() == 12, "ifactors.size() should be 12");

    VeryLong v80(v79);
    for (auto& f: ifactors)
    {
        v80 *= f;
    }
    t.check(v80 == v23, "v80 should equal v23");

    VeryLong v81;
    std::vector<long int> iifactors;
    t.check(v23.factorise_trial_division(&iifactors, &v81, 10000), "factorise_trial_division should succeed");
    t.check(iifactors.size() == 10, "iifactors.size() should be 10");

    VeryLong v82(v81);
    for (auto& iif: iifactors)
    {
        v82 *= iif;
    }
    t.check(v82 == v23, "v82 should equal v23");

    VeryLong v83;
    VeryLong v84;
    t.check(v23.factorise_trial_division(&v83, &v84), "factorise_trial_division should succeed");
    t.check(v83 * v84 == v23, "Should have v83 * v84 == v23");

    VeryLong v85;
    VeryLong v86;
    t.check(v73.factorise_ecm(&v85, &v86), "factorise_ecm should succeed");
    t.check(v85 * v86 == v73, "Should have v85 * v86 == v73");

    VeryLong v87;
    VeryLong v88;
    t.check(v73.factorise_ecm(&v87, &v88), "factorise_fermat should succeed");
    t.check(v87 * v88 == v73, "Should have v87 * v88 == v73");

    std::vector<VeryLong> factors89;
    t.check(v23.factorise(&factors89), "factorise should succeed");
    t.check(factors89.size() == 13, "factors89.size() should be 13");
    VeryLong v90(1L);
    for (auto& f: factors89)
    {
        v90 *= f;
    }
    t.check(v90 == v23, "v90 should equal v23");

    std::vector<VeryLong> factors91;
    t.check(v73.factorise_no_trial(&factors91), "factorise_no_trial should succeed");
    t.check(factors91.size() == 2, "factors91.size() should be 2");
    if (factors91.size() == 2)
    {
        t.check(factors91[0] * factors91[1] == v73, "product should be v73");
    }

    long long int lli92 = v73.get_long_long();
    std::vector<long int> ifactors93;
    t.check(VeryLong::factorise_no_trial_ll(lli92, &ifactors93), "factorise_no_trial_ll should succeed");
    t.check(ifactors93.size() == 2, "ifactors93.size() should be 2");
    if (ifactors93.size() == 2)
    {
        t.check(ifactors93[0] * (long long int)ifactors93[1] == lli92, "product should be lli92");
    }

    long int li94 = 1237L * 2351L;
    std::vector<long int> ifactors95;
    t.check(VeryLong::factorise_no_trial_l(li94, &ifactors95), "factorise_no_trial_l should succeed");
    t.check(ifactors95.size() == 2, "ifactors95.size() should be 2");
    if (ifactors95.size() == 2)
    {
        t.check(ifactors95[0] * ifactors95[1] == li94, "product should be li94");
    }

    long long int lli96(0x7FFFFFFFFFFFFFFFLL);
    const VeryLong v97(lli96);
    t.check(lli96 == v97.get_long_long(), "v97.get_long_long() should be lli96");

    long long int lli96m(0x8000000000000000LL);
    const VeryLong v97m(lli96m);
    t.check(lli96m == v97m.get_long_long(), "v97m.get_long_long() should be lli96m");

    long int li98(0x7FFFFFFFL);
    const VeryLong v99(li98);
    t.check(li98 == v99.get_long(), "v99.get_long() should be li98");

    long int li98m(0x80000000L);
    const VeryLong v99m(li98m);
    t.check(li98m == v99m.get_long(), "v99m.get_long() should be li98m");

    const VeryLong v100("123346620910319027308203983058475023948203498720394");
    const VeryLong v101("504930582034509586048345983759385739857394520570157");

    VeryLong v102 = VeryLong::random(v100, v101);
    t.check(v102 >= v100 && v102 <= v101, "v102 should be in range v100 - v101");

    const VeryLong v103(-v102);
    VeryLong v104 = abs(v103);
    t.check(v104 == v102, "v104 should equal v102");

    VeryLong v105 = v103.abs();
    t.check(v105 == v102, "v105 should equal v102");

    VeryLong::clear_prime_table();
    VeryLong::generate_prime_table();
    long int p = VeryLong::firstPrime();
    t.check(p == 2L, "firstPrime should be 2");

    t.check(VeryLong::nextPrime() == 3L, "nextPrime should be 3");
    t.check(VeryLong::nextPrime() == 5L, "nextPrime should be 5");
    t.check(VeryLong::nextPrime() == 7L, "nextPrime should be 7");
    t.check(VeryLong::nextPrime() == 11L, "nextPrime should be 11");
    t.check(VeryLong::nextPrime() == 13L, "nextPrime should be 13");
    t.check(VeryLong::nextPrime() == 17L, "nextPrime should be 17");

    VeryLong::set_max_prime(97L);
    t.check(VeryLong::get_max_prime() == 97L, "get_max_prime() should return 97");

    // FastVeryLong tests
    FastVeryLong fvl1;
    t.check(fvl1.get_long() == 0L, "fvl1 should be 0");

    const VeryLong v106(0x7FFFFFFFL);
    const FastVeryLong fvl2(v106);
    t.check(fvl2.is_faster(), "fvl2 should be faster");
    t.check(fvl2.is_fast(), "fvl2 should be fast");

    const FastVeryLong fvl3(fvl2);
    t.check(fvl2 == fvl3, "fvl2 should be equal to fvl3");

    t.check(!(fvl2 != v106), "fvl2 should be equal to v106");

    const VeryLong v2000("24304511608315849");
    std::vector<VeryLong> v2000_factors;
    VeryLong v2000_unfactored;
    t.check(v2000.factorise_no_trial(&v2000_factors, v2000_unfactored), "v2000 should factor");

    VeryLong::setDebug();
    VeryLong v2001("68454066430929935312919949146284654064784560463655945766546376196659857165127911441901256275347366429914574181910564173262954907243781843526369433234400129532781555616626107827873263556651889521214284201244347354065119176141583589197364828686588669194047973922985132474269362219900353120931603959821770405789500417740915047315436720600435558247495565875952517451686427355414486482150319474507487");
    VeryLong v2001_factor;
    VeryLong v2001_unfactored;
    t.check(v2001.factorise_ecm(&v2001_factor, &v2001_unfactored), "v2001 should factor with ECM");

    const VeryLong v200("5544779380905324760346515880849056979247549397556131607090256471929448430375360826794001758303136680823080508734755698034299347486746329325635924091986410492155306004946714734057734348088803051218357020300792135679274653267468270724986551123613682204717885887761795730415818339811928602795459920745563402868949533837014118832550374368635280218047140835952153913586600615788573405054175877435106447");
    std::vector<VeryLong> factors200;

    VeryLong v200_unfactored;
    t.check(v200.factorise(&factors200, v200_unfactored), "v200 should factor");
    VeryLong v200_product(1L);

    for (size_t i = 0; i < factors200.size(); ++i)
    {
        v200_product *= factors200[i];
        //std::cout << factors200[i] << std::endl;
    }
    //std::cout << "Unfactored : " << v200_unfactored << std::endl;
    t.check(v200_unfactored == 1L, "Should be no unfactored part in v200");
    t.check(v200 == v200_product, "v200 should equal v200_product");

    long double alpha = -13501.183274494388;
    const VeryLong a("-434022535");
    const VeryLong b("32147");
    long double delta = 1.0;

    long double aa = a.get_long_double();
    long double bb = b.get_long_double();
    long double mod2 = aa * aa + bb * bb * alpha * alpha - 2.0 * aa * bb * alpha;
    long double answer = (long double)log((double)mod2) / 2.0 + (long double)log(delta);

    std::cout << "mod2 = " << mod2 << ", answer = " << answer << std::endl;

    std::cout << "a - b alpha = " << aa - bb*alpha << std::endl;
    t.test_summary();
#endif
    return 0;
}
