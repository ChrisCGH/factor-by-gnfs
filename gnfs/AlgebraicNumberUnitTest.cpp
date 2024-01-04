#include <iostream>
#include <iomanip>
#include <sstream>
#include "NumberField.h"
#include "Polynomial.inl"
#include "AlgebraicNumber.h"
#include "AlgebraicNumber_in_O_pO.h"
#include "Ideal.h"
#include <limits>
#include <complex>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/ui/text/TestRunner.h>

template<> Matrix<Quotient<VeryLong> > AlgebraicNumber_in_O_pO_<VeryLong, VeryLong, VeryLongModular>::W_mult_;
template<> Matrix<Quotient<VeryLong> > AlgebraicNumber_in_O_pO_<long, VeryLong, LongModular>::W_mult_;
template<> Matrix<VeryLongModular> AlgebraicNumber_in_O_pO_<VeryLong, VeryLong, VeryLongModular>::M_;
template<> Matrix<LongModular> AlgebraicNumber_in_O_pO_<long, VeryLong, LongModular>::M_;
template<> VeryLong AlgebraicNumber_in_O_pO_<VeryLong, VeryLong, VeryLongModular>::p_;
template<> long AlgebraicNumber_in_O_pO_<long, VeryLong, LongModular>::p_;
template<> VeryLongModular AlgebraicNumber_in_O_pO_<VeryLong, VeryLong, VeryLongModular>::w01_;
template<> LongModular AlgebraicNumber_in_O_pO_<long, VeryLong, LongModular>::w01_;
template<> VeryLongModular AlgebraicNumber_in_O_pO_<VeryLong, VeryLong, VeryLongModular>::w11_;
template<> LongModular AlgebraicNumber_in_O_pO_<long, VeryLong, LongModular>::w11_;
template<> bool AlgebraicNumber_in_O_pO_<VeryLong, VeryLong, VeryLongModular>::optimisation_ok_;
template<> bool AlgebraicNumber_in_O_pO_<long, VeryLong, LongModular>::optimisation_ok_;
class AlgebraicNumberTest : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE( AlgebraicNumberTest );
    CPPUNIT_TEST( testExceptions );
    CPPUNIT_TEST( testNumberField );
    CPPUNIT_TEST( testNumberFieldThrow );
    CPPUNIT_TEST( testAlgebraicNumber );
    CPPUNIT_TEST( testAlgebraicNumber_in_O_pO );
    CPPUNIT_TEST( testAlgebraicNumber_in_O_pO_1 );
    CPPUNIT_TEST( testIdeal );
    CPPUNIT_TEST( testPrimeIdeal );
    CPPUNIT_TEST_SUITE_END();
private:
    Polynomial<VeryLong> f;
public:
    void setUp()
    {
        // f1 = 9401401242022932575419604204700
        //      + 623872110646578368801362410 X
        //      - 8591401659640532521423 X^2
        //      - 222405543007291322 X^3
        //      + 1249913278668 X^4
        //      + 2691780 X^5
        //
        std::vector<VeryLong> c;
        c.resize(6);
        c[0] = "9401401242022932575419604204700";
        c[1] = "623872110646578368801362410";
        c[2] = "-8591401659640532521423";
        c[3] = "-222405543007291322";
        c[4] = "1249913278668";
        c[5] = "2691780";
        f = Polynomial<VeryLong>(c);
    }

    void tearDown()
    {
    }

    void testExceptions()
    {
        CPPUNIT_ASSERT_THROW_MESSAGE("", const VeryLong c_d = AlgebraicNumber::c_d(), std::string);
        CPPUNIT_ASSERT_THROW_MESSAGE("", const VeryLong index = AlgebraicNumber::index(), std::string);
        CPPUNIT_ASSERT_THROW_MESSAGE("", AlgebraicNumber::degree(), std::string);
        CPPUNIT_ASSERT_THROW_MESSAGE("", AlgebraicNumber::nf(), std::string);
        CPPUNIT_ASSERT_THROW_MESSAGE("", AlgebraicNumber().make_ibc(), std::string);
        CPPUNIT_ASSERT_THROW_MESSAGE("", AlgebraicNumber(), std::string);
    }

    void testNumberFieldThrow()
    {
        Polynomial<VeryLong> f11 = Polynomial<VeryLong>::read_polynomial("X^11 - 1");
        CPPUNIT_ASSERT_THROW_MESSAGE("", NumberField nf11(f11, "nf11.fb.dat"), std::string);
    }

    void testNumberField()
    {
        NumberField nf(f, "nf.fb.dat");
        CPPUNIT_ASSERT(nf.conjugates() == 5);
        // nf.conjugate(0) = (-13501.1832744944,0)
        // nf.conjugate(1) = (49973.9880026806,0)
        // nf.conjugate(2) = (150553.590506638,0)
        // nf.conjugate(3) = (-57939.61757667,0)
        // nf.conjugate(4) = (-593431.292688357,0)
        CPPUNIT_ASSERT_DOUBLES_EQUAL(nf.conjugate(0).real(), -13501.1832744944, std::numeric_limits<float>::epsilon());
        CPPUNIT_ASSERT_DOUBLES_EQUAL(nf.conjugate(1).real(), 49973.9880026806, std::numeric_limits<float>::epsilon());
        CPPUNIT_ASSERT_DOUBLES_EQUAL(nf.conjugate(2).real(), 150553.590506638, std::numeric_limits<float>::epsilon());
        CPPUNIT_ASSERT_DOUBLES_EQUAL(nf.conjugate(3).real(), -57939.61757667, std::numeric_limits<float>::epsilon());
        CPPUNIT_ASSERT_DOUBLES_EQUAL(nf.conjugate(4).real(), -593431.292688357, std::numeric_limits<float>::epsilon());
        CPPUNIT_ASSERT_DOUBLES_EQUAL(nf.conjugate(0).imag(), 0.0, std::numeric_limits<float>::epsilon());
        CPPUNIT_ASSERT_DOUBLES_EQUAL(nf.conjugate(1).imag(), 0.0, std::numeric_limits<float>::epsilon());
        CPPUNIT_ASSERT_DOUBLES_EQUAL(nf.conjugate(2).imag(), 0.0, std::numeric_limits<float>::epsilon());
        CPPUNIT_ASSERT_DOUBLES_EQUAL(nf.conjugate(3).imag(), 0.0, std::numeric_limits<float>::epsilon());
        CPPUNIT_ASSERT_DOUBLES_EQUAL(nf.conjugate(4).imag(), 0.0, std::numeric_limits<float>::epsilon());

        CPPUNIT_ASSERT_THROW_MESSAGE("", nf.conjugate(-1), std::string);
        CPPUNIT_ASSERT_THROW_MESSAGE("", nf.conjugate(6), std::string);

        CPPUNIT_ASSERT(nf.degree() == 5);
        CPPUNIT_ASSERT(nf.c_d() == VeryLong("2691780"));
        const VeryLong index("798834230984453137044709875604358400000000");
        CPPUNIT_ASSERT(nf.index() == index);
        const Polynomial<VeryLong> mp = Polynomial<VeryLong>::read_polynomial("9401401242022932575419604204700 + 623872110646578368801362410 X - 8591401659640532521423 X^2 - 222405543007291322 X^3 + 1249913278668 X^4 + 2691780 X^5");
        CPPUNIT_ASSERT(nf.min_poly() == mp);

        const Polynomial<VeryLong> mmp = Polynomial<VeryLong>::read_polynomial("493572365661951165444447584389504458344834234261232000000 + 12167861492325100495180675345143661535602320000 X - 62250543469175257378970366393833200 X^2 - 598666792556166634733160 X^3 + 1249913278668 X^4 + 1 X^5");

        CPPUNIT_ASSERT(nf.monic_min_poly() == mmp);

        const VeryLong disc("66792821331670213274787635927696441943633124232181335670826527788207293595077180699627767499828407596673634024356167317150449784300992672779002195101757800000");
        CPPUNIT_ASSERT(nf.discriminant() == disc);
        const VeryLong field_disc("15145764474301635663217151004012798626674177830426606728078577729752220769858771133702441609938414421014429483980990321349308341111336206979365577120580");
        CPPUNIT_ASSERT(nf.fieldDiscriminant() == field_disc);

        Matrix<Quotient<VeryLong> > checkIntegralBasisAlpha(5,5);
        // 1       0       0       0           1/7
        // 0 1345890   40458 4638574/5  32730265/42
        // 0       0 2691780 4078128/5 275165447/105
        // 0       0       0  538356     7409138/35
        // 0       0       0       0       12818
        checkIntegralBasisAlpha(0,0) = Quotient<VeryLong>(1L);

        checkIntegralBasisAlpha(1,1) = Quotient<VeryLong>(1345890L);

        checkIntegralBasisAlpha(1,2) = Quotient<VeryLong>(40458L);
        checkIntegralBasisAlpha(2,2) = Quotient<VeryLong>(2691780L);

        checkIntegralBasisAlpha(1,3) = Quotient<VeryLong>(4638574L, 5L);
        checkIntegralBasisAlpha(2,3) = Quotient<VeryLong>(4078128L, 5L);
        checkIntegralBasisAlpha(3,3) = Quotient<VeryLong>(538356L);

        checkIntegralBasisAlpha(0,4) = Quotient<VeryLong>(1L, 7L);
        checkIntegralBasisAlpha(1,4) = Quotient<VeryLong>(32730265L, 42L);
        checkIntegralBasisAlpha(2,4) = Quotient<VeryLong>(275165447L, 105L);
        checkIntegralBasisAlpha(3,4) = Quotient<VeryLong>(7409138L, 35L);
        checkIntegralBasisAlpha(4,4) = Quotient<VeryLong>(12818L);

        CPPUNIT_ASSERT(nf.w() == checkIntegralBasisAlpha);

        Matrix<Quotient<VeryLong> > checkIntegralBasisAlphaInv(5,5);
        // 1 0             0                          0                                     -1/89726
        // 0 1/1345890 -6743/603806630700 -513376159129/406328653096411500 -501825203419378151/22786402954830386405625
        // 0 0             1/2691780             -84961/150951657675           -54169883960389/812657306192823000
        // 0 0             0                          1/538356                        -3704569/120761326140
        // 0 0             0                          0                                      1/12818

        checkIntegralBasisAlphaInv(0,0) = Quotient<VeryLong>(1L);

        checkIntegralBasisAlphaInv(1,1) = Quotient<VeryLong>(1L, 1345890L);

        checkIntegralBasisAlphaInv(1,2) = Quotient<VeryLong>(-6743L, "603806630700");
        checkIntegralBasisAlphaInv(2,2) = Quotient<VeryLong>(1L, 2691780L);

        checkIntegralBasisAlphaInv(1,3) = Quotient<VeryLong>("-513376159129", "406328653096411500");
        checkIntegralBasisAlphaInv(2,3) = Quotient<VeryLong>(-84961L, "150951657675");
        checkIntegralBasisAlphaInv(3,3) = Quotient<VeryLong>(1L, 538356L);

        checkIntegralBasisAlphaInv(0,4) = Quotient<VeryLong>(-1L, 89726L);
        checkIntegralBasisAlphaInv(1,4) = Quotient<VeryLong>("-501825203419378151", "22786402954830386405625");
        checkIntegralBasisAlphaInv(2,4) = Quotient<VeryLong>("-54169883960389", "812657306192823000");
        checkIntegralBasisAlphaInv(3,4) = Quotient<VeryLong>(-3704569L, "120761326140");
        checkIntegralBasisAlphaInv(4,4) = Quotient<VeryLong>(1L, 12818L);

        CPPUNIT_ASSERT(nf.winv() == checkIntegralBasisAlphaInv);

        Matrix<Quotient<VeryLong> > checkIntegralBasisTheta(5,5);
        // 1 0      0              0                        1/7
        // 0 1/2 6743/448630 2319287/6729450          6546053/22610952
        // 0 0      1/2691780  84961/754758288375   275165447/760796354682000
        // 0 0      0              1/36228397842000   3704569/341316068600985660000
        // 0 0      0              0                        1/4095792823211827920000

        checkIntegralBasisTheta(0,0) = Quotient<VeryLong>(1L);

        checkIntegralBasisTheta(1,1) = Quotient<VeryLong>(1L, 2L);

        checkIntegralBasisTheta(1,2) = Quotient<VeryLong>(6743L, 448630L);
        checkIntegralBasisTheta(2,2) = Quotient<VeryLong>(1L, 2691780L);

        checkIntegralBasisTheta(1,3) = Quotient<VeryLong>(2319287L, 6729450L);
        checkIntegralBasisTheta(2,3) = Quotient<VeryLong>(84961L, "754758288375");
        checkIntegralBasisTheta(3,3) = Quotient<VeryLong>(1L, "36228397842000");

        checkIntegralBasisTheta(0,4) = Quotient<VeryLong>(1L, 7L);
        checkIntegralBasisTheta(1,4) = Quotient<VeryLong>(6546053L, 22610952L);
        checkIntegralBasisTheta(2,4) = Quotient<VeryLong>(275165447L, "760796354682000");
        checkIntegralBasisTheta(3,4) = Quotient<VeryLong>(3704569L, "341316068600985660000");
        checkIntegralBasisTheta(4,4) = Quotient<VeryLong>(1L, "4095792823211827920000");

        CPPUNIT_ASSERT(nf.W() == checkIntegralBasisTheta);

        Matrix<Quotient<VeryLong> > checkIntegralBasisThetaInv(5,5);
        // 1 0       0               0  -585113260458832560000
        // 0 2  -80916 -24642055638192 -1156205268678247259904
        // 0 0 2691780 -10977423387840 -3499521845925501658080
        // 0 0       0  36228397842000 -1610527194781681176000
        // 0 0       0               0  4095792823211827920000

        checkIntegralBasisThetaInv(0,0) = Quotient<VeryLong>(1L);

        checkIntegralBasisThetaInv(1,1) = Quotient<VeryLong>(2L);

        checkIntegralBasisThetaInv(1,2) = Quotient<VeryLong>(-80916L);
        checkIntegralBasisThetaInv(2,2) = Quotient<VeryLong>(2691780L);

        checkIntegralBasisThetaInv(1,3) = Quotient<VeryLong>("-24642055638192");
        checkIntegralBasisThetaInv(2,3) = Quotient<VeryLong>("-10977423387840");
        checkIntegralBasisThetaInv(3,3) = Quotient<VeryLong>("36228397842000");

        checkIntegralBasisThetaInv(0,4) = Quotient<VeryLong>("-585113260458832560000");
        checkIntegralBasisThetaInv(1,4) = Quotient<VeryLong>("-1156205268678247259904");
        checkIntegralBasisThetaInv(2,4) = Quotient<VeryLong>("-3499521845925501658080");
        checkIntegralBasisThetaInv(3,4) = Quotient<VeryLong>("-1610527194781681176000");
        checkIntegralBasisThetaInv(4,4) = Quotient<VeryLong>("4095792823211827920000");

        CPPUNIT_ASSERT(nf.Winv() == checkIntegralBasisThetaInv);

        std::vector<std::pair<Polynomial<VeryLongModular>, int> > factors;
        VeryLong p(11L);
        nf.factorise_monic_min_poly_over_p(p, factors);
        CPPUNIT_ASSERT(factors.size() == 4);

        // (1 + 1 X, 1)
        Polynomial<VeryLong> p0 = Polynomial<VeryLong>::read_polynomial("1 + 1 X");
        Polynomial<VeryLongModular> pm0 = convert_to_F_p<VeryLong, VeryLong, VeryLongModular>(p0, p);
        // (4 + 5 X, 1)
        Polynomial<VeryLong> p1 = Polynomial<VeryLong>::read_polynomial("4 + 5 X");
        Polynomial<VeryLongModular> pm1 = convert_to_F_p<VeryLong, VeryLong, VeryLongModular>(p1, p);
        // (4 + 9 X, 1)
        Polynomial<VeryLong> p2 = Polynomial<VeryLong>::read_polynomial("4 + 9 X");
        Polynomial<VeryLongModular> pm2 = convert_to_F_p<VeryLong, VeryLong, VeryLongModular>(p2, p);
        // (9 + 8 X + 1 X^2, 1)
        Polynomial<VeryLong> p3 = Polynomial<VeryLong>::read_polynomial("9 + 8 X + 1 X^2");
        Polynomial<VeryLongModular> pm3 = convert_to_F_p<VeryLong, VeryLong, VeryLongModular>(p3, p);

        CPPUNIT_ASSERT(factors[0].first == pm0);
        CPPUNIT_ASSERT(factors[0].second == 1);
        CPPUNIT_ASSERT(factors[1].first == pm1);
        CPPUNIT_ASSERT(factors[1].second == 1);
        CPPUNIT_ASSERT(factors[2].first == pm2);
        CPPUNIT_ASSERT(factors[2].second == 1);
        CPPUNIT_ASSERT(factors[3].first == pm3);
        CPPUNIT_ASSERT(factors[3].second == 1);

        // gcc support long double, which leads to a slightly different result for the ideal bound,
        // but this is not too important
#ifdef USING_GCC
        CPPUNIT_ASSERT(nf.idealBound() == VeryLong("37740227251498267579591702209487087561625622897229031140623634963687187973188114373313552295994367747328291494511444396891172432802474305883555275067000760486701338287134007615985800132410941340233778924675727101257944007906118073104515418265035743306541510911214588725646956953600000"));
#else
        CPPUNIT_ASSERT(nf.idealBound() == VeryLong("37740227251498528623991403793957042009614011788307957572937127172913054839465973958706327060082047725629179516419458695364014374281607269631063265784193060887943255560142212230236643239239982592477863724308010398517597149967524376244250767628584897477988084305625053774479360000000000"));
#endif
        const Matrix<Quotient<VeryLong> >& sm = nf.structureMatrix();
        // 1 0 0 0 0
        // 0 1 0 0 0
        // 0 0 1 0 0
        // 0 0 0 1 0
        // 0 0 0 0 1
        // -156690020700382209590326736745/44863 -1223278648326624252551691/5278 8591401659640532521423/2691780 15886110214806523/192270 -104159439889/224315
        // 3264148958469525287829727410372071004261/2012688769 123280969191591102760858690997101749/1183934570 -1034819457229270111310786182241197/603806630700 -1011355072598385651906063953/28752696700 90039583197781502441141/301903315350
        // -940553610340944653585779867927434154671401197628403/902952562436470 -35853090966551487099091348618272519800238028577/531148566139100 858186762567649786586236481160677803589484743/812657306192823000 332724093463615655169759468744438222287/14511737610586125 -23521049482747918031001748948107793/135442884365470500
        // 245700915356432369323094884625744872343268708758596452279596919/405091608085873536100 9342720441724197272689552186352838088786465418254655257621/238289181226984433000 -226688485029063262277660101367617942280993513378649999739/364582447277286182490000 -14423224568515460308347916819655111237547808539346/1085066807372875543125 9439598159881575411373406390522350411132500271/91145611819321545622500
        Matrix<Quotient<VeryLong> > checksm(9, 5);
        checksm(0, 0) = Quotient<VeryLong>(1L);
        checksm(1, 1) = Quotient<VeryLong>(1L);
        checksm(2, 2) = Quotient<VeryLong>(1L);
        checksm(3, 3) = Quotient<VeryLong>(1L);
        checksm(4, 4) = Quotient<VeryLong>(1L);
        checksm(5, 0) = Quotient<VeryLong>("-156690020700382209590326736745", "44863");
        checksm(5, 1) = Quotient<VeryLong>("-1223278648326624252551691", "5278");
        checksm(5, 2) = Quotient<VeryLong>("8591401659640532521423", "2691780");
        checksm(5, 3) = Quotient<VeryLong>("15886110214806523", "192270");
        checksm(5, 4) = Quotient<VeryLong>("-104159439889", "224315");
        checksm(6, 0) = Quotient<VeryLong>("3264148958469525287829727410372071004261", "2012688769");
        checksm(6, 1) = Quotient<VeryLong>("123280969191591102760858690997101749", "1183934570");
        checksm(6, 2) = Quotient<VeryLong>("-1034819457229270111310786182241197", "603806630700");
        checksm(6, 3) = Quotient<VeryLong>("-1011355072598385651906063953", "28752696700");
        checksm(6, 4) = Quotient<VeryLong>("90039583197781502441141", "301903315350");
        checksm(7, 0) = Quotient<VeryLong>("-940553610340944653585779867927434154671401197628403", "902952562436470");
        checksm(7, 1) = Quotient<VeryLong>("-35853090966551487099091348618272519800238028577", "531148566139100");
        checksm(7, 2) = Quotient<VeryLong>("858186762567649786586236481160677803589484743", "812657306192823000");
        checksm(7, 3) = Quotient<VeryLong>("332724093463615655169759468744438222287", "14511737610586125");
        checksm(7, 4) = Quotient<VeryLong>("-23521049482747918031001748948107793", "135442884365470500");
        checksm(8, 0) = Quotient<VeryLong>("245700915356432369323094884625744872343268708758596452279596919", "405091608085873536100");
        checksm(8, 1) = Quotient<VeryLong>("9342720441724197272689552186352838088786465418254655257621", "238289181226984433000");
        checksm(8, 2) = Quotient<VeryLong>("-226688485029063262277660101367617942280993513378649999739", "364582447277286182490000");
        checksm(8, 3) = Quotient<VeryLong>("-14423224568515460308347916819655111237547808539346", "1085066807372875543125");
        checksm(8, 4) = Quotient<VeryLong>("9439598159881575411373406390522350411132500271", "91145611819321545622500");
        CPPUNIT_ASSERT(sm == checksm);

    }

    void testAlgebraicNumber()
    {
        NumberField nf(f, "nf.fb.dat");
        AlgebraicNumber::setNumberField(nf);
        CPPUNIT_ASSERT(AlgebraicNumber::index() == nf.index());
        CPPUNIT_ASSERT(AlgebraicNumber::c_d() == nf.c_d());
        CPPUNIT_ASSERT(AlgebraicNumber::degree() == 5);

        const AlgebraicNumber an1;
        const AlgebraicNumber an2(0L);
        CPPUNIT_ASSERT(an1 == an2);

        const AlgebraicNumber an3("11");
        std::ostringstream oss;
        oss << an3;
        CPPUNIT_ASSERT(oss.str() == "11");
        CPPUNIT_ASSERT(an3.norm() == VeryLong("161051"));
        CPPUNIT_ASSERT(an3.trace() == VeryLong("55"));

        const AlgebraicNumber an4("230945019420469103498230593480394580149820458");
        oss.str("");
        oss << an4;
        CPPUNIT_ASSERT(oss.str() == "230945019420469103498230593480394580149820458");
        CPPUNIT_ASSERT(an4.norm() == VeryLong("656966165172963101825137420658750722510626461989988189958767745390542918683071809791023805605390785573925651934709383759825823701608822934545941180892066143302661885222651628487564300597697388860498332016998098970302148768"));
        CPPUNIT_ASSERT(an4.trace() == VeryLong("1154725097102345517491152967401972900749102290"));

        // 6924553 49395 3 5 5 7 43 3209 6079 9173 48479 101839 201809 207433 534913 1041253 6000047
        const AlgebraicNumber an5(6924553L, 49395L);
        const Quotient<VeryLong> norm_an5 = an5.norm();
        Quotient<VeryLong> check_norm_an5(3L);
        check_norm_an5 *= (5L);
        check_norm_an5 *= (5L);
        check_norm_an5 *= (7L);
        check_norm_an5 *= (43L);
        check_norm_an5 *= (3209L);
        check_norm_an5 *= (6079L);
        check_norm_an5 *= (9173L);
        check_norm_an5 *= (48479L);
        check_norm_an5 *= (101839L);
        check_norm_an5 *= (201809L);
        check_norm_an5 *= (207433L);
        check_norm_an5 *= (534913L);
        check_norm_an5 *= (1041253L);
        check_norm_an5 *= (6000047L);
        /*
        N(a - b alpha) = prod(a - b alpha_i)
                       = b^d prod (a / b - alpha_i)
                       = b^d f(a / b) / c_d
                       = F(a, b) / c_d
        */
        check_norm_an5 /= AlgebraicNumber::c_d();
        oss.str("");
        oss << "an5.norm() should be " << check_norm_an5;
        CPPUNIT_ASSERT(norm_an5 == check_norm_an5);
        oss.str("");
        oss << an5;
        CPPUNIT_ASSERT(oss.str() == "6924553 - 49395 alpha");

        // -2250353386 14449 2 2 7 11 19 83 431 27983 370919 2453779 3242429 10189651 26925421 137133001
        const AlgebraicNumber an6(VeryLong("-2250353386"), 14449L);
        const Quotient<VeryLong> norm_an6 = an6.norm();
        Quotient<VeryLong> check_norm_an6(2L);
        check_norm_an6 *= (2L);
        check_norm_an6 *= (7L);
        check_norm_an6 *= (11L);
        check_norm_an6 *= (19L);
        check_norm_an6 *= (83L);
        check_norm_an6 *= (431L);
        check_norm_an6 *= (27983L);
        check_norm_an6 *= (370919L);
        check_norm_an6 *= (2453779L);
        check_norm_an6 *= (3242429L);
        check_norm_an6 *= (10189651L);
        check_norm_an6 *= (26925421L);
        check_norm_an6 *= (137133001L);
        check_norm_an6 /= AlgebraicNumber::c_d();
        oss.str("");
        oss << "an6.norm() should be " << check_norm_an6;
        CPPUNIT_ASSERT(norm_an6 == check_norm_an6);
        oss.str("");
        oss << an6;
        CPPUNIT_ASSERT(oss.str() == "-2250353386 - 14449 alpha");

        const std::vector<AlgebraicNumber>& ib = AlgebraicNumber::integralBasis();
        CPPUNIT_ASSERT(ib.size() == 5);

        // Integral basis:
        std::string ib1("1");
        AlgebraicNumber an_ib1;
        read(ib1, an_ib1);
        CPPUNIT_ASSERT(an_ib1 == ib[0]);
        std::string ib2("1345890 alpha");
        AlgebraicNumber an_ib2;
        read(ib2, an_ib2);
        CPPUNIT_ASSERT(an_ib2 == ib[1]);
        std::string ib3("40458 alpha + 2691780 alpha^2");
        AlgebraicNumber an_ib3;
        read(ib3, an_ib3);
        CPPUNIT_ASSERT(an_ib3 == ib[2]);
        std::string ib4("4638574/5 alpha + 4078128/5 alpha^2 + 538356 alpha^3");
        AlgebraicNumber an_ib4;
        read(ib4, an_ib4);
        CPPUNIT_ASSERT(an_ib4 == ib[3]);
        std::string ib5("1/7 + 32730265/42 alpha + 275165447/105 alpha^2 + 7409138/35 alpha^3 + 12818 alpha^4");
        AlgebraicNumber an_ib5;
        read(ib5, an_ib5);
        CPPUNIT_ASSERT(an_ib5 == ib[4]);
        Polynomial<Quotient<VeryLong> > p_ib5 = Polynomial<Quotient<VeryLong> >::read_polynomial("1/7 + 32730265/42 a + 275165447/105 a^2 + 7409138/35 a^3 + 12818 a^4");
        AlgebraicNumber pan_ib5(p_ib5);
        CPPUNIT_ASSERT(pan_ib5 == ib[4]);

        std::string an7_str("2474902038376743239461332322533453486017145088647579917686/7 + 848687138292690641752893485837966789178774334641756986027067185671/210 alpha + 21211646115309304218501624599244752000939259665291200926111702457/105 alpha^2 + 40188048938846254718786675529434797027834368184592286864002025034/35 alpha^3 - 5073943593046678921562417863134788148090475853480933891772146 alpha^4");
        AlgebraicNumber an7;
        read(an7_str, an7);
        oss.str("");
        oss << an7;
        AlgebraicNumber an8;
        read(oss.str(), an8);
        CPPUNIT_ASSERT(an7 == an8);

        std::vector<Quotient<VeryLong> > coeff;
        coeff.resize(5);
        coeff[0] = Quotient<VeryLong>(VeryLong("2474902038376743239461332322533453486017145088647579917686"), 7L);
        coeff[1] = Quotient<VeryLong>(VeryLong("848687138292690641752893485837966789178774334641756986027067185671"), 210L);
        coeff[2] = Quotient<VeryLong>(VeryLong("21211646115309304218501624599244752000939259665291200926111702457"), 105L);
        coeff[3] = Quotient<VeryLong>(VeryLong("40188048938846254718786675529434797027834368184592286864002025034"), 35L);
        coeff[4] = Quotient<VeryLong>(VeryLong("-5073943593046678921562417863134788148090475853480933891772146"));

        AlgebraicNumber an9(coeff);
        CPPUNIT_ASSERT(an7 == an9);

        AlgebraicNumber an10(an9);
        CPPUNIT_ASSERT(an9 == an10);

        Matrix<Quotient<VeryLong> > M(5, 2);
        M(0, 0) = Quotient<VeryLong>(VeryLong(1L));

        M(0, 1) = Quotient<VeryLong>(VeryLong("2474902038376743239461332322533453486017145088647579917686"), 7L);
        M(1, 1) = Quotient<VeryLong>(VeryLong("848687138292690641752893485837966789178774334641756986027067185671"), 210L);
        M(2, 1) = Quotient<VeryLong>(VeryLong("21211646115309304218501624599244752000939259665291200926111702457"), 105L);
        M(3, 1) = Quotient<VeryLong>(VeryLong("40188048938846254718786675529434797027834368184592286864002025034"), 35L);
        M(4, 1) = Quotient<VeryLong>(VeryLong("-5073943593046678921562417863134788148090475853480933891772146"));

        AlgebraicNumber an11(M, 0);
        CPPUNIT_ASSERT(an11 == ib[0]);

        AlgebraicNumber an12(M, 1);
        CPPUNIT_ASSERT(an12 == an10);

        Matrix<VeryLong> M1(5, 1);
        M1(0, 0) = VeryLong("2474902038376743239461332322533453486017145088647579917686") * VeryLong(30L);
        M1(1, 0) = VeryLong("848687138292690641752893485837966789178774334641756986027067185671");
        M1(2, 0) = VeryLong("21211646115309304218501624599244752000939259665291200926111702457") * VeryLong(2L);
        M1(3, 0) = VeryLong("40188048938846254718786675529434797027834368184592286864002025034") * VeryLong(6L);
        M1(4, 0) = VeryLong("-5073943593046678921562417863134788148090475853480933891772146") * VeryLong(210L);

        AlgebraicNumber an13(M1, 210L, 0);
        CPPUNIT_ASSERT(an13 == an10);

        AlgebraicNumber an14 = AlgebraicNumber::alpha();
        std::vector<VeryLong> poly_coeff;
        poly_coeff.resize(5);
        poly_coeff[0] = VeryLong("2474902038376743239461332322533453486017145088647579917686") * 30L;
        poly_coeff[1] = VeryLong("848687138292690641752893485837966789178774334641756986027067185671");
        poly_coeff[2] = VeryLong("21211646115309304218501624599244752000939259665291200926111702457") * 2L;
        poly_coeff[3] = VeryLong("40188048938846254718786675529434797027834368184592286864002025034") * 6L;
        poly_coeff[4] = VeryLong("-5073943593046678921562417863134788148090475853480933891772146") * 210L;
        Polynomial<VeryLong> poly(poly_coeff);
        AlgebraicNumber an15(poly, an14);
        an15 = an15 / AlgebraicNumber(210L);
        CPPUNIT_ASSERT(an15 == an10);

        AlgebraicNumber an16;
        an16 = an15;

        CPPUNIT_ASSERT(an16 == an10);

        CPPUNIT_ASSERT(!(an16 != an10));

        AlgebraicNumber an17 = -an16;
        AlgebraicNumber an18 = an17 + an16;

        AlgebraicNumber an19(0L);

        CPPUNIT_ASSERT(an18 == an19);

        AlgebraicNumber an20;
        read(std::string("-99613514790170994695510715502130670121157763057692218497/7 + 808581101117173264056387701675422447113900167454304054501974769/7 alpha + 3255298074118373942254947442634818164093435585383176298316007126/35 alpha^2 + 621448775887439848968335739705960511347337985032517146939084332/35 alpha^3 - 517674952379884810385442290443176235151383151399503587124548 alpha^4"), an20);
        const std::vector<Quotient<VeryLong> >& an20_ibc = an20.ib_coefficients();

        Quotient<VeryLong> an20_ibc0(VeryLong("-8460993248339689717825536197569653104580802169649770073"));
        Quotient<VeryLong> an20_ibc1(VeryLong("73754095890121180984190321559121911906185653017087414758"));
        Quotient<VeryLong> an20_ibc2(VeryLong("59066292933216668405020150301351584598530802211166504080"));
        Quotient<VeryLong> an20_ibc3(VeryLong("48861899533008883568747892157258880690721880440177070740"));
        Quotient<VeryLong> an20_ibc4(VeryLong("-40386562051793166670731962119143098389092147870143827986"));

        CPPUNIT_ASSERT(an20_ibc0 == an20_ibc[0]);
        CPPUNIT_ASSERT(an20_ibc1 == an20_ibc[1]);
        CPPUNIT_ASSERT(an20_ibc2 == an20_ibc[2]);
        CPPUNIT_ASSERT(an20_ibc3 == an20_ibc[3]);
        CPPUNIT_ASSERT(an20_ibc4 == an20_ibc[4]);

        AlgebraicNumber an21;
        read(std::string("216554214229275955286817457004043512484881047489100707614 + 1792309608863510291413033290984827934383779445445860206653001236/5 alpha - 1685051724190470296600982455239705985643198761771037705017087232/5 alpha^2 - 63065500470735307822694010526860900389991748964929361059965108/5 alpha^3 - 2352342375004986357099452396494467149726814909117873361294656 alpha^4"), an21);

        AlgebraicNumber an22 = an20 * an21;

        // These coefficients were calculated independently using pari
        Quotient<VeryLong> q4(VeryLong("30369051015411934990424505037185631224112309249181711246704929265585137346524028269272966051369374824551808526380675288085934753891542877"),VeryLong("12641335874110580000"));
        Quotient<VeryLong> q3(VeryLong("-5246015098668913636704850596397092984636491292222759162200340195930863161351751730480251319498783561857408491082349322632780634482826972023071396277"),VeryLong("6320667937055290000"));
        Quotient<VeryLong> q2(VeryLong("-736170271261009472022571698391436427350466950378199437989484744362712484986976326858033610468015815162947318976576875311222447104709933978730736733452137"),VeryLong("7044410691500"));
        Quotient<VeryLong> q1(VeryLong("278525615298216050967095784054039975024797806400290351695333832743718972394104163166174606750631243742417426065220570423171680004449076609499767514379676981463"),VeryLong("15702050"));
        Quotient<VeryLong> q0(VeryLong("5170479660301472920727431940233008764522036335425604415612231183240001750655417036922757924046493144156461094897437512873471883172730368084093617038313355183796882"),VeryLong("7"));

        AlgebraicNumber theta = AlgebraicNumber::alpha() * AlgebraicNumber::c_d();
        AlgebraicNumber an23 = q0 * AlgebraicNumber(1L) + q1 * theta + q2 * theta * theta + q3 * theta * theta * theta +
                               q4 * theta * theta * theta * theta;

        CPPUNIT_ASSERT(an23 == an22);

        long int pp = 15998663L;

        an22 *= pp;

        AlgebraicNumber an24 = an22 * Quotient<VeryLong>(1L, pp);
        CPPUNIT_ASSERT(an24 == an23);

        VeryLong ppp(pp);

        an24 *= ppp;
        AlgebraicNumber an25 = an24 * Quotient<VeryLong>(1L, ppp);
        CPPUNIT_ASSERT(an25 == an23);

        an25 *= an22;
        AlgebraicNumber an26 = an25 / an22;
        CPPUNIT_ASSERT(an26 == an23);

        std::vector<AlgebraicNumber> basis;
        const VeryLong pppp(2L);
        AlgebraicNumber::createSpecialBasis(pppp, basis);
        // 2
        // 1249913278668 + 2691780 alpha
        // -222405543007291322 + 1249913278668 alpha + 2691780 alpha^2
        // -8591401659640532521423 - 222405543007291322 alpha + 1249913278668 alpha^2 + 2691780 alpha^3
        // 623872110646578368801362410 - 8591401659640532521423 alpha - 222405543007291322 alpha^2 + 1249913278668 alpha^3 + 2691780 alpha^4

        CPPUNIT_ASSERT(basis.size() == 5);
        AlgebraicNumber an30(2L);
        CPPUNIT_ASSERT(basis[0] == an30);
        AlgebraicNumber an31;
        read(std::string("1249913278668 + 2691780 alpha"), an31);
        CPPUNIT_ASSERT(basis[1] == an31);
        AlgebraicNumber an32;
        read(std::string("-222405543007291322 + 1249913278668 alpha + 2691780 alpha^2"), an32);
        CPPUNIT_ASSERT(basis[2] == an32);
        AlgebraicNumber an33;
        read(std::string("-8591401659640532521423 - 222405543007291322 alpha + 1249913278668 alpha^2 + 2691780 alpha^3"), an33);
        CPPUNIT_ASSERT(basis[3] == an33);
        AlgebraicNumber an34;
        read(std::string("623872110646578368801362410 - 8591401659640532521423 alpha - 222405543007291322 alpha^2 + 1249913278668 alpha^3 + 2691780 alpha^4"), an34);
        CPPUNIT_ASSERT(basis[4] == an34);

        std::string s40("623872110646578368801362410 - 8591401659640532521423 alpha - 222405543007291322 alpha^2 + 1249913278668 alpha^3 + 2691780 alpha^4");
        std::istringstream iss(s40);
        AlgebraicNumber an40;
        iss >> an40;
        std::ostringstream oss1;
        oss1 << an40;
        CPPUNIT_ASSERT(s40 == oss1.str());
        std::istringstream iss1(oss1.str());
        AlgebraicNumber an41;
        iss1 >> an41;
        CPPUNIT_ASSERT(an40 == an41);

        CPPUNIT_ASSERT(an41.coefficients().size() == 5);
        // 623872110646578368801362410
        // -8591401659640532521423
        // -222405543007291322
        // 1249913278668
        // 2691780
        CPPUNIT_ASSERT(an41.coefficient(0) == Quotient<VeryLong>("623872110646578368801362410"));
        CPPUNIT_ASSERT(an41.coefficient(1) == Quotient<VeryLong>("-8591401659640532521423"));
        CPPUNIT_ASSERT(an41.coefficient(2) == Quotient<VeryLong>("-222405543007291322"));
        CPPUNIT_ASSERT(an41.coefficient(3) == Quotient<VeryLong>("1249913278668"));
        CPPUNIT_ASSERT(an41.coefficient(4) == Quotient<VeryLong>("2691780"));

        CPPUNIT_ASSERT_THROW_MESSAGE("", an41.coefficient(6), std::string);
        CPPUNIT_ASSERT_THROW_MESSAGE("", an41.coefficient(-1), std::string);

        an41.set_coefficient(0, Quotient<VeryLong>("623872110646578368801362411"));
        std::string s41("623872110646578368801362411 - 8591401659640532521423 alpha - 222405543007291322 alpha^2 + 1249913278668 alpha^3 + 2691780 alpha^4");
        std::ostringstream oss2;
        oss2 << an41;
        CPPUNIT_ASSERT(s41 == oss2.str());

        CPPUNIT_ASSERT_THROW_MESSAGE("", an41.set_coefficient(6, Quotient<VeryLong>("1")), std::string);

        long double ln_re;
        long int re_sign;
        long double ln_im;
        long int im_sign;
        // 61.80787892574, 1, -inf, 1
        an41.ln_sigma(0, ln_re, re_sign, ln_im, im_sign);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(ln_re, 61.80787892574, std::numeric_limits<float>::epsilon());
        CPPUNIT_ASSERT_DOUBLES_EQUAL(an41.ln_sigma(0), 61.80787892574, std::numeric_limits<float>::epsilon());
        CPPUNIT_ASSERT(re_sign == 1);
        // 60.4991536271959, -1, -inf, -1
        an41.ln_sigma(1, ln_re, re_sign, ln_im, im_sign);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(ln_re, 60.4991536271959, std::numeric_limits<float>::epsilon());
        CPPUNIT_ASSERT_DOUBLES_EQUAL(an41.ln_sigma(1), 60.4991536271959, std::numeric_limits<float>::epsilon());
        CPPUNIT_ASSERT(re_sign == -1);
        // 59.396337153399, -1, -inf, -1
        an41.ln_sigma(2, ln_re, re_sign, ln_im, im_sign);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(ln_re, 59.396337153399, std::numeric_limits<float>::epsilon());
        CPPUNIT_ASSERT_DOUBLES_EQUAL(an41.ln_sigma(2), 59.396337153399, std::numeric_limits<float>::epsilon());
        CPPUNIT_ASSERT(re_sign == -1);
        // 60.3512548653201, 1, -inf, 1
        an41.ln_sigma(3, ln_re, re_sign, ln_im, im_sign);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(ln_re, 60.3512548653201, std::numeric_limits<float>::epsilon());
        CPPUNIT_ASSERT_DOUBLES_EQUAL(an41.ln_sigma(3), 60.3512548653201, std::numeric_limits<float>::epsilon());
        CPPUNIT_ASSERT(re_sign == 1);
        // 58.0247348162748, 1, -inf, 1
        an41.ln_sigma(4, ln_re, re_sign, ln_im, im_sign);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(ln_re, 58.0247348162748, std::numeric_limits<float>::epsilon());
        CPPUNIT_ASSERT_DOUBLES_EQUAL(an41.ln_sigma(4), 58.0247348162748, std::numeric_limits<float>::epsilon());
        CPPUNIT_ASSERT(re_sign == 1);

        CPPUNIT_ASSERT_THROW_MESSAGE("", an41.ln_sigma(-1, ln_re, re_sign, ln_im, im_sign), std::string);
        CPPUNIT_ASSERT_THROW_MESSAGE("", an41.ln_sigma(-1), std::string);
        CPPUNIT_ASSERT_THROW_MESSAGE("", an41.ln_sigma(5, ln_re, re_sign, ln_im, im_sign), std::string);
        CPPUNIT_ASSERT_THROW_MESSAGE("", an41.ln_sigma(6), std::string);

        // 4.84888083570085e+53
        CPPUNIT_ASSERT_DOUBLES_EQUAL(an41.mod_sigma_2(0), 4.84888083570085e+53, 1e39);
        // 3.53913525374421e+52
        CPPUNIT_ASSERT_DOUBLES_EQUAL(an41.mod_sigma_2(1), 3.53913525374421e+52, 1e38);
        // 3.89944628207922e+51
        CPPUNIT_ASSERT_DOUBLES_EQUAL(an41.mod_sigma_2(2), 3.89944628207922e+51, 1e37);
        // 2.63289735301011e+52
        CPPUNIT_ASSERT_DOUBLES_EQUAL(an41.mod_sigma_2(3), 2.63289735301011e+52, 1e38);
        // 2.50982990232009e+50
        CPPUNIT_ASSERT_DOUBLES_EQUAL(an41.mod_sigma_2(4), 2.50982990232009e+50, 1e36);
        CPPUNIT_ASSERT_THROW_MESSAGE("", an41.mod_sigma_2(-1), std::string);
        CPPUNIT_ASSERT_THROW_MESSAGE("", an41.mod_sigma_2(7), std::string);

    }

    void testAlgebraicNumber_in_O_pO()
    {
        NumberField nf(f, "nf.fb.dat");
        AlgebraicNumber::setNumberField(nf);
        CPPUNIT_ASSERT(AlgebraicNumber::index() == nf.index());
        CPPUNIT_ASSERT(AlgebraicNumber::c_d() == nf.c_d());
        CPPUNIT_ASSERT(AlgebraicNumber::degree() == 5);

        long int p(23L);
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", AlgebraicNumber_in_O_pO::set_basis(p));

        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", AlgebraicNumber_in_O_pO());
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", AlgebraicNumber_in_O_pO(0L));
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", AlgebraicNumber_in_O_pO(1L));
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", AlgebraicNumber_in_O_pO(11L));

        AlgebraicNumber_in_O_pO an1;
        AlgebraicNumber_in_O_pO an2(0L);

        CPPUNIT_ASSERT(an1 == an2);
        CPPUNIT_ASSERT(!(an1 != an2));

        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", AlgebraicNumber_in_O_pO(0, 0));
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", AlgebraicNumber_in_O_pO(0, 1));
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", AlgebraicNumber_in_O_pO(0, 2));
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", AlgebraicNumber_in_O_pO(0, 3));
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", AlgebraicNumber_in_O_pO(0, 4));
        CPPUNIT_ASSERT_THROW_MESSAGE("", AlgebraicNumber_in_O_pO(0, 5), std::string);
        CPPUNIT_ASSERT_THROW_MESSAGE("", AlgebraicNumber_in_O_pO(0, NumberField::MAX_DEGREE), std::string);
        std::vector<VeryLongModular> c;
        c.push_back(VeryLongModular(1L));
        c.push_back(VeryLongModular(1L));
        c.push_back(VeryLongModular(1L));
        c.push_back(VeryLongModular(1L));
        c.push_back(VeryLongModular(1L));
        c.push_back(VeryLongModular(1L));
        c.push_back(VeryLongModular(1L));
        CPPUNIT_ASSERT_THROW_MESSAGE("", {AlgebraicNumber_in_O_pO a(c);}, std::string);

        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", AlgebraicNumber::read_algebraic_number("1/7 + 32730265/42 alpha + 275165447/105 alpha^2 + 7409138/35 alpha^3 + 12818 alpha^4"));
        AlgebraicNumber an3 = AlgebraicNumber::read_algebraic_number("1/7 + 32730265/42 alpha + 275165447/105 alpha^2 + 7409138/35 alpha^3 + 12818 alpha^4");
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", (AlgebraicNumber_in_O_pO(an3)));
        AlgebraicNumber_in_O_pO an4(an3);
        AlgebraicNumber_in_O_pO an5(0, 4);
        CPPUNIT_ASSERT(an4 == an5);

        AlgebraicNumber an6 = AlgebraicNumber::read_algebraic_number("1 + alpha + alpha^2 + alpha^3 + alpha^4");
        AlgebraicNumber_in_O_pO an7(an6);
        c.clear();
        c.push_back(VeryLongModular(16L));
        c.push_back(VeryLongModular(1L));
        c.push_back(VeryLongModular(14L));
        c.push_back(VeryLongModular(3L));
        c.push_back(VeryLongModular(10L));
        AlgebraicNumber_in_O_pO an8(c);
        CPPUNIT_ASSERT(an8 == an7);

        AlgebraicNumber_in_O_pO an9(10LL, 7L);
        AlgebraicNumber_in_O_pO an10(AlgebraicNumber::read_algebraic_number("10 - 7 alpha"));
        CPPUNIT_ASSERT(an9 == an10);
        AlgebraicNumber_in_O_pO an11(10LL + 23LL*134023580LL, 7L + 23L*1048012L);
        CPPUNIT_ASSERT(an9 == an11);

        AlgebraicNumber an12(AlgebraicNumber::read_algebraic_number("1241 + 323 alpha - 230923029 alpha^2 + 2342 alpha^3 + 989225208 alpha^4 "));
        AlgebraicNumber an13(AlgebraicNumber::read_algebraic_number("-82738 + 33 alpha + 3003629 alpha^2 + alpha^3 + 985908 alpha^4 "));
        AlgebraicNumber an14 = an12 * an13;
        AlgebraicNumber_in_O_pO an15(an12);
        AlgebraicNumber_in_O_pO an16(an13);
        AlgebraicNumber_in_O_pO an17 = an15 * an16;
        CPPUNIT_ASSERT(an17 == AlgebraicNumber_in_O_pO(an14));

        AlgebraicNumber_in_O_pO an18(an15);
        an18 *= an16;
        CPPUNIT_ASSERT(an18 == AlgebraicNumber_in_O_pO(an14));

        AlgebraicNumber_in_O_pO an19(an15);
        an19 *= an13;
        CPPUNIT_ASSERT(an19 == AlgebraicNumber_in_O_pO(an14));

        AlgebraicNumber_in_O_pO an20 = an19 / an13;
        CPPUNIT_ASSERT(an20 == an15);

        CPPUNIT_ASSERT_THROW_MESSAGE("", an20.coefficient(-1), std::string);
        CPPUNIT_ASSERT_THROW_MESSAGE("", an20.coefficient(10), std::string);
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", an20.coefficient(0));
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", an20.coefficient(1));
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", an20.coefficient(2));
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", an20.coefficient(3));
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", an20.coefficient(4));

        CPPUNIT_ASSERT(an20.coefficient(0) == VeryLongModular(17L));
        CPPUNIT_ASSERT(an20.coefficient(1) == VeryLongModular(6L));
        CPPUNIT_ASSERT(an20.coefficient(2) == VeryLongModular(9L));
        CPPUNIT_ASSERT(an20.coefficient(3) == VeryLongModular(12L));
        CPPUNIT_ASSERT(an20.coefficient(4) == VeryLongModular(12L));

        const VeryLongModular* b = an20.basis();
        CPPUNIT_ASSERT(an20.coefficient(3) == b[3]);

        const Matrix<Quotient<VeryLong> >& w = AlgebraicNumber_in_O_pO::W_mult();
        CPPUNIT_ASSERT(w.rows() == 5);
        CPPUNIT_ASSERT(w.columns() == 25);

        CPPUNIT_ASSERT(VeryLong("1") == w(0,0).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(0,1).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(0,2).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(0,3).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(0,4).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(0,5).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(0,6).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(0,7).numerator());
        CPPUNIT_ASSERT(VeryLong("-8075340") == w(0,8).numerator());
        CPPUNIT_ASSERT(VeryLong("-60253580560124974875864154071577740") == w(0,9).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(0,10).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(0,11).numerator());
        CPPUNIT_ASSERT(VeryLong("-80753400") == w(0,12).numerator());
        CPPUNIT_ASSERT(VeryLong("-5061300767050497889572588941770512708") == w(0,13).numerator());
        CPPUNIT_ASSERT(VeryLong("55954847297396316526588059915345865764012") == w(0,14).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(0,15).numerator());
        CPPUNIT_ASSERT(VeryLong("-8075340") == w(0,16).numerator());
        CPPUNIT_ASSERT(VeryLong("-5061300767050497889572588941770512708") == w(0,17).numerator());
        CPPUNIT_ASSERT(VeryLong("470034382810200095167641549279048628733028") == w(0,18).numerator());
        CPPUNIT_ASSERT(VeryLong("-7187804455625270788403051126189423952358060800") == w(0,19).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(0,20).numerator());
        CPPUNIT_ASSERT(VeryLong("-60253580560124974875864154071577740") == w(0,21).numerator());
        CPPUNIT_ASSERT(VeryLong("55954847297396316526588059915345865764012") == w(0,22).numerator());
        CPPUNIT_ASSERT(VeryLong("-7187804455625270788403051126189423952358060800") == w(0,23).numerator());
        CPPUNIT_ASSERT(VeryLong("99648192905093100537294590751735723982880115924947") == w(0,24).numerator());

        CPPUNIT_ASSERT(VeryLong("0") == w(1,0).numerator());
        CPPUNIT_ASSERT(VeryLong("1") == w(1,1).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(1,2).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(1,3).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(1,4).numerator());
        CPPUNIT_ASSERT(VeryLong("1") == w(1,5).numerator());
        CPPUNIT_ASSERT(VeryLong("-20229") == w(1,6).numerator());
        CPPUNIT_ASSERT(VeryLong("-4577887") == w(1,7).numerator());
        CPPUNIT_ASSERT(VeryLong("-17358055") == w(1,8).numerator());
        CPPUNIT_ASSERT(VeryLong("-2970820191215456207309219") == w(1,9).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(1,10).numerator());
        CPPUNIT_ASSERT(VeryLong("-4577887") == w(1,11).numerator());
        CPPUNIT_ASSERT(VeryLong("-159846889") == w(1,12).numerator());
        CPPUNIT_ASSERT(VeryLong("-249548896062098320563624894") == w(1,13).numerator());
        CPPUNIT_ASSERT(VeryLong("2669332840668229707871499869061") == w(1,14).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(1,15).numerator());
        CPPUNIT_ASSERT(VeryLong("-17358055") == w(1,16).numerator());
        CPPUNIT_ASSERT(VeryLong("-249548896062098320563624894") == w(1,17).numerator());
        CPPUNIT_ASSERT(VeryLong("22423069643632497211587271584927") == w(1,18).numerator());
        CPPUNIT_ASSERT(VeryLong("-346081888125177508551699363688172746") == w(1,19).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(1,20).numerator());
        CPPUNIT_ASSERT(VeryLong("-2970820191215456207309219") == w(1,21).numerator());
        CPPUNIT_ASSERT(VeryLong("2669332840668229707871499869061") == w(1,22).numerator());
        CPPUNIT_ASSERT(VeryLong("-346081888125177508551699363688172746") == w(1,23).numerator());
        CPPUNIT_ASSERT(VeryLong("4786029865631083395935292174929311650534") == w(1,24).numerator());

        CPPUNIT_ASSERT(VeryLong("0") == w(2,0).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(2,1).numerator());
        CPPUNIT_ASSERT(VeryLong("1") == w(2,2).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(2,3).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(2,4).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(2,5).numerator());
        CPPUNIT_ASSERT(VeryLong("672945") == w(2,6).numerator());
        CPPUNIT_ASSERT(VeryLong("-2018835") == w(2,7).numerator());
        CPPUNIT_ASSERT(VeryLong("-48452040") == w(2,8).numerator());
        CPPUNIT_ASSERT(VeryLong("20454916506633907935") == w(2,9).numerator());
        CPPUNIT_ASSERT(VeryLong("1") == w(2,10).numerator());
        CPPUNIT_ASSERT(VeryLong("-2018835") == w(2,11).numerator());
        CPPUNIT_ASSERT(VeryLong("-483102469") == w(2,12).numerator());
        CPPUNIT_ASSERT(VeryLong("1718212986558828369317") == w(2,13).numerator());
        CPPUNIT_ASSERT(VeryLong("-21966462075117512935209077") == w(2,14).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(2,15).numerator());
        CPPUNIT_ASSERT(VeryLong("-48452040") == w(2,16).numerator());
        CPPUNIT_ASSERT(VeryLong("1718212986558828369317") == w(2,17).numerator());
        CPPUNIT_ASSERT(VeryLong("-184522920606050814041768997") == w(2,18).numerator());
        CPPUNIT_ASSERT(VeryLong("2707058551024412718403172372547") == w(2,19).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(2,20).numerator());
        CPPUNIT_ASSERT(VeryLong("20454916506633907935") == w(2,21).numerator());
        CPPUNIT_ASSERT(VeryLong("-21966462075117512935209077") == w(2,22).numerator());
        CPPUNIT_ASSERT(VeryLong("2707058551024412718403172372547") == w(2,23).numerator());
        CPPUNIT_ASSERT(VeryLong("-37948643836859776737649671244042563") == w(2,24).numerator());

        CPPUNIT_ASSERT(VeryLong("0") == w(3,0).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(3,1).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(3,2).numerator());
        CPPUNIT_ASSERT(VeryLong("1") == w(3,3).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(3,4).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(3,5).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(3,6).numerator());
        CPPUNIT_ASSERT(VeryLong("6729450") == w(3,7).numerator());
        CPPUNIT_ASSERT(VeryLong("-20188350") == w(3,8).numerator());
        CPPUNIT_ASSERT(VeryLong("2647930775940450") == w(3,9).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(3,10).numerator());
        CPPUNIT_ASSERT(VeryLong("6729450") == w(3,11).numerator());
        CPPUNIT_ASSERT(VeryLong("-221869560") == w(3,12).numerator());
        CPPUNIT_ASSERT(VeryLong("222426185299521030") == w(3,13).numerator());
        CPPUNIT_ASSERT(VeryLong("-2254546995920819152860") == w(3,14).numerator());
        CPPUNIT_ASSERT(VeryLong("1") == w(3,15).numerator());
        CPPUNIT_ASSERT(VeryLong("-20188350") == w(3,16).numerator());
        CPPUNIT_ASSERT(VeryLong("222426185299521030") == w(3,17).numerator());
        CPPUNIT_ASSERT(VeryLong("-18938795316433623221689") == w(3,18).numerator());
        CPPUNIT_ASSERT(VeryLong("293918829423809106677588968") == w(3,19).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(3,20).numerator());
        CPPUNIT_ASSERT(VeryLong("2647930775940450") == w(3,21).numerator());
        CPPUNIT_ASSERT(VeryLong("-2254546995920819152860") == w(3,22).numerator());
        CPPUNIT_ASSERT(VeryLong("293918829423809106677588968") == w(3,23).numerator());
        CPPUNIT_ASSERT(VeryLong("-4057027952193313585761437354175") == w(3,24).numerator());

        CPPUNIT_ASSERT(VeryLong("0") == w(4,0).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(4,1).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(4,2).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(4,3).numerator());
        CPPUNIT_ASSERT(VeryLong("1") == w(4,4).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(4,5).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(4,6).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(4,7).numerator());
        CPPUNIT_ASSERT(VeryLong("56527380") == w(4,8).numerator());
        CPPUNIT_ASSERT(VeryLong("-624934411920") == w(4,9).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(4,10).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(4,11).numerator());
        CPPUNIT_ASSERT(VeryLong("565273800") == w(4,12).numerator());
        CPPUNIT_ASSERT(VeryLong("-52496184723444") == w(4,13).numerator());
        CPPUNIT_ASSERT(VeryLong("802775257629195096") == w(4,14).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(4,15).numerator());
        CPPUNIT_ASSERT(VeryLong("56527380") == w(4,16).numerator());
        CPPUNIT_ASSERT(VeryLong("-52496184723444") == w(4,17).numerator());
        CPPUNIT_ASSERT(VeryLong("6743453899752474996") == w(4,18).numerator());
        CPPUNIT_ASSERT(VeryLong("-93488159898472500455580") == w(4,19).numerator());
        CPPUNIT_ASSERT(VeryLong("1") == w(4,20).numerator());
        CPPUNIT_ASSERT(VeryLong("-624934411920") == w(4,21).numerator());
        CPPUNIT_ASSERT(VeryLong("802775257629195096") == w(4,22).numerator());
        CPPUNIT_ASSERT(VeryLong("-93488159898472500455580") == w(4,23).numerator());
        CPPUNIT_ASSERT(VeryLong("1327437115515222519031202729") == w(4,24).numerator());

        for (size_t i = 0; i < 5; ++i)
        {
            for (size_t j = 0; j < 25; ++j)
            {
                CPPUNIT_ASSERT(VeryLong("1") == w(i,j).denominator());
            }
        }

        AlgebraicNumber_in_O_pO::set_basis(29L);
        Matrix<VeryLongModular> beta1 = AlgebraicNumber_in_O_pO::make_beta(VeryLong(29L));
        CPPUNIT_ASSERT(beta1.rows() == 5);
        CPPUNIT_ASSERT(beta1.columns() == 0);

        AlgebraicNumber_in_O_pO::set_basis(2L);
        Matrix<VeryLongModular> beta2 = AlgebraicNumber_in_O_pO::make_beta(VeryLong(8L));
        CPPUNIT_ASSERT(beta2(1,0) == VeryLongModular(1L));
        CPPUNIT_ASSERT(beta2(2,1) == VeryLongModular(1L));

        AlgebraicNumber_in_O_pO::set_basis(5L);
        Matrix<VeryLongModular> beta3 = AlgebraicNumber_in_O_pO::make_beta(VeryLong(5L));
        CPPUNIT_ASSERT(beta3.rows() == 5);
        CPPUNIT_ASSERT(beta3.columns() == 1);
        CPPUNIT_ASSERT(beta3(0,0) == VeryLongModular(0L));
        CPPUNIT_ASSERT(beta3(1,0) == VeryLongModular(2L));
        CPPUNIT_ASSERT(beta3(2,0) == VeryLongModular(4L));
        CPPUNIT_ASSERT(beta3(3,0) == VeryLongModular(4L));
        CPPUNIT_ASSERT(beta3(4,0) == VeryLongModular(1L));

        CPPUNIT_ASSERT_THROW_MESSAGE("", AlgebraicNumber_in_O_pO(1LL, 2L), std::string);
        AlgebraicNumber_in_O_pO::set_basis(113L);
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", AlgebraicNumber_in_O_pO(1LL, 2L));
    }

    void testAlgebraicNumber_in_O_pO_1()
    {
        NumberField nf(f, "nf.fb.dat");
        AlgebraicNumber::setNumberField(nf);
        CPPUNIT_ASSERT(AlgebraicNumber::index() == nf.index());
        CPPUNIT_ASSERT(AlgebraicNumber::c_d() == nf.c_d());
        CPPUNIT_ASSERT(AlgebraicNumber::degree() == 5);

        long int p(23L);
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", AlgebraicNumber_in_O_pO_1::set_basis(p));

        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", AlgebraicNumber_in_O_pO_1());
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", AlgebraicNumber_in_O_pO_1(0L));
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", AlgebraicNumber_in_O_pO_1(1L));
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", AlgebraicNumber_in_O_pO_1(11L));

        AlgebraicNumber_in_O_pO_1 an1;
        AlgebraicNumber_in_O_pO_1 an2(0L);

        CPPUNIT_ASSERT(an1 == an2);
        CPPUNIT_ASSERT(!(an1 != an2));

        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", AlgebraicNumber_in_O_pO_1(0, 0));
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", AlgebraicNumber_in_O_pO_1(0, 1));
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", AlgebraicNumber_in_O_pO_1(0, 2));
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", AlgebraicNumber_in_O_pO_1(0, 3));
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", AlgebraicNumber_in_O_pO_1(0, 4));
        CPPUNIT_ASSERT_THROW_MESSAGE("", AlgebraicNumber_in_O_pO_1(0, 5), std::string);
        CPPUNIT_ASSERT_THROW_MESSAGE("", AlgebraicNumber_in_O_pO_1(0, NumberField::MAX_DEGREE), std::string);
        std::vector<LongModular> c;
        c.push_back(LongModular(1L));
        c.push_back(LongModular(1L));
        c.push_back(LongModular(1L));
        c.push_back(LongModular(1L));
        c.push_back(LongModular(1L));
        c.push_back(LongModular(1L));
        c.push_back(LongModular(1L));
        CPPUNIT_ASSERT_THROW_MESSAGE("", {AlgebraicNumber_in_O_pO_1 a(c);}, std::string);

        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", AlgebraicNumber::read_algebraic_number("1/7 + 32730265/42 alpha + 275165447/105 alpha^2 + 7409138/35 alpha^3 + 12818 alpha^4"));
        AlgebraicNumber an3 = AlgebraicNumber::read_algebraic_number("1/7 + 32730265/42 alpha + 275165447/105 alpha^2 + 7409138/35 alpha^3 + 12818 alpha^4");
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", (AlgebraicNumber_in_O_pO_1(an3)));
        AlgebraicNumber_in_O_pO_1 an4(an3);
        AlgebraicNumber_in_O_pO_1 an5(0, 4);
        CPPUNIT_ASSERT(an4 == an5);

        AlgebraicNumber an6 = AlgebraicNumber::read_algebraic_number("1 + alpha + alpha^2 + alpha^3 + alpha^4");
        AlgebraicNumber_in_O_pO_1 an7(an6);
        c.clear();
        c.push_back(LongModular(16L));
        c.push_back(LongModular(1L));
        c.push_back(LongModular(14L));
        c.push_back(LongModular(3L));
        c.push_back(LongModular(10L));
        AlgebraicNumber_in_O_pO_1 an8(c);
        CPPUNIT_ASSERT(an8 == an7);

        AlgebraicNumber_in_O_pO_1 an9(10LL, 7L);
        AlgebraicNumber_in_O_pO_1 an10(AlgebraicNumber::read_algebraic_number("10 - 7 alpha"));
        CPPUNIT_ASSERT(an9 == an10);
        AlgebraicNumber_in_O_pO_1 an11(10LL + 23LL*180LL, 7L + 23L*18012L);
        CPPUNIT_ASSERT(an9 == an11);

        AlgebraicNumber an12(AlgebraicNumber::read_algebraic_number("1241 + 323 alpha - 230923029 alpha^2 + 2342 alpha^3 + 989225208 alpha^4 "));
        AlgebraicNumber an13(AlgebraicNumber::read_algebraic_number("-82738 + 33 alpha + 3003629 alpha^2 + alpha^3 + 985908 alpha^4 "));
        AlgebraicNumber an14 = an12 * an13;
        AlgebraicNumber_in_O_pO_1 an15(an12);
        AlgebraicNumber_in_O_pO_1 an16(an13);
        AlgebraicNumber_in_O_pO_1 an17 = an15 * an16;
        CPPUNIT_ASSERT(an17 == AlgebraicNumber_in_O_pO_1(an14));

        AlgebraicNumber_in_O_pO_1 an18(an15);
        an18 *= an16;
        CPPUNIT_ASSERT(an18 == AlgebraicNumber_in_O_pO_1(an14));

        AlgebraicNumber_in_O_pO_1 an19(an15);
        an19 *= an13;
        CPPUNIT_ASSERT(an19 == AlgebraicNumber_in_O_pO_1(an14));

        AlgebraicNumber_in_O_pO_1 an20 = an19 / an13;
        CPPUNIT_ASSERT(an20 == an15);

        CPPUNIT_ASSERT_THROW_MESSAGE("", an20.coefficient(-1), std::string);
        CPPUNIT_ASSERT_THROW_MESSAGE("", an20.coefficient(10), std::string);
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", an20.coefficient(0));
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", an20.coefficient(1));
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", an20.coefficient(2));
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", an20.coefficient(3));
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", an20.coefficient(4));

        CPPUNIT_ASSERT(an20.coefficient(0) == LongModular(17L));
        CPPUNIT_ASSERT(an20.coefficient(1) == LongModular(6L));
        CPPUNIT_ASSERT(an20.coefficient(2) == LongModular(9L));
        CPPUNIT_ASSERT(an20.coefficient(3) == LongModular(12L));
        CPPUNIT_ASSERT(an20.coefficient(4) == LongModular(12L));

        const LongModular* b = an20.basis();
        CPPUNIT_ASSERT(an20.coefficient(3) == b[3]);

        const Matrix<Quotient<VeryLong> >& w = AlgebraicNumber_in_O_pO_1::W_mult();
        CPPUNIT_ASSERT(w.rows() == 5);
        CPPUNIT_ASSERT(w.columns() == 25);

        CPPUNIT_ASSERT(VeryLong("1") == w(0,0).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(0,1).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(0,2).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(0,3).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(0,4).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(0,5).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(0,6).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(0,7).numerator());
        CPPUNIT_ASSERT(VeryLong("-8075340") == w(0,8).numerator());
        CPPUNIT_ASSERT(VeryLong("-60253580560124974875864154071577740") == w(0,9).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(0,10).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(0,11).numerator());
        CPPUNIT_ASSERT(VeryLong("-80753400") == w(0,12).numerator());
        CPPUNIT_ASSERT(VeryLong("-5061300767050497889572588941770512708") == w(0,13).numerator());
        CPPUNIT_ASSERT(VeryLong("55954847297396316526588059915345865764012") == w(0,14).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(0,15).numerator());
        CPPUNIT_ASSERT(VeryLong("-8075340") == w(0,16).numerator());
        CPPUNIT_ASSERT(VeryLong("-5061300767050497889572588941770512708") == w(0,17).numerator());
        CPPUNIT_ASSERT(VeryLong("470034382810200095167641549279048628733028") == w(0,18).numerator());
        CPPUNIT_ASSERT(VeryLong("-7187804455625270788403051126189423952358060800") == w(0,19).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(0,20).numerator());
        CPPUNIT_ASSERT(VeryLong("-60253580560124974875864154071577740") == w(0,21).numerator());
        CPPUNIT_ASSERT(VeryLong("55954847297396316526588059915345865764012") == w(0,22).numerator());
        CPPUNIT_ASSERT(VeryLong("-7187804455625270788403051126189423952358060800") == w(0,23).numerator());
        CPPUNIT_ASSERT(VeryLong("99648192905093100537294590751735723982880115924947") == w(0,24).numerator());

        CPPUNIT_ASSERT(VeryLong("0") == w(1,0).numerator());
        CPPUNIT_ASSERT(VeryLong("1") == w(1,1).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(1,2).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(1,3).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(1,4).numerator());
        CPPUNIT_ASSERT(VeryLong("1") == w(1,5).numerator());
        CPPUNIT_ASSERT(VeryLong("-20229") == w(1,6).numerator());
        CPPUNIT_ASSERT(VeryLong("-4577887") == w(1,7).numerator());
        CPPUNIT_ASSERT(VeryLong("-17358055") == w(1,8).numerator());
        CPPUNIT_ASSERT(VeryLong("-2970820191215456207309219") == w(1,9).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(1,10).numerator());
        CPPUNIT_ASSERT(VeryLong("-4577887") == w(1,11).numerator());
        CPPUNIT_ASSERT(VeryLong("-159846889") == w(1,12).numerator());
        CPPUNIT_ASSERT(VeryLong("-249548896062098320563624894") == w(1,13).numerator());
        CPPUNIT_ASSERT(VeryLong("2669332840668229707871499869061") == w(1,14).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(1,15).numerator());
        CPPUNIT_ASSERT(VeryLong("-17358055") == w(1,16).numerator());
        CPPUNIT_ASSERT(VeryLong("-249548896062098320563624894") == w(1,17).numerator());
        CPPUNIT_ASSERT(VeryLong("22423069643632497211587271584927") == w(1,18).numerator());
        CPPUNIT_ASSERT(VeryLong("-346081888125177508551699363688172746") == w(1,19).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(1,20).numerator());
        CPPUNIT_ASSERT(VeryLong("-2970820191215456207309219") == w(1,21).numerator());
        CPPUNIT_ASSERT(VeryLong("2669332840668229707871499869061") == w(1,22).numerator());
        CPPUNIT_ASSERT(VeryLong("-346081888125177508551699363688172746") == w(1,23).numerator());
        CPPUNIT_ASSERT(VeryLong("4786029865631083395935292174929311650534") == w(1,24).numerator());

        CPPUNIT_ASSERT(VeryLong("0") == w(2,0).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(2,1).numerator());
        CPPUNIT_ASSERT(VeryLong("1") == w(2,2).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(2,3).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(2,4).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(2,5).numerator());
        CPPUNIT_ASSERT(VeryLong("672945") == w(2,6).numerator());
        CPPUNIT_ASSERT(VeryLong("-2018835") == w(2,7).numerator());
        CPPUNIT_ASSERT(VeryLong("-48452040") == w(2,8).numerator());
        CPPUNIT_ASSERT(VeryLong("20454916506633907935") == w(2,9).numerator());
        CPPUNIT_ASSERT(VeryLong("1") == w(2,10).numerator());
        CPPUNIT_ASSERT(VeryLong("-2018835") == w(2,11).numerator());
        CPPUNIT_ASSERT(VeryLong("-483102469") == w(2,12).numerator());
        CPPUNIT_ASSERT(VeryLong("1718212986558828369317") == w(2,13).numerator());
        CPPUNIT_ASSERT(VeryLong("-21966462075117512935209077") == w(2,14).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(2,15).numerator());
        CPPUNIT_ASSERT(VeryLong("-48452040") == w(2,16).numerator());
        CPPUNIT_ASSERT(VeryLong("1718212986558828369317") == w(2,17).numerator());
        CPPUNIT_ASSERT(VeryLong("-184522920606050814041768997") == w(2,18).numerator());
        CPPUNIT_ASSERT(VeryLong("2707058551024412718403172372547") == w(2,19).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(2,20).numerator());
        CPPUNIT_ASSERT(VeryLong("20454916506633907935") == w(2,21).numerator());
        CPPUNIT_ASSERT(VeryLong("-21966462075117512935209077") == w(2,22).numerator());
        CPPUNIT_ASSERT(VeryLong("2707058551024412718403172372547") == w(2,23).numerator());
        CPPUNIT_ASSERT(VeryLong("-37948643836859776737649671244042563") == w(2,24).numerator());

        CPPUNIT_ASSERT(VeryLong("0") == w(3,0).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(3,1).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(3,2).numerator());
        CPPUNIT_ASSERT(VeryLong("1") == w(3,3).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(3,4).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(3,5).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(3,6).numerator());
        CPPUNIT_ASSERT(VeryLong("6729450") == w(3,7).numerator());
        CPPUNIT_ASSERT(VeryLong("-20188350") == w(3,8).numerator());
        CPPUNIT_ASSERT(VeryLong("2647930775940450") == w(3,9).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(3,10).numerator());
        CPPUNIT_ASSERT(VeryLong("6729450") == w(3,11).numerator());
        CPPUNIT_ASSERT(VeryLong("-221869560") == w(3,12).numerator());
        CPPUNIT_ASSERT(VeryLong("222426185299521030") == w(3,13).numerator());
        CPPUNIT_ASSERT(VeryLong("-2254546995920819152860") == w(3,14).numerator());
        CPPUNIT_ASSERT(VeryLong("1") == w(3,15).numerator());
        CPPUNIT_ASSERT(VeryLong("-20188350") == w(3,16).numerator());
        CPPUNIT_ASSERT(VeryLong("222426185299521030") == w(3,17).numerator());
        CPPUNIT_ASSERT(VeryLong("-18938795316433623221689") == w(3,18).numerator());
        CPPUNIT_ASSERT(VeryLong("293918829423809106677588968") == w(3,19).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(3,20).numerator());
        CPPUNIT_ASSERT(VeryLong("2647930775940450") == w(3,21).numerator());
        CPPUNIT_ASSERT(VeryLong("-2254546995920819152860") == w(3,22).numerator());
        CPPUNIT_ASSERT(VeryLong("293918829423809106677588968") == w(3,23).numerator());
        CPPUNIT_ASSERT(VeryLong("-4057027952193313585761437354175") == w(3,24).numerator());

        CPPUNIT_ASSERT(VeryLong("0") == w(4,0).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(4,1).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(4,2).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(4,3).numerator());
        CPPUNIT_ASSERT(VeryLong("1") == w(4,4).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(4,5).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(4,6).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(4,7).numerator());
        CPPUNIT_ASSERT(VeryLong("56527380") == w(4,8).numerator());
        CPPUNIT_ASSERT(VeryLong("-624934411920") == w(4,9).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(4,10).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(4,11).numerator());
        CPPUNIT_ASSERT(VeryLong("565273800") == w(4,12).numerator());
        CPPUNIT_ASSERT(VeryLong("-52496184723444") == w(4,13).numerator());
        CPPUNIT_ASSERT(VeryLong("802775257629195096") == w(4,14).numerator());
        CPPUNIT_ASSERT(VeryLong("0") == w(4,15).numerator());
        CPPUNIT_ASSERT(VeryLong("56527380") == w(4,16).numerator());
        CPPUNIT_ASSERT(VeryLong("-52496184723444") == w(4,17).numerator());
        CPPUNIT_ASSERT(VeryLong("6743453899752474996") == w(4,18).numerator());
        CPPUNIT_ASSERT(VeryLong("-93488159898472500455580") == w(4,19).numerator());
        CPPUNIT_ASSERT(VeryLong("1") == w(4,20).numerator());
        CPPUNIT_ASSERT(VeryLong("-624934411920") == w(4,21).numerator());
        CPPUNIT_ASSERT(VeryLong("802775257629195096") == w(4,22).numerator());
        CPPUNIT_ASSERT(VeryLong("-93488159898472500455580") == w(4,23).numerator());
        CPPUNIT_ASSERT(VeryLong("1327437115515222519031202729") == w(4,24).numerator());

        for (size_t i = 0; i < 5; ++i)
        {
            for (size_t j = 0; j < 25; ++j)
            {
                CPPUNIT_ASSERT(VeryLong("1") == w(i,j).denominator());
            }
        }

        AlgebraicNumber_in_O_pO_1::set_basis(29L);
        Matrix<LongModular> beta1 = AlgebraicNumber_in_O_pO_1::make_beta(VeryLong(29L));
        CPPUNIT_ASSERT(beta1.rows() == 5);
        CPPUNIT_ASSERT(beta1.columns() == 0);

        AlgebraicNumber_in_O_pO_1::set_basis(2L);
        Matrix<LongModular> beta2 = AlgebraicNumber_in_O_pO_1::make_beta(VeryLong(8L));
        CPPUNIT_ASSERT(beta2(1,0) == LongModular(1L));
        CPPUNIT_ASSERT(beta2(2,1) == LongModular(1L));

        AlgebraicNumber_in_O_pO_1::set_basis(5L);
        Matrix<LongModular> beta3 = AlgebraicNumber_in_O_pO_1::make_beta(VeryLong(5L));
        CPPUNIT_ASSERT(beta3.rows() == 5);
        CPPUNIT_ASSERT(beta3.columns() == 1);
        CPPUNIT_ASSERT(beta3(0,0) == LongModular(0L));
        CPPUNIT_ASSERT(beta3(1,0) == LongModular(2L));
        CPPUNIT_ASSERT(beta3(2,0) == LongModular(4L));
        CPPUNIT_ASSERT(beta3(3,0) == LongModular(4L));
        CPPUNIT_ASSERT(beta3(4,0) == LongModular(1L));

        CPPUNIT_ASSERT_THROW_MESSAGE("", AlgebraicNumber_in_O_pO_1(1LL, 2L), std::string);
        AlgebraicNumber_in_O_pO_1::set_basis(113L);
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", AlgebraicNumber_in_O_pO_1(1LL, 2L));
    }

    void testIdeal()
    {
        AlgebraicNumber::clearNumberField();
        Matrix<VeryLong> mvl;
        std::vector<AlgebraicNumber> van;
        CPPUNIT_ASSERT_THROW_MESSAGE("", Ideal(), std::string);
        CPPUNIT_ASSERT_THROW_MESSAGE("", Ideal(VeryLong(1L), VeryLong(2L)), std::string);
        CPPUNIT_ASSERT_THROW_MESSAGE("", Ideal(VeryLong(1L), AlgebraicNumber(2L)), std::string);
        CPPUNIT_ASSERT_THROW_MESSAGE("", Ideal(AlgebraicNumber(2L)), std::string);
        CPPUNIT_ASSERT_THROW_MESSAGE("", (Ideal(mvl)), std::string);
        CPPUNIT_ASSERT_THROW_MESSAGE("", (Ideal(van)), std::string);
        CPPUNIT_ASSERT_THROW_MESSAGE("", Ideal(mvl, VeryLong(1L)), std::string);

        NumberField nf(f, "nf.fb.dat");
        AlgebraicNumber::setNumberField(nf);
        CPPUNIT_ASSERT(AlgebraicNumber::index() == nf.index());
        CPPUNIT_ASSERT(AlgebraicNumber::c_d() == nf.c_d());
        CPPUNIT_ASSERT(AlgebraicNumber::degree() == 5);
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", Ideal());
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", Ideal(VeryLong(1L), VeryLong(2L)));
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", Ideal(VeryLong(1L), AlgebraicNumber(2L)));
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", Ideal(AlgebraicNumber(2L)));
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", Ideal(Matrix<VeryLong>(5, 5)));
        std::vector<AlgebraicNumber> van1;
        van1.push_back(AlgebraicNumber(2L));
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", (Ideal(van1)));
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", Ideal(Matrix<VeryLong>(5, 5), VeryLong(1L)));

        Ideal i1(VeryLong(23L), VeryLong(19L));
        // denominator = 210, HNF basis =
        // 210 0 0 0 30
        // 0 282636900 8496180 194820108 163651325
        // 0 0 565273800 171281376 550330894
        // 0 0 0 113054760 44454828
        // 0 0 0 0 2691780
        Matrix<VeryLong> hnf1(5, 5);
        hnf1(0,0) = VeryLong("210");
        hnf1(0,4) = VeryLong("30");
        hnf1(1,1) = VeryLong("282636900");
        hnf1(1,2) = VeryLong("8496180");
        hnf1(1,3) = VeryLong("194820108");
        hnf1(1,4) = VeryLong("163651325");
        hnf1(2,2) = VeryLong("565273800");
        hnf1(2,3) = VeryLong("171281376");
        hnf1(2,4) = VeryLong("550330894");
        hnf1(3,3) = VeryLong("113054760");
        hnf1(3,4) = VeryLong("44454828");
        hnf1(4,4) = VeryLong("2691780");

        Ideal i2(hnf1, VeryLong("210"));
        CPPUNIT_ASSERT(i1 == i2);

        Ideal i3(VeryLong(23L), AlgebraicNumber::read_algebraic_number("1 + 2 alpha - 3 alpha^3"));
        // 210 0 0 0 30
        // 0 105 63 0 5
        // 0 0 126 0 64
        // 0 0 0 630 138
        // 0 0 0 0 2691780
        Matrix<VeryLong> hnf2(5, 5);
        hnf2(0,0) = VeryLong("210");
        hnf2(0,4) = VeryLong("30");
        hnf2(1,1) = VeryLong("105");
        hnf2(1,2) = VeryLong("63");
        hnf2(1,4) = VeryLong("5");
        hnf2(2,2) = VeryLong("126");
        hnf2(2,4) = VeryLong("64");
        hnf2(3,3) = VeryLong("630");
        hnf2(3,4) = VeryLong("138");
        hnf2(4,4) = VeryLong("2691780");

        Ideal i4(hnf2, VeryLong("210"));
        CPPUNIT_ASSERT(i3 == i4);

        AlgebraicNumber an5(AlgebraicNumber::read_algebraic_number("1 + 2 alpha - 3 alpha^3"));
        Ideal i5(an5);

        // 4711517646197236099128792239074896465353634419503848839621552048694300863215344347756895299738190 4465738774691251984589339740974944354690036501222430362951667300992990535220445877275401845386220 2888456736368779725166881730379466308254293175914232339218528381566960882341781748241207999961700 3728402160173299640970982246675088022699242746378174932942013057889059551235750465830921482330100 2289464754271834762134198923494447536616817151415578629526368972205424951074558203522629853749050
        // 0 105 63 0 5
        // 0 0 126 0 64
        // 0 0 0 630 138
        // 0 0 0 0 2691780
        Matrix<VeryLong> hnf3(5, 5);
        hnf3(0,0) = VeryLong("4711517646197236099128792239074896465353634419503848839621552048694300863215344347756895299738190");
        hnf3(0,1) = VeryLong("4465738774691251984589339740974944354690036501222430362951667300992990535220445877275401845386220");
        hnf3(0,2) = VeryLong("2888456736368779725166881730379466308254293175914232339218528381566960882341781748241207999961700");
        hnf3(0,3) = VeryLong("3728402160173299640970982246675088022699242746378174932942013057889059551235750465830921482330100");
        hnf3(0,4) = VeryLong("2289464754271834762134198923494447536616817151415578629526368972205424951074558203522629853749050");
        hnf3(1,1) = VeryLong("105");
        hnf3(1,2) = VeryLong("63");
        hnf3(1,4) = VeryLong("5");
        hnf3(2,2) = VeryLong("126");
        hnf3(2,4) = VeryLong("64");
        hnf3(3,3) = VeryLong("630");
        hnf3(3,4) = VeryLong("138");
        hnf3(4,4) = VeryLong("2691780");

        Ideal i6(hnf3, VeryLong("210"));
        CPPUNIT_ASSERT(i5 == i6);

        Ideal i7(hnf3);
        Ideal i8(hnf3, VeryLong(1L));
        CPPUNIT_ASSERT(i7 == i8);

        Ideal i9(23L, 29L);
        i9 = i7;
        CPPUNIT_ASSERT(i9 == i8);

        Ideal i10 = i5 * i6;
        Matrix<VeryLong> hnf4(5, 5);
        // Ideal I : denominator = 603806630700, HNF basis =
        // 303935152462438283884540564339184716664921807694174140512188807845388505307046562141018653873102309779281162019395599446498362497700262652216962928765073716694394495175696990197831509031753516935094700 155391732017519170042973412432766172592456137322624869785467928571087558363957911327990881250657194384788565048661966863617071967041912012476487387575940786199386993240418204610693715501602770628099900 90015451235931039872765129033586300821440472080934164618839805161854434983572425206621071858640966916016466531366896793659048127285301707115850612975169741093519103078810561671915904245671029678570200 174638139897914539955401703756827879739668877493073244076682483348002948090684525136881874958403473805872076252733919592537686677093800773196303054420313008560793469629603339622921658600726899269011400 234308281863861870519707739201626744047018070257406888849878552945650794890710408030570187350861357014229430738603948663286618523754095418459720351836894554402676496369820478384184619642507201507048400
        // 0 301903315350 241522652280 103509708120 60980886320
        // 0 0 30190331535 8625809010 9331744129
        // 0 0 0 12938713515 1579389231
        // 0 0 0 0 6
        hnf4(0,0) = VeryLong("303935152462438283884540564339184716664921807694174140512188807845388505307046562141018653873102309779281162019395599446498362497700262652216962928765073716694394495175696990197831509031753516935094700");
        hnf4(0,1) = VeryLong("155391732017519170042973412432766172592456137322624869785467928571087558363957911327990881250657194384788565048661966863617071967041912012476487387575940786199386993240418204610693715501602770628099900");
        hnf4(0,2) = VeryLong("90015451235931039872765129033586300821440472080934164618839805161854434983572425206621071858640966916016466531366896793659048127285301707115850612975169741093519103078810561671915904245671029678570200");
        hnf4(0,3) = VeryLong("174638139897914539955401703756827879739668877493073244076682483348002948090684525136881874958403473805872076252733919592537686677093800773196303054420313008560793469629603339622921658600726899269011400");
        hnf4(0,4) = VeryLong("234308281863861870519707739201626744047018070257406888849878552945650794890710408030570187350861357014229430738603948663286618523754095418459720351836894554402676496369820478384184619642507201507048400");
        hnf4(1,1) = VeryLong("301903315350");
        hnf4(1,2) = VeryLong("241522652280");
        hnf4(1,3) = VeryLong("103509708120");
        hnf4(1,4) = VeryLong("60980886320");
        hnf4(2,2) = VeryLong("30190331535");
        hnf4(2,3) = VeryLong("8625809010");
        hnf4(2,4) = VeryLong("9331744129");
        hnf4(3,3) = VeryLong("12938713515");
        hnf4(3,4) = VeryLong("1579389231");
        hnf4(4,4) = 6L;
        Ideal i11(hnf4, VeryLong("603806630700"));
        CPPUNIT_ASSERT(i10 == i11);

        {
            Ideal I(AlgebraicNumber(VeryLong(1L)));
            Ideal J = I.invert();
            CPPUNIT_ASSERT(I == J);
        }

        {
            VeryLong lcm;
            Matrix<VeryLong> AA;
            Ideal::integralPart(i6, lcm, AA);
            Ideal J1(i6);
            J1 *= lcm;
            Ideal J1inv(J1.invert());
            J1inv *= lcm;
            Ideal i14(AlgebraicNumber(VeryLong(1L)));
            CPPUNIT_ASSERT(i6 * J1inv == i14);
            CPPUNIT_ASSERT(J1inv * i6 == i14);
            Ideal i15 = i10 * J1inv;
            Ideal i16 = i15 * i6;
            Ideal i17 = i10 * i14;

            CPPUNIT_ASSERT(i17 == i10);
            CPPUNIT_ASSERT(i16 == i10);
        }

        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", (i10 / i6));
        CPPUNIT_ASSERT(i10 / i6 == i5);
        CPPUNIT_ASSERT(i10 == (i10 / i6) * i6);

        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", i7.invert());

        CPPUNIT_ASSERT(!i7.isIntegral());

        VeryLong lcm7;
        Matrix<VeryLong> AA7;
        Ideal::integralPart(i7, lcm7, AA7);

        Ideal i12 = i7;
        i12 *= lcm7;
        CPPUNIT_ASSERT(i12.isIntegral());
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", i12.invert());
        Ideal i13 = i12.invert();
        Ideal i14(AlgebraicNumber(VeryLong(1L)));
        CPPUNIT_ASSERT(i14 == i13 * i12);

        Ideal i30 = intersection(i5, i5);
        CPPUNIT_ASSERT(i30 == i5);
        AlgebraicNumber an6(AlgebraicNumber::read_algebraic_number("1 + alpha"));
        Ideal i35(an6);
        CPPUNIT_ASSERT(i35.isPrincipal());
        Ideal i40 = intersection(i5, i6);
        // denominator = 210, HNF basis =
        // 4711517646197236099128792239074896465353634419503848839621552048694300863215344347756895299738190 4465738774691251984589339740974944354690036501222430362951667300992990535220445877275401845386220 2888456736368779725166881730379466308254293175914232339218528381566960882341781748241207999961700 3728402160173299640970982246675088022699242746378174932942013057889059551235750465830921482330100 2289464754271834762134198923494447536616817151415578629526368972205424951074558203522629853749050
        // 0 105 63 0 5
        // 0 0 126 0 64
        // 0 0 0 630 138
        // 0 0 0 0 2691780

        Matrix<VeryLong> hnf40(5, 5);
        hnf40(0,0) = VeryLong("4711517646197236099128792239074896465353634419503848839621552048694300863215344347756895299738190");
        hnf40(0,1) = VeryLong("4465738774691251984589339740974944354690036501222430362951667300992990535220445877275401845386220");
        hnf40(0,2) = VeryLong("2888456736368779725166881730379466308254293175914232339218528381566960882341781748241207999961700");
        hnf40(0,3) = VeryLong("3728402160173299640970982246675088022699242746378174932942013057889059551235750465830921482330100");
        hnf40(0,4) = VeryLong("2289464754271834762134198923494447536616817151415578629526368972205424951074558203522629853749050");
        hnf40(1,1) = 105L;
        hnf40(1,2) = 63L;
        hnf40(1,4) = 5L;
        hnf40(2,2) = 126L;
        hnf40(2,4) = 64L;
        hnf40(3,3) = 630L;
        hnf40(3,4) = 138L;
        hnf40(4,4) = 2691780L;
        Ideal i41(hnf40, 210L);
        CPPUNIT_ASSERT(i40 == i41);

        CPPUNIT_ASSERT_THROW_MESSAGE("", i41.reducedBasisOmega(), std::string);

        VeryLong lcm42;
        Matrix<VeryLong> AA42;
        Ideal::integralPart(i41, lcm42, AA42);
        Ideal i42(i41);
        i42 *= lcm42;

        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", i42.reducedBasisOmega());

        Matrix<VeryLong> m43 = i42.reducedBasisOmega();

        Matrix<Quotient<VeryLong> > m44 = nf.w() * m43;
        //Matrix<Quotient<VeryLong> > m44 = m43 * nf.w();

        VeryLong lcm44(1L);
        for (size_t i = 0; i < m44.rows(); ++i)
        {
            for (size_t j = 0; j < m44.columns(); ++j)
            {
                lcm44 = lcm44 * m44(i,j).denominator() / gcd(lcm44, m44(i,j).denominator());
            }
        }
        Matrix<VeryLong> m45(m44.rows(), m44.columns());
        for (size_t i = 0; i < m44.rows(); ++i)
        {
            for (size_t j = 0; j < m44.columns(); ++j)
            {
                m45(i,j) = m44(i,j).numerator() * (lcm44 / m44(i,j).denominator());
            }
        }

        Ideal i45(HNF1(m45), lcm44);
        CPPUNIT_ASSERT(i42 == i45);

        CPPUNIT_ASSERT(i45.norm() == Quotient<VeryLong>(VeryLong("471896026035834322178876176087679575306370954866538858653814155957908622510215097006005310090297773290789903620950041511769639913259859455513867891898746093750000")));


        Ideal i46(AlgebraicNumber(VeryLong(1L)));
        CPPUNIT_ASSERT(i46.norm() == Quotient<VeryLong>(VeryLong(1L)));

        CPPUNIT_ASSERT(!i45.isPrincipal());
        CPPUNIT_ASSERT(i46.isPrincipal());
        CPPUNIT_ASSERT(Ideal(an5).isPrincipal());

        Ideal_mod_pO::set_basis(VeryLong(5L));
        Matrix<VeryLongModular> beta = AlgebraicNumber_in_O_pO::make_beta(VeryLong(5L));
        CPPUNIT_ASSERT(beta.rows() == 5);
        CPPUNIT_ASSERT(beta.columns() == 1);
        CPPUNIT_ASSERT(beta(0,0) == VeryLongModular(0L));
        CPPUNIT_ASSERT(beta(0,1) == VeryLongModular(3L));
        CPPUNIT_ASSERT(beta(0,2) == VeryLongModular(1L));
        CPPUNIT_ASSERT(beta(0,3) == VeryLongModular(0L));
        CPPUNIT_ASSERT(beta(0,4) == VeryLongModular(0L));

        Ideal_mod_pO i50(beta);
        CPPUNIT_ASSERT(i50.basis().rows() == 5);
        CPPUNIT_ASSERT(i50.basis().columns() == 1);
        CPPUNIT_ASSERT(i50.rank() == 1);
        CPPUNIT_ASSERT(i50.basis()(0,0) == VeryLongModular(0L));
        CPPUNIT_ASSERT(i50.basis()(0,1) == VeryLongModular(3L));
        CPPUNIT_ASSERT(i50.basis()(0,2) == VeryLongModular(1L));
        CPPUNIT_ASSERT(i50.basis()(0,3) == VeryLongModular(0L));
        CPPUNIT_ASSERT(i50.basis()(0,4) == VeryLongModular(0L));
        Ideal i500 = i50.makeIdeal();
        // i500 = Ideal I : denominator = 210, HNF basis =
        // 1050 0 0 0 450
        // 0 282636900 8496180 117693660 163651325
        // 0 0 565273800 291133080 550330894
        // 0 0 0 565273800 44454828
        // 0 0 0 0 2691780
        Matrix<VeryLong> hnf500(5, 5);
        hnf500(0,0) = VeryLong("1050");
        hnf500(0,4) = VeryLong("450");
        hnf500(1,1) = VeryLong("282636900");
        hnf500(1,2) = VeryLong("8496180");
        hnf500(1,3) = VeryLong("117693660");
        hnf500(1,4) = VeryLong("163651325");
        hnf500(2,2) = VeryLong("565273800");
        hnf500(2,3) = VeryLong("291133080");
        hnf500(2,4) = VeryLong("550330894");
        hnf500(3,3) = VeryLong("565273800");
        hnf500(3,4) = VeryLong("44454828");
        hnf500(4,4) = VeryLong("2691780");
        CPPUNIT_ASSERT(i500 == Ideal(hnf500, VeryLong("210")));

        Matrix<VeryLongModular> m200(5,1);
        m200(1,0) = VeryLongModular(2L);
        m200(2,0) = VeryLongModular(4L);
        m200(3,0) = VeryLongModular(4L);
        m200(4,0) = VeryLongModular(1L);
        Ideal_mod_pO i200(m200);
        Ideal_mod_pO i201(i200);
        Ideal i2010 = i201.makeIdeal();

        // i2010 = Ideal I : denominator = 210, HNF basis =
        // 1050 0 0 0 30
        // 0 1413184500 42480900 974100540 86524877
        // 0 0 2826369000 856406880 670182598
        // 0 0 0 565273800 496673868
        // 0 0 0 0 2691780
        Matrix<VeryLong> hnf2010(5,5);
        hnf2010(0,0) = VeryLong("1050");
        hnf2010(0,4) = VeryLong("30");
        hnf2010(1,1) = VeryLong("1413184500");
        hnf2010(1,2) = VeryLong("42480900");
        hnf2010(1,3) = VeryLong("974100540");
        hnf2010(1,4) = VeryLong("86524877");
        hnf2010(2,2) = VeryLong("2826369000");
        hnf2010(2,3) = VeryLong("856406880");
        hnf2010(2,4) = VeryLong("670182598");
        hnf2010(3,3) = VeryLong("565273800");
        hnf2010(3,4) = VeryLong("496673868");
        hnf2010(4,4) = VeryLong("2691780");
        CPPUNIT_ASSERT(i2010 == Ideal(hnf2010, VeryLong("210")));

        Ideal_mod_pO i202 = i200 * i201;
        CPPUNIT_ASSERT(i202.basis().rows() == 5);
        CPPUNIT_ASSERT(i202.basis().columns() == 0);
        CPPUNIT_ASSERT(i202.rank() == 0);

        Ideal_mod_pO::set_basis(VeryLong(3L));
        Matrix<VeryLongModular> m100(5,2);
        m100(0,1) = VeryLongModular(2L);
        m100(1,0) = VeryLongModular(1L);
        m100(4,1) = VeryLongModular(1L);
        Ideal_mod_pO i100(m100);
        Ideal_mod_pO i101(m100);
        Ideal_mod_pO i102 = i100 * i101;
        CPPUNIT_ASSERT(i102.basis().rows() == 5);
        CPPUNIT_ASSERT(i102.basis().columns() == 0);
        CPPUNIT_ASSERT(i102.rank() == 0);

        Ideal_mod_pO i103 = i102 / i100;
        CPPUNIT_ASSERT(i103.basis().rows() == 5);
        CPPUNIT_ASSERT(i103.basis().columns() == 3);
        CPPUNIT_ASSERT(i103.rank() == 3);
        CPPUNIT_ASSERT(i103.basis()(0,2) == VeryLongModular(2L));
        CPPUNIT_ASSERT(i103.basis()(1,0) == VeryLongModular(1L));
        CPPUNIT_ASSERT(i103.basis()(2,1) == VeryLongModular(2L));
        CPPUNIT_ASSERT(i103.basis()(3,1) == VeryLongModular(1L));
        CPPUNIT_ASSERT(i103.basis()(4,2) == VeryLongModular(1L));

        Ideal_mod_pO i104 = i100 / i103;
        CPPUNIT_ASSERT(i104.basis().rows() == 5);
        CPPUNIT_ASSERT(i104.basis().columns() == 4);
        CPPUNIT_ASSERT(i104.rank() == 4);
        CPPUNIT_ASSERT(i104.basis()(0,2) == VeryLongModular(1L));
        CPPUNIT_ASSERT(i104.basis()(0,3) == VeryLongModular(2L));
        CPPUNIT_ASSERT(i104.basis()(1,0) == VeryLongModular(1L));
        CPPUNIT_ASSERT(i104.basis()(2,1) == VeryLongModular(1L));
        CPPUNIT_ASSERT(i104.basis()(3,2) == VeryLongModular(1L));
        CPPUNIT_ASSERT(i104.basis()(4,3) == VeryLongModular(1L));
    }

    void testPrimeIdeal()
    {
        AlgebraicNumber::clearNumberField();
        NumberField nf(f, "nf.fb.dat");
        AlgebraicNumber::setNumberField(nf);

        std::vector<std::pair<PrimeIdeal, int> > primeIdeals;
        PrimeIdeal::primeDecomposition(VeryLong(113L), primeIdeals);
        CPPUNIT_ASSERT(primeIdeals.size() == 2);

        const PrimeIdeal& pi1 = primeIdeals[0].first;
        int e1 = primeIdeals[0].second;
        CPPUNIT_ASSERT(pi1.norm() == VeryLong(12769L));
        CPPUNIT_ASSERT(pi1.ramificationIndex() == 1);
        CPPUNIT_ASSERT(pi1.residualIndex() == 2);
        CPPUNIT_ASSERT(e1 == 1);

        const PrimeIdeal& pi2 = primeIdeals[1].first;
        int e2 = primeIdeals[1].second;
        CPPUNIT_ASSERT(pi2.norm() == VeryLong(1442897L));
        CPPUNIT_ASSERT(pi2.ramificationIndex() == 1);
        CPPUNIT_ASSERT(pi2.residualIndex() == 3);
        CPPUNIT_ASSERT(e2 == 1);

        primeIdeals.clear();
        PrimeIdeal::primeDecomposition(VeryLong(29L), primeIdeals);
        CPPUNIT_ASSERT(primeIdeals.size() == 3);
        CPPUNIT_ASSERT(primeIdeals[0].second == 1);
        CPPUNIT_ASSERT(primeIdeals[0].first.norm() == VeryLong(29L));
        CPPUNIT_ASSERT(primeIdeals[0].first.ramificationIndex() == 1);
        CPPUNIT_ASSERT(primeIdeals[0].first.residualIndex() == 1);

        CPPUNIT_ASSERT(primeIdeals[1].second == 1);
        CPPUNIT_ASSERT(primeIdeals[1].first.norm() == VeryLong(29L));
        CPPUNIT_ASSERT(primeIdeals[1].first.ramificationIndex() == 1);
        CPPUNIT_ASSERT(primeIdeals[1].first.residualIndex() == 1);

        CPPUNIT_ASSERT(primeIdeals[2].second == 1);
        CPPUNIT_ASSERT(primeIdeals[2].first.norm() == VeryLong(24389L));
        CPPUNIT_ASSERT(primeIdeals[2].first.ramificationIndex() == 1);
        CPPUNIT_ASSERT(primeIdeals[2].first.residualIndex() == 3);

        primeIdeals.clear();
        PrimeIdeal::primeDecomposition(VeryLong(31L), primeIdeals);
        CPPUNIT_ASSERT(primeIdeals.size() == 4);
        CPPUNIT_ASSERT(primeIdeals[0].second == 1);
        CPPUNIT_ASSERT(primeIdeals[0].first.norm() == VeryLong(31L));
        CPPUNIT_ASSERT(primeIdeals[0].first.ramificationIndex() == 1);
        CPPUNIT_ASSERT(primeIdeals[0].first.residualIndex() == 1);

        CPPUNIT_ASSERT(primeIdeals[1].second == 1);
        CPPUNIT_ASSERT(primeIdeals[1].first.norm() == VeryLong(31L));
        CPPUNIT_ASSERT(primeIdeals[1].first.ramificationIndex() == 1);
        CPPUNIT_ASSERT(primeIdeals[1].first.residualIndex() == 1);

        CPPUNIT_ASSERT(primeIdeals[2].second == 1);
        CPPUNIT_ASSERT(primeIdeals[2].first.norm() == VeryLong(31L));
        CPPUNIT_ASSERT(primeIdeals[2].first.ramificationIndex() == 1);
        CPPUNIT_ASSERT(primeIdeals[2].first.residualIndex() == 1);

        CPPUNIT_ASSERT(primeIdeals[3].second == 1);
        CPPUNIT_ASSERT(primeIdeals[3].first.norm() == VeryLong(961L));
        CPPUNIT_ASSERT(primeIdeals[3].first.ramificationIndex() == 1);
        CPPUNIT_ASSERT(primeIdeals[3].first.residualIndex() == 2);

        primeIdeals.clear();
        PrimeIdeal::primeDecomposition(VeryLong(3L), primeIdeals);
        CPPUNIT_ASSERT(primeIdeals.size() == 3);
        CPPUNIT_ASSERT(primeIdeals[0].second == 1);
        CPPUNIT_ASSERT(primeIdeals[0].first.norm() == VeryLong(3L));
        CPPUNIT_ASSERT(primeIdeals[0].first.ramificationIndex() == 1);
        CPPUNIT_ASSERT(primeIdeals[0].first.residualIndex() == 1);

        CPPUNIT_ASSERT(primeIdeals[1].second == 2);
        CPPUNIT_ASSERT(primeIdeals[1].first.norm() == VeryLong(3L));
        CPPUNIT_ASSERT(primeIdeals[1].first.ramificationIndex() == 2);
        CPPUNIT_ASSERT(primeIdeals[1].first.residualIndex() == 1);

        CPPUNIT_ASSERT(primeIdeals[2].second == 2);
        CPPUNIT_ASSERT(primeIdeals[2].first.norm() == VeryLong(3L));
        CPPUNIT_ASSERT(primeIdeals[2].first.ramificationIndex() == 2);
        CPPUNIT_ASSERT(primeIdeals[2].first.residualIndex() == 1);

        primeIdeals.clear();
        PrimeIdeal::primeDecomposition(VeryLong(5L), primeIdeals);
        CPPUNIT_ASSERT(primeIdeals.size() == 3);
        CPPUNIT_ASSERT(primeIdeals[0].second == 1);
        CPPUNIT_ASSERT(primeIdeals[0].first.norm() == VeryLong(25L));
        CPPUNIT_ASSERT(primeIdeals[0].first.ramificationIndex() == 1);
        CPPUNIT_ASSERT(primeIdeals[0].first.residualIndex() == 2);

        CPPUNIT_ASSERT(primeIdeals[1].second == 1);
        CPPUNIT_ASSERT(primeIdeals[1].first.norm() == VeryLong(5L));
        CPPUNIT_ASSERT(primeIdeals[1].first.ramificationIndex() == 1);
        CPPUNIT_ASSERT(primeIdeals[1].first.residualIndex() == 1);

        CPPUNIT_ASSERT(primeIdeals[2].second == 2);
        CPPUNIT_ASSERT(primeIdeals[2].first.norm() == VeryLong(5L));
        CPPUNIT_ASSERT(primeIdeals[2].first.ramificationIndex() == 2);
        CPPUNIT_ASSERT(primeIdeals[2].first.residualIndex() == 1);

        primeIdeals.clear();
        PrimeIdeal::primeDecomposition(VeryLong(2L), primeIdeals);
        CPPUNIT_ASSERT(primeIdeals.size() == 2);
        CPPUNIT_ASSERT(primeIdeals[0].second == 1);
        CPPUNIT_ASSERT(primeIdeals[0].first.norm() == VeryLong(4L));
        CPPUNIT_ASSERT(primeIdeals[0].first.ramificationIndex() == 1);
        CPPUNIT_ASSERT(primeIdeals[0].first.residualIndex() == 2);

        CPPUNIT_ASSERT(primeIdeals[1].second == 3);
        CPPUNIT_ASSERT(primeIdeals[1].first.norm() == VeryLong(2L));
        CPPUNIT_ASSERT(primeIdeals[1].first.ramificationIndex() == 3);
        CPPUNIT_ASSERT(primeIdeals[1].first.residualIndex() == 1);

        primeIdeals.clear();
        PrimeIdeal::primeDecomposition(VeryLong(19L), primeIdeals);
        CPPUNIT_ASSERT(primeIdeals.size() == 4);
        CPPUNIT_ASSERT(primeIdeals[0].second == 1);
        CPPUNIT_ASSERT(primeIdeals[0].first.norm() == VeryLong(19L));
        CPPUNIT_ASSERT(primeIdeals[0].first.ramificationIndex() == 1);
        CPPUNIT_ASSERT(primeIdeals[0].first.residualIndex() == 1);

        CPPUNIT_ASSERT(primeIdeals[1].second == 1);
        CPPUNIT_ASSERT(primeIdeals[1].first.norm() == VeryLong(19L));
        CPPUNIT_ASSERT(primeIdeals[1].first.ramificationIndex() == 1);
        CPPUNIT_ASSERT(primeIdeals[1].first.residualIndex() == 1);

        CPPUNIT_ASSERT(primeIdeals[2].second == 1);
        CPPUNIT_ASSERT(primeIdeals[2].first.norm() == VeryLong(19L));
        CPPUNIT_ASSERT(primeIdeals[2].first.ramificationIndex() == 1);
        CPPUNIT_ASSERT(primeIdeals[2].first.residualIndex() == 1);

        CPPUNIT_ASSERT(primeIdeals[3].second == 1);
        CPPUNIT_ASSERT(primeIdeals[3].first.norm() == VeryLong(361L));
        CPPUNIT_ASSERT(primeIdeals[3].first.ramificationIndex() == 1);
        CPPUNIT_ASSERT(primeIdeals[3].first.residualIndex() == 2);

        primeIdeals.clear();
        PrimeIdeal::primeDecomposition(VeryLong(11L), primeIdeals);
        CPPUNIT_ASSERT(primeIdeals.size() == 4);
        CPPUNIT_ASSERT(primeIdeals[0].second == 1);
        CPPUNIT_ASSERT(primeIdeals[0].first.norm() == VeryLong(11L));
        CPPUNIT_ASSERT(primeIdeals[0].first.ramificationIndex() == 1);
        CPPUNIT_ASSERT(primeIdeals[0].first.residualIndex() == 1);

        CPPUNIT_ASSERT(primeIdeals[1].second == 1);
        CPPUNIT_ASSERT(primeIdeals[1].first.norm() == VeryLong(11L));
        CPPUNIT_ASSERT(primeIdeals[1].first.ramificationIndex() == 1);
        CPPUNIT_ASSERT(primeIdeals[1].first.residualIndex() == 1);

        CPPUNIT_ASSERT(primeIdeals[2].second == 1);
        CPPUNIT_ASSERT(primeIdeals[2].first.norm() == VeryLong(11L));
        CPPUNIT_ASSERT(primeIdeals[2].first.ramificationIndex() == 1);
        CPPUNIT_ASSERT(primeIdeals[2].first.residualIndex() == 1);

        CPPUNIT_ASSERT(primeIdeals[3].second == 1);
        CPPUNIT_ASSERT(primeIdeals[3].first.norm() == VeryLong(121L));
        CPPUNIT_ASSERT(primeIdeals[3].first.ramificationIndex() == 1);
        CPPUNIT_ASSERT(primeIdeals[3].first.residualIndex() == 2);

    }
};
int main()
{
    CppUnit::TextUi::TestRunner runner;
    runner.addTest(AlgebraicNumberTest::suite());
    runner.run();
    return 0;
}

