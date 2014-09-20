#include <iostream>
#include "Matrix.h"
#ifdef USING_CPPUNIT
#include <limits>
#include <complex>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/ui/text/TestRunner.h>
#endif
#include "UnitTest.h"
bool compare_double_poly(const Polynomial<double>& p1, const Polynomial<double>& p2)
{
   if (p1.deg() != p2.deg()) return false;

   for (size_t i = 0; i <= static_cast<size_t>(p1.deg()); ++i)
   {
      if (!UnitTest::compare_double(p1.coefficient(i), p2.coefficient(i))) return false;
   }
   return true;
}

#ifdef USING_CPPUNIT
class MatrixTest : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE(MatrixTest);
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
        Matrix<long int> m1;
        CPPUNIT_ASSERT(m1.rows() == 0 && m1.columns() == 0);
     
        Matrix<long int> m2(5UL, 4UL);
        CPPUNIT_ASSERT(m2.rows() == 5 && m2.columns() == 4);
        CPPUNIT_ASSERT(m2(1,2) == 0 && m2(4,3) == 0);
     
        Matrix<long int> m3(m2);
        CPPUNIT_ASSERT(m3.rows() == 5 && m3.columns() == 4);
        CPPUNIT_ASSERT(m3(1,2) == 0 && m3(4,3) == 0);
     
        Matrix<long int> m4(7UL, 1L, 0L);
        CPPUNIT_ASSERT(m4.rows() == 7 && m4.columns() == 7);
        CPPUNIT_ASSERT(m4(1,2) == 0 && m4(4,4) == 1);
     
        m1.set_size(4UL, 5UL);
        CPPUNIT_ASSERT(m1.rows() == 4 && m1.columns() == 5);
        CPPUNIT_ASSERT(m1(1,2) == 0 && m1(3,4) == 0);
    
        CPPUNIT_ASSERT_THROW_MESSAGE("", m1.at(3, 5), std::string);
     
        m1.add_column();
        CPPUNIT_ASSERT(m1.rows() == 4 && m1.columns() == 6);
     
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", m1.at(3, 5));
     
     
        m2 = m1;
        CPPUNIT_ASSERT(m2.rows() == 4 && m2.columns() == 6);
     
        for (size_t i = 0; i < m1.rows(); ++i)
        {
           for (size_t j = 0; j < m1.columns(); ++j)
           {
              m1(i,j) = i + j;
              m2(i,j) = i - j;
           }
        }
        for (size_t i = 0; i < m1.rows(); ++i)
        {
           for (size_t j = 0; j < m1.columns(); ++j)
           {
              CPPUNIT_ASSERT(m1(i,j) == static_cast<long int>(i + j));
              CPPUNIT_ASSERT(m2(i,j) == static_cast<long int>(i - j));
           }
        }
     
        Matrix<long int> m5 = m1 + m2;
        for (size_t i = 0; i < m1.rows(); ++i)
        {
           for (size_t j = 0; j < m1.columns(); ++j)
           {
              CPPUNIT_ASSERT(m5(i,j) == static_cast<long int>(i + i));
           }
        }
     
        Matrix<long int> m6 = m1 - m2;
        for (size_t i = 0; i < m1.rows(); ++i)
        {
           for (size_t j = 0; j < m1.columns(); ++j)
           {
              CPPUNIT_ASSERT(m6(i,j) == static_cast<long int>(j + j));
           }
        }
     
        Matrix<long int> m7 = 24L * m4;
        CPPUNIT_ASSERT(m7.rows() == 7 && m7.columns() == 7);
        CPPUNIT_ASSERT(m7(1,2) == 0 && m7(4,4) == 24);
     
        Matrix<long int> m8 = m7 * 2L;
        CPPUNIT_ASSERT(m8(1,2) == 0 && m8(4,4) == 48);
     
        m8 /= 6L;
        CPPUNIT_ASSERT(m8(1,2) == 0 && m8(4,4) == 8);
     
        Matrix<long int> m9 = m8 / 4L;
        CPPUNIT_ASSERT(m9(1,2) == 0 && m9(4,4) == 2);
     
        std::vector<long int> v;
        v.resize(5L);
        CPPUNIT_ASSERT_THROW_MESSAGE("", std::vector<long int> v1 = m1 * v, std::string);

        v.resize(6L);
        std::vector<long int> v2;
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", v2 = m1 * v);
        for (size_t i = 0; i < v2.size(); ++i)
        {
            CPPUNIT_ASSERT(v2[i] == 0);
        }

        Matrix<long int> m10(5L, 7L);
        CPPUNIT_ASSERT(m10.rows() == 5 && m10.columns() == 7);
        for (size_t i = 0; i < 5L; ++i)
        {
           for (size_t j = 0; j < 7L; ++j)
           {
              m10(i, j) = i * j;
           }
        }
        for (size_t i = 0; i < m10.rows(); ++i)
        {
           for (size_t j = 0; j < m10.columns(); ++j)
           {
              CPPUNIT_ASSERT(m10(i, j) == static_cast<long int>(i * j));
           }
        }
     
        Matrix<long int> m11(7L, 2L, 0L);
        CPPUNIT_ASSERT(m11.rows() == 7 && m11.columns() == 7);
        CPPUNIT_ASSERT(m11(1,1) == 2 && m11(1,2) == 0);
     
        Matrix<long int> m12 = m10 * m11;
        for (size_t i = 0; i < m10.rows(); ++i)
        {
           for (size_t j = 0; j < m10.columns(); ++j)
           {
              CPPUNIT_ASSERT(m12(i, j) == static_cast<long int>(2 * i * j));
           }
        }
     
        m12 *= 3L;
        for (size_t i = 0; i < m10.rows(); ++i)
        {
           for (size_t j = 0; j < m10.columns(); ++j)
           {
              CPPUNIT_ASSERT(m12(i, j) == static_cast<long int>(6 * i * j));
           }
        }
     
        m12 *= m11;
        for (size_t i = 0; i < m10.rows(); ++i)
        {
           for (size_t j = 0; j < m10.columns(); ++j)
           {
              CPPUNIT_ASSERT(m12(i, j) == static_cast<long int>(12 * i * j));
           }
        }
     
        Matrix<long int> m13(5L, 5L);
        for (size_t i = 0; i < 5L; ++i)
        {
           for (size_t j = 0; j < 5L; ++j)
           {
              m13.at(i, j) = i*i - j*j;
           }
        }
        for (size_t i = 0; i < 5L; ++i)
        {
           for (size_t j = 0; j < 5L; ++j)
           {
              CPPUNIT_ASSERT(m13.at(i, j) == static_cast<long int>(i*i - j*j));
           }
        }
        long int tr = m13.trace();
        CPPUNIT_ASSERT(tr == 0);
     
        Matrix<long int> m14(m13);
        CPPUNIT_ASSERT(m13 == m14);
     
        for (size_t i = 0; i < m14.rows(); ++i)
        {
           m14(i,i) = 1L;
        }
        CPPUNIT_ASSERT(m13 != m14);

        Matrix<double> m13d(m13.rows(), m13.columns());
        Matrix<double> m14d(m14.rows(), m14.columns());
        for (size_t i = 0; i < m14.rows(); ++i)
        {
           for (size_t j = 0; j < m14.rows(); ++j)
           {
              m13d(i, j) = m13(i, j);
              m14d(i, j) = m14(i, j);
           }
        }
        std::vector<double> b;
        std::vector<double> x;
        b.resize(m14d.rows());
        for (size_t i = 0; i < b.size(); ++i)
        {
           b[i] = i;
        }
        bool solved = ::solve(m14d, b, x);
        CPPUNIT_ASSERT(solved);
        std::vector<double> b1 = m14d * x;
        CPPUNIT_ASSERT(b1.size() == b.size());
        for (size_t i = 0; i < b.size(); ++i)
        {
           CPPUNIT_ASSERT(b1[i] == b1[i]);
        }
     
        Matrix<double> m14dinv(m14d.rows(), m14d.columns());
        CPPUNIT_ASSERT(::invert(m14d, m14dinv));
     
        Matrix<double> m15d = m14d * m14dinv;
        Matrix<double> m16d(m14d.rows(), 1.0, 0);
	const double epsilon = std::numeric_limits<float>::epsilon();
        for (size_t i = 0; i < m15d.rows(); ++i)
        {
           for (size_t j = 0; j < m15d.columns(); ++j)
           {
                CPPUNIT_ASSERT_DOUBLES_EQUAL(m15d(i,j), m16d(i,j), epsilon);
           }
        }
     
        double det = ::determinant(m14d);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(det, 871.0, epsilon);
     
        long int det1 = ::determinant_in_integral_domain(m14);
        CPPUNIT_ASSERT(det1 == 871);
     
        Matrix<double> m14d_adj;
        Polynomial<double> p1 = ::characteristic_polynomial(m14d, m14d_adj);
        // p1 = -871 + 2615 X - 2620 X^2 + 880 X^3 - 5 X^4 + 1 X^5
        std::vector<double> c;
        c.resize(6);
        c[0] = -871;
        c[1] = 2615;
        c[2] = -2620;
        c[3] = 880;
        c[4] = -5;
        c[5] = 1;
        Polynomial<double> p2(c);
        CPPUNIT_ASSERT(compare_double_poly(p1, p2));
     
        Polynomial<double> p3 = ::characteristic_polynomial(m14d);
        CPPUNIT_ASSERT(compare_double_poly(p3, p2));
     
        CPPUNIT_ASSERT_THROW_MESSAGE("", ::solve(m13d, b, x), std::string);
     
        Matrix<double> m17d = ::kernel(m13d);
        Matrix<double> m18d = m13d * m17d;
        for (size_t i = 0; i < m18d.rows(); ++i)
        {
           for (size_t j = 0; j < m18d.columns(); ++j)
           {
              CPPUNIT_ASSERT_DOUBLES_EQUAL(m18d(i, j), 0, epsilon);
           }
        }
     
        Matrix<double> m19d = ::image(m13d);
        Matrix<double> m20d(5L, 2L);
        m20d(0,0) = -1;
        m20d(1,1) = 1;
        m20d(2,0) = 3;
        m20d(2,1) = 4;
        m20d(3,0) = 8;
        m20d(3,1) = 9;
        m20d(4,0) = 15;
        m20d(4,1) = 16;
     
        CPPUNIT_ASSERT(m19d == m20d);
     
        std::vector<double> v3;
        v3.resize(m19d.rows());
        Matrix<double> m21d = ::inverse_image(m19d, v3);
        Matrix<double> m22d = m19d * m21d;
        CPPUNIT_ASSERT(m22d.rows() == v3.size() && m22d.columns() == 1);
        for (size_t i = 0; i < v3.size(); ++i)
        {
           CPPUNIT_ASSERT_DOUBLES_EQUAL(m22d(i, 0), v3[i], epsilon);
        }
     
        Matrix<double> m23d(m19d.rows(), 2L);
        m20d(0,0) = -1;
        m20d(1,0) = 1;
        m20d(2,0) = 7;
        m20d(3,0) = 17;
        m20d(4,0) = 31;
        m20d(0,1) = -1;
        m20d(1,1) = -1;
        m20d(2,1) = -1;
        m20d(3,1) = -1;
        m20d(4,1) = -1;
     
        Matrix<double> m24d = ::inverse_image_matrix(m19d, m23d);
        Matrix<double> m25d = m19d * m24d;
        for (size_t i = 0; i < m25d.rows(); ++i)
        {
           for (size_t j = 0; j < m25d.columns(); ++j)
           {
              CPPUNIT_ASSERT_DOUBLES_EQUAL(m23d(i, j), m25d(i, j), epsilon);
           }
        }
     
        Matrix<double> m26d = ::supplement(m19d);
        Matrix<double> m26dinv(m26d.rows(), m26d.columns());
        CPPUNIT_ASSERT(::invert(m26d, m26dinv));
        for (size_t i = 0; i < m19d.rows(); ++i)
        {
           for (size_t j = 0; j < m19d.columns(); ++j)
           {
              CPPUNIT_ASSERT_DOUBLES_EQUAL(m19d(i, j), m26d(i, j), epsilon);
           }
        }
     
        Matrix<double> m27d(m26d.rows(), m26d.columns() - 1);
        for (size_t i = 0; i < m27d.rows(); ++i)
        {
           for (size_t j = 0; j < m27d.columns(); ++j)
           {
              m27d(i, j) = m26d(i, j);
           }
        }
        Matrix<double> m28d = ::supplement_subspace(m19d, m27d);
        CPPUNIT_ASSERT(m28d.rows() == m19d.rows() && m28d.columns() == m27d.columns() - m19d.columns());
        Matrix<double> m29d(m19d.rows(), m27d.columns() - m19d.columns());
        m29d(2,0) = 1;
        m29d(3,1) = 1;
     
        for (size_t i = 0; i < m28d.rows(); ++i)
        {
           for (size_t j = 0; j < m28d.columns(); ++j)
           {
              CPPUNIT_ASSERT_DOUBLES_EQUAL(m28d(i, j), m29d(i, j), epsilon);
           }
        }
     
        Matrix<long int> m15 = HNF(m14);
     
        Matrix<long int> m16(m14.rows(), 1L, 0);
        m16(0,0) = 871;
        m16(0,1) = 720;
        m16(0,2) = 270;
        m16(0,3) = 391;
        m16(0,4) = 212;
     
        CPPUNIT_ASSERT(m15 == m16);
     
        Matrix<long int> m17 = HNF1(m14);
        CPPUNIT_ASSERT(m17 == m16);
     
        Matrix<long int> m18 = ::HNF_mod_D(m14, det1);
        CPPUNIT_ASSERT(m18 == m16);
     
        Matrix<long int> m19(3L, 4L);
        Matrix<long int> m20(3L, 4L);
        Matrix<long int> m21(3L, 4L);
        for (size_t i = 0; i < m19.rows(); ++i)
        {
           for (size_t j = 0; j < m19.columns(); ++j)
           {
              m19(i, j) = i + j;
              m20(i, j) = i + j;
              m21(i, j) = i + j;
              if (i == 0) m20(i, j) = 1 + j;
              if (i == 1) m20(i, j) = j;
              if (i == 0) m21(i, j) = 1 + j;
              if (i == 1) m21(i, j) = 2 + j;
              if (i == 2) m21(i, j) = j;
           }
        }
        m19.swap(0, 1);
        CPPUNIT_ASSERT(m19 == m20);
        m19.swap(1, 2);
        CPPUNIT_ASSERT(m19 == m21);
     
        m19.transform_row(0, 2, 2L, 3L, 4L, 5L);
     
        Matrix<long int> m22(3L, 4L);
        m22(0,0) = 2;
        m22(0,1) = 7;
        m22(0,2) = 12;
        m22(0,3) = 17;
        m22(1,0) = 2;
        m22(1,1) = 3;
        m22(1,2) = 4;
        m22(1,3) = 5;
        m22(2,0) = 4;
        m22(2,1) = 13;
        m22(2,2) = 22;
        m22(2,3) = 31;
     
        CPPUNIT_ASSERT(m19 == m22);
     
        long int x1 = m22.dot(0, 1);
        long int y1 = m22.dot(0, 2);
        long int z1 = m22.dot(1, 2);
        CPPUNIT_ASSERT(x1 == 158);
        CPPUNIT_ASSERT(y1 == 890);
        CPPUNIT_ASSERT(z1 == 290);
    }
};


#endif
int main()
{
#ifdef USING_CPPUNIT
   CppUnit::TextUi::TestRunner runner;
   runner.addTest(MatrixTest::suite());
   runner.run();
#else
   std::cout << "Testing Matrix class ..." << std::endl;

   UnitTest t;

   Matrix<long int> m1;
   t.check(m1.rows() == 0 && m1.columns() == 0, "m1 should be empty");

   Matrix<long int> m2(5UL, 4UL);
   t.check(m2.rows() == 5 && m2.columns() == 4, "m2 should be 5x4");
   t.check(m2(1,2) == 0 && m2(4,3) == 0, "m2 should have all zeros");

   Matrix<long int> m3(m2);
   t.check(m3.rows() == 5 && m3.columns() == 4, "m3 should be 5x4");
   t.check(m3(1,2) == 0 && m3(4,3) == 0, "m3 should have all zeros");

   Matrix<long int> m4(7UL, 1L, 0L);
   t.check(m4.rows() == 7 && m4.columns() == 7, "m4 should be 7x7");
   t.check(m4(1,2) == 0 && m4(4,4) == 1, "m4 should have ones on the diagonal and zeros everywhere else");

   m1.set_size(4UL, 5UL);
   t.check(m1.rows() == 4 && m1.columns() == 5, "m1 should be 4x5");
   t.check(m1(1,2) == 0 && m1(3,4) == 0, "m1 should have all zeros");

   try
   {
      m1.at(3, 5);
      t.check(false, "m1.at(3, 5) should throw");
   }
   catch (const char* c)
   {
      t.check(true, "m1.at(3, 5) should throw");
   }

   m1.add_column();
   t.check(m1.rows() == 4 && m1.columns() == 6, "m1 should be 4x6");

   try
   {
      m1.at(3, 5);
      t.check(true, "m1.at(3, 5) should not throw");
   }
   catch (const char* c)
   {
      t.check(false, "m1.at(3, 5) should not throw");
   }


   m2 = m1;
   t.check(m2.rows() == 4 && m2.columns() == 6, "m2 should be 4x6");

   for (size_t i = 0; i < m1.rows(); ++i)
   {
      for (size_t j = 0; j < m1.columns(); ++j)
      {
         m1(i,j) = i + j;
         m2(i,j) = i - j;
      }
   }
   for (size_t i = 0; i < m1.rows(); ++i)
   {
      for (size_t j = 0; j < m1.columns(); ++j)
      {
         t.check(m1(i,j) == static_cast<long int>(i + j), "m1(i,j) should be i + j");
         t.check(m2(i,j) == static_cast<long int>(i - j), "m2(i,j) should be i - j");
      }
   }

   Matrix<long int> m5 = m1 + m2;
   for (size_t i = 0; i < m1.rows(); ++i)
   {
      for (size_t j = 0; j < m1.columns(); ++j)
      {
         t.check(m5(i,j) == static_cast<long int>(i + i), "m2(i,j) should be 2i");
      }
   }

   Matrix<long int> m6 = m1 - m2;
   for (size_t i = 0; i < m1.rows(); ++i)
   {
      for (size_t j = 0; j < m1.columns(); ++j)
      {
         t.check(m6(i,j) == static_cast<long int>(j + j), "m2(i,j) should be 2j");
      }
   }

   Matrix<long int> m7 = 24L * m4;
   t.check(m7.rows() == 7 && m7.columns() == 7, "m7 should be 7x7");
   t.check(m7(1,2) == 0 && m7(4,4) == 24, "m7 should have 24 on the diagonal and zeros everywhere else");

   Matrix<long int> m8 = m7 * 2L;
   t.check(m8(1,2) == 0 && m8(4,4) == 48, "m8 should have 48 on the diagonal and zeros everywhere else");

   m8 /= 6L;
   t.check(m8(1,2) == 0 && m8(4,4) == 8, "m8 should have 8 on the diagonal and zeros everywhere else");

   Matrix<long int> m9 = m8 / 4L;
   t.check(m9(1,2) == 0 && m9(4,4) == 2, "m9 should have 8 on the diagonal and zeros everywhere else");

   std::vector<long int> v;
   v.resize(5L);
   try
   {
      std::vector<long int> v1 = m1 * v;
      t.check(false, "m1 * v should throw");
   }
   catch (const char* c)
   {
      t.check(true, "m1 * v should throw");
   }

   v.resize(6L);
   try
   {
      std::vector<long int> v2 = m1 * v;
      t.check(true, "m1 * v should not throw");
      for (size_t i = 0; i < v2.size(); ++i)
      {
         t.check(v2[i] == 0, "v2[i] should be 0");
      }
   }
   catch (const char* c)
   {
      t.check(false, "m1 * v should not throw");
   }

   Matrix<long int> m10(5L, 7L);
   t.check(m10.rows() == 5 && m10.columns() == 7, "m10 should be 5x7");
   for (size_t i = 0; i < 5L; ++i)
   {
      for (size_t j = 0; j < 7L; ++j)
      {
         m10(i, j) = i * j;
      }
   }
   for (size_t i = 0; i < m10.rows(); ++i)
   {
      for (size_t j = 0; j < m10.columns(); ++j)
      {
         t.check(m10(i, j) == static_cast<long int>(i * j), "m10(i, j) should be i*j");
      }
   }

   Matrix<long int> m11(7L, 2L, 0L);
   t.check(m11.rows() == 7 && m11.columns() == 7, "m11 should be 7x7");
   t.check(m11(1,1) == 2 && m11(1,2) == 0, "m11 should have 2 on the diagonal and zero everywhere else");

   Matrix<long int> m12 = m10 * m11;
   for (size_t i = 0; i < m10.rows(); ++i)
   {
      for (size_t j = 0; j < m10.columns(); ++j)
      {
         t.check(m12(i, j) == static_cast<long int>(2 * i * j), "m12(i, j) should be 2*i*j");
      }
   }

   m12 *= 3L;
   for (size_t i = 0; i < m10.rows(); ++i)
   {
      for (size_t j = 0; j < m10.columns(); ++j)
      {
         t.check(m12(i, j) == static_cast<long int>(6 * i * j), "m12(i, j) should be 6*i*j");
      }
   }

   m12 *= m11;
   for (size_t i = 0; i < m10.rows(); ++i)
   {
      for (size_t j = 0; j < m10.columns(); ++j)
      {
         t.check(m12(i, j) == static_cast<long int>(12 * i * j), "m12(i, j) should be 12*i*j");
      }
   }

   Matrix<long int> m13(5L, 5L);
   for (size_t i = 0; i < 5L; ++i)
   {
      for (size_t j = 0; j < 5L; ++j)
      {
         m13.at(i, j) = i*i - j*j;
      }
   }
   for (size_t i = 0; i < 5L; ++i)
   {
      for (size_t j = 0; j < 5L; ++j)
      {
         t.check(m13.at(i, j) == static_cast<long int>(i*i - j*j), "m13(i,j) should be i^2 - j^2");
      }
   }
   long int tr = m13.trace();
   t.check(tr == 0, "trace(m13) should be zero");

   Matrix<long int> m14(m13);
   t.check(m13 == m14, "m13 should be identical to m14");

   for (size_t i = 0; i < m14.rows(); ++i)
   {
      m14(i,i) = 1L;
   }
   t.check(m13 != m14, "m13 should not be identical to m14");

   Matrix<double> m13d(m13.rows(), m13.columns());
   Matrix<double> m14d(m14.rows(), m14.columns());
   for (size_t i = 0; i < m14.rows(); ++i)
   {
      for (size_t j = 0; j < m14.rows(); ++j)
      {
         m13d(i, j) = m13(i, j);
         m14d(i, j) = m14(i, j);
      }
   }
   std::vector<double> b;
   std::vector<double> x;
   b.resize(m14d.rows());
   for (size_t i = 0; i < b.size(); ++i)
   {
      b[i] = i;
   }
   bool solved = ::solve(m14d, b, x);
   t.check(solved, "m14d should be invertible");
   std::vector<double> b1 = m14d * x;
   t.check(b1.size() == b.size(), "b1 should be same size as b");
   for (size_t i = 0; i < b.size(); ++i)
   {
      t.check(b1[i] == b1[i], "b1 should be same as b");
   }

   Matrix<double> m14dinv(m14d.rows(), m14d.columns());
   t.check(::invert(m14d, m14dinv), "m14d should be invertible");

   Matrix<double> m15d = m14d * m14dinv;
   Matrix<double> m16d(m14d.rows(), 1.0, 0);
   for (size_t i = 0; i < m15d.rows(); ++i)
   {
      for (size_t j = 0; j < m15d.columns(); ++j)
      {
         t.check(UnitTest::compare_double(m15d(i,j), m16d(i,j)), "m15d should be identity");
      }
   }

   double det = ::determinant(m14d);
   t.check(UnitTest::compare_double(det, 871.0), "determinant of m14d should be 871");

   long int det1 = ::determinant_in_integral_domain(m14);
   t.check(det1 == 871, "determinant of m14 should be 871");

   Matrix<double> m14d_adj;
   Polynomial<double> p1 = ::characteristic_polynomial(m14d, m14d_adj);
   // p1 = -871 + 2615 X - 2620 X^2 + 880 X^3 - 5 X^4 + 1 X^5
   std::vector<double> c;
   c.resize(6);
   c[0] = -871;
   c[1] = 2615;
   c[2] = -2620;
   c[3] = 880;
   c[4] = -5;
   c[5] = 1;
   Polynomial<double> p2(c);
   t.check(compare_double_poly(p1, p2), "characteristic polynomial of m14d should be -871 + 2615 X - 2620 X^2 + 880 X^3 - 5 X^4 + 1 X^5");

   Polynomial<double> p3 = ::characteristic_polynomial(m14d);
   t.check(compare_double_poly(p3, p2), "characteristic polynomial of m14d should be -871 + 2615 X - 2620 X^2 + 880 X^3 - 5 X^4 + 1 X^5");

   t.check(!::solve(m13d, b, x), "m13d should not be invertible");

   Matrix<double> m17d = ::kernel(m13d);
   Matrix<double> m18d = m13d * m17d;
   for (size_t i = 0; i < m18d.rows(); ++i)
   {
      for (size_t j = 0; j < m18d.columns(); ++j)
      {
         t.check(UnitTest::compare_double(m18d(i, j), 0), "m18d should be zero matrix");
      }
   }

   Matrix<double> m19d = ::image(m13d);
   Matrix<double> m20d(5L, 2L);
   m20d(0,0) = -1;
   m20d(1,1) = 1;
   m20d(2,0) = 3;
   m20d(2,1) = 4;
   m20d(3,0) = 8;
   m20d(3,1) = 9;
   m20d(4,0) = 15;
   m20d(4,1) = 16;

   t.check(m19d == m20d, "image(m13d) should be m20d");

   std::vector<double> v3;
   v3.resize(m19d.rows());
   Matrix<double> m21d = ::inverse_image(m19d, v3);
   Matrix<double> m22d = m19d * m21d;
   t.check(m22d.rows() == v3.size() && m22d.columns() == 1, "inverse_image of m19 should be a column vector");
   for (size_t i = 0; i < v3.size(); ++i)
   {
      t.check(UnitTest::compare_double(m22d(i, 0), v3[i]), "m19d * m21d should be v3");
   }

   Matrix<double> m23d(m19d.rows(), 2L);
   m20d(0,0) = -1;
   m20d(1,0) = 1;
   m20d(2,0) = 7;
   m20d(3,0) = 17;
   m20d(4,0) = 31;
   m20d(0,1) = -1;
   m20d(1,1) = -1;
   m20d(2,1) = -1;
   m20d(3,1) = -1;
   m20d(4,1) = -1;

   Matrix<double> m24d = ::inverse_image_matrix(m19d, m23d);
   Matrix<double> m25d = m19d * m24d;
   for (size_t i = 0; i < m25d.rows(); ++i)
   {
      for (size_t j = 0; j < m25d.columns(); ++j)
      {
         t.check(UnitTest::compare_double(m23d(i, j), m25d(i, j)), "m23d should be same as m25d");
      }
   }

   Matrix<double> m26d = ::supplement(m19d);
   Matrix<double> m26dinv(m26d.rows(), m26d.columns());
   t.check(::invert(m26d, m26dinv), "supplement(m19d) should be invertible");
   for (size_t i = 0; i < m19d.rows(); ++i)
   {
      for (size_t j = 0; j < m19d.columns(); ++j)
      {
         t.check(UnitTest::compare_double(m19d(i, j), m26d(i, j)), "first columns of supplement(m19d) should be m19d");
      }
   }

   Matrix<double> m27d(m26d.rows(), m26d.columns() - 1);
   for (size_t i = 0; i < m27d.rows(); ++i)
   {
      for (size_t j = 0; j < m27d.columns(); ++j)
      {
         m27d(i, j) = m26d(i, j);
      }
   }
   Matrix<double> m28d = ::supplement_subspace(m19d, m27d);
   t.check(m28d.rows() == m19d.rows() && m28d.columns() == m27d.columns() - m19d.columns(), "m28d should be 5x2");
   Matrix<double> m29d(m19d.rows(), m27d.columns() - m19d.columns());
   m29d(2,0) = 1;
   m29d(3,1) = 1;

   for (size_t i = 0; i < m28d.rows(); ++i)
   {
      for (size_t j = 0; j < m28d.columns(); ++j)
      {
         t.check(UnitTest::compare_double(m28d(i, j), m29d(i, j)), "m28d should be same as m29d");
      }
   }

   Matrix<long int> m15 = HNF(m14);

   Matrix<long int> m16(m14.rows(), 1L, 0);
   m16(0,0) = 871;
   m16(0,1) = 720;
   m16(0,2) = 270;
   m16(0,3) = 391;
   m16(0,4) = 212;

   t.check(m15 == m16, "HNF(m14) should equal m16");

   Matrix<long int> m17 = HNF1(m14);
   t.check(m17 == m16, "HNF1(m14) should equal m16");

   Matrix<long int> m18 = ::HNF_mod_D(m14, det1);
   t.check(m18 == m16, "HNF_mod_D(m14, 1000L) should equal m16");

   Matrix<long int> m19(3L, 4L);
   Matrix<long int> m20(3L, 4L);
   Matrix<long int> m21(3L, 4L);
   for (size_t i = 0; i < m19.rows(); ++i)
   {
      for (size_t j = 0; j < m19.columns(); ++j)
      {
         m19(i, j) = i + j;
         m20(i, j) = i + j;
         m21(i, j) = i + j;
         if (i == 0) m20(i, j) = 1 + j;
         if (i == 1) m20(i, j) = j;
         if (i == 0) m21(i, j) = 1 + j;
         if (i == 1) m21(i, j) = 2 + j;
         if (i == 2) m21(i, j) = j;
      }
   }
   m19.swap(0, 1);
   t.check(m19 == m20, "m19 should be equal m20");
   m19.swap(1, 2);
   t.check(m19 == m21, "m19 should be equal m21");

   m19.transform_row(0, 2, 2L, 3L, 4L, 5L);

   Matrix<long int> m22(3L, 4L);
   m22(0,0) = 2;
   m22(0,1) = 7;
   m22(0,2) = 12;
   m22(0,3) = 17;
   m22(1,0) = 2;
   m22(1,1) = 3;
   m22(1,2) = 4;
   m22(1,3) = 5;
   m22(2,0) = 4;
   m22(2,1) = 13;
   m22(2,2) = 22;
   m22(2,3) = 31;

   t.check(m19 == m22, "m19 should be equal m22");

   long int x1 = m22.dot(0, 1);
   long int y1 = m22.dot(0, 2);
   long int z1 = m22.dot(1, 2);
   t.check(x1 == 158, "m22.dot(0,1) should be 158");
   t.check(y1 == 890, "m22.dot(0,2) should be 890");
   t.check(z1 == 290, "m22.dot(1,2) should be 290");

   t.test_summary();
#endif
   return 0;
}
