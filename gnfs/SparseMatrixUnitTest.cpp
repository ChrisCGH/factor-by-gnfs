#include <iostream>
#include "SparseMatrix.h"
#include <limits>
#include <complex>
#include <sstream>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/ui/text/TestRunner.h>
#include "UnitTest.h"

class SparseMatrixTest : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE(SparseMatrixTest);
    CPPUNIT_TEST(testFileBasedSparseRowManager);
    CPPUNIT_TEST(testSparseMatrix);
    CPPUNIT_TEST_SUITE_END();

    public:
    void setUp()
    {
    }

    void tearDown()
    {
    }

    void check_row_operations(FileBasedSparseRow& row)
    {
        
        CPPUNIT_ASSERT(row.highest_column() == -1);
        CPPUNIT_ASSERT(row.xor(2) == ISparseRow::XOR_ADDED);
        CPPUNIT_ASSERT(row.highest_column() == 2);
        CPPUNIT_ASSERT(row.xor(9) == ISparseRow::XOR_ADDED);
        CPPUNIT_ASSERT(row.highest_column() == 9);
        CPPUNIT_ASSERT(row.xor(3) == ISparseRow::XOR_ADDED);
        CPPUNIT_ASSERT(row.highest_column() == 9);
        CPPUNIT_ASSERT(row.xor(10) == ISparseRow::XOR_ADDED);
        CPPUNIT_ASSERT(row.highest_column() == 10);
        CPPUNIT_ASSERT(row.xor(4) == ISparseRow::XOR_ADDED);
        CPPUNIT_ASSERT(row.highest_column() == 10);
        CPPUNIT_ASSERT(row.size() == 5);
        FileBasedSparseRow::const_iterator it = row.begin();
        CPPUNIT_ASSERT(*it == 2); ++it;
        CPPUNIT_ASSERT(*it == 3); ++it;
        CPPUNIT_ASSERT(*it == 4); ++it;
        CPPUNIT_ASSERT(*it == 9); ++it;
        CPPUNIT_ASSERT(*it == 10); ++it;
        CPPUNIT_ASSERT(it == row.end());
        CPPUNIT_ASSERT(row.xor(4) == ISparseRow::XOR_REMOVED);
        CPPUNIT_ASSERT(row.size() == 4);
        CPPUNIT_ASSERT(row.highest_column() == 10);
        CPPUNIT_ASSERT(row.xor(10) == ISparseRow::XOR_REMOVED);
        CPPUNIT_ASSERT(row.size() == 3);
        CPPUNIT_ASSERT(row.highest_column() == 9);
 
        for (size_t i = 0; i < 17; ++i)
        {
            CPPUNIT_ASSERT(row.xor(i + 1000) == ISparseRow::XOR_ADDED);
        }        
        CPPUNIT_ASSERT(row.size() == 20);
        CPPUNIT_ASSERT(row.highest_column() == 1016);
        static size_t expected_row_data[] = { 2, 3, 9, 1000, 1001, 1002, 1003, 1004, 1005, 1006, 1007, 1008, 1009, 1010, 1011, 1012, 1013, 1014, 1015, 1016 };
        size_t j = 0;
        for (FileBasedSparseRow::const_iterator it = row.begin();
            it != row.end();
            ++it, ++j)
        {
            CPPUNIT_ASSERT(*it == expected_row_data[j]);
        }

        CPPUNIT_ASSERT(row.xor(1017) == ISparseRow::XOR_ADDED);
        CPPUNIT_ASSERT(row.size() == 21);
        CPPUNIT_ASSERT(row.highest_column() == 1017);
        static size_t expected_row_data1[] = { 2, 3, 9, 1000, 1001, 1002, 1003, 1004, 1005, 1006, 1007, 1008, 1009, 1010, 1011, 1012, 1013, 1014, 1015, 1016, 1017 };
        j = 0;
        for (FileBasedSparseRow::const_iterator it = row.begin();
            it != row.end();
            ++it, ++j)
        {
            CPPUNIT_ASSERT(*it == expected_row_data1[j]);
        }

        SparseRow sr(row.size());
        row.copy(sr);
        j = 0;
        for (ISparseRow::const_iterator it = sr.begin();
            it != sr.end();
            ++it, ++j)
        {
            CPPUNIT_ASSERT(*it == expected_row_data1[j]);
        }

        CPPUNIT_ASSERT(row.add_next(1020) == ISparseRow::XOR_ADDED);
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", row.clear());
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", row.compress());
        CPPUNIT_ASSERT(row.memory_usage() == 0);
    }

    void testFileBasedSparseRowManager()
    {
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", FileBasedSparseRowManager("testfile", 10000, 100000));

        FileBasedSparseRowManager fbsrm1("testfile", 10000, 100000);
        CPPUNIT_ASSERT_THROW_MESSAGE("", fbsrm1.row(10), std::string);
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", fbsrm1.row(10000));
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", fbsrm1.row(50000));
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", fbsrm1.row(99999));
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", fbsrm1.row(100000));
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", fbsrm1.row(100001));
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", fbsrm1.row(110001));
        FileBasedSparseRow row1 = fbsrm1.row(10000);
        check_row_operations(row1);

        FileBasedSparseRow row2 = fbsrm1.row(110000);
        check_row_operations(row2);

        FileBasedSparseRow row3 = fbsrm1.row(210000);
        check_row_operations(row3);

#ifdef WIN32
        /***
        std::cout << "        CPPUNIT_ASSERT(fbsrm1.extend(210001, 20) == " << fbsrm1.extend(210001, 20) << "ull);" << std::endl;
        std::cout << "        CPPUNIT_ASSERT(fbsrm1.extend(210001, 50) == " << fbsrm1.extend(210001, 50) << "ull);" << std::endl;
        std::cout << "        CPPUNIT_ASSERT(fbsrm1.extend(210001, 1000) == " << fbsrm1.extend(210001, 1000) << "ull);" << std::endl;
        std::cout << "        CPPUNIT_ASSERT(fbsrm1.extend(110000, 50) == " << fbsrm1.extend(110000, 50) << "ull);" << std::endl;
        std::cout << "        CPPUNIT_ASSERT(fbsrm1.extend(109999, 50) == " << fbsrm1.extend(109999, 50) << "ull);" << std::endl;
        std::cout << "        CPPUNIT_ASSERT(fbsrm1.extend(109998, 50) == " << fbsrm1.extend(109998, 50) << "ull);" << std::endl;
        std::cout << "        CPPUNIT_ASSERT(fbsrm1.extend(109995, 50) == " << fbsrm1.extend(109995, 50) << "ull);" << std::endl;
        std::cout << "        CPPUNIT_ASSERT(fbsrm1.extend(109997, 50) == " << fbsrm1.extend(109997, 50) << "ull);" << std::endl;
        std::cout << "        CPPUNIT_ASSERT(fbsrm1.extend(109996, 50) == " << fbsrm1.extend(109996, 50) << "ull);" << std::endl;
        ***/
        CPPUNIT_ASSERT(fbsrm1.extend(210001, 20) == 17600400ull);
        CPPUNIT_ASSERT(fbsrm1.extend(210001, 50) == 17600776ull);
        CPPUNIT_ASSERT(fbsrm1.extend(210001, 1000) == 17600984ull);
        CPPUNIT_ASSERT(fbsrm1.extend(110000, 50) == 8800464ull);
        CPPUNIT_ASSERT(fbsrm1.extend(109999, 50) == 17600776ull);
        CPPUNIT_ASSERT(fbsrm1.extend(109998, 50) == 17604992ull);
        CPPUNIT_ASSERT(fbsrm1.extend(109995, 50) == 17605200ull);
        CPPUNIT_ASSERT(fbsrm1.extend(109997, 50) == 8799736ull);
        CPPUNIT_ASSERT(fbsrm1.extend(109996, 50) == 17605408ull);

#else
        /*std::cout << "        CPPUNIT_ASSERT(fbsrm1.extend(210001, 20) == " << fbsrm1.extend(210001, 20) << "ull);" << std::endl;
        std::cout << "        CPPUNIT_ASSERT(fbsrm1.extend(210001, 50) == " << fbsrm1.extend(210001, 50) << "ull);" << std::endl;
        std::cout << "        CPPUNIT_ASSERT(fbsrm1.extend(210001, 1000) == " << fbsrm1.extend(210001, 1000) << "ull);" << std::endl;
        std::cout << "        CPPUNIT_ASSERT(fbsrm1.extend(110000, 50) == " << fbsrm1.extend(110000, 50) << "ull);" << std::endl;
        std::cout << "        CPPUNIT_ASSERT(fbsrm1.extend(109999, 50) == " << fbsrm1.extend(109999, 50) << "ull);" << std::endl;
        std::cout << "        CPPUNIT_ASSERT(fbsrm1.extend(109998, 50) == " << fbsrm1.extend(109998, 50) << "ull);" << std::endl;
        std::cout << "        CPPUNIT_ASSERT(fbsrm1.extend(109995, 50) == " << fbsrm1.extend(109995, 50) << "ull);" << std::endl;
        std::cout << "        CPPUNIT_ASSERT(fbsrm1.extend(109997, 50) == " << fbsrm1.extend(109997, 50) << "ull);" << std::endl;
        std::cout << "        CPPUNIT_ASSERT(fbsrm1.extend(109996, 50) == " << fbsrm1.extend(109996, 50) << "ull);" << std::endl;
        CPPUNIT_ASSERT(fbsrm1.extend(210001, 20) == 19199712ull);
        CPPUNIT_ASSERT(fbsrm1.extend(210001, 50) == 19200792ull);
        CPPUNIT_ASSERT(fbsrm1.extend(210001, 1000) == 19201008ull);
        CPPUNIT_ASSERT(fbsrm1.extend(110000, 50) == 19200792ull);
        CPPUNIT_ASSERT(fbsrm1.extend(109999, 50) == 19205024ull);
        CPPUNIT_ASSERT(fbsrm1.extend(109998, 50) == 19205240ull);
        CPPUNIT_ASSERT(fbsrm1.extend(109995, 50) == 19205456ull);
        CPPUNIT_ASSERT(fbsrm1.extend(109997, 50) == 9599712ull);
        CPPUNIT_ASSERT(fbsrm1.extend(109996, 50) == 19205672ull);*/
        CPPUNIT_ASSERT(fbsrm1.extend(210001, 20) == 17600400ull);
        CPPUNIT_ASSERT(fbsrm1.extend(210001, 50) == 17600776ull);
        CPPUNIT_ASSERT(fbsrm1.extend(210001, 1000) == 17600984ull);
        CPPUNIT_ASSERT(fbsrm1.extend(110000, 50) == 8800464ull);
        CPPUNIT_ASSERT(fbsrm1.extend(109999, 50) == 17600776ull);
        CPPUNIT_ASSERT(fbsrm1.extend(109998, 50) == 17604992ull);
        CPPUNIT_ASSERT(fbsrm1.extend(109995, 50) == 17605200ull);
        CPPUNIT_ASSERT(fbsrm1.extend(109997, 50) == 8799736ull);
        CPPUNIT_ASSERT(fbsrm1.extend(109996, 50) == 17605408ull);
#endif

        //std::cerr << std::endl << "        CPPUNIT_ASSERT(fbsrm1.extend(109996, 50) == " << fbsrm1.extend(109996, 50) << "ull);" << std::endl;
        FileBasedSparseRow row4 = fbsrm1.row(109995);
        check_row_operations(row4);

        std::vector<size_t> columns;
        columns.push_back(1);
        columns.push_back(2);
        columns.push_back(4);
        columns.push_back(5);
        columns.push_back(100);
        columns.push_back(10000);
        columns.push_back(100000);

        fbsrm1.set_row(123456, columns);
        FileBasedSparseRow row5 = fbsrm1.row(123456);
        CPPUNIT_ASSERT(row5.highest_column() == 100000);
        FileBasedSparseRow::const_iterator it5 = row5.begin();
        CPPUNIT_ASSERT(*it5 == 1);
        ++it5;
        CPPUNIT_ASSERT(*it5 == 2);
        ++it5;
        CPPUNIT_ASSERT(*it5 == 4);
        ++it5;
        CPPUNIT_ASSERT(*it5 == 5);
        ++it5;
        CPPUNIT_ASSERT(*it5 == 100);
        ++it5;
        CPPUNIT_ASSERT(*it5 == 10000);
        ++it5;
        CPPUNIT_ASSERT(*it5 == 100000);

        int columns1[7] = { 1, 4, 10000, 100000, 100, 5, 2 };
        fbsrm1.set_row(123457, columns1, 7);
        FileBasedSparseRow row6 = fbsrm1.row(123457);
        CPPUNIT_ASSERT(row6.highest_column() == 100000);
        FileBasedSparseRow::const_iterator it6 = row6.begin();
        CPPUNIT_ASSERT(*it6 == 1);
        ++it6;
        CPPUNIT_ASSERT(*it6 == 2);
        ++it6;
        CPPUNIT_ASSERT(*it6 == 4);
        ++it6;
        CPPUNIT_ASSERT(*it6 == 5);
        ++it6;
        CPPUNIT_ASSERT(*it6 == 100);
        ++it6;
        CPPUNIT_ASSERT(*it6 == 10000);
        ++it6;
        CPPUNIT_ASSERT(*it6 == 100000);

    }

    void testSparseMatrix()
    {
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", SparseMatrix());
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", SparseMatrix(1000));
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", SparseMatrix(2000, 1000));
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", SparseMatrix(20000, 1000));

        SparseMatrix sm1;
        CPPUNIT_ASSERT(sm1.rows() == 0);
        CPPUNIT_ASSERT(sm1.cols() == 0);

        SparseMatrix sm2(2000);
        CPPUNIT_ASSERT(sm2.rows() == 2000);
        CPPUNIT_ASSERT(sm2.cols() == 0);

        SparseMatrix sm3(2000, 1500);
        CPPUNIT_ASSERT(sm3.rows() == 2000);
        CPPUNIT_ASSERT(sm3.cols() == 1500);

        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", sm3.clear());
        CPPUNIT_ASSERT(sm3.rows() == 0);
        CPPUNIT_ASSERT(sm3.cols() == 0);

        sm1.set_size(1000, 1000);
        CPPUNIT_ASSERT(sm1.rows() == 1000);
        CPPUNIT_ASSERT(sm1.cols() == 1000);

        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", sm1.clear());
        CPPUNIT_ASSERT(sm1.rows() == 0);
        CPPUNIT_ASSERT(sm1.cols() == 0);

        CPPUNIT_ASSERT(sm1.xor(10, 5) == ISparseRow::XOR_ADDED);
        CPPUNIT_ASSERT(sm1.rows() == 11);
        CPPUNIT_ASSERT(sm1.cols() == 6);

        for (size_t row = 0; row < 10; ++row)
        {
            CPPUNIT_ASSERT(sm1.row_size(row) == 0);
            CPPUNIT_ASSERT(sm1.row_size(row) == 0);
        }
        CPPUNIT_ASSERT(sm1.row_size(10) == 1);
        CPPUNIT_ASSERT(sm1.row_size(10) == 1);
        CPPUNIT_ASSERT(*sm1.begin(10) == 5);

        CPPUNIT_ASSERT(sm1.xor(10, 5) == ISparseRow::XOR_REMOVED);
        CPPUNIT_ASSERT(sm1.rows() == 11);
        CPPUNIT_ASSERT(sm1.cols() == 6);
        for (size_t row = 0; row < 11; ++row)
        {
            CPPUNIT_ASSERT(sm1.row_size(row) == 0);
        }

        sm1.removeEmptyRows();
        CPPUNIT_ASSERT(sm1.rows() == 0);
        CPPUNIT_ASSERT(sm1.cols() == 6);

        sm1.set_cols();
        CPPUNIT_ASSERT(sm1.rows() == 0);
        CPPUNIT_ASSERT(sm1.cols() == 0);

        CPPUNIT_ASSERT(sm1.xor(10, 5) == ISparseRow::XOR_ADDED);
        CPPUNIT_ASSERT(sm1.rows() == 11);
        CPPUNIT_ASSERT(sm1.cols() == 6);

        CPPUNIT_ASSERT(sm1.xor(11, 3) == ISparseRow::XOR_ADDED);
        CPPUNIT_ASSERT(sm1.rows() == 12);
        CPPUNIT_ASSERT(sm1.cols() == 6);

        CPPUNIT_ASSERT(sm1.xor(6, 20) == ISparseRow::XOR_ADDED);
        CPPUNIT_ASSERT(sm1.rows() == 12);
        CPPUNIT_ASSERT(sm1.cols() == 21);

        sm1.clear_row(1);
        CPPUNIT_ASSERT(sm1.rows() == 12);
        CPPUNIT_ASSERT(sm1.cols() == 21);

        sm1.clear_row(6);
        CPPUNIT_ASSERT(sm1.rows() == 12);
        CPPUNIT_ASSERT(sm1.cols() == 21);

        sm1.removeEmptyRows();
        CPPUNIT_ASSERT(sm1.rows() == 2);
        CPPUNIT_ASSERT(sm1.cols() == 21);
        
        sm1.set_cols();
        CPPUNIT_ASSERT(sm1.rows() == 2);
        CPPUNIT_ASSERT(sm1.cols() == 6);

        std::ostringstream oss;
        oss << sm1;
        std::ostringstream oss1;
        oss1 << "2" << std::endl;
        oss1 << "1 5" << std::endl;
        oss1 << "1 3" << std::endl;
        CPPUNIT_ASSERT(oss.str() == oss1.str());

        CPPUNIT_ASSERT(sm1.xor(6, 20) == ISparseRow::XOR_ADDED);
        CPPUNIT_ASSERT(sm1.rows() == 7);
        CPPUNIT_ASSERT(sm1.cols() == 21);

        oss.str("");
        oss << sm1;
        oss1.str("");
        oss1 << "7" << std::endl;
        oss1 << "1 5" << std::endl;
        oss1 << "1 3" << std::endl;
        oss1 << "0" << std::endl;
        oss1 << "0" << std::endl;
        oss1 << "0" << std::endl;
        oss1 << "0" << std::endl;
        oss1 << "1 20" << std::endl;
        CPPUNIT_ASSERT(oss.str() == oss1.str());

        CPPUNIT_ASSERT(sm1.xor(4, 2) == ISparseRow::XOR_ADDED);
        CPPUNIT_ASSERT(sm1.xor(4, 3) == ISparseRow::XOR_ADDED);
        CPPUNIT_ASSERT(sm1.xor(4, 20) == ISparseRow::XOR_ADDED);
        CPPUNIT_ASSERT(sm1.xor(4, 21) == ISparseRow::XOR_ADDED);
        CPPUNIT_ASSERT(sm1.xor(4, 22) == ISparseRow::XOR_ADDED);
        CPPUNIT_ASSERT(sm1.rows() == 7);
        CPPUNIT_ASSERT(sm1.cols() == 23);

        oss.str("");
        oss << sm1;
        oss1.str("");
        oss1 << "7" << std::endl;
        oss1 << "1 5" << std::endl;
        oss1 << "1 3" << std::endl;
        oss1 << "0" << std::endl;
        oss1 << "0" << std::endl;
        oss1 << "5 2 3 20 21 22" << std::endl;
        oss1 << "0" << std::endl;
        oss1 << "1 20" << std::endl;
        CPPUNIT_ASSERT(oss.str() == oss1.str());

        sm1.removeEmptyRows();
        CPPUNIT_ASSERT(sm1.rows() == 4);
        CPPUNIT_ASSERT(sm1.cols() == 23);

        oss.str("");
        oss << sm1;
        oss1.str("");
        oss1 << "4" << std::endl;
        oss1 << "1 5" << std::endl;
        oss1 << "1 3" << std::endl;
        oss1 << "5 2 3 20 21 22" << std::endl;
        oss1 << "1 20" << std::endl;
        CPPUNIT_ASSERT(oss.str() == oss1.str());

        std::istringstream iss;
        iss.str(oss1.str());
        SparseMatrix sm10;
        iss >> sm10;
        std::ostringstream oss2;
        oss2 << sm10;
        CPPUNIT_ASSERT(oss.str() == oss2.str());

        sm1.set_write_row_count(false);

        oss.str("");
        oss << sm1;
        oss1.str("");
        oss1 << "1 5" << std::endl;
        oss1 << "1 3" << std::endl;
        oss1 << "5 2 3 20 21 22" << std::endl;
        oss1 << "1 20" << std::endl;
        CPPUNIT_ASSERT(oss.str() == oss1.str());

        SparseMatrix sm1t;
        transpose(sm1, sm1t, 0, false);
        CPPUNIT_ASSERT(sm1t.rows() == 23);
        CPPUNIT_ASSERT(sm1t.cols() == 4);

        oss1.str("");
        oss1 << "5" << std::endl;
        oss1 << "1 5" << std::endl;
        oss1 << "1 3" << std::endl;
        oss1 << "5 2 3 20 21 12000" << std::endl;
        oss1 << "5 2 3 20 21 120000" << std::endl;
        oss1 << "1 20" << std::endl;

        std::istringstream iss2;
        iss2.str(oss1.str());
        SparseMatrix sm30;
        iss2 >> sm30;

        SparseMatrix sm30t;
        transpose(sm30, sm30t, 0, false);
        CPPUNIT_ASSERT(sm30t.rows() == 120001);
        CPPUNIT_ASSERT(sm30t.cols() == 5);

        SparseMatrix sm30tt;
        transpose(sm30t, sm30tt, 0, false);
        CPPUNIT_ASSERT(sm30tt.rows() == 5);
        CPPUNIT_ASSERT(sm30tt.cols() == 120001);

        SparseMatrix sm40(10);
        std::vector<size_t> columns;
        columns.push_back(1);
        columns.push_back(100);
        columns.push_back(2);
        columns.push_back(5);
        columns.push_back(100000);
        columns.push_back(10000);
        columns.push_back(4);
        sm40.set_row(0, columns);
        CPPUNIT_ASSERT(sm40.rows() == 10);
        CPPUNIT_ASSERT(sm40.cols() == 100001);
        columns.clear();
        columns.push_back(100000);
        columns.push_back(100);
        columns.push_back(10000000);
        columns.push_back(10000);
        columns.push_back(1000000);
        sm40.set_row(9, columns);
        CPPUNIT_ASSERT(sm40.rows() == 10);
        CPPUNIT_ASSERT(sm40.cols() == 10000001);

        int columns1[7] = { 1, 100, 10000, 20000000, 2, 5, 4 };
        sm40.set_row(5, columns1, 7);
        CPPUNIT_ASSERT(sm40.rows() == 10);
        CPPUNIT_ASSERT(sm40.cols() == 20000001);
    }
};

int main()
{
   CppUnit::TextUi::TestRunner runner;
   runner.addTest(SparseMatrixTest::suite());
   runner.run();
   return 0;
}
