#include <iostream>
#include "SparseMatrix.h"
#include <cstdio>
#include <fstream>
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
    CPPUNIT_TEST(testBitMatrix);
    CPPUNIT_TEST(testBitMatrix64);
    CPPUNIT_TEST(testSparseMatrix3MediumRowStripes);
    CPPUNIT_TEST(testSparseMatrix3MultiplyMatchesSparseMatrix4);
    CPPUNIT_TEST(testSparseMatrix3VeryDenseRowsInMemory);
    CPPUNIT_TEST_SUITE_END();

public:
    void setUp()
    {
    }

    void tearDown()
    {
        std::remove("SparseMatrix3MediumRowStripes.dat");
        std::remove("SparseMatrix3MixedRows.dat");
        std::remove("SparseMatrix3VeryDenseRows.dat");
        std::remove("dense_rows.txt");
    }

    static std::string sparseMatrix3MediumRowContent()
    {
        std::ostringstream oss;
        oss << "10000\n";
        oss << "1 4\n";
        oss << "2 1 65540\n";
        oss << "2 32769 32770\n";
        oss << "3 2 32772 65543\n";
        return oss.str();
    }

    static void writeSparseMatrix3MediumRowFile()
    {
        std::ofstream os("SparseMatrix3MediumRowStripes.dat");
        os << sparseMatrix3MediumRowContent();
    }

    static std::string sparseMatrix3VeryDenseRowContent()
    {
        std::ostringstream oss;
        oss << "10\n";
        oss << "5 0 1 2 3 4\n";
        for (int i = 0; i < 9; ++i)
        {
            oss << "1 " << (i + 5) << "\n";
        }
        return oss.str();
    }

    static void writeSparseMatrix3VeryDenseRowFile()
    {
        std::ofstream os("SparseMatrix3VeryDenseRows.dat");
        os << sparseMatrix3VeryDenseRowContent();
    }

    static void writeSparseMatrix3MixedRowFile()
    {
        const size_t sparse_rows = 20000;
        const size_t medium_rows = 10;
        std::ofstream os("SparseMatrix3MixedRows.dat");
        os << (sparse_rows + medium_rows) << "\n";
        for (size_t row = 0; row < sparse_rows; ++row)
        {
            const size_t col0 = (row * 37) % 98306;
            if (row % 5 == 0)
            {
                const size_t col1 = 65536 + (row % 257);
                os << "2 " << col0 << " " << col1 << "\n";
            }
            else
            {
                os << "1 " << col0 << "\n";
            }
        }

        os << "4 0 32768 65536 98305\n";
        os << "4 1 32769 65537 98304\n";
        os << "4 0 32769 65536 98304\n";
        os << "4 2 32770 65538 98303\n";
        os << "4 2 32768 65538 98303\n";
        os << "4 3 32771 65539 98302\n";
        os << "4 3 32769 65539 98302\n";
        os << "4 4 32772 65540 98301\n";
        os << "4 4 32768 65540 98301\n";
        os << "4 5 32773 65541 98300\n";
    }

    void check_row_operations(FileBasedSparseRow& row)
    {

        CPPUNIT_ASSERT(row.highest_column() == -1);
        CPPUNIT_ASSERT(row.do_xor(2) == ISparseRow::XOR_ADDED);
        CPPUNIT_ASSERT(row.highest_column() == 2);
        CPPUNIT_ASSERT(row.do_xor(9) == ISparseRow::XOR_ADDED);
        CPPUNIT_ASSERT(row.highest_column() == 9);
        CPPUNIT_ASSERT(row.do_xor(3) == ISparseRow::XOR_ADDED);
        CPPUNIT_ASSERT(row.highest_column() == 9);
        CPPUNIT_ASSERT(row.do_xor(10) == ISparseRow::XOR_ADDED);
        CPPUNIT_ASSERT(row.highest_column() == 10);
        CPPUNIT_ASSERT(row.do_xor(4) == ISparseRow::XOR_ADDED);
        CPPUNIT_ASSERT(row.highest_column() == 10);
        CPPUNIT_ASSERT(row.size() == 5);
        auto it = row.begin();
        CPPUNIT_ASSERT(*it == 2);
        ++it;
        CPPUNIT_ASSERT(*it == 3);
        ++it;
        CPPUNIT_ASSERT(*it == 4);
        ++it;
        CPPUNIT_ASSERT(*it == 9);
        ++it;
        CPPUNIT_ASSERT(*it == 10);
        ++it;
        CPPUNIT_ASSERT(it == row.end());
        CPPUNIT_ASSERT(row.do_xor(4) == ISparseRow::XOR_REMOVED);
        CPPUNIT_ASSERT(row.size() == 4);
        CPPUNIT_ASSERT(row.highest_column() == 10);
        CPPUNIT_ASSERT(row.do_xor(10) == ISparseRow::XOR_REMOVED);
        CPPUNIT_ASSERT(row.size() == 3);
        CPPUNIT_ASSERT(row.highest_column() == 9);

        for (size_t i = 0; i < 17; ++i)
        {
            CPPUNIT_ASSERT(row.do_xor(i + 1000) == ISparseRow::XOR_ADDED);
        }
        CPPUNIT_ASSERT(row.size() == 20);
        CPPUNIT_ASSERT(row.highest_column() == 1016);
        static size_t expected_row_data[] = { 2, 3, 9, 1000, 1001, 1002, 1003, 1004, 1005, 1006, 1007, 1008, 1009, 1010, 1011, 1012, 1013, 1014, 1015, 1016 };
        size_t j = 0;
        for (auto it = row.begin();
                it != row.end();
                ++it, ++j)
        {
            CPPUNIT_ASSERT(*it == expected_row_data[j]);
        }

        CPPUNIT_ASSERT(row.do_xor(1017) == ISparseRow::XOR_ADDED);
        CPPUNIT_ASSERT(row.size() == 21);
        CPPUNIT_ASSERT(row.highest_column() == 1017);
        static size_t expected_row_data1[] = { 2, 3, 9, 1000, 1001, 1002, 1003, 1004, 1005, 1006, 1007, 1008, 1009, 1010, 1011, 1012, 1013, 1014, 1015, 1016, 1017 };
        j = 0;
        for (auto it = row.begin();
                it != row.end();
                ++it, ++j)
        {
            CPPUNIT_ASSERT(*it == expected_row_data1[j]);
        }

        SparseRow sr(row.size());
        row.copy(sr);
        j = 0;
        for (auto it = sr.begin();
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
        auto row1 = fbsrm1.row(10000);
        check_row_operations(row1);

        auto row2 = fbsrm1.row(110000);
        check_row_operations(row2);

        auto row3 = fbsrm1.row(210000);
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
#if 0
        std::cout << "        CPPUNIT_ASSERT(fbsrm1.extend(210001, 20) == " << fbsrm1.extend(210001, 20) << "ull);" << std::endl;
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
        CPPUNIT_ASSERT(fbsrm1.extend(109996, 50) == 19205672ull);
#endif

        CPPUNIT_ASSERT(fbsrm1.extend(210001, 20) == 19200400ull);
        CPPUNIT_ASSERT(fbsrm1.extend(210001, 50) == 19200792ull);
        CPPUNIT_ASSERT(fbsrm1.extend(210001, 1000) == 19201008ull);
        CPPUNIT_ASSERT(fbsrm1.extend(110000, 50) == 9600488ull);
        CPPUNIT_ASSERT(fbsrm1.extend(109999, 50) == 19200792ull);
        CPPUNIT_ASSERT(fbsrm1.extend(109998, 50) == 19205024ull);
        CPPUNIT_ASSERT(fbsrm1.extend(109995, 50) == 19205240ull);
        CPPUNIT_ASSERT(fbsrm1.extend(109997, 50) == 9599712ull);
        CPPUNIT_ASSERT(fbsrm1.extend(109996, 50) == 19205456ull);
#endif

        //std::cerr << std::endl << "        CPPUNIT_ASSERT(fbsrm1.extend(109996, 50) == " << fbsrm1.extend(109996, 50) << "ull);" << std::endl;
        auto row4 = fbsrm1.row(109995);
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
        auto row5 = fbsrm1.row(123456);
        CPPUNIT_ASSERT(row5.highest_column() == 100000);
        auto it5 = row5.begin();
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
        auto row6 = fbsrm1.row(123457);
        CPPUNIT_ASSERT(row6.highest_column() == 100000);
        auto it6 = row6.begin();
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

        CPPUNIT_ASSERT(sm1.do_xor(10, 5) == ISparseRow::XOR_ADDED);
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

        CPPUNIT_ASSERT(sm1.do_xor(10, 5) == ISparseRow::XOR_REMOVED);
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

        CPPUNIT_ASSERT(sm1.do_xor(10, 5) == ISparseRow::XOR_ADDED);
        CPPUNIT_ASSERT(sm1.rows() == 11);
        CPPUNIT_ASSERT(sm1.cols() == 6);

        CPPUNIT_ASSERT(sm1.do_xor(11, 3) == ISparseRow::XOR_ADDED);
        CPPUNIT_ASSERT(sm1.rows() == 12);
        CPPUNIT_ASSERT(sm1.cols() == 6);

        CPPUNIT_ASSERT(sm1.do_xor(6, 20) == ISparseRow::XOR_ADDED);
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

        CPPUNIT_ASSERT(sm1.do_xor(6, 20) == ISparseRow::XOR_ADDED);
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

        CPPUNIT_ASSERT(sm1.do_xor(4, 2) == ISparseRow::XOR_ADDED);
        CPPUNIT_ASSERT(sm1.do_xor(4, 3) == ISparseRow::XOR_ADDED);
        CPPUNIT_ASSERT(sm1.do_xor(4, 20) == ISparseRow::XOR_ADDED);
        CPPUNIT_ASSERT(sm1.do_xor(4, 21) == ISparseRow::XOR_ADDED);
        CPPUNIT_ASSERT(sm1.do_xor(4, 22) == ISparseRow::XOR_ADDED);
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

    void testBitMatrix()
    {
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", BitMatrix());
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", BitMatrix(100, 30));
        CPPUNIT_ASSERT_THROW_MESSAGE("", BitMatrix(2000, 1000), char*);
        std::ostringstream oss;
        oss << "32" << std::endl;
        oss << "32" << std::endl;
        oss << "11101100011001010000111101110010" << std::endl;
        oss << "11101110111000000001010010000001" << std::endl;
        oss << "11100011100111101100111011101010" << std::endl;
        oss << "00000011010010011110100101010101" << std::endl;
        oss << "11000011011100001101001110101011" << std::endl;
        oss << "11000010010101110011110111011011" << std::endl;
        oss << "01111100101010100110011110100101" << std::endl;
        oss << "00111000100101011010000111010110" << std::endl;
        oss << "01100011010110100111000001001001" << std::endl;
        oss << "11011100111101011010001101101111" << std::endl;
        oss << "11001010010010110000111101010011" << std::endl;
        oss << "00101101110111011100110010000110" << std::endl;
        oss << "00110010101110000010010110101000" << std::endl;
        oss << "10100101010100110110111010101001" << std::endl;
        oss << "00100110101001110000010111111001" << std::endl;
        oss << "10010101011101110000001111101101" << std::endl;
        oss << "00111001010100000010011000000101" << std::endl;
        oss << "00111010100101000101000101000001" << std::endl;
        oss << "00010111110011001000001101110101" << std::endl;
        oss << "01001100100000000100000111110111" << std::endl;
        oss << "10110100001101000000011110100001" << std::endl;
        oss << "11100110001111101000101010011011" << std::endl;
        oss << "10101010011001011010111100000110" << std::endl;
        oss << "10011111011010110111101000010111" << std::endl;
        oss << "01101111000111110001110011011100" << std::endl;
        oss << "10110101111000110111000011001011" << std::endl;
        oss << "10101010010011110011100000100010" << std::endl;
        oss << "10010101001000100011010110000000" << std::endl;
        oss << "00101100110011110000010011001101" << std::endl;
        oss << "00010011010100011011001110001101" << std::endl;
        oss << "10101101011100000001011101100000" << std::endl;
        oss << "01011110111001111111110101001101" << std::endl;
        BitMatrix bm;
        CPPUNIT_ASSERT(bm.isZero());
        std::istringstream iss(oss.str());
        bm.read(iss);
        CPPUNIT_ASSERT(bm.rows() == 32);
        CPPUNIT_ASSERT(bm.cols() == 32);
        BitMatrix bm1(bm);

        CPPUNIT_ASSERT(bm == bm1);
        BitMatrix kerbm;
        kernel(bm, kerbm);
        CPPUNIT_ASSERT(kerbm.isZero());
    }

    void testBitMatrix64()
    {
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", BitMatrix64());
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", BitMatrix64(100, 30));
        CPPUNIT_ASSERT_THROW_MESSAGE("", BitMatrix64(2000, 1000), char*);
        std::ostringstream oss;
        oss << "64" << std::endl;
        oss << "64" << std::endl;
        oss << "1110110001100101000011110111001011101100011001010000111101110010" << std::endl;
        oss << "1110111011100000000101001000000111101110111000000001010010000001" << std::endl;
        oss << "1110001110011110110011101110101011100011100111101100111011101010" << std::endl;
        oss << "0000001101001001111010010101010100000011010010011110100101010101" << std::endl;
        oss << "1100001101110000110100111010101111000011011100001101001110101011" << std::endl;
        oss << "1100001001010111001111011101101111000010010101110011110111011011" << std::endl;
        oss << "0111110010101010011001111010010101111100101010100110011110100101" << std::endl;
        oss << "0011100010010101101000011101011000111000100101011010000111010110" << std::endl;
        oss << "0110001101011010011100000100100101100011010110100111000001001001" << std::endl;
        oss << "1101110011110101101000110110111111011100111101011010001101101111" << std::endl;
        oss << "1100101001001011000011110101001111001010010010110000111101010011" << std::endl;
        oss << "0010110111011101110011001000011000101101110111011100110010000110" << std::endl;
        oss << "0011001010111000001001011010100000110010101110000010010110101000" << std::endl;
        oss << "1010010101010011011011101010100110100101010100110110111010101001" << std::endl;
        oss << "0010011010100111000001011111100100100110101001110000010111111001" << std::endl;
        oss << "1001010101110111000000111110110110010101011101110000001111101101" << std::endl;
        oss << "0011100101010000001001100000010100111001010100000010011000000101" << std::endl;
        oss << "0011101010010100010100010100000100111010100101000101000101000001" << std::endl;
        oss << "0001011111001100100000110111010100010111110011001000001101110101" << std::endl;
        oss << "0100110010000000010000011111011101001100100000000100000111110111" << std::endl;
        oss << "1011010000110100000001111010000110110100001101000000011110100001" << std::endl;
        oss << "1110011000111110100010101001101111100110001111101000101010011011" << std::endl;
        oss << "1010101001100101101011110000011010101010011001011010111100000110" << std::endl;
        oss << "1001111101101011011110100001011110011111011010110111101000010111" << std::endl;
        oss << "0110111100011111000111001101110001101111000111110001110011011100" << std::endl;
        oss << "1011010111100011011100001100101110110101111000110111000011001011" << std::endl;
        oss << "1010101001001111001110000010001010101010010011110011100000100010" << std::endl;
        oss << "1001010100100010001101011000000010010101001000100011010110000000" << std::endl;
        oss << "0010110011001111000001001100110100101100110011110000010011001101" << std::endl;
        oss << "0001001101010001101100111000110100010011010100011011001110001101" << std::endl;
        oss << "1010110101110000000101110110000010101101011100000001011101100000" << std::endl;
        oss << "0101111011100111111111010100110101011110111001111111110101001101" << std::endl;
        oss << "0001001110011010111100001000110100010011100110101111000010001101" << std::endl;
        oss << "0001000100011111111010110111111000010001000111111110101101111110" << std::endl;
        oss << "0001110001100001001100010001010100011100011000010011000100010101" << std::endl;
        oss << "1111110010110110000101101010101011111100101101100001011010101010" << std::endl;
        oss << "0011110010001111001011000101010000111100100011110010110001010100" << std::endl;
        oss << "0011110110101000110000100010010000111101101010001100001000100100" << std::endl;
        oss << "1000001101010101100110000101101010000011010101011001100001011010" << std::endl;
        oss << "1100011101101010010111100010100111000111011010100101111000101001" << std::endl;
        oss << "1001110010100101100011111011011010011100101001011000111110110110" << std::endl;
        oss << "0010001100001010010111001001000000100011000010100101110010010000" << std::endl;
        oss << "0011010110110100111100001010110000110101101101001111000010101100" << std::endl;
        oss << "1101001000100010001100110111100111010010001000100011001101111001" << std::endl;
        oss << "1100110101000111110110100101011111001101010001111101101001010111" << std::endl;
        oss << "0101101010101100100100010101011001011010101011001001000101010110" << std::endl;
        oss << "1101100101011000111110100000011011011001010110001111101000000110" << std::endl;
        oss << "0110101010001000111111000001001001101010100010001111110000010010" << std::endl;
        oss << "1100011010101111110110011111101011000110101011111101100111111010" << std::endl;
        oss << "1100010101101011101011101011111011000101011010111010111010111110" << std::endl;
        oss << "1110100000110011011111001000101011101000001100110111110010001010" << std::endl;
        oss << "1011001101111111101111100000100010110011011111111011111000001000" << std::endl;
        oss << "0100101111001011111110000101111001001011110010111111100001011110" << std::endl;
        oss << "0001100111000001011101010110010000011001110000010111010101100100" << std::endl;
        oss << "0101010110011010010100001111100101010101100110100101000011111001" << std::endl;
        oss << "0110000010010100100001011110100001100000100101001000010111101000" << std::endl;
        oss << "1001000011100000111000110010001110010000111000001110001100100011" << std::endl;
        oss << "0100101000011100100011110011010001001010000111001000111100110100" << std::endl;
        oss << "0101010110110000110001111101110101010101101100001100011111011101" << std::endl;
        oss << "0110101011011101110010100111111101101010110111011100101001111111" << std::endl;
        oss << "1101001100110000111110110011001011010011001100001111101100110010" << std::endl;
        oss << "1110110010101110010011000111001011101100101011100100110001110010" << std::endl;
        oss << "0101001010001111111010001001111101010010100011111110100010011111" << std::endl;
        oss << "1010000100011000000000101011001010100001000110000000001010110010" << std::endl;
        BitMatrix64 bm;
        CPPUNIT_ASSERT(bm.isZero());
        std::istringstream iss(oss.str());
        bm.read(iss);
        CPPUNIT_ASSERT(bm.rows() == 64);
        CPPUNIT_ASSERT(bm.cols() == 64);
        BitMatrix64 bm1(bm);

        CPPUNIT_ASSERT(bm == bm1);
        BitMatrix64 kerbm;
        kernel(bm, kerbm);
        CPPUNIT_ASSERT(!kerbm.isZero());
    }

    void testSparseMatrix3MediumRowStripes()
    {
        writeSparseMatrix3MediumRowFile();

        SparseMatrix3 sm("SparseMatrix3MediumRowStripes.dat");
        CPPUNIT_ASSERT_EQUAL(static_cast<size_t>(4), sm.rows());
        CPPUNIT_ASSERT_EQUAL(static_cast<size_t>(65544), sm.cols());

        std::ostringstream serialized;
        serialized << sm;
        CPPUNIT_ASSERT_EQUAL(std::string("4\n1 4\n2 1 65540\n2 32769 32770\n3 2 32772 65543\n"), serialized.str());

        BitMatrix X(sm.cols(), 3);
        X.row_[4] = 0x1UL;
        X.row_[1] = 0x2UL;
        X.row_[65540] = 0x4UL;
        X.row_[32769] = 0x3UL;
        X.row_[32770] = 0x5UL;
        X.row_[2] = 0x6UL;
        X.row_[32772] = 0x7UL;
        X.row_[65543] = 0x1UL;

        BitMatrix AX;
        multiply(sm, X, AX);
        CPPUNIT_ASSERT_EQUAL(static_cast<size_t>(4), AX.rows());
        CPPUNIT_ASSERT_EQUAL(static_cast<size_t>(3), AX.cols());
        CPPUNIT_ASSERT_EQUAL(static_cast<uint32_t>(0x1UL), AX.row_[0]);
        CPPUNIT_ASSERT_EQUAL(static_cast<uint32_t>(0x6UL), AX.row_[1]);
        CPPUNIT_ASSERT_EQUAL(static_cast<uint32_t>(0x6UL), AX.row_[2]);
        CPPUNIT_ASSERT_EQUAL(static_cast<uint32_t>(0x0UL), AX.row_[3]);

        BitMatrix Xt(sm.rows(), 3);
        Xt.row_[0] = 0x1UL;
        Xt.row_[1] = 0x2UL;
        Xt.row_[2] = 0x4UL;
        Xt.row_[3] = 0x7UL;

        BitMatrix AtX;
        multiplyt(sm, Xt, AtX);
        CPPUNIT_ASSERT_EQUAL(sm.cols(), AtX.rows());
        CPPUNIT_ASSERT_EQUAL(static_cast<size_t>(3), AtX.cols());
        CPPUNIT_ASSERT_EQUAL(static_cast<uint32_t>(0x0UL), AtX.row_[0]);
        CPPUNIT_ASSERT_EQUAL(static_cast<uint32_t>(0x2UL), AtX.row_[1]);
        CPPUNIT_ASSERT_EQUAL(static_cast<uint32_t>(0x7UL), AtX.row_[2]);
        CPPUNIT_ASSERT_EQUAL(static_cast<uint32_t>(0x1UL), AtX.row_[4]);
        CPPUNIT_ASSERT_EQUAL(static_cast<uint32_t>(0x4UL), AtX.row_[32769]);
        CPPUNIT_ASSERT_EQUAL(static_cast<uint32_t>(0x4UL), AtX.row_[32770]);
        CPPUNIT_ASSERT_EQUAL(static_cast<uint32_t>(0x7UL), AtX.row_[32772]);
        CPPUNIT_ASSERT_EQUAL(static_cast<uint32_t>(0x2UL), AtX.row_[65540]);
        CPPUNIT_ASSERT_EQUAL(static_cast<uint32_t>(0x7UL), AtX.row_[65543]);

        BitMatrix64 X64(sm.cols(), 3);
        X64.row_[4] = 0x1ULL;
        X64.row_[1] = 0x2ULL;
        X64.row_[65540] = 0x4ULL;
        X64.row_[32769] = 0x3ULL;
        X64.row_[32770] = 0x5ULL;
        X64.row_[2] = 0x6ULL;
        X64.row_[32772] = 0x7ULL;
        X64.row_[65543] = 0x1ULL;

        BitMatrix64 AX64;
        multiply(sm, X64, AX64);
        CPPUNIT_ASSERT_EQUAL(static_cast<size_t>(4), AX64.rows());
        CPPUNIT_ASSERT_EQUAL(static_cast<size_t>(3), AX64.cols());
        CPPUNIT_ASSERT_EQUAL(0x1ULL, AX64.row_[0]);
        CPPUNIT_ASSERT_EQUAL(0x6ULL, AX64.row_[1]);
        CPPUNIT_ASSERT_EQUAL(0x6ULL, AX64.row_[2]);
        CPPUNIT_ASSERT_EQUAL(0x0ULL, AX64.row_[3]);

        BitMatrix64 Xt64(sm.rows(), 3);
        Xt64.row_[0] = 0x1ULL;
        Xt64.row_[1] = 0x2ULL;
        Xt64.row_[2] = 0x4ULL;
        Xt64.row_[3] = 0x7ULL;

        BitMatrix64 AtX64;
        multiplyt(sm, Xt64, AtX64);
        CPPUNIT_ASSERT_EQUAL(sm.cols(), AtX64.rows());
        CPPUNIT_ASSERT_EQUAL(static_cast<size_t>(3), AtX64.cols());
        CPPUNIT_ASSERT_EQUAL(0x0ULL, AtX64.row_[0]);
        CPPUNIT_ASSERT_EQUAL(0x2ULL, AtX64.row_[1]);
        CPPUNIT_ASSERT_EQUAL(0x7ULL, AtX64.row_[2]);
        CPPUNIT_ASSERT_EQUAL(0x1ULL, AtX64.row_[4]);
        CPPUNIT_ASSERT_EQUAL(0x4ULL, AtX64.row_[32769]);
        CPPUNIT_ASSERT_EQUAL(0x4ULL, AtX64.row_[32770]);
        CPPUNIT_ASSERT_EQUAL(0x7ULL, AtX64.row_[32772]);
        CPPUNIT_ASSERT_EQUAL(0x2ULL, AtX64.row_[65540]);
        CPPUNIT_ASSERT_EQUAL(0x7ULL, AtX64.row_[65543]);
    }

    void testSparseMatrix3MultiplyMatchesSparseMatrix4()
    {
        writeSparseMatrix3MixedRowFile();

        SparseMatrix3 sm3("SparseMatrix3MixedRows.dat");
        SparseMatrix4 sm4("SparseMatrix3MixedRows.dat");
        CPPUNIT_ASSERT_EQUAL(static_cast<size_t>(sm4.rows()), sm3.rows());
        CPPUNIT_ASSERT_EQUAL(static_cast<size_t>(sm4.cols()), sm3.cols());

        BitMatrix X(sm3.cols(), 31);
        for (size_t row = 0; row < X.rows(); ++row)
        {
            X.row_[row] = static_cast<uint32_t>((row * 2654435761UL) ^ (row << 7) ^ (row >> 3));
        }

        BitMatrix AX3;
        BitMatrix AX4;
        multiply(sm3, X, AX3);
        multiply(sm4, X, AX4);
        CPPUNIT_ASSERT(AX3 == AX4);

        BitMatrix Xt(sm3.rows(), 29);
        for (size_t row = 0; row < Xt.rows(); ++row)
        {
            Xt.row_[row] = static_cast<uint32_t>((row * 2246822519UL) ^ (row << 5) ^ 0x5a5a5a5aUL);
        }

        BitMatrix AtX3;
        BitMatrix AtX4;
        multiplyt(sm3, Xt, AtX3);
        multiplyt(sm4, Xt, AtX4);
        CPPUNIT_ASSERT(AtX3 == AtX4);

        BitMatrix64 X64(sm3.cols(), 61);
        for (size_t row = 0; row < X64.rows(); ++row)
        {
            X64.row_[row] = (static_cast<unsigned long long int>(row) << 32)
                ^ (static_cast<unsigned long long int>(row * 11400714819323198485ULL) >> 3)
                ^ 0x13579bdf2468ace0ULL;
        }

        BitMatrix64 AX364;
        BitMatrix64 AX464;
        multiply(sm3, X64, AX364);
        multiply(sm4, X64, AX464);
        CPPUNIT_ASSERT(AX364 == AX464);

        BitMatrix64 Xt64(sm3.rows(), 63);
        for (size_t row = 0; row < Xt64.rows(); ++row)
        {
            Xt64.row_[row] = (static_cast<unsigned long long int>(row) * 6364136223846793005ULL)
                ^ (static_cast<unsigned long long int>(row) << 11)
                ^ 0xfedcba9876543210ULL;
        }

        BitMatrix64 AtX364;
        BitMatrix64 AtX464;
        multiplyt(sm3, Xt64, AtX364);
        multiplyt(sm4, Xt64, AtX464);
        CPPUNIT_ASSERT(AtX364 == AtX464);
    }

    void testSparseMatrix3VeryDenseRowsInMemory()
    {
        std::remove("dense_rows.txt");
        writeSparseMatrix3VeryDenseRowFile();

        SparseMatrix3 sm("SparseMatrix3VeryDenseRows.dat");
        CPPUNIT_ASSERT_EQUAL(static_cast<size_t>(9), sm.rows());
        CPPUNIT_ASSERT_EQUAL(static_cast<size_t>(14), sm.cols());

        std::ifstream dense_rows_file("dense_rows.txt");
        CPPUNIT_ASSERT(!dense_rows_file.good());

        std::ostringstream serialized;
        serialized << sm;
        CPPUNIT_ASSERT_EQUAL(std::string("9\n1 5\n1 6\n1 7\n1 8\n1 9\n1 10\n1 11\n1 12\n1 13\n5 0 1 2 3 4\n"), serialized.str());

        BitMatrix L(sm.cols(), 3);
        L.row_[0] = 0x1UL;
        L.row_[1] = 0x2UL;
        L.row_[2] = 0x4UL;
        L.row_[3] = 0x7UL;
        L.row_[4] = 0x3UL;
        L.row_[5] = 0x6UL;

        BitMatrix R(sm.cols(), 3);
        R.row_[0] = 0x7UL;
        R.row_[1] = 0x4UL;
        R.row_[2] = 0x2UL;
        R.row_[3] = 0x1UL;
        R.row_[4] = 0x3UL;

        BitMatrix BL;
        sm.multiply_dense_part_by_bit_matrix(L, BL);
        CPPUNIT_ASSERT_EQUAL(static_cast<size_t>(1), BL.rows());
        CPPUNIT_ASSERT_EQUAL(static_cast<size_t>(3), BL.cols());
        CPPUNIT_ASSERT_EQUAL(static_cast<uint32_t>(0x3UL), BL.row_[0]);

        BitMatrix BL2;
        BitMatrix BR2;
        sm.multiply_dense_part_by_bit_matrix(L, R, BL2, BR2);
        CPPUNIT_ASSERT_EQUAL(static_cast<size_t>(1), BL2.rows());
        CPPUNIT_ASSERT_EQUAL(static_cast<size_t>(1), BR2.rows());
        CPPUNIT_ASSERT_EQUAL(static_cast<uint32_t>(0x3UL), BL2.row_[0]);
        CPPUNIT_ASSERT_EQUAL(static_cast<uint32_t>(0x3UL), BR2.row_[0]);

        BitMatrix64 L64(sm.cols(), 3);
        L64.row_[0] = 0x1ULL;
        L64.row_[1] = 0x2ULL;
        L64.row_[2] = 0x4ULL;
        L64.row_[3] = 0x7ULL;
        L64.row_[4] = 0x3ULL;

        BitMatrix64 BL64;
        sm.multiply_dense_part_by_bit_matrix(L64, BL64);
        CPPUNIT_ASSERT_EQUAL(static_cast<size_t>(1), BL64.rows());
        CPPUNIT_ASSERT_EQUAL(static_cast<size_t>(3), BL64.cols());
        CPPUNIT_ASSERT_EQUAL(0x3ULL, BL64.row_[0]);
    }
};

int main()
{
    CppUnit::TextUi::TestRunner runner;
    runner.addTest(SparseMatrixTest::suite());
    runner.run();
    return 0;
}
