#include "blockLanczos.h"
#include "SparseMatrix.h"
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/ui/text/TestRunner.h>
#include <fstream>
#include <sstream>
#include <string>
#include <unistd.h>

// Helper: write a sparse matrix file in SparseMatrix3 format.
// Format: first line = row count, then one line per row:
//   "count col0 col1 ..." or "0" for an empty row.
static void writeMatrixFile(const char* path, const std::string& content)
{
    std::ofstream f(path);
    f << content;
}

class BlockLanczosTest : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE(BlockLanczosTest);
    CPPUNIT_TEST(testKernelNonTrivial);
    CPPUNIT_TEST(testKernelFullRank);
    CPPUNIT_TEST_SUITE_END();

    static const char* const MATRIX_FILE;

public:
    void setUp()
    {
    }

    void tearDown()
    {
        ::unlink(MATRIX_FILE);
    }

    // Build a 100-row, 100-column matrix over GF(2) whose null space has
    // dimension 4, spanned by:
    //   v0 = e_0 + e_1 + e_96
    //   v1 = e_2 + e_3 + e_97
    //   v2 = e_4 + e_5 + e_98
    //   v3 = e_6 + e_7 + e_99
    //
    // Row layout:
    //   row 0 : columns 0, 96   => constraint x[0] = x[96]
    //   row 1 : columns 1, 96   => constraint x[1] = x[96]
    //   row 2 : columns 2, 97   => constraint x[2] = x[97]
    //   row 3 : columns 3, 97   => constraint x[3] = x[97]
    //   row 4 : columns 4, 98   => constraint x[4] = x[98]
    //   row 5 : columns 5, 98   => constraint x[5] = x[98]
    //   row 6 : columns 6, 99   => constraint x[6] = x[99]
    //   row 7 : columns 7, 99   => constraint x[7] = x[99]
    //   rows 8-95: column i     => constraint x[i] = 0
    //   rows 96-99: empty
    static std::string nonTrivialMatrixContent()
    {
        std::ostringstream oss;
        oss << "100\n";
        oss << "2 0 96\n";
        oss << "2 1 96\n";
        oss << "2 2 97\n";
        oss << "2 3 97\n";
        oss << "2 4 98\n";
        oss << "2 5 98\n";
        oss << "2 6 99\n";
        oss << "2 7 99\n";
        for (int i = 8; i <= 95; ++i)
        {
            oss << "1 " << i << "\n";
        }
        // rows 96-99: zero
        oss << "0\n0\n0\n0\n";
        return oss.str();
    }

    // The Block Lanczos algorithm finds vectors x such that B*x = 0 over GF(2).
    // For this matrix the null space has dimension 4.
    void testKernelNonTrivial()
    {
        writeMatrixFile(MATRIX_FILE, nonTrivialMatrixContent());

        BlockLanczos bl(MATRIX_FILE);
        BITMATRIX kerL, kerR;
        bl.kernel(kerL, kerR);

        // The kernel should be non-empty (dim 4 null space).
        CPPUNIT_ASSERT(!kerL.isZero());

        // Primary correctness check: B * kerL = 0 over GF(2).
        SparseMatrix3 B(MATRIX_FILE);
        BITMATRIX BkerL;
        multiply(B, kerL, BkerL);
        CPPUNIT_ASSERT(BkerL.isZero());
    }

    // A 100x100 identity matrix has trivial null space.
    // The correctness invariant B * kerL = 0 should still hold.
    void testKernelFullRank()
    {
        std::ostringstream oss;
        oss << "100\n";
        for (int i = 0; i < 100; ++i)
        {
            oss << "1 " << i << "\n";
        }
        writeMatrixFile(MATRIX_FILE, oss.str());

        BlockLanczos bl(MATRIX_FILE);
        BITMATRIX kerL, kerR;
        bl.kernel(kerL, kerR);

        // Verify B * kerL = 0 (trivially true when kerL is zero / has no columns).
        SparseMatrix3 B(MATRIX_FILE);
        BITMATRIX BkerL;
        multiply(B, kerL, BkerL);
        CPPUNIT_ASSERT(BkerL.isZero());
    }
};

const char* const BlockLanczosTest::MATRIX_FILE = "/tmp/blut_matrix.dat";

int main()
{
    CppUnit::TextUi::TestRunner runner;
    runner.addTest(BlockLanczosTest::suite());
    runner.run();
    return 0;
}
