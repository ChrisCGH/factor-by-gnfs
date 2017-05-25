#include <iostream>
#include "BitOperations.h"
#ifdef USING_CPPUNIT
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/ui/text/TestRunner.h>
#else
#include "UnitTest.h"
#endif

#ifdef USING_CPPUNIT
class BitOperationsTest : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE( BitOperationsTest );
    CPPUNIT_TEST( test_colMask_32 );
    CPPUNIT_TEST( test_bitSet_32 );
    CPPUNIT_TEST( test_bitClear_32 );
    CPPUNIT_TEST( test_setBit_32 );
    CPPUNIT_TEST( test_clearBit_32 );
    CPPUNIT_TEST( test_copyBit_32 );
    CPPUNIT_TEST( test_toggleBit_32 );
    CPPUNIT_TEST( test_highestBitSet_32 );
    CPPUNIT_TEST( test_bitCount_32 );
    CPPUNIT_TEST( test_colMask_64 );
    CPPUNIT_TEST( test_bitSet_64 );
    CPPUNIT_TEST( test_bitClear_64 );
    CPPUNIT_TEST( test_setBit_64 );
    CPPUNIT_TEST( test_clearBit_64 );
    CPPUNIT_TEST( test_copyBit_64 );
    CPPUNIT_TEST( test_toggleBit_64 );
    CPPUNIT_TEST( test_highestBitSet_64 );
    CPPUNIT_TEST( test_bitCount_64 );
    CPPUNIT_TEST_SUITE_END();
private:
    uint32_t w;
    unsigned long long int w64;
public:
    void setUp()
    {
        w = 0UL;
        w64 = 0UL;
    }

    void tearDown()
    {
    }

    void test_colMask_32()
    {
        CPPUNIT_ASSERT(BitOperations::colMask(0) == 0);
        CPPUNIT_ASSERT(BitOperations::colMask(17) == 0x0001FFFF);
        CPPUNIT_ASSERT(BitOperations::colMask(35) == 0xFFFFFFFF);
    }

    void test_bitSet_32()
    {
        CPPUNIT_ASSERT(!BitOperations::bitSet(0, 2UL));
        CPPUNIT_ASSERT(BitOperations::bitSet(5, 34UL));
        CPPUNIT_ASSERT(!BitOperations::bitSet(32, 0xFFFFFFFFUL));
    }

    void test_bitClear_32()
    {
        CPPUNIT_ASSERT(BitOperations::bitClear(0, 2UL));
        CPPUNIT_ASSERT(!BitOperations::bitClear(5, 34UL));
        CPPUNIT_ASSERT(BitOperations::bitClear(32, 0xFFFFFFFFUL));
    }

    void test_setBit_32()
    {
        BitOperations::setBit(0, w);
        CPPUNIT_ASSERT(w == 1UL);
        BitOperations::setBit(5, w);
        CPPUNIT_ASSERT(w == 33UL);
        BitOperations::setBit(17, w);
        CPPUNIT_ASSERT(w == 0x20021UL);
        BitOperations::setBit(32, w);
        CPPUNIT_ASSERT(w == 0x20021UL);
    }

    void test_clearBit_32()
    {
        w = 0x20021UL;
        BitOperations::clearBit(5, w);
        CPPUNIT_ASSERT(w == 0x20001UL);
        BitOperations::clearBit(4, w);
        CPPUNIT_ASSERT(w == 0x20001UL);
        BitOperations::clearBit(32, w);
        CPPUNIT_ASSERT(w == 0x20001UL);
        BitOperations::clearBit(17, w);
        CPPUNIT_ASSERT(w == 0x1UL);
    }

    void test_copyBit_32()
    {
        w = 0x1UL;
        BitOperations::copyBit(1, 6, w);
        CPPUNIT_ASSERT(w == 0x41UL);
        BitOperations::copyBit(0, 0, w);
        CPPUNIT_ASSERT(w == 0x40UL);
        BitOperations::copyBit(0, 32, w);
        CPPUNIT_ASSERT(w == 0x40UL);
        BitOperations::copyBit(1, 35, w);
        CPPUNIT_ASSERT(w == 0x40UL);
    }

    void test_toggleBit_32()
    {
        w = 0x40UL;
        BitOperations::toggleBit(0, w);
        CPPUNIT_ASSERT(w == 0x41UL);
        BitOperations::toggleBit(6, w);
        CPPUNIT_ASSERT(w == 0x1UL);
        BitOperations::toggleBit(100, w);
        CPPUNIT_ASSERT(w == 0x1UL);
    }

    void test_highestBitSet_32()
    {
        CPPUNIT_ASSERT(BitOperations::highestBitSetInByte(23) == 4);
        CPPUNIT_ASSERT(BitOperations::highestBitSetInByte(255) == 7);
        CPPUNIT_ASSERT(BitOperations::highestBitSetInByte(0) == -1);

        CPPUNIT_ASSERT(BitOperations::highestBitSet(0x23F14UL) == 17);
        CPPUNIT_ASSERT(BitOperations::highestBitSet(0xFFFFFFFFUL) == 31);
        CPPUNIT_ASSERT(BitOperations::highestBitSet(0) == -1);
    }

    void test_bitCount_32()
    {
        CPPUNIT_ASSERT(BitOperations::bitCountInByte(0) == 0);
        CPPUNIT_ASSERT(BitOperations::bitCountInByte(0xFF) == 8);
        CPPUNIT_ASSERT(BitOperations::bitCountInByte(0x15) == 3);

        CPPUNIT_ASSERT(BitOperations::bitCount(0) == 0);
        CPPUNIT_ASSERT(BitOperations::bitCount(0xFFFFFFFF) == 32);
        CPPUNIT_ASSERT(BitOperations::bitCount(0x15236E14) == 13);
    }

    void test_colMask_64()
    {
        CPPUNIT_ASSERT(BitOperations64::colMask(0) == 0);
        CPPUNIT_ASSERT(BitOperations64::colMask(17) == 0x0001FFFF);
        CPPUNIT_ASSERT(BitOperations64::colMask(67) == 0xFFFFFFFFFFFFFFFFULL);
    }

    void test_bitSet_64()
    {
        CPPUNIT_ASSERT(!BitOperations64::bitSet(0, 2UL));
        CPPUNIT_ASSERT(BitOperations64::bitSet(5, 34UL));
        CPPUNIT_ASSERT(!BitOperations64::bitSet(64, 0xFFFFFFFFFFFFFFFFULL));
    }

    void test_bitClear_64()
    {
        CPPUNIT_ASSERT(BitOperations64::bitClear(0, 2UL));
        CPPUNIT_ASSERT(!BitOperations64::bitClear(5, 34UL));
        CPPUNIT_ASSERT(BitOperations64::bitClear(32, 0xFFFFFFFFUL));
    }

    void test_setBit_64()
    {
        BitOperations64::setBit(0, w64);
        CPPUNIT_ASSERT(w64 == 1UL);
        BitOperations64::setBit(5, w64);
        CPPUNIT_ASSERT(w64 == 33UL);
        BitOperations64::setBit(17, w64);
        CPPUNIT_ASSERT(w64 == 0x20021UL);
        BitOperations64::setBit(64, w64);
        CPPUNIT_ASSERT(w64 == 0x20021UL);
    }

    void test_clearBit_64()
    {
        w64 = 0x20021ULL;
        BitOperations64::clearBit(5, w64);
        CPPUNIT_ASSERT(w64 == 0x20001UL);
        BitOperations64::clearBit(4, w64);
        CPPUNIT_ASSERT(w64 == 0x20001UL);
        BitOperations64::clearBit(64, w64);
        CPPUNIT_ASSERT(w64 == 0x20001UL);
        BitOperations64::clearBit(17, w64);
        CPPUNIT_ASSERT(w64 == 0x1UL);
    }

    void test_copyBit_64()
    {
        w64 = 0x1ULL;
        BitOperations64::copyBit(1, 6, w64);
        CPPUNIT_ASSERT(w64 == 0x41UL);
        BitOperations64::copyBit(0, 0, w64);
        CPPUNIT_ASSERT(w64 == 0x40UL);
        BitOperations64::copyBit(0, 64, w64);
        CPPUNIT_ASSERT(w64 == 0x40UL);
        BitOperations64::copyBit(1, 67, w64);
        CPPUNIT_ASSERT(w64 == 0x40UL);
    }

    void test_toggleBit_64()
    {
        w64 = 0x40ULL;
        BitOperations64::toggleBit(0, w64);
        CPPUNIT_ASSERT(w64 == 0x41UL);
        BitOperations64::toggleBit(6, w64);
        CPPUNIT_ASSERT(w64 == 0x1UL);
        BitOperations64::toggleBit(100, w64);
        CPPUNIT_ASSERT(w64 == 0x1UL);
    }

    void test_highestBitSet_64()
    {
        CPPUNIT_ASSERT(BitOperations64::highestBitSet(0x23F14UL) == 17);
        CPPUNIT_ASSERT(BitOperations64::highestBitSet(0xFFFFFFFFUL) == 31);
        CPPUNIT_ASSERT(BitOperations64::highestBitSet(0) == -1);
    }

    void test_bitCount_64()
    {
        CPPUNIT_ASSERT(BitOperations64::bitCount(0) == 0);
        CPPUNIT_ASSERT(BitOperations64::bitCount(0xFFFFFFFF) == 32);
        CPPUNIT_ASSERT(BitOperations64::bitCount(0x15236E14) == 13);
        CPPUNIT_ASSERT(BitOperations64::bitCount(0xFFFFFFFFFFFFFFFFULL) == 64);
    }
};
#endif

int main()
{
#ifdef USING_CPPUNIT
    CppUnit::TextUi::TestRunner runner;
    runner.addTest(BitOperationsTest::suite());
    runner.run();
#else
    std::cout << "Testing BitOperations ..." << std::endl;
    UnitTest t;

    t.check(BitOperations::colMask(0) == 0, "BitOperations::colMask(0) should be 0");
    t.check(BitOperations::colMask(17) == 0x0001FFFF, "BitOperations::colMask(0) should be 0x0001FFFF");
    t.check(BitOperations::colMask(35) == 0xFFFFFFFF, "BitOperations::colMask(35) should be 0xFFFFFFFF");

    t.check(!BitOperations::bitSet(0, 2UL), "BitOperations::bitSet(0, 2UL) should be false");
    t.check(BitOperations::bitSet(5, 34UL), "BitOperations::bitSet(5, 18UL) should be true");
    t.check(!BitOperations::bitSet(32, 0xFFFFFFFFUL), "BitOperations::bitSet(32, 0xFFFFFFFFUL) should be false");

    t.check(BitOperations::bitClear(0, 2UL), "BitOperations::bitClear(0, 2UL) should be true");
    t.check(!BitOperations::bitClear(5, 34UL), "BitOperations::bitClear(5, 18UL) should be false");
    t.check(BitOperations::bitClear(32, 0xFFFFFFFFUL), "BitOperations::bitClear(32, 0xFFFFFFFFUL) should be true");

    uint32_t w(0UL);
    BitOperations::setBit(0, w);
    t.check(w == 1UL, "w should be 1UL");
    BitOperations::setBit(5, w);
    t.check(w == 33UL, "w should be 33UL");
    BitOperations::setBit(17, w);
    t.check(w == 0x20021UL, "w should be 0x20021UL");
    BitOperations::setBit(32, w);
    t.check(w == 0x20021UL, "w should be 0x20021UL");

    BitOperations::clearBit(5, w);
    t.check(w == 0x20001UL, "w should be 0x20001UL");
    BitOperations::clearBit(4, w);
    t.check(w == 0x20001UL, "w should be 0x20001UL");
    BitOperations::clearBit(32, w);
    t.check(w == 0x20001UL, "w should be 0x20001UL");
    BitOperations::clearBit(17, w);
    t.check(w == 0x1UL, "w should be 0x1UL");

    BitOperations::copyBit(1, 6, w);
    t.check(w == 0x41UL, "w should be 0x41UL");
    BitOperations::copyBit(0, 0, w);
    t.check(w == 0x40UL, "w should be 0x40UL");
    BitOperations::copyBit(0, 32, w);
    t.check(w == 0x40UL, "w should be 0x40UL");
    BitOperations::copyBit(1, 35, w);
    t.check(w == 0x40UL, "w should be 0x40UL");

    BitOperations::toggleBit(0, w);
    t.check(w == 0x41UL, "w should be 0x41UL");
    BitOperations::toggleBit(6, w);
    t.check(w == 0x1UL, "w should be 0x1UL");
    BitOperations::toggleBit(100, w);
    t.check(w == 0x1UL, "w should be 0x1UL");

    t.check(BitOperations::highestBitSetInByte(23) == 4, "BitOperations::highestBitSetInByte(23) should be 4");
    t.check(BitOperations::highestBitSetInByte(255) == 7, "BitOperations::highestBitSetInByte(255) should be 7");
    t.check(BitOperations::highestBitSetInByte(0) == -1, "BitOperations::highestBitSetInByte(0) should be -1");

    t.check(BitOperations::highestBitSet(0x23F14UL) == 17, "BitOperations::highestBitSet(0x23F14UL) should be 17");
    t.check(BitOperations::highestBitSet(0xFFFFFFFFUL) == 31, "BitOperations::highestBitSet(0xFFFFFFFFUL) should be 31");
    t.check(BitOperations::highestBitSet(0) == -1, "BitOperations::highestBitSet(0) should be -1");

    t.check(BitOperations::bitCountInByte(0) == 0, "BitOperations::bitCountInByte(0) should be 0");
    t.check(BitOperations::bitCountInByte(0xFF) == 8, "BitOperations::bitCountInByte(0xFF) should be 8");
    t.check(BitOperations::bitCountInByte(0x15) == 3, "BitOperations::bitCountInByte(0x15) should be 3");

    t.check(BitOperations::bitCount(0) == 0, "BitOperations::bitCount(0) should be 0");
    t.check(BitOperations::bitCount(0xFFFFFFFF) == 32, "BitOperations::bitCount(0xFFFFFFFF) should be 32");
    t.check(BitOperations::bitCount(0x15236E14) == 13, "BitOperations::bitCount(0x15236E14) should be 13");

    t.check(BitOperations64::colMask(0) == 0, "BitOperations64::colMask(0) should be 0");
    t.check(BitOperations64::colMask(17) == 0x0001FFFF, "BitOperations64::colMask(0) should be 0x0001FFFF");
    t.check(BitOperations64::colMask(67) == 0xFFFFFFFFFFFFFFFFULL, "BitOperations64::colMask(67) should be 0xFFFFFFFFFFFFFFFFULL");

    t.check(!BitOperations64::bitSet(0, 2UL), "BitOperations64::bitSet(0, 2UL) should be false");
    t.check(BitOperations64::bitSet(5, 34UL), "BitOperations64::bitSet(5, 18UL) should be true");
    t.check(!BitOperations64::bitSet(64, 0xFFFFFFFFFFFFFFFFULL), "BitOperations64::bitSet(64, 0xFFFFFFFFFFFFFFFFUL) should be false");

    t.check(BitOperations64::bitClear(0, 2UL), "BitOperations64::bitClear(0, 2UL) should be true");
    t.check(!BitOperations64::bitClear(5, 34UL), "BitOperations64::bitClear(5, 18UL) should be false");
    t.check(BitOperations64::bitClear(32, 0xFFFFFFFFUL), "BitOperations64::bitClear(32, 0xFFFFFFFFUL) should be true");

    unsigned long long int w64(0UL);
    BitOperations64::setBit(0, w64);
    t.check(w64 == 1UL, "w64 should be 1UL");
    BitOperations64::setBit(5, w64);
    t.check(w64 == 33UL, "w64 should be 33UL");
    BitOperations64::setBit(17, w64);
    t.check(w64 == 0x20021UL, "w64 should be 0x20021UL");
    BitOperations64::setBit(64, w64);
    t.check(w64 == 0x20021UL, "w64 should be 0x20021UL");

    BitOperations64::clearBit(5, w64);
    t.check(w64 == 0x20001UL, "w64 should be 0x20001UL");
    BitOperations64::clearBit(4, w64);
    t.check(w64 == 0x20001UL, "w64 should be 0x20001UL");
    BitOperations64::clearBit(64, w64);
    t.check(w64 == 0x20001UL, "w64 should be 0x20001UL");
    BitOperations64::clearBit(17, w64);
    t.check(w64 == 0x1UL, "w64 should be 0x1UL");

    BitOperations64::copyBit(1, 6, w64);
    t.check(w64 == 0x41UL, "w64 should be 0x41UL");
    BitOperations64::copyBit(0, 0, w64);
    t.check(w64 == 0x40UL, "w64 should be 0x40UL");
    BitOperations64::copyBit(0, 64, w64);
    t.check(w64 == 0x40UL, "w64 should be 0x40UL");
    BitOperations64::copyBit(1, 67, w64);
    t.check(w64 == 0x40UL, "w64 should be 0x40UL");

    BitOperations64::toggleBit(0, w64);
    t.check(w64 == 0x41UL, "w64 should be 0x41UL");
    BitOperations64::toggleBit(6, w64);
    t.check(w64 == 0x1UL, "w64 should be 0x1UL");
    BitOperations64::toggleBit(100, w64);
    t.check(w64 == 0x1UL, "w64 should be 0x1UL");
    t.check(BitOperations64::highestBitSet(0x23F14UL) == 17, "BitOperations64::highestBitSet(0x23F14UL) should be 17");
    t.check(BitOperations64::highestBitSet(0xFFFFFFFFUL) == 31, "BitOperations64::highestBitSet(0xFFFFFFFFUL) should be 31");
    t.check(BitOperations64::highestBitSet(0) == -1, "BitOperations64::highestBitSet(0) should be -1");
    t.check(BitOperations64::bitCount(0) == 0, "BitOperations64::bitCount(0) should be 0");
    t.check(BitOperations64::bitCount(0xFFFFFFFF) == 32, "BitOperations64::bitCount(0xFFFFFFFF) should be 32");
    t.check(BitOperations64::bitCount(0x15236E14) == 13, "BitOperations64::bitCount(0x15236E14) should be 13");
    t.check(BitOperations64::bitCount(0xFFFFFFFFFFFFFFFFULL) == 64, "BitOperations64::bitCount(0xFFFFFFFFFFFFFFFFULL) should be 64");

    t.test_summary();
#endif
}
