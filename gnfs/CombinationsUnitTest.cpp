#include "Combinations.h"
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/ui/text/TestRunner.h>

class CombinationsTest : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE(CombinationsTest);
    CPPUNIT_TEST(test);
    CPPUNIT_TEST(test_knuth);
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
        CPPUNIT_ASSERT_THROW_MESSAGE("", Combinations(2, 1), std::string);
        CPPUNIT_ASSERT_THROW_MESSAGE("", Combinations(0, 1), std::string);
        CPPUNIT_ASSERT_THROW_MESSAGE("", Combinations(1, 0), std::string);
        CPPUNIT_ASSERT_THROW_MESSAGE("", Combinations(0, 0), std::string);

        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", Combinations(2, 2));
        Combinations c0(2, 2);
        CPPUNIT_ASSERT_THROW_MESSAGE("", c0(2), std::string);
        CPPUNIT_ASSERT(c0.size() == 2);
        CPPUNIT_ASSERT(c0(0) == 0);
        CPPUNIT_ASSERT(c0(1) == 1);
        CPPUNIT_ASSERT(!c0.next());

        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", Combinations(1, 2));
        Combinations c1(1, 5);
        CPPUNIT_ASSERT_THROW_MESSAGE("", c1(1), std::string);
        CPPUNIT_ASSERT(c1.size() == 1);
        CPPUNIT_ASSERT(c1(0) == 0);
        CPPUNIT_ASSERT(c1.next());
        CPPUNIT_ASSERT(!c1.done());
        CPPUNIT_ASSERT(c1(0) == 1);
        CPPUNIT_ASSERT(c1.next());
        CPPUNIT_ASSERT(!c1.done());
        CPPUNIT_ASSERT(c1(0) == 2);
        CPPUNIT_ASSERT(c1.next());
        CPPUNIT_ASSERT(!c1.done());
        CPPUNIT_ASSERT(c1(0) == 3);
        CPPUNIT_ASSERT(c1.next());
        CPPUNIT_ASSERT(!c1.done());
        CPPUNIT_ASSERT(c1(0) == 4);
        CPPUNIT_ASSERT(!c1.next());
        CPPUNIT_ASSERT(c1.done());

        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", Combinations c2(4, 10));

        Combinations c3(4, 10);
        // 10 C 4 = 10! / (4! 6!) = 10.9.8.7 / 4.3.2 = 10.3.7 = 210
        CPPUNIT_ASSERT(c3.size() == 4);
        CPPUNIT_ASSERT(c3(0) == 0);
        CPPUNIT_ASSERT(c3(1) == 1);
        CPPUNIT_ASSERT(c3(2) == 2);
        CPPUNIT_ASSERT(c3(3) == 3);
        CPPUNIT_ASSERT_THROW_MESSAGE("", c1(4), std::string);

        CPPUNIT_ASSERT(c3.next());
        CPPUNIT_ASSERT(c3(0) == 0);
        CPPUNIT_ASSERT(c3(1) == 1);
        CPPUNIT_ASSERT(c3(2) == 2);
        CPPUNIT_ASSERT(c3(3) == 4);
        CPPUNIT_ASSERT(c3.next());
        CPPUNIT_ASSERT(c3.next());
        CPPUNIT_ASSERT(c3.next());
        CPPUNIT_ASSERT(c3.next());
        CPPUNIT_ASSERT(c3.next());
        CPPUNIT_ASSERT(c3(0) == 0);
        CPPUNIT_ASSERT(c3(1) == 1);
        CPPUNIT_ASSERT(c3(2) == 2);
        CPPUNIT_ASSERT(c3(3) == 9);
        CPPUNIT_ASSERT(c3.next());
        CPPUNIT_ASSERT(c3(0) == 0);
        CPPUNIT_ASSERT(c3(1) == 1);
        CPPUNIT_ASSERT(c3(2) == 3);
        CPPUNIT_ASSERT(c3(3) == 4);
        CPPUNIT_ASSERT(c3.next());
        CPPUNIT_ASSERT(c3(0) == 0);
        CPPUNIT_ASSERT(c3(1) == 1);
        CPPUNIT_ASSERT(c3(2) == 3);
        CPPUNIT_ASSERT(c3(3) == 5);
        CPPUNIT_ASSERT(c3.next());
        CPPUNIT_ASSERT(c3.next());
        CPPUNIT_ASSERT(c3.next());
        CPPUNIT_ASSERT(c3.next());
        CPPUNIT_ASSERT(c3(0) == 0);
        CPPUNIT_ASSERT(c3(1) == 1);
        CPPUNIT_ASSERT(c3(2) == 3);
        CPPUNIT_ASSERT(c3(3) == 9);
        CPPUNIT_ASSERT(c3.next());
        CPPUNIT_ASSERT(c3(0) == 0);
        CPPUNIT_ASSERT(c3(1) == 1);
        CPPUNIT_ASSERT(c3(2) == 4);
        CPPUNIT_ASSERT(c3(3) == 5);
        CPPUNIT_ASSERT(c3.next());
        CPPUNIT_ASSERT(c3.next());
        CPPUNIT_ASSERT(c3.next());
        CPPUNIT_ASSERT(c3.next());
        CPPUNIT_ASSERT(c3.next());
        CPPUNIT_ASSERT(c3.next());
        CPPUNIT_ASSERT(c3.next());
        CPPUNIT_ASSERT(c3.next());
        CPPUNIT_ASSERT(c3.next());
        CPPUNIT_ASSERT(c3(0) == 0);
        CPPUNIT_ASSERT(c3(1) == 1);
        CPPUNIT_ASSERT(c3(2) == 6);
        CPPUNIT_ASSERT(c3(3) == 7);
        CPPUNIT_ASSERT(c3.next());
        CPPUNIT_ASSERT(c3.next());
        CPPUNIT_ASSERT(c3.next());
        CPPUNIT_ASSERT(c3.next());
        CPPUNIT_ASSERT(c3.next());
        CPPUNIT_ASSERT(c3(0) == 0);
        CPPUNIT_ASSERT(c3(1) == 1);
        CPPUNIT_ASSERT(c3(2) == 8);
        CPPUNIT_ASSERT(c3(3) == 9);
        CPPUNIT_ASSERT(c3.next());
        CPPUNIT_ASSERT(c3(0) == 0);
        CPPUNIT_ASSERT(c3(1) == 2);
        CPPUNIT_ASSERT(c3(2) == 3);
        CPPUNIT_ASSERT(c3(3) == 4);
        CPPUNIT_ASSERT(c3.next());
        CPPUNIT_ASSERT(c3.next());
        CPPUNIT_ASSERT(c3.next());
        CPPUNIT_ASSERT(c3.next());
        CPPUNIT_ASSERT(c3.next());
        for (size_t i = 34; i < 210; ++i)
        {
            CPPUNIT_ASSERT(c3.next());
        }
        CPPUNIT_ASSERT(c3(0) == 6);
        CPPUNIT_ASSERT(c3(1) == 7);
        CPPUNIT_ASSERT(c3(2) == 8);
        CPPUNIT_ASSERT(c3(3) == 9);
        CPPUNIT_ASSERT(!c3.next());
    }

    void test_knuth()
    {
        CPPUNIT_ASSERT_THROW_MESSAGE("", KnuthCombinations(2, 1), std::string);
        CPPUNIT_ASSERT_THROW_MESSAGE("", KnuthCombinations(0, 1), std::string);
        CPPUNIT_ASSERT_THROW_MESSAGE("", KnuthCombinations(1, 0), std::string);
        CPPUNIT_ASSERT_THROW_MESSAGE("", KnuthCombinations(0, 0), std::string);

        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", KnuthCombinations(2, 2));
        KnuthCombinations c0(2, 2);
        CPPUNIT_ASSERT_THROW_MESSAGE("", c0(2), std::string);
        CPPUNIT_ASSERT(c0.size() == 2);
        CPPUNIT_ASSERT(c0(0) == 0);
        CPPUNIT_ASSERT(c0(1) == 1);
        CPPUNIT_ASSERT(!c0.next());
        CPPUNIT_ASSERT(c0.done());

        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", KnuthCombinations(1, 2));
        KnuthCombinations c1(1, 5);
        CPPUNIT_ASSERT_THROW_MESSAGE("", c1(1), std::string);
        CPPUNIT_ASSERT(c1.size() == 1);
        CPPUNIT_ASSERT(c1(0) == 0);
        CPPUNIT_ASSERT(c1.next());
        CPPUNIT_ASSERT(!c1.done());
        CPPUNIT_ASSERT(c1(0) == 1);
        CPPUNIT_ASSERT(c1.next());
        CPPUNIT_ASSERT(!c1.done());
        CPPUNIT_ASSERT(c1(0) == 2);
        CPPUNIT_ASSERT(c1.next());
        CPPUNIT_ASSERT(!c1.done());
        CPPUNIT_ASSERT(c1(0) == 3);
        CPPUNIT_ASSERT(c1.next());
        CPPUNIT_ASSERT(!c1.done());
        CPPUNIT_ASSERT(c1(0) == 4);
        CPPUNIT_ASSERT(!c1.next());
        CPPUNIT_ASSERT(c1.done());

        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", KnuthCombinations c2(4, 10));

        KnuthCombinations c3(4, 10);
        // 10 C 4 = 10! / (4! 6!) = 10.9.8.7 / 4.3.2 = 10.3.7 = 210
        CPPUNIT_ASSERT(c3.size() == 4);
        CPPUNIT_ASSERT(c3(0) == 0);
        CPPUNIT_ASSERT(c3(1) == 1);
        CPPUNIT_ASSERT(c3(2) == 2);
        CPPUNIT_ASSERT(c3(3) == 3);
        CPPUNIT_ASSERT_THROW_MESSAGE("", c1(4), std::string);

        CPPUNIT_ASSERT(c3.next());
        CPPUNIT_ASSERT(c3(0) == 0);
        CPPUNIT_ASSERT(c3(1) == 1);
        CPPUNIT_ASSERT(c3(2) == 2);
        CPPUNIT_ASSERT(c3(3) == 4);
        CPPUNIT_ASSERT(c3.next());
        CPPUNIT_ASSERT(c3.next());
        CPPUNIT_ASSERT(c3.next());
        CPPUNIT_ASSERT(c3.next());
        CPPUNIT_ASSERT(c3.next());
        CPPUNIT_ASSERT(c3(0) == 0);
        CPPUNIT_ASSERT(c3(1) == 1);
        CPPUNIT_ASSERT(c3(2) == 3);
        CPPUNIT_ASSERT(c3(3) == 5);
        CPPUNIT_ASSERT(c3.next());
        CPPUNIT_ASSERT(c3(0) == 0);
        CPPUNIT_ASSERT(c3(1) == 2);
        CPPUNIT_ASSERT(c3(2) == 3);
        CPPUNIT_ASSERT(c3(3) == 5);
        CPPUNIT_ASSERT(c3.next());
        CPPUNIT_ASSERT(c3(0) == 1);
        CPPUNIT_ASSERT(c3(1) == 2);
        CPPUNIT_ASSERT(c3(2) == 3);
        CPPUNIT_ASSERT(c3(3) == 5);
        CPPUNIT_ASSERT(c3.next());
        CPPUNIT_ASSERT(c3.next());
        CPPUNIT_ASSERT(c3.next());
        CPPUNIT_ASSERT(c3.next());
        CPPUNIT_ASSERT(c3(0) == 0);
        CPPUNIT_ASSERT(c3(1) == 3);
        CPPUNIT_ASSERT(c3(2) == 4);
        CPPUNIT_ASSERT(c3(3) == 5);
        CPPUNIT_ASSERT(c3.next());
        CPPUNIT_ASSERT(c3(0) == 1);
        CPPUNIT_ASSERT(c3(1) == 3);
        CPPUNIT_ASSERT(c3(2) == 4);
        CPPUNIT_ASSERT(c3(3) == 5);
        CPPUNIT_ASSERT(c3.next()); // 2 3 4 5
        CPPUNIT_ASSERT(c3.next()); // 0 1 2 6
        CPPUNIT_ASSERT(c3.next()); // 0 1 3 6
        CPPUNIT_ASSERT(c3.next()); // 0 2 3 6
        CPPUNIT_ASSERT(c3.next()); // 1 2 3 6
        CPPUNIT_ASSERT(c3.next()); // 0 1 4 6
        CPPUNIT_ASSERT(c3.next()); // 0 2 4 6
        CPPUNIT_ASSERT(c3.next()); // 1 2 4 6
        CPPUNIT_ASSERT(c3.next()); // 0 3 4 6
        CPPUNIT_ASSERT(c3(0) == 0);
        CPPUNIT_ASSERT(c3(1) == 3);
        CPPUNIT_ASSERT(c3(2) == 4);
        CPPUNIT_ASSERT(c3(3) == 6);
        CPPUNIT_ASSERT(c3.next()); // 1 3 4 6
        CPPUNIT_ASSERT(c3.next()); // 2 3 4 6
        CPPUNIT_ASSERT(c3.next()); // 0 1 5 6
        CPPUNIT_ASSERT(c3.next()); // 0 2 5 6
        CPPUNIT_ASSERT(c3.next()); // 1 2 5 6
        CPPUNIT_ASSERT(c3(0) == 1);
        CPPUNIT_ASSERT(c3(1) == 2);
        CPPUNIT_ASSERT(c3(2) == 5);
        CPPUNIT_ASSERT(c3(3) == 6);
        CPPUNIT_ASSERT(c3.next());
        CPPUNIT_ASSERT(c3(0) == 0);
        CPPUNIT_ASSERT(c3(1) == 3);
        CPPUNIT_ASSERT(c3(2) == 5);
        CPPUNIT_ASSERT(c3(3) == 6);
        CPPUNIT_ASSERT(c3.next()); // 1 3 5 6
        CPPUNIT_ASSERT(c3.next()); // 2 3 5 6
        CPPUNIT_ASSERT(c3.next()); // 0 4 5 6
        CPPUNIT_ASSERT(c3.next()); // 1 4 5 6
        CPPUNIT_ASSERT(c3.next()); // 2 4 5 6
        for (size_t i = 34; i < 210; ++i)
        {
            CPPUNIT_ASSERT(c3.next());
        }
        CPPUNIT_ASSERT(c3(0) == 6);
        CPPUNIT_ASSERT(c3(1) == 7);
        CPPUNIT_ASSERT(c3(2) == 8);
        CPPUNIT_ASSERT(c3(3) == 9);
        CPPUNIT_ASSERT(!c3.next());
    }
};

int main()
{
    CppUnit::TextUi::TestRunner runner;
    runner.addTest(CombinationsTest::suite());
    runner.run();

    return 0;
}

