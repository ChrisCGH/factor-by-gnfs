#include "MemoryMappedFile.h"
#ifndef WIN32
#include <unistd.h>
#endif
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/ui/text/TestRunner.h>

class MemoryMappedFileTest : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE(MemoryMappedFileTest);
    CPPUNIT_TEST(test);
    CPPUNIT_TEST_SUITE_END();
public:
    void setUp()
    {
    }

    void tearDown()
    {
        ::unlink("newfile");
    }

    void test()
    {
        CPPUNIT_ASSERT_THROW_MESSAGE("", MemoryMappedFile("nonexistentfile"), std::string);
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", MemoryMappedFile("newfile", 102400));
        MemoryMappedFile mmf1("newfile", 1024*1024*300);
        char* pos = mmf1.at_offset(1024*1024*33, 10);
        for (size_t i = 0; i < 9; ++i)
        {
            pos[i] = 'a' + i;
        }
        pos[9] = '\0';
        pos = mmf1.at_offset(1024, 100);
        for (size_t i = 0; i < 99; ++i)
        {
            pos[i] = 'a' + (i % 26);
        }
        pos[99] = '\0';

        pos = mmf1.at_offset(1024*1024*33, 10);
        CPPUNIT_ASSERT(std::string(pos) == "abcdefghi");
        pos = mmf1.at_offset(1024, 100);
        CPPUNIT_ASSERT(std::string(pos) == "abcdefghijklmnopqrstuvwxyzabcdefghijklmnopqrstuvwxyzabcdefghijklmnopqrstuvwxyzabcdefghijklmnopqrstu");

        mmf1.set_size(1024*1024*400);
        pos = mmf1.at_offset(1024*1024*33, 10);
        CPPUNIT_ASSERT(std::string(pos) == "abcdefghi");
        pos = mmf1.at_offset(1024, 100);
        CPPUNIT_ASSERT(std::string(pos) == "abcdefghijklmnopqrstuvwxyzabcdefghijklmnopqrstuvwxyzabcdefghijklmnopqrstuvwxyzabcdefghijklmnopqrstu");

        mmf1.set_size(1024*1024*100);
        pos = mmf1.at_offset(1024*1024*33, 10);
        CPPUNIT_ASSERT(std::string(pos) == "abcdefghi");
        pos = mmf1.at_offset(1024, 100);
        CPPUNIT_ASSERT(std::string(pos) == "abcdefghijklmnopqrstuvwxyzabcdefghijklmnopqrstuvwxyzabcdefghijklmnopqrstuvwxyzabcdefghijklmnopqrstu");
    }

};

int main()
{
    CppUnit::TextUi::TestRunner runner;
    runner.addTest(MemoryMappedFileTest::suite());
    runner.run();

    return 0;
}

