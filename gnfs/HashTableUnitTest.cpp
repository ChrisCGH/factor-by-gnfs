#include <cstdlib>
#include "HashTable.h"
#include <string>
#include <sstream>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/ui/text/TestRunner.h>

namespace
{
struct StringHasher
{
    static unsigned long int hash(const std::string& key)
    {
        unsigned long int h = 0;
        const char* c = key.c_str();
        while (c && *c)
        {
            h += static_cast<unsigned int>(*c);
            h <<= 1;
            ++c;
        }
        return h;
    }
};
struct IntHasher
{
    static unsigned long int hash(int key)
    {
        return key;
    }
};
};

class HashTableTest : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE(HashTableTest);
    CPPUNIT_TEST(test1);
    CPPUNIT_TEST_SUITE_END();

public:
    void setUp()
    {
    }

    void tearDown()
    {
    }

    void test1()
    {
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", (HashTable<std::string, std::string, StringHasher>()));
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", (HashTable<std::string, int, StringHasher>()));
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", (HashTable<int, int, IntHasher>()));
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", (HashTable<std::string, std::string, StringHasher, 1009>()));
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", (HashTable<std::string, int, StringHasher, 1009>()));
        HashTable<std::string, int, StringHasher, 1009> h1;
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", h1["key1"]);
        CPPUNIT_ASSERT(h1["key1"] == 0);
        HashTable<std::string, std::string, StringHasher, 1009> h2;
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", h2["key2"]);
        CPPUNIT_ASSERT(h2["key2"] == "");
        HashTable<std::string, std::string, StringHasher, 1009> h3;
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", h3["key3"] = "value3");
        CPPUNIT_ASSERT(h3["key3"] == "value3");
        HashTable<std::string, int, StringHasher, 1009> h4;
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", h4["key4"] = 10001);
        CPPUNIT_ASSERT(h4["key4"] == 10001);
        HashTable<int, int, IntHasher, 1009> h5;
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", h5[1231234] = 10001);
        CPPUNIT_ASSERT(h5[1231234] == 10001);
        HashTable<int, int, IntHasher, 5> h6;
        for (int i = 0; i < 1000; ++i)
        {
            CPPUNIT_ASSERT_NO_THROW_MESSAGE("", h6[i] = i);
        }
        for (int i = 0; i < 1000; ++i)
        {
            CPPUNIT_ASSERT(h6[i] == i);
        }
        CPPUNIT_ASSERT(h6.find(0));
        CPPUNIT_ASSERT(h6.find(100));
        CPPUNIT_ASSERT(h6.find(874));
        CPPUNIT_ASSERT(!h6.find(-1));
        CPPUNIT_ASSERT(!h6.find(1000));
        CPPUNIT_ASSERT(!h6.find(1000000));

        HashTable<int, std::string, IntHasher, 7> h7;
        for (int i = 0; i < 1000; ++i)
        {
            std::ostringstream oss;
            oss << i;
            CPPUNIT_ASSERT_NO_THROW_MESSAGE("", h7[i] = oss.str());
        }
        for (int i = 0; i < 1000; ++i)
        {
            std::ostringstream oss;
            oss << i;
            CPPUNIT_ASSERT(h7[i] == oss.str());
        }
        CPPUNIT_ASSERT(h7.find(0));
        CPPUNIT_ASSERT(h7.find(100));
        CPPUNIT_ASSERT(h7.find(874));
        CPPUNIT_ASSERT(!h7.find(-1));
        CPPUNIT_ASSERT(!h7.find(1000));
        CPPUNIT_ASSERT(!h7.find(1000000));
        CPPUNIT_ASSERT_NO_THROW_MESSAGE("", h7.clear());
        CPPUNIT_ASSERT(!h7.find(0));
        CPPUNIT_ASSERT(!h7.find(100));
        CPPUNIT_ASSERT(!h7.find(874));

        HashTable<std::string, std::string, StringHasher, 7> h8;
        for (int i = 0; i < 1000; ++i)
        {
            std::ostringstream oss;
            oss << i;
            CPPUNIT_ASSERT_NO_THROW_MESSAGE("", h8[oss.str()] = oss.str());
        }
        for (int i = 0; i < 1000; ++i)
        {
            std::ostringstream oss;
            oss << i;
            CPPUNIT_ASSERT(h8[oss.str()] == oss.str());
        }
        CPPUNIT_ASSERT(h8.find("0"));
        CPPUNIT_ASSERT(h8.find("100"));
        CPPUNIT_ASSERT(h8.find("874"));
        CPPUNIT_ASSERT(!h8.find("-1"));
        CPPUNIT_ASSERT(!h8.find("1000"));
        CPPUNIT_ASSERT(!h8.find("1000000"));
        CPPUNIT_ASSERT(!h8.find("A random string"));
        CPPUNIT_ASSERT(h8.find("A random string") == h8.end());

        for (int i = 0; i < 1000; ++i)
        {
            std::ostringstream oss;
            oss << i;
            CPPUNIT_ASSERT(h8.find(oss.str()));
            CPPUNIT_ASSERT_NO_THROW_MESSAGE("", h8.remove(oss.str()));
            CPPUNIT_ASSERT(!h8.find(oss.str()));
        }
    }
};

int main()
{
    CppUnit::TextUi::TestRunner runner;
    runner.addTest(HashTableTest::suite());
    runner.run();

    return 0;
}
