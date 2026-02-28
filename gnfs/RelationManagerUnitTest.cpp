#include "RelationManager.h"
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/ui/text/TestRunner.h>

namespace
{
// Helper to populate a RelationTable with a single relation containing the given prime indices
void addRelation(RelationTable& rt, int id, const std::vector<int>& primes)
{
    Relation r(id);
    r.primes_index_ = rt.next_prime_index();
    r.prime_count_ = static_cast<unsigned char>(primes.size());
    for (size_t i = 0; i < primes.size(); ++i)
    {
        rt.set_prime(r.primes_index_, static_cast<unsigned char>(i), primes[i]);
    }
    rt[id] = r;
    rt.increment_next_prime(primes.size());
}
} // namespace

class RelationManagerTest : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE(RelationManagerTest);
    CPPUNIT_TEST(testPrimeFrequencyTable);
    CPPUNIT_TEST(testPrimeFrequencyTableDecrementAndReset);
    CPPUNIT_TEST(testRelationManagerConstruction);
    CPPUNIT_TEST(testRemoveSingletons);
    CPPUNIT_TEST(testBuildPrimeRelationSetMap);
    CPPUNIT_TEST(testMerge);
    CPPUNIT_TEST_SUITE_END();

public:
    void setUp()
    {
    }

    void tearDown()
    {
    }

    void testPrimeFrequencyTable()
    {
        PrimeFrequencyTable ft(10);
        CPPUNIT_ASSERT(ft.prime_count() == 0);
        CPPUNIT_ASSERT(ft.capacity() == 0);
        CPPUNIT_ASSERT(ft.check());

        ft.add_new_prime(0);
        ft.add_new_prime(1);
        CPPUNIT_ASSERT(ft.capacity() == 2);
        CPPUNIT_ASSERT(ft.prime_count() == 0);
        CPPUNIT_ASSERT(ft.frequency(0) == 0);
        CPPUNIT_ASSERT(ft.frequency(1) == 0);

        ft.increment_prime(0);
        CPPUNIT_ASSERT(ft.prime_count() == 1);
        CPPUNIT_ASSERT(ft.frequency(0) == 1);

        ft.increment_prime(0);
        CPPUNIT_ASSERT(ft.prime_count() == 1);
        CPPUNIT_ASSERT(ft.frequency(0) == 2);

        ft.increment_prime(1);
        CPPUNIT_ASSERT(ft.prime_count() == 2);
        CPPUNIT_ASSERT(ft.frequency(1) == 1);
        CPPUNIT_ASSERT(ft.check());

        ft.decrement_prime(1);
        CPPUNIT_ASSERT(ft.prime_count() == 1);
        CPPUNIT_ASSERT(ft.frequency(1) == 0);
        CPPUNIT_ASSERT(ft.check());
    }

    void testPrimeFrequencyTableDecrementAndReset()
    {
        PrimeFrequencyTable ft(10);
        ft.add_new_prime(0);
        ft.add_new_prime(1);
        ft.add_new_prime(2);
        ft.increment_prime(0);
        ft.increment_prime(1);
        ft.increment_prime(2);
        CPPUNIT_ASSERT(ft.prime_count() == 3);
        CPPUNIT_ASSERT(ft.check());

        ft.decrement_prime(0);
        CPPUNIT_ASSERT(ft.prime_count() == 2);
        CPPUNIT_ASSERT(ft.frequency(0) == 0);
        CPPUNIT_ASSERT(ft.check());

        // Test reset: resets to prime_index_ = 3, all frequencies zero
        ft.reset(3);
        CPPUNIT_ASSERT(ft.prime_count() == 0);
        CPPUNIT_ASSERT(ft.capacity() == 3);
        CPPUNIT_ASSERT(ft.check());

        // Add frequency counts after reset
        ft.increment_prime(0);
        ft.increment_prime(1);
        CPPUNIT_ASSERT(ft.prime_count() == 2);
        CPPUNIT_ASSERT(ft.check());
    }

    void testRelationManagerConstruction()
    {
        // Set up 3 relations with 4 distinct primes (indices 0..3):
        //   Relation 0: primes [0, 1]
        //   Relation 1: primes [1, 2]
        //   Relation 2: primes [0, 2, 3]
        // Prime frequencies: 0->2, 1->2, 2->2, 3->1

        PrimeFrequencyTable ft(10);
        ft.add_new_prime(0);
        ft.add_new_prime(1);
        ft.add_new_prime(2);
        ft.add_new_prime(3);
        ft.increment_prime(0); // relation 0
        ft.increment_prime(1); // relation 0
        ft.increment_prime(1); // relation 1
        ft.increment_prime(2); // relation 1
        ft.increment_prime(0); // relation 2
        ft.increment_prime(2); // relation 2
        ft.increment_prime(3); // relation 2
        CPPUNIT_ASSERT(ft.check());

        RelationTable rt(10, 100);
        addRelation(rt, 0, {0, 1});
        addRelation(rt, 1, {1, 2});
        addRelation(rt, 2, {0, 2, 3});
        rt.set_size(3);

        RelationManager rm(rt, ft, 1);

        CPPUNIT_ASSERT(rm.prime_weight(0) == 2);
        CPPUNIT_ASSERT(rm.prime_weight(1) == 2);
        CPPUNIT_ASSERT(rm.prime_weight(2) == 3);

        // combined_weight = weight(rs1) + weight(rs2) - 1
        CPPUNIT_ASSERT(rm.combined_weight(0, 1) == 3); // 2 + 2 - 1
        CPPUNIT_ASSERT(rm.combined_weight(0, 2) == 4); // 2 + 3 - 1
        CPPUNIT_ASSERT(rm.combined_weight(1, 2) == 4); // 2 + 3 - 1
    }

    void testRemoveSingletons()
    {
        // 3 relations, 4 primes, designed to cascade-remove as singletons:
        //   Relation 0: primes [0, 1]  (prime 0 freq=2, prime 1 freq=2)
        //   Relation 1: primes [1, 2]  (prime 2 freq=2)
        //   Relation 2: primes [2, 3]  (prime 3 freq=1 -> singleton)
        //
        // Removing R2 drops prime 2 to freq 1 -> R1 becomes singleton.
        // Removing R1 drops prime 1 to freq 1 -> R0 becomes singleton.
        // Removing R0 -> all removed.

        PrimeFrequencyTable ft(10);
        ft.add_new_prime(0);
        ft.add_new_prime(1);
        ft.add_new_prime(2);
        ft.add_new_prime(3);
        ft.increment_prime(0); // R0
        ft.increment_prime(1); // R0
        ft.increment_prime(1); // R1
        ft.increment_prime(2); // R1
        ft.increment_prime(2); // R2
        ft.increment_prime(3); // R2
        CPPUNIT_ASSERT(ft.check());
        CPPUNIT_ASSERT(ft.prime_count() == 4);

        RelationTable rt(10, 100);
        addRelation(rt, 0, {0, 1});
        addRelation(rt, 1, {1, 2});
        addRelation(rt, 2, {2, 3});
        rt.set_size(3);

        RelationManager rm(rt, ft, 1);

        CPPUNIT_ASSERT(rm.prime_weight(0) == 2);
        CPPUNIT_ASSERT(rm.prime_weight(1) == 2);
        CPPUNIT_ASSERT(rm.prime_weight(2) == 2);

        rm.remove_singletons();

        // All relation sets should be removed (weight 0)
        CPPUNIT_ASSERT(rm.prime_weight(0) == 0);
        CPPUNIT_ASSERT(rm.prime_weight(1) == 0);
        CPPUNIT_ASSERT(rm.prime_weight(2) == 0);

        // All prime frequencies should be 0
        CPPUNIT_ASSERT(ft.frequency(0) == 0);
        CPPUNIT_ASSERT(ft.frequency(1) == 0);
        CPPUNIT_ASSERT(ft.frequency(2) == 0);
        CPPUNIT_ASSERT(ft.frequency(3) == 0);
        CPPUNIT_ASSERT(ft.prime_count() == 0);
    }

    void testBuildPrimeRelationSetMap()
    {
        // 3 relations, 3 primes (indices 0..2), no singletons:
        //   Relation 0: primes [0, 1]
        //   Relation 1: primes [1, 2]
        //   Relation 2: primes [0, 2]
        // After build_prime_relation_set_map:
        //   prime 0 -> relation sets {0, 2}
        //   prime 1 -> relation sets {0, 1}
        //   prime 2 -> relation sets {1, 2}

        PrimeFrequencyTable ft(10);
        ft.add_new_prime(0);
        ft.add_new_prime(1);
        ft.add_new_prime(2);
        ft.increment_prime(0); // R0
        ft.increment_prime(1); // R0
        ft.increment_prime(1); // R1
        ft.increment_prime(2); // R1
        ft.increment_prime(0); // R2
        ft.increment_prime(2); // R2
        CPPUNIT_ASSERT(ft.check());

        RelationTable rt(10, 100);
        addRelation(rt, 0, {0, 1});
        addRelation(rt, 1, {1, 2});
        addRelation(rt, 2, {0, 2});
        rt.set_size(3);

        RelationManager rm(rt, ft, 1);
        rm.build_prime_relation_set_map(0);

        std::vector<long int> rs;

        CPPUNIT_ASSERT(rm.sets_including_prime(0, rs) == 2);
        CPPUNIT_ASSERT(rs.size() == 2);
        CPPUNIT_ASSERT(std::find(rs.begin(), rs.end(), 0L) != rs.end());
        CPPUNIT_ASSERT(std::find(rs.begin(), rs.end(), 2L) != rs.end());

        CPPUNIT_ASSERT(rm.sets_including_prime(1, rs) == 2);
        CPPUNIT_ASSERT(rs.size() == 2);
        CPPUNIT_ASSERT(std::find(rs.begin(), rs.end(), 0L) != rs.end());
        CPPUNIT_ASSERT(std::find(rs.begin(), rs.end(), 1L) != rs.end());

        CPPUNIT_ASSERT(rm.sets_including_prime(2, rs) == 2);
        CPPUNIT_ASSERT(rs.size() == 2);
        CPPUNIT_ASSERT(std::find(rs.begin(), rs.end(), 1L) != rs.end());
        CPPUNIT_ASSERT(std::find(rs.begin(), rs.end(), 2L) != rs.end());
    }

    void testMerge()
    {
        // 4 relations, 5 primes (indices 0..4):
        //   Relation 0: primes [0, 1, 2]
        //   Relation 1: primes [1, 2, 3]
        //   Relation 2: primes [2, 3, 4]
        //   Relation 3: primes [0, 4]
        // Merge relation sets 0 and 1 (shared primes 1,2 cancel out, result has [0,3])

        PrimeFrequencyTable ft(10);
        ft.add_new_prime(0);
        ft.add_new_prime(1);
        ft.add_new_prime(2);
        ft.add_new_prime(3);
        ft.add_new_prime(4);
        ft.increment_prime(0); // R0
        ft.increment_prime(1); // R0
        ft.increment_prime(2); // R0
        ft.increment_prime(1); // R1
        ft.increment_prime(2); // R1
        ft.increment_prime(3); // R1
        ft.increment_prime(2); // R2
        ft.increment_prime(3); // R2
        ft.increment_prime(4); // R2
        ft.increment_prime(0); // R3
        ft.increment_prime(4); // R3
        CPPUNIT_ASSERT(ft.check());

        RelationTable rt(10, 100);
        addRelation(rt, 0, {0, 1, 2});
        addRelation(rt, 1, {1, 2, 3});
        addRelation(rt, 2, {2, 3, 4});
        addRelation(rt, 3, {0, 4});
        rt.set_size(4);

        RelationManager rm(rt, ft, 1);
        rm.build_prime_relation_set_map(0);

        // Verify initial prime weights
        CPPUNIT_ASSERT(rm.prime_weight(0) == 3);
        CPPUNIT_ASSERT(rm.prime_weight(1) == 3);
        CPPUNIT_ASSERT(rm.prime_weight(2) == 3);
        CPPUNIT_ASSERT(rm.prime_weight(3) == 2);

        // Merge rs0 and rs1 via prime 1 (symmetric diff: {0,1,2} XOR {1,2,3} = {0,3})
        long int new_rs = rm.merge(0, 1, 1);
        CPPUNIT_ASSERT(new_rs >= 0);

        // Merged set should have primes [0, 3] (symmetric difference)
        CPPUNIT_ASSERT(rm.prime_weight(new_rs) == 2);

        // After merge, frequency table incremented for primes in merged set
        // primes 0 and 3 each gained one count from the new relation set
        CPPUNIT_ASSERT(ft.frequency(0) == 3); // was 2 (R0, R3), +1 from merged set
        CPPUNIT_ASSERT(ft.frequency(3) == 3); // was 2 (R1, R2), +1 from merged set
    }
};

int main()
{
    CppUnit::TextUi::TestRunner runner;
    runner.addTest(RelationManagerTest::suite());
    runner.run();
    return 0;
}
