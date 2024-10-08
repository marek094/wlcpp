#include <gtest/gtest.h>
#include "wl.hpp"


class TrieTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Code here will be called immediately after the constructor (right before each test).
    }

    void TearDown() override {
        // Code here will be called immediately after each test (right before the destructor).
    }

    // Objects declared here can be used by all tests in the test case for trie.
    wl::trie<int> trie_;
};


TEST_F(TrieTest, InsertTest) {
    std::vector<int> vec = {10, 20, 30};
    trie_.insert(vec);

    EXPECT_EQ(trie_.size(), 1);
    EXPECT_EQ(trie_.root.next.size(), 1);

    {
        auto it = trie_.root.next.find(10);
        EXPECT_NE(it, trie_.root.next.end());
    }

    auto&& [node, count] = trie_.root.next[10];

    EXPECT_EQ(count, 1);
    EXPECT_EQ(node.next.size(), 1);

    {
        auto it = node.next.find(20);
        EXPECT_NE(it, node.next.end());
    }
    
    auto&& [node2, count2] = node.next[20];

    EXPECT_EQ(count2, 1);
    EXPECT_EQ(node2.next.size(), 1);

    {
        auto it = node2.next.find(30);
        EXPECT_NE(it, node2.next.end());
    }

    auto&& [node3, count3] = node2.next[30];

    EXPECT_EQ(count3, 1);
    EXPECT_EQ(node3.next.size(), 0);


    // update
    trie_.insert({11});

    EXPECT_EQ(trie_.size(), 2);
    EXPECT_EQ(trie_.root.next.size(), 2);
    
    trie_.insert({11});
    trie_.insert({11});

    EXPECT_EQ(trie_.size(), 4);
    EXPECT_EQ(trie_.root.next.size(), 2);

    {
        auto it = trie_.root.next.find(11);
        EXPECT_NE(it, trie_.root.next.end());
    }

    auto&& [node4, count4] = trie_.root.next[11];

    EXPECT_EQ(count4, 3);
    EXPECT_EQ(node4.next.size(), 0);


    // check number of legs
    EXPECT_EQ(trie_.root.next[11].second, 3);
    EXPECT_EQ(trie_.root.next[11].first.size(), 0);
    EXPECT_EQ(trie_.count_legs(), 3);
}


TEST_F(TrieTest, MergeTest) {
    std::vector<int> vec = {10, 20, 30};
    trie_.insert(vec);

    auto trie2 = wl::trie<int>{};
    trie2.insert({10, 20, 30});
    trie2.insert({10, 20, 40});
    trie2.insert({10, 20, 50});

    trie_.merge(trie2);

    EXPECT_EQ(trie_.size(), 4);
    EXPECT_EQ(trie_.root.next.size(), 1);

    EXPECT_EQ(trie_.root.next[10].second, 4);
    EXPECT_EQ(trie_.root.next[10].first.size(), 4);
    EXPECT_EQ(trie_.count_legs(), 0);

    {
        auto it = trie_.root.next.find(10);
        EXPECT_NE(it, trie_.root.next.end());
    }

    auto&& [node, count] = trie_.root.next[10];

    EXPECT_EQ(count, 4);
    EXPECT_EQ(node.next.size(), 1);

    {
        auto it = node.next.find(20);
        EXPECT_NE(it, node.next.end());
    }
    
    auto&& [node2, count2] = node.next[20];

    EXPECT_EQ(count2, 4);
    EXPECT_EQ(node2.next.size(), 3);

}


TEST_F(TrieTest, IncTest) {
    trie_.insert({2, 1, 0, 0});
    trie_.insert({2, 1, 1, 0});
    trie_.insert({2, 1, 0, 0});

    auto trie2 = wl::trie<int>{};
    trie2.insert({3, 0});

    trie_.merge(trie2);
    trie_.inc();

    EXPECT_EQ(trie_.size(), 4);

    {
        auto it = trie_.root.next.find(3);
        EXPECT_NE(it, trie_.root.next.end());
    }

    auto&& [node, count] = trie_.root.next[3];

    EXPECT_EQ(count, 4);
    EXPECT_EQ(node.size(), 4);

    {
        auto it = node.next.find(2);
        EXPECT_NE(it, node.next.end());
    }

}


TEST_F(TrieTest, IncTestEmpty) {
    EXPECT_EQ(trie_.size(), 0);

    trie_.inc();

    EXPECT_EQ(trie_.size(), 1);
    // EXPECT_EQ(trie_.max_depth(), 1);

    trie_.inc();

    EXPECT_EQ(trie_.size(), 1);
    EXPECT_EQ(trie_.size() - trie_.count_legs() - 1, 0);
    EXPECT_EQ(trie_.root.next.size(), 1);

    {
        // std::cout << trie_.root.next.begin()->first << "\n";
        auto it = trie_.root.next.find(0);
        EXPECT_NE(it, trie_.root.next.end());
    }

    // EXPECT_EQ(trie_.max_depth(), 2);
    {
        auto trie2 = wl::trie<int>{};
        trie2.insert({0, 0});
        EXPECT_EQ(trie_, trie2);

        // EXPECT_EQ(trie2.max_depth(), 1);
    }

    
    trie_.inc();

    {
        auto trie3 = wl::trie<int>{};
        trie3.insert({0, 0, 0});
        EXPECT_EQ(trie_, trie3);
    }


    trie_.inc();

    {
        auto trie4 = wl::trie<int>{};
        trie4.insert({0, 0, 0, 0});
        EXPECT_EQ(trie_, trie4);
    }

}


// TEST_F(TrieTest, Hash) {

//     trie_.insert({2, 1, 0, 0});
//     trie_.insert({2, 1, 1, 0});

// }



int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}