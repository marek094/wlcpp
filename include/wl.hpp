#pragma once


#include "unlabelled_graph.hpp"

#include <sstream>
#include <string>
#include <unordered_map>

namespace wl {

using rehash_t = std::unordered_map<std::string, unsigned long long>;
using colvec_t = std::vector<unsigned long long>;
// auto format = [](auto&&... elems) { return ((std::ostringstream{} << ... << elems)).str(); };

template <typename T>
auto stringyfy_vector(std::vector<T> const& vec) -> std::string {
    std::ostringstream ss;
    for (auto &&value : vec) {
        ss << value << ",";
    }
    return std::move(ss).str();
}

auto relabel(std::string const& col, rehash_t& rehash) -> int {
    if (rehash.find(col) == rehash.end()) {
        rehash[col] = rehash.size();
    }
    return rehash[col];
}

auto relabel(colvec_t const& col, rehash_t& rehash) -> int {
    return relabel(stringyfy_vector(col), rehash);
}

auto is_refined(colvec_t const& vec1, colvec_t const& vec2) -> bool {
    std::unordered_map<unsigned long long, int> rho1, rho2;
    colvec_t sizes1, sizes2;
    
    for (auto &&value : vec1) {
        rho1[value] += 1;
    }

    for (auto &&value : vec2) {
        rho2[value] += 1;
    }

    for (auto &&[key, value] : rho1) {
        sizes1.push_back(value);
    }

    for (auto &&[key, value] : rho2) {
        sizes2.push_back(value);
    }

    if (sizes1.size() != sizes2.size()) {
        return true;
    }

    std::sort(sizes1.begin(), sizes1.end());
    std::sort(sizes2.begin(), sizes2.end());

    if (sizes1 != sizes2) {
        return true;
    }

    return false;
}

auto is_refined(std::vector<colvec_t> const& vec2d1, std::vector<colvec_t> const& vec2d2) -> bool {
    colvec_t vec1, vec2;
    auto rehash = rehash_t{};

    vec1.reserve(vec2d1.size());
    for (auto &&vec : vec2d1) {
        vec1.emplace_back(relabel(stringyfy_vector(vec), rehash));
    }
    vec2.reserve(vec2d2.size());
    for (auto &&vec : vec2d2) {
        vec2.emplace_back(relabel(stringyfy_vector(vec), rehash));
    }

    return is_refined(vec1, vec2);
}




auto colors_1(SmallGraph const& graph, rehash_t& rehash, colvec_t const& labels = {}) -> colvec_t {
    auto A = graph.to_adjacency_matrix();
    auto n = A.rows();
    
    auto colvec = (labels.empty()) ? colvec_t(n, -1) : labels;

    for (int t = 0;; ++t) {
        auto colvec2 = colvec_t{};
        colvec2.reserve(n);

        // std::cout << "t=" << t << ": ";
        for (int i = 0; i < n; ++i) {
            std::vector<std::string> col;
            col.reserve(n);
            for (int j = 0; j < n; ++j) {
                std::ostringstream ss;
                ss << (i==j) << A(i, j) << "&";
                ss << colvec[j] << "&" << colvec[i];
                col.emplace_back(std::move(ss).str());
            }
            std::sort(col.begin(), col.end());
            colvec2.push_back(relabel(stringyfy_vector(col), rehash));
        }
        // std::cout << "" << stringyfy_vector(colvec2) << " \n";

        if (!is_refined(colvec, colvec2)) {
            std::sort(colvec2.begin(), colvec2.end());
            // std::cout << "t=" << t;
            return colvec2;
        }
        colvec = std::move(colvec2);
    }

    return {};
}



auto colors_p1_03(SmallGraph const& graph, rehash_t& rehash, colvec_t const& labels = {}) -> colvec_t {
    auto A = graph.to_adjacency_matrix();
    auto n = A.rows();

    auto labels_ = (labels.empty() ? colvec_t(n, -1) : labels);
    auto colvec = labels_;

    while (true) {
        auto colvec2 = colvec_t{};
        colvec2.reserve(n);
        for (int i = 0; i < n; ++i) {
            std::vector<std::string> col;
            col.reserve(n);
            for (int j = 0; j < n; ++j) {
                std::ostringstream ss;
                ss << (i==j) << A(i, j) << " & ";
                ss << colvec[j];
                col.emplace_back(std::move(ss).str());
            }
            std::sort(col.begin(), col.end());
            colvec2.push_back(relabel(stringyfy_vector(col), rehash));
        }

        if (!is_refined(colvec, colvec2)) {
            std::sort(colvec2.begin(), colvec2.end());
            return colvec2;
        }

        colvec = std::move(colvec2);
    }
}


auto colors_p1_04(SmallGraph const& graph, rehash_t& rehash, colvec_t const& labels = {}) -> colvec_t {
    auto A = graph.to_adjacency_matrix();
    auto n = A.rows();

    auto labels_ = (labels.empty() ? colvec_t(n, -1) : labels);


    // every vertex has set of formulas that models 
    auto colvec2d = std::vector<colvec_t>{};
    colvec2d.reserve(n);
    for (auto label : labels_) {
        colvec2d.emplace_back(colvec_t{label});
    }

    for (int t = 0;; ++t) {
        auto colvec2d2 = std::vector<colvec_t>{};
        colvec2d2.reserve(n);
        
        // std::cout << "t=" << t << ": ";

        for (int i = 0; i < n; ++i) {
            auto sort_counter = std::unordered_map<std::string, unsigned long long>{};
            for (int j = 0; j < n; ++j) {
                std::ostringstream atp;
                atp << (i==j) << A(i, j) << -labels_[i] << -labels_[j] << " &";
                assert (colvec2d.size() > j);
                for (auto &&elem : colvec2d[j]) {
                    auto atp_color = atp.str() + std::to_string(elem);
                    sort_counter[atp_color] += 1;
                }
            }

            auto col = colvec_t{};
            col.reserve(n);
            for (auto &&[key, value] : sort_counter) {
                auto grouped_key = key + "X" + std::to_string(value);
                col.emplace_back(relabel(grouped_key, rehash));
            }

            std::sort(col.begin(), col.end());
            // std::cout << "(" << stringyfy_vector(col) << ") ";
            colvec2d2.emplace_back(std::move(col));
        }

        // std::cout << "\n";

        if (!is_refined(colvec2d, colvec2d2)) {
            auto colvec2 = colvec_t{};
            colvec2.reserve(n);
            for (auto &&vec : colvec2d2) {
                colvec2.emplace_back(relabel(stringyfy_vector(vec), rehash));
            }
            std::sort(colvec2.begin(), colvec2.end());
            // std::cout << "t=" << t;
            return colvec2;
        }

        colvec2d = std::move(colvec2d2);
    }
}


auto colors_p1_05(SmallGraph const& graph, rehash_t& rehash, colvec_t const& labels = {}) -> colvec_t {
    auto A = graph.to_adjacency_matrix();
    auto n = A.rows();

    auto labels_ = (labels.empty() ? colvec_t(n, -1) : labels);

    // every vertex has set of formulas that models 
    auto colvec2d = std::vector<colvec_t>{};
    colvec2d.reserve(n);
    for (auto label : colvec_t(n, -1)) {
        colvec2d.emplace_back(colvec_t{label});
    }

    for (int t = 0;; ++t) {
        auto colvec2d2 = std::vector<colvec_t>{};
        colvec2d2.reserve(n);
        
        // std::cout << "t=" << t << ": ";

        


        if (!is_refined(colvec2d, colvec2d2)) {
            auto colvec2 = colvec_t{};
            colvec2.reserve(n);
            for (auto &&vec : colvec2d2) {
                colvec2.emplace_back(relabel(stringyfy_vector(vec), rehash));
            }
            std::sort(colvec2.begin(), colvec2.end());
            // std::cout << "t=" << t;
            return colvec2;
        }

        colvec2d = std::move(colvec2d2);
    }
}


template <typename Int=uint64_t>   
class trie {
    struct node {
        std::map<Int, std::pair<node, size_t>> next;

        auto is_leaf() -> bool {
            return next.empty();
        }

        auto operator<=>(node const& other) const {
            auto it = next.begin();
            auto it_other = other.next.begin();

            while (it != next.end() && it_other != other.next.end()) {
                if (it->first != it_other->first) {
                    return it->first <=> it_other->first;
                }
                if (it->second != it_other->second) {
                    return it->second <=> it_other->second;
                }
                ++it;
                ++it_other;
            }
            
            if (it == next.end() && it_other == other.next.end()) {
                return std::strong_ordering::equal;
            } else if (it == next.end()) {
                return std::strong_ordering::less;
            } else {
                return std::strong_ordering::greater;
            }
        }

        auto operator==(node const& other) const {
            return next == other.next;
        }
    };

public:
    trie() : root{} {}

    auto insert(std::vector<Int> const& vec) -> void {
        auto current = &root;
        for (Int elem : vec) {
            if (auto it = current->next.find(elem); it != current->next.end()) {
                auto& [node, count] = it->second;
                count += 1;
            } else {
                current->next[elem] = {node{}, 1};
            }
            current = &current->next[elem].first;
        }
    }

    friend auto merge_helper(node* current, node const& other) -> void {
        for (auto &&[key, value] : other.next) {
            if (auto it = current->next.find(key); it != current->next.end()) {
                it->second.second += value.second;
                merge_helper(&it->second.first, value.first);
            } else {
                current->next[key] = value;
            }
        }
    }

    auto merge(trie const& other) -> void {
        merge_helper(&root, other.root);
    }

    // friend auto inc(node* current) -> void {
    //     auto root = node{};
    //     root.next[0] = {*current, 1};
    //     return root;
    // }

    auto inc() -> void {
        auto was_root = std::move(root);
        root = node{};
        root.next[0] = {std::move(was_root), 1};
    }

    auto operator<=>(trie const& other) const {
        return root <=> other.root;
    }

    auto operator==(trie const& other) const -> bool {
        return root == other.root;
    }

    auto size() -> size_t {
        int sum = 0;
        for (auto &&[key, value] : root.next) {
            sum += value.second;
        }
        return sum;
    }

    node root;
};




auto colors_path(SmallGraph const& graph) -> trie<uint64_t> {

    auto n = graph.number_of_vertices();

    auto coltries = std::vector<trie<uint64_t>>{};
    coltries.resize(n);
    
    for (int t = 0;; ++t) {
        auto coltries2 = std::vector<trie<uint64_t>>{};
        coltries2.resize(n);

        for (int i = 0; i < n; ++i) {
            if (i >= graph.adj_list.size()) {
                continue;
            }

            auto&& list = graph.adj_list[i];
            auto sum = trie<uint64_t>{};
            for (int j : list) {
                auto trie_inc = coltries[j];
                trie_inc.inc();
                sum.merge(trie_inc);
            }
        }

        if (coltries == coltries2) {
            auto sum = trie<uint64_t>{};
            for (auto &&trie : coltries) {
                sum.merge(trie);
            }
            return sum;
        }

        coltries = std::move(coltries2);
    }

};




// auto colors_1(SmallGraph const& graph, rehash_t& rehash, colvec_t const& labels = {}) -> colvec_t {
//     auto A = graph.to_adjacency_matrix();
//     auto n = A.rows();
    
//     auto colvec = (labels.empty()) ? colvec_t(n, -1) : labels;

//     for (int t = 0;; ++t) {
//         auto colvec2 = colvec_t{};
//         colvec2.reserve(n);

//         // std::cout << "t=" << t << ": ";
//         for (int i = 0; i < n; ++i) {
//             std::vector<std::string> col;
//             col.reserve(n);
//             for (int j = 0; j < n; ++j) {
//                 std::ostringstream ss;
//                 ss << (i==j) << A(i, j) << "&";
//                 ss << colvec[j] << "&" << colvec[i];
//                 col.emplace_back(std::move(ss).str());
//             }
//             std::sort(col.begin(), col.end());
//             colvec2.push_back(relabel(stringyfy_vector(col), rehash));
//         }
//         // std::cout << "" << stringyfy_vector(colvec2) << " \n";

//         if (!is_refined(colvec, colvec2)) {
//             std::sort(colvec2.begin(), colvec2.end());
//             // std::cout << "t=" << t;
//             return colvec2;
//         }
//         colvec = std::move(colvec2);
//     }

//     return {};
// }





} // namespace wl
