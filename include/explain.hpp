#pragma once 

#include <string>
#include <vector>
#include <unordered_map>
#include <sstream>
#include <set>

#include "unlabelled_graph.hpp"
#include "wl.hpp"

namespace wl {

auto get_colvec(std::vector<std::string> const& col_to_str, unsigned long long col) -> colvec_t {
    auto colvec_str = col_to_str[col];
    std::istringstream ss{colvec_str};
    auto colvec = colvec_t{};
    
    unsigned long long value;
    char c; 
    while (ss >> value >> c) {
        assert(c == ',');
        colvec.push_back(value);
    }

    return colvec;
}


auto rec_(std::vector<std::string> const& col_to_str, unsigned long long col) -> std::string {
    if (col >= col_to_str.size()) {
        // return std::to_string((int)col);
        return ".";
    }

    auto colvec_str = col_to_str[col];
    std::istringstream ss{colvec_str};
    auto colvec = colvec_t{};
    std::string atp;
    char sep, sep2, sep3;
    unsigned long long count;
    unsigned long long next_col;

    std::stringstream fla;

    if (ss >> atp >> sep >> next_col >> sep2 >> count) {
        assert(sep == '&');
        assert(sep2 == 'X');
        // assert(sep3 == 'x' || sep3 == 'y');
        assert(ss.eof());

        fla << "E[" << count << "]" << atp << "& ";
        fla << rec_(col_to_str, next_col);
    } else {
        fla << "ERR " << colvec_str;
    }

    return std::move(fla).str();
}


auto build_pw_formula(rehash_t const& rehash, unsigned long long col1, unsigned long long col2) -> std::string {
    auto col_to_str = std::vector<std::string>{rehash.size()};
    for (auto &&[key, value] : rehash) {
        if (value < col_to_str.size()) {
            col_to_str[value] = key;
        }
    }

    auto colvec1 = get_colvec(col_to_str, col1);
    auto colvec2 = get_colvec(col_to_str, col2);
    
    // std::cout << col_to_str[col1] << '\n';
    // std::cout << col_to_str[col2] << '\n';

    auto colvec1_set = std::set<int>{colvec1.begin(), colvec1.end()};
    auto colvec2_set = std::set<int>{colvec2.begin(), colvec2.end()};
    
    // get element that is not in colvec2:
    unsigned long long colvec_diff = ~0;

    std::string result = "";

    for (auto &&elem : colvec1) {
        auto count = std::count(colvec1.begin(), colvec1.end(), elem);
        if (colvec2_set.find(elem) == colvec2_set.end()) {
            colvec_diff = elem;
            // break;
            auto x1 = rec_(col_to_str, colvec_diff);
            // count number of occurences of colvec_diff in colvec1
            result += "[" + std::to_string(count) + "] G_1 models " + x1 + " \t\tG_2 does not model " + x1 + "\n";
        }
    }

    for (auto &&elem : colvec2) {
        auto count = std::count(colvec2.begin(), colvec2.end(), elem);
        if (colvec1_set.find(elem) == colvec1_set.end()) {
            colvec_diff = elem;
            // break;
            auto x1 = rec_(col_to_str, colvec_diff);
            result += "[" + std::to_string(count) + "] G_1 does not model " + x1 + " \t\tG_2 models " + x1 + "\n";
        }
    }

    if (colvec_diff == ~0) {
        return "A:ERR";
    }

    return result;
}

auto explain(rehash_t& rehash, colvec_t const& colvec1, colvec_t const& colvec2) -> std::string {
    std::ostringstream ss;

    auto colvec2_set = std::set<int>{colvec2.begin(), colvec2.end()};
    auto colvec1_set = std::set<int>{colvec1.begin(), colvec1.end()};

    for (auto &&elem : colvec1_set) {
        auto count = std::count(colvec1.begin(), colvec1.end(), elem);
        if (colvec2_set.find(elem) == colvec2_set.end()) {
            ss << "[" << count << "]\n" << build_pw_formula(rehash, elem, colvec2[0]) << '\n';
        }
    }

    for (auto &&elem : colvec2_set) {
        auto count = std::count(colvec2.begin(), colvec2.end(), elem);
        if (colvec1_set.find(elem) == colvec1_set.end()) {
            ss << "[" << count << "]\n" << build_pw_formula(rehash, colvec1[0], elem) << '\n';
        }
    }

    return std::move(ss).str();
}


} // namespace wl