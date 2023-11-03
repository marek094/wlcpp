#include <iostream>
#include <vector>
#include <algorithm>
#include <ranges>
#include <string>
#include <chrono>
#include <filesystem>
#include <tuple>
#include <set>
#include <map>
#include <coroutine>
#include <compare>
#include <array>
#include <functional>
// #include <generator>
#include "ref_generator.hpp"

#include <omp.h>

#include "read_g6.hpp"
#include "read_labeled_txt.hpp"
#include "read_tree_txt.hpp"
#include "homcounts.hpp"
#include "wl.hpp"
#include "unlabelled_graph.hpp"
#include "logic.hpp"



int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <path_to_graph6_file>\n";
        return 1;
    }
    
    if (argc >= 3) {
        omp_set_num_threads(std::atoi(argv[2]));
    }
    
    int q_depth = 1;
    if (argc >= 4) {
        q_depth = std::atoi(argv[3]);
    }


    std::filesystem::path file_path = argv[1];

    auto graph_list = [&]() {
        auto limit = ~0ULL;
        auto skip = 0;
        if (file_path.extension() == ".txt" && file_path.stem().string().substr(0, 4) == "tree") {    
            return wl::read_graph_from_tree_txt_file(file_path, limit, skip);
        } else
        if (file_path.extension() == ".txt") {
            return wl::read_graph_from_labeled_txt_file(file_path, limit, skip);
        }
        assert(file_path.extension() == ".g6");
        return wl::read_graph_from_graph6_file(file_path, limit, skip);
    }();

    

    std::cout << "Read " << graph_list.size() << " graphs.\n";

    auto start = std::chrono::steady_clock::now();
    
    assert (graph_list.size() == 2);

    using QuantifierT = wl::LogicQuantifierBound;

    auto b1 = wl::check_formulas_on_graph<QuantifierT>(graph_list[0], q_depth);
    std::cout << "Graph 1 done.\n";

    // for (auto &&[args, result] : b1) {
    //     std::cout << std::get<0>(args) << " " << std::get<1>(args) << " " << result << "\n";
    // }

    auto b2 = wl::check_formulas_on_graph<QuantifierT>(graph_list[1], q_depth);
    std::cout << "Graph 2 done.\n";


    // types
    auto logic_type1 = std::vector<std::vector<wl::LogicFormula<QuantifierT>>>{}; // vector of formulas
    logic_type1.resize(graph_list[0].number_of_vertices());
    
    for (int i = 0; i < graph_list[0].number_of_vertices(); ++i) {
        for (auto &&[args, result] : b1) {
            if (! result) {
                continue;
            }
            auto [node_i, formula] = args;
            if (node_i == i) {
                logic_type1[i].push_back(formula);
            }
        }
    }

    auto logic_type2 = std::vector<std::vector<wl::LogicFormula<QuantifierT>>>{}; // vector of formulas
    logic_type2.resize(graph_list[1].number_of_vertices());

    for (int i = 0; i < graph_list[1].number_of_vertices(); ++i) {
        for (auto &&[args, result] : b2) {
            if (! result) {
                continue;
            }
            auto [node_i, formula] = args;
            if (node_i == i) {
                logic_type2[i].push_back(formula);
            }
        }
    }

    // auto sort_function = [](auto&& a, auto&& b) {
    //     if (a.quantfiers.size() < b.size()) {
    //         return true;
    //     }
    //     if (a.size() > b.size()) {
    //         return false;
    //     }

    //     for (int i = 0; i < a.size(); ++i) {
    //         if (a[i] < b[i]) {
    //             return true;
    //         }
    //         if (b[i] < a[i]) {
    //             return false;
    //         }
    //     }

    //     return false;
    // };

    int i = 0;
    for (auto&& formulas1 : logic_type1) {
        std::sort(formulas1.begin(), formulas1.end());
        int n_found = 0;
        for (auto&& formulas2 : logic_type2) {
            std::sort(formulas2.begin(), formulas2.end());
            if (formulas1 == formulas2) {
                n_found += 1;
            }
        }
        std::cout << "#" << i << " " << n_found << "\n";
        i += 1;
    }


    int diff_count = 0;
    for (auto &&[args1, result1] : b1) {
        auto [model1, formula] = args1;
        auto model2 = model1;
        auto args2 = std::make_tuple(model2, formula);
        assert (b2.find(args2) != b2.end());
        auto&& result2 = b2[args2];
        if (result1 != result2) {
            diff_count += 1;
            if (model1 == ~0ULL)
                std::cout << "DIFF: " << model1 << " " << model2 << " [] " << formula << " " << result1 << " " << result2 << "\n";
        }
    }

    std::cout << "DIFF COUNT: " << diff_count << "\n";


    std::cout << "\ndebug eval\n";


    std::cout << std::get<1>(b1.begin()->first) << "\n";



    auto end = std::chrono::steady_clock::now();
    std::cout << "Time elapsed: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms\n";

    return 0;
}