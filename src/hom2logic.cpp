#include "read_g6.hpp"
#include "read_labeled_txt.hpp"
#include "read_tree_txt.hpp"
#include "homcounts.hpp"
#include "wl.hpp"
#include "explain.hpp"
#include "logic.hpp"

#include <omp.h>

#include <iostream>
#include <vector>
#include <algorithm>
#include <ranges>
#include <chrono>
#include <filesystem>
#include <unordered_map>
#include <unordered_set>
#include <ranges>
#include <string>
#include <random>



// hash std::vector<ull> to ull

// struct hash_vec {
//     size_t operator()(std::vector<unsigned long long> const& vec) const {
//         size_t seed = vec.size();
//         for (auto &&i : vec) {
//             seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2);
//         }
//         return seed;
//     }
// };


int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <path_to_graph6_file>\n";
        return 1;
    }
    
    if (argc >= 3) {
        omp_set_num_threads(std::atoi(argv[2]));
    }

    int rank = 0;
    if (argc >= 4) {
        rank = std::atoi(argv[3]);
    }
    
    size_t limit = -1;
    if (argc >= 5) {
        limit = std::atoi(argv[3]);
    }

    size_t skip = 0;
    if (argc >= 6) {
        skip = std::atoi(argv[4]);
    }


    std::filesystem::path file_path = argv[1];

    auto graph_list = [&]() {
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


    auto results = std::vector<std::vector<unsigned long long>>{};
    results.resize(graph_list.size());
    {
        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(graph_list.size()); i++) {
            results[i] = wl::compute_path_homvec(graph_list[i], graph_list[i].number_of_vertices()+1);
        }
    }


    // homvec -> list of graphs
    auto eq_classes_vec = std::vector<std::pair<std::vector<unsigned long long>, std::vector<int>>>{};
    {
        auto eq_classes = std::unordered_map<std::vector<unsigned long long>, std::vector<int>, wl::hash_vec>{};
        int i = 0;
        for (auto&& hom_counts : results) {
            eq_classes[hom_counts].push_back(i);
            i += 1;
        }

        std::cout << "Found " << eq_classes.size() << " equivalence classes.\n";

        eq_classes_vec.reserve(eq_classes.size());
        for (auto &&[homvec, graph_idcs] : eq_classes) {
            eq_classes_vec.emplace_back(std::move(homvec), std::move(graph_idcs));
        }
    }



    using QuantifierT = wl::LogicQuantifierBound;
    
    auto graph_types = std::vector< std::unordered_set< wl::LogicFormula<QuantifierT>>>{};
    {
        graph_types.resize(graph_list.size());

        #pragma omp parallel for
        for (size_t i = 0; i < graph_list.size(); i++) {
            graph_types[i] = wl::compute_type_set<QuantifierT>(graph_list[i], rank);
        }
    }


    {
        bool all_equal = true;

        for (auto &&[homvec, graph_idcs] : eq_classes_vec) {
            
            for (int i = 0; i < graph_idcs.size() && all_equal; ++i) {
                auto type1 = graph_types[graph_idcs[i]];
                for (int j = i+1; j < graph_idcs.size() && all_equal; ++j) {
                    auto type2 = graph_types[graph_idcs[j]];
                    if (type1 != type2) {
                        all_equal = false;

                        // std::cout << "Found non-equal types for graphs " << graph_idcs[i] << " and " << graph_idcs[j] << ".\n";
                        // graph_list[graph_idcs[i]].print_as_txt_line(std::cout);
                        // graph_list[graph_idcs[j]].print_as_txt_line(std::cout);

                        // for (auto &&fla : type1) {
                        //     if (type2.find(fla) == type2.end()) {
                        //         std::cout << fla << "\n";
                        //     }
                        // }

                        // std::cout << "-\n";
                        
                        // for (auto &&fla : type2) {
                        //     if (type1.find(fla) == type1.end()) {
                        //         std::cout << fla << "\n";
                        //     }
                        // }

                        // std::cout << "\n";
                    }
                }
            }

            if (! all_equal) break;
        }
        // std::cout << "\n";

        if (all_equal) {
            std::cout << "All rank-" << rank << " 0-types of graphs " << file_path.filename() << " in hom(P, -) eq. classes are equal.\n";
        } else {
            std::cout << "Not all types are equal.\n";
        }
    }


    // return 0;

    
    size_t sqrt_m = static_cast<size_t>(std::sqrt(eq_classes_vec.size()));

    sqrt_m = std::min(eq_classes_vec.size()/2, sqrt_m);
    
    // std::random_device rd;
    // std::mt19937 g(rd());
    std::mt19937 g(0);

    auto indices1 = std::vector<size_t>(sqrt_m);
    auto indices2 = std::vector<size_t>(sqrt_m);
    
    {
        auto range = std::ranges::iota_view(0UL, eq_classes_vec.size()) | std::views::common;
        auto indices = std::vector<size_t>(range.begin(), range.end());
        std::ranges::shuffle(indices, g);
        std::copy(indices.begin(), indices.begin()+sqrt_m, indices1.begin());
        std::copy(indices.begin()+sqrt_m, indices.begin()+2*sqrt_m, indices2.begin());
    }

    size_t all_count = 0;
    size_t all_nonhypothesis = 0;

    #pragma omp parallel for
    for (int cl1 : indices1) {

        size_t count = 0;
        size_t nonhypothesis = 0;
        
        const auto class1 = eq_classes_vec[cl1];
        for (auto&& idx1 : class1.second) {
            auto&& graph1 = graph_list[idx1];
            auto logic_tp1 = graph_types[idx1];

            for (int cl2 : indices2) {
                assert(cl1 != cl2);
                const auto class2 = eq_classes_vec[cl2];
                count += class1.second.size() * class2.second.size();

                    for (auto&& idx2 : class2.second) {
                        auto&& graph2 = graph_list[idx2];
                        auto logic_tp2 = graph_types[idx2];
                        if (logic_tp1 == logic_tp2) {
                            nonhypothesis += 1;
                            std::cout << "Found non-hypothesis pair: " << idx1 << " and " << idx2 << "\n";
                            // graph1.print_as_txt_line(std::cout);
                            // graph2.print_as_txt_line(std::cout);
                            // std::cout << "\n";
                            for (auto &&fla : logic_tp1) {
                                if (logic_tp2.find(fla) == logic_tp2.end()) {
                                    std::cout << fla << "\n";
                                }
                            }
                            std::cout << "-\n";
                        }
                    }
            }
        }

        all_count += count;
        all_nonhypothesis += nonhypothesis;
    }


    std::cout << "Found " << all_nonhypothesis << " non-hypothesis pairs out of " << all_count << " pairs.\n";




    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::cout << "Elapsed time: " << elapsed_seconds.count() << "s\n";


    return 0;
}