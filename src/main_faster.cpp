#include "read_g6.hpp"
#include "read_labeled_txt.hpp"
#include "read_tree_txt.hpp"
#include "homcounts.hpp"
#include "wl.hpp"
#include "explain.hpp"
#include "crtu64t.hpp"

#include <ska/unordered_map.hpp>
#include <omp.h>

#include <iostream>
#include <vector>
#include <algorithm>
#include <ranges>
#include <chrono>
#include <filesystem>
#include <complex>
#include <map>
#include <functional>
#include <concepts>


using map_hashmap_t = std::map<std::string, ska::unordered_map<std::string, int>>;
using vecmap_hashmap_t = std::vector<map_hashmap_t>;

using namespace std::string_literals;

struct config_t {
    int plus = 0;
};



template<std::three_way_comparable T>
auto unique_size(std::vector< T> vec) -> size_t {
    std::sort(vec.begin(), vec.end());
    auto proposed_end = std::unique(vec.begin(), vec.end());
    return std::distance(vec.begin(), proposed_end);
}



template<std::three_way_comparable T>
struct data_with_index {
    T data;
    size_t index;

    auto operator<=>(data_with_index const& other) const {
        return data <=> other.data;
    }

    auto operator==(data_with_index const& other) const -> bool {
        return data == other.data;
    }
};


auto partitions_to_map(std::vector<std::vector<size_t>> const& partitions) {
    auto result = std::unordered_map<size_t, size_t>{};
    for (size_t i = 0; i < partitions.size(); ++i) {
        for (auto&& elem : partitions[i]) {
            result[elem] = i;
        }
    }
    return result;
}


template<typename GetGraphs, typename GetColors>
auto compute_equivalence(std::string name, GetGraphs&& get_graphs, int threads, GetColors&& get_colors, bool quiet=false) {
    
    auto start = std::chrono::steady_clock::now();

    auto graph_list = get_graphs();

    if (!quiet) {
        std::cout << "read " << graph_list.size() << " graphs\n";
    }
    
    using return_t = decltype(get_colors(graph_list[0]));
    auto colvecs = std::vector<data_with_index<return_t>>{};
    colvecs.resize(graph_list.size());

    #pragma omp parallel for if (threads > 1)
    for (size_t i = 0; i < graph_list.size(); ++i) {
        colvecs[i].index = i;
        colvecs[i].data = get_colors(graph_list[i]);
        std::sort(colvecs[i].data.begin(), colvecs[i].data.end());
    }

    if (!quiet) {
        std::cout << "color classes " << unique_size(colvecs) << " (" << name << ") \n";

        std::cout << "\ttime elapsed ";
        std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::steady_clock::now() - start).count() << "ms\n\n";
    }
    
    std::sort(colvecs.begin(), colvecs.end());

    auto result = std::vector<std::vector<size_t>>{};
    for (size_t i = 0; i < colvecs.size(); ++i) {
        if (i == 0 || colvecs[i].data != colvecs[i-1].data) {
            result.emplace_back();
        }
        result.back().push_back(colvecs[i].index);
    }


    return result;
}


auto print_subclasses(
    std::string name,
    std::vector<std::vector<size_t>> const& eq_parts, 
    std::unordered_map<size_t, size_t> const& part2graph,
    std::vector<wl::SmallGraph> const& graph_list
) {
    for (auto const& part : eq_parts) {
        auto subclasses = std::unordered_map<size_t, std::vector<size_t>>{};
        for (auto idx : part) {
            subclasses[part2graph.at(idx)].emplace_back(idx);
        }

        if (subclasses.size() > 1) {
            std::cout << "# Subclasses (" << name << ") " << subclasses.size() << ":\n";
            // std::cout << "# ";
            for (auto const& [part, idxs] : subclasses) {
                // for (auto idx : idxs) {
                //     std::cout << idx << " ";
                // }
                for (auto idx : idxs) {
                    std::cout << "g" << idx << " = " << graph_list[idx] << "\n";
                }
                std::cout << "#-\n";
            }
            std::cout << "\n";
        }
    }
}



int main(int argc, char* argv[]) {
    using std::to_string;

    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <path_to_graph6_file>\n";
        return 1;
    }
    
    size_t threads = 1;
    if (argc >= 3) {
        threads = std::atoi(argv[2]);
        omp_set_num_threads(threads);
    }
    
    size_t limit = -1;
    if (argc >= 4) {
        limit = std::atoi(argv[3]);
    }

    size_t skip = 0;
    if (argc >= 5) {
        skip = std::atoi(argv[4]);
    }

    bool do_pathwith = false;
    if (argc >= 6) {
        do_pathwith = std::atoi(argv[5]);
    }

    std::string file_path_str = argv[1];

    auto get_graph_list = [=]() -> std::vector<wl::SmallGraph> {
        std::istringstream iss(file_path_str);
        auto graph_list = std::vector<wl::SmallGraph>{};
        for (std::string token; std::getline(iss, token, ' '); ) {
            std::filesystem::path file_path = token;
            auto r = std::vector<wl::SmallGraph>{};
            if (file_path.extension() == ".txt" && file_path.stem().string().substr(0, 4) == "tree") {    
                r = wl::read_graph_from_tree_txt_file(file_path, limit, skip);
            } else if (file_path.extension() == ".txt") {
                r = wl::read_graph_from_labeled_txt_file(file_path, limit, skip);
            } else if (file_path.extension() == ".g6") {
                r = wl::read_graph_from_graph6_file(file_path, limit, skip);
            } else if (file_path.extension() == ".is6") {
                r = wl::read_graph_from_sparse6_file(file_path, limit, skip);
            } else {
                std::cerr << "Unknown file extension: " << file_path.extension() << std::endl;
            }

            if (graph_list.empty()) {
                graph_list = std::move(r);
            } else {
                graph_list.insert(graph_list.end(), r.begin(), r.end());
            }
        }
        return graph_list;
    };


    std::cout << "threads: " << threads << "\n";


    // auto eq_parts_tree = [&](){
    //     auto rehash = wl::rehash_t{};
    //     return compute_equivalence("colors_1", get_graph_list, 1, [&rehash](auto&& graph) {
    //         return wl::colors_1(graph, rehash);
    //     });
    // }();


    // auto eq_parts_cat = compute_equivalence("colors_caterpillar", get_graph_list, threads, [](auto&& graph) {
    //     return wl::colors_caterpillar(graph);
    // });

    auto eq_parts_path = compute_equivalence("colors_path", get_graph_list, threads, [](auto&& graph) {
        return wl::colors_path(graph);
    });


    auto graph_list = get_graph_list();

    for (auto const& part : eq_parts_path) {
        if (part.size() <= 1) continue;

        auto subglgen = [&]() { 
            auto result = std::vector<wl::SmallGraph>{};
            for (auto idx : part) {
                result.emplace_back(graph_list[idx]);
            }
            return result;
        };

        auto eq_parts_cat = compute_equivalence("colors_caterpillar", subglgen, threads, [](auto&& graph) {
            return wl::colors_caterpillar(graph);
        }, true);

        auto rehash = wl::rehash_t{};
        {
            // speedup (?)
            auto eq_parts_tree = compute_equivalence("colors_1", subglgen, 1, [&rehash](auto&& graph) {
                return wl::colors_1(graph, rehash);
            }, true);

            if (eq_parts_tree.size() <= 1) continue;
        }
        auto const& crehash = rehash;


        #pragma omp parallel for if (threads > 1)
        for (auto const& part2 : eq_parts_cat) {
            if (part2.size() <= 1) continue;


            auto graph_list2 = std::vector<wl::SmallGraph>{};
            for (auto idx : part2) {
                graph_list2.emplace_back(graph_list[part[idx]]);
            }

            auto subglgen2 = [&]() {
                return graph_list2;
            };

            auto eq_parts_tree = compute_equivalence("colors_1", subglgen2, 1, [&crehash](auto&& graph) {
                return wl::colors_1(graph, crehash);
            }, true);

            if (eq_parts_tree.size() <= 1) continue;
            
            std::cerr << "eq_parts_tree.size() " << eq_parts_tree.size() << "\n";
            for (auto idx : part2) {
                std::cerr << "idx " << idx << " part[idx] " << part[idx] << "\n";
            }
            std::cerr << std::endl;


            auto const part2graph_tree = partitions_to_map(eq_parts_tree);  

            print_subclasses("cat < tree", eq_parts_cat, part2graph_tree, graph_list2);
        }
    }




    // auto const part2graph_tree = partitions_to_map(eq_parts_tree);
    // auto const part2graph_cat = partitions_to_map(eq_parts_cat);
    // auto const part2graph_path = partitions_to_map(eq_parts_path);    

    // auto graph_list = get_graph_list();

    // print_subclasses("cat < tree", eq_parts_cat, part2graph_tree, graph_list);
    // print_subclasses("path < cat", eq_parts_path, part2graph_cat, graph_list);



    return 0;
}