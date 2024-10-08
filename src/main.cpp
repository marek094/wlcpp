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



using map_hashmap_t = std::map<std::string, ska::unordered_map<std::string, int>>;
using vecmap_hashmap_t = std::vector<map_hashmap_t>;

using namespace std::string_literals;

struct config_t {
    int plus = 0;
};


template<typename InpFunc, typename ThreadFunc>
void test_run(std::string name, InpFunc get_graph_list, ThreadFunc&& thread_func) {
    // timer 
    auto start = std::chrono::high_resolution_clock::now();

    auto par_classes = vecmap_hashmap_t{};

    // get omp parallel threads count 
    int num_threads = omp_get_max_threads();
    par_classes.resize(num_threads);

    {
        auto graph_list = get_graph_list();
        std::cout << "Read " << graph_list.size() << " graphs.\n";
        int percent = std::max(graph_list.size() / 100, 1UL);

        #pragma omp parallel for
        for (int i = 0; i < graph_list.size(); i++) {
            // std::cout << "X" << std::endl;
            thread_func(
                graph_list[i],
                par_classes[omp_get_thread_num()]
            );

            if (i % percent == percent-1) {
                std::cout << "|" << std::flush;
            }
        }
        
        std::cout << "\nDone with computation" << std::endl;
    }
    
    // aggregate classes 
    auto& classes = par_classes.front();
    {
        for (int i = 1; i < par_classes.size(); ++i) {
            for (auto &&[key, value] : par_classes[i]) {
                for (auto &&[key2, value2] : value) {
                    classes[key][key2] += value2;
                }
            }
            std::cout << "/" << std::flush;

            // make the space available in par_classes[i]
            par_classes[i] = map_hashmap_t{};
        }

        std::cout << "\nDone with aggregation" << std::endl;
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "\n\t " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms";
    // for (auto &&[key, value] : classes) {
    //     std::cout << "\n\t" << name << " classes " << key << ": " << value.size() << "";
    // }
    std::cout << "\n\t" << "Classes: " << classes.size() << "\n";
    std::cout << "\n\n\n";

}


// void trie_test() {
//     auto trie = trie<uint64_t>{};
//     trie.root.next[0] = {node{}, 1};
//     trie.root.next[1] = {node{}, 1};
//     trie.root.next[0].first.next[0] = {node{}, 1};
//     trie.root.next[0].first.next[1] = {node{}, 1};
//     trie.root.next[1].first.next[0] = {node{}, 1};
//     trie.root.next[1].first.next[1] = {node{}, 1};

//     auto trie2 = trie<uint64_t>{};
//     trie2.root.next[0] = {node{}, 1};
//     trie2.root.next[1] = {node{}, 1};
//     trie2.root.next[0].first.next[0] = {node{}, 1};
//     trie2.root.next[0].first.next[1] = {node{}, 1};
//     trie2.root.next[1].first.next[0] = {node{}, 1};
//     trie2.root.next[1].first.next[1] = {node{}, 1};
    
//     trie.merge(trie2);
// }



int main(int argc, char* argv[]) {
    using std::to_string;

    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <path_to_graph6_file>\n";
        return 1;
    }
    
    size_t threads = 1;
    if (argc >= 3) {
        omp_set_num_threads(std::atoi(argv[2]));
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


    {
        auto rehash = wl::rehash_t{};
        auto graph_list = get_graph_list();
        std::cout << "Read " << graph_list.size() << " graphs.\n";

        for (auto &&graph : graph_list) {
            auto colvec = wl::colors_1(graph, rehash);
        }

        std::cout << "colors_1 color classes " << rehash.size() << "\n";
    }


    {
        // auto path_rehash = wl::rehash_t{};
        auto graph_list = get_graph_list();
        std::cout << "Read " << graph_list.size() << " graphs.\n";

        auto colvecs = std::vector<wl::trie<uint64_t>>{};
        for (auto &&graph : graph_list) {
            auto colvec = wl::colors_path(graph);
            colvecs.emplace_back(std::move(colvec));
        }

        std::sort(colvecs.begin(), colvecs.end());
        colvecs.erase(std::unique(colvecs.begin(), colvecs.end()), colvecs.end());

        std::cout << "colors_path color classes " << colvecs.size() << "\n";
    }



    std::vector<std::pair<std::string, std::function<void(wl::SmallGraph&, map_hashmap_t&)>>> runs;


    // list of all runs
    for (auto&& [name, func] : runs) {
        std::cout << name << std::endl;
    }

    for (auto&& [name, func] : runs) {
        test_run(name, get_graph_list, func);
    }
    





    return 0;
}