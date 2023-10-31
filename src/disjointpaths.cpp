#include <iostream>
#include <vector>
#include <algorithm>
#include <ranges>
#include <chrono>
#include <filesystem>

#include "read_g6.hpp"
#include "read_labeled_txt.hpp"
#include "read_tree_txt.hpp"
#include "homcounts.hpp"
#include "wl.hpp"
#include "explain.hpp"

#include <omp.h>


auto graph_check(wl::SmallGraph const& graph) {
    
    // for all pairs check number of disjoint paths
    for (auto i = 0; i < graph.n; ++i) {
        for (auto j = i + 1; j < graph.n; ++j) {
            
        }
    }


}






int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <path_to_graph6_file>\n";
        return 1;
    }
    
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
    
    for (auto &&graph : graph_list) {
        if (!graph_check(graph)) {
            std::cerr << "Graph is not binary!\n";
            return 1;
        }
    }

    auto end = std::chrono::steady_clock::now();
    std::cout << "Time elapsed: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms\n";

    return 0;
}