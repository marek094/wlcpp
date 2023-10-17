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

    // std::string s = "";
    // // s += (char)126;
    // // s += (char)66;
    // // s += (char)63;
    // // s += (char)120;
    // s += (char)68;

    // s += (char)81;
    // s += (char)99;

    // std::cout << s << '\n';

    // auto g = wl::graph6_to_graph(s);
    
    // std::cout << "vertices: " << (int)g.number_of_vertices() << '\n';

    // for (auto &&edge : g.edges) {
    //     std::cout << "(" << (int)edge.first << " " << (int)edge.second << "), ";
    // }
    // std::cout << '\n';

    // return 0;

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


    int percent = std::max(graph_list.size() / 100, 1UL);


    // {
    //     // timer 
    //     auto start = std::chrono::high_resolution_clock::now();
    //     classes.clear();


    //     std::vector<std::vector<unsigned long long>> homvec_list(graph_list.size());
    //     #pragma omp parallel for
    //     for (int i = 0; i < graph_list.size(); i++) {
    //         homvec_list[i] = wl::compute_pathwidth_one_homvec(graph_list[i], graph_list[i].number_of_vertices()+4);
    //         if (i % percent == percent-1) {
    //             std::cout << "|";
    //             std::cout.flush();
    //         }
    //     }

    //     std::cout << "\nDone with computation";
    //     for (int i = 0; i < graph_list.size(); i++) {
    //         classes[wl::stringyfy_vector(homvec_list[i])].push_back(i);
    //     }

    //     auto end = std::chrono::high_resolution_clock::now();
    //     std::cout << "\n\t " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms";
    //     std::cout << "\n\t #HomCount classes: " << classes.size() << "\n\n\n";
    // }
    
    {
        auto rehash = std::unordered_map<std::string, unsigned long long>{};
        auto classes = std::unordered_map<std::string, std::vector<int>>{};
        auto start = std::chrono::high_resolution_clock::now();
        classes.clear();
        int i = 0;
        for (auto &&graph : graph_list) {
            auto vec = wl::colors_1(graph, rehash, graph.labels);
            classes[wl::stringyfy_vector(vec)].push_back(i);
            if (i % percent == percent-1) {
                std::cout << "|";
                std::cout.flush();
            }
            i += 1;
        }
        auto end = std::chrono::high_resolution_clock::now();
        std::cout << "\n\t " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms";
        std::cout << "\n\t #TW classes: " << classes.size() << "\n\n\n";
    }

    {
        auto rehash = wl::rehash_t{};
        auto start = std::chrono::high_resolution_clock::now();
        auto classes = std::unordered_map<std::string, std::vector<std::pair<int, wl::colvec_t>>>{};
        int i = 0;

        auto vecs = std::vector<wl::colvec_t>{};
        for (auto &&graph : graph_list) {
            auto vec = wl::colors_p1_05(graph, rehash, graph.labels);
            classes[wl::stringyfy_vector(vec)].emplace_back(i, vec);
            if (i % percent == percent-1) {
                std::cout << "|";
                std::cout.flush();
            }
            i += 1;
            vecs.push_back(vec);

            // std::cout << "\n\t";
            // std::cout << graph.labels.size() <<" -- "<< (int)graph.number_of_vertices() << "\t";
            // for (auto &&label : graph.labels) {
                // std::cout << (int)label << " ";
            // }

            // std::cout << "\n\t" << wl::stringyfy_vector(vec);
        }

        auto end = std::chrono::high_resolution_clock::now();
        std::cout << "\n\t " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms";
        std::cout << "\n\t #PW classes: " << classes.size() << "\n\n\n";

        // if (vecs.size() >= 2) {
        //     std::cout << "Explain:\n";
        //     std::cout << wl::explain(rehash, vecs[0], vecs[1]) << '\n';
        // }

    }

    return 0;
}