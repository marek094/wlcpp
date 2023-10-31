#include <iostream>
#include <vector>
#include <algorithm>
#include <ranges>
#include <chrono>
#include <filesystem>
#include <unordered_map>
#include <unordered_set>

#include "read_g6.hpp"
#include "read_labeled_txt.hpp"
#include "read_tree_txt.hpp"
#include "homcounts.hpp"
#include "wl.hpp"
#include "explain.hpp"

#include <omp.h>


// hash std::vector<ull> to ull

struct hash_vec {
    size_t operator()(std::vector<unsigned long long> const& vec) const {
        size_t seed = vec.size();
        for (auto &&i : vec) {
            seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};


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


    auto results = std::vector<std::vector<unsigned long long>>{};
    results.resize(graph_list.size());
    {
        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(graph_list.size()); i++) {
            results[i] = wl::compute_path_homvec(graph_list[i], graph_list[i].number_of_vertices());
        }
    }

    // homvec -> list of graphs
    auto eq_classes = std::unordered_map<std::vector<unsigned long long>, std::vector<int>, hash_vec>{};
    int i = 0;
    for (auto&& hom_counts : results) {
        eq_classes[hom_counts].push_back(i);
        i += 1;
    }

    std::cout << "Found " << eq_classes.size() << " equivalence classes.\n";

    for (auto &&[homvec, graphs] : eq_classes) {
        if (graphs.size() > 1) {
            auto graph_degree_classes = std::unordered_map<std::vector<unsigned long long>, std::vector<int>, hash_vec>{};
            for (int i : graphs) {
                auto dvec = graph_list[i].degree_vector();
                std::sort(dvec.begin(), dvec.end());
                graph_degree_classes[dvec].push_back(i);
            }

            for (auto&& [dvec, degree_class] : graph_degree_classes) {
                assert(degree_class.size() > 0);
                if (degree_class.size() <= 1) {
                    continue;
                }

                auto oss_name = std::ostringstream{};
                oss_name << "n" << graph_list.size() << "_c";
                for (int i : degree_class) {
                    oss_name << i << "+";
                }
                // del last +
                oss_name.seekp(-1, std::ios_base::end);
                oss_name << ".txt";
                auto out_path = std::filesystem::path("wild_outputs") / oss_name.str();

                
                std::cout << "Equivalence class of size " << degree_class.size() << ": " << out_path.string() << "\n";

                std::filesystem::create_directories(out_path.parent_path());
                {
                    auto outss = std::ofstream(out_path);
                    for (int i : degree_class) {
                        std::cout << "\t";
                        graph_list[i].print_python();
                        graph_list[i].print_as_txt_line(outss);
                    }
                }

                std::cout << '\n';
                for (int i : degree_class) {
                    std::cout << "\t";
                    for (auto deg : graph_list[i].degree_vector()) {
                        std::cout << deg << ' ';
                    }
                    std::cout << '\n';
                }
                std::cout << '\n';

            }
        }   
    }




    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::cout << "Elapsed time: " << elapsed_seconds.count() << "s\n";


    return 0;
}