#include "read_g6.hpp"
#include "read_labeled_txt.hpp"
#include "read_tree_txt.hpp"
#include "homcounts.hpp"
#include "wl.hpp"
#include "explain.hpp"

#include <omp.h>

#include <iostream>
#include <vector>
#include <algorithm>
#include <ranges>
#include <chrono>
#include <filesystem>
#include <complex>
#include <map>

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

    auto get_graph_list = [=]() {
        if (file_path.extension() == ".txt" && file_path.stem().string().substr(0, 4) == "tree") {    
            return wl::read_graph_from_tree_txt_file(file_path, limit, skip);
        } else
        if (file_path.extension() == ".txt") {
            return wl::read_graph_from_labeled_txt_file(file_path, limit, skip);
        }
        assert(file_path.extension() == ".g6");
        return wl::read_graph_from_graph6_file(file_path, limit, skip);
    };

    int percent = 100000;
    int plus = 9;

    // {
    //     auto classes = std::unordered_map<std::string, std::vector<int>>{};
    //     auto classes_N1 = std::unordered_map<std::string, std::vector<int>>{};
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
    //         uint64_t n4 = graph_list[i].number_of_vertices()+4;
    //         classes_N1[std::to_string(homvec_list[i][n4])].push_back(i);
    //     }

    //     auto end = std::chrono::high_resolution_clock::now();
    //     std::cout << "\n\t " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms";
    //     std::cout << "\n\t #HomCount classes: " << classes.size() << "\n";
    //     std::cout << "\n\t #HomCount classes N+4: " << classes_N1.size() << "\n\n\n";
    // }
    
    {
        auto classes = std::map<std::string, std::unordered_map<std::string, int>>{};
        // timer 
        auto start = std::chrono::high_resolution_clock::now();
        classes.clear();

        std::vector<std::vector<unsigned long long>> homvec_list{};
        {
            auto graph_list = get_graph_list();
            std::cout << "Read " << graph_list.size() << " graphs.\n";
            percent = std::max(graph_list.size() / 100, 1UL);

            homvec_list.resize(graph_list.size());

            #pragma omp parallel for
            for (int i = 0; i < graph_list.size(); i++) {
                homvec_list[i] = wl::compute_path_homvec(graph_list[i], graph_list[i].number_of_vertices()+plus);
                if (i % percent == percent-1) {
                    std::cout << "|";
                    std::cout.flush();
                }
            }
            
            std::cout << "\nDone with computation";
        }
        
        
        for (int i = 0; i < homvec_list.size(); i++) {
            classes["   "][wl::stringyfy_vector(homvec_list[i])] += 1;
            uint64_t n = homvec_list[i].size() -plus-1;
            for (int j = 0; j <= plus; ++j) {
                classes["N+"+std::to_string(j)][std::to_string(homvec_list[i][n+j])] += 1;
            }
        }

        auto end = std::chrono::high_resolution_clock::now();
        std::cout << "\n\t " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms";
        for (auto &&[key, value] : classes) {
            std::cout << "\n\tcompute_path_homvec classes " << key << ": " << value.size() << "";
        }
        std::cout << "\n\n\n";
    }


    {
        auto classes = std::map<std::string, std::unordered_map<std::string, int>>{};
        // timer 
        auto start = std::chrono::high_resolution_clock::now();
        classes.clear();
        std::vector<std::vector<std::complex<long long>>> homvec_list{};

        {
            auto graph_list = get_graph_list();
            std::cout << "Read " << graph_list.size() << " graphs.\n";

            homvec_list.resize(graph_list.size());
            
            #pragma omp parallel for
            for (int i = 0; i < graph_list.size(); i++) {
                homvec_list[i] = wl::compute_complex_path_homvec(graph_list[i], graph_list[i].number_of_vertices()+plus, 1);
                if (i % percent == percent-1) {
                    std::cout << "|";
                    std::cout.flush();
                }
            }

            std::cout << "\nDone with computation";
        }

        for (int i = 0; i < homvec_list.size(); i++) {
            classes["   "][wl::stringyfy_vector(homvec_list[i])] += 1;
            uint64_t n = homvec_list[i].size() -plus-1;
            for (int j = 0; j <= plus; ++j) {
                classes["N+"+std::to_string(j)][std::to_string(homvec_list[i][n+j].real()) + "|" + std::to_string(homvec_list[i][n+j].imag())] += 1;
            }
        }
    

        auto end = std::chrono::high_resolution_clock::now();
        std::cout << "\n\t " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms";
        for (auto &&[key, value] : classes) {
            std::cout << "\n\tcompute_complex_path_homvec(+) classes " << key << ": " << value.size() << "";
        }
        std::cout << "\n\n\n";
    }
    

    {
        auto classes = std::map<std::string, std::unordered_map<std::string, int>>{};
        // timer 
        auto start = std::chrono::high_resolution_clock::now();
        classes.clear();

        std::vector<std::vector<std::complex<long long>>> homvec_list{};

        {
            auto graph_list = get_graph_list();
            std::cout << "Read " << graph_list.size() << " graphs.\n";

            homvec_list.resize(graph_list.size());
            
            #pragma omp parallel for
            for (int i = 0; i < graph_list.size(); i++) {
                homvec_list[i] = wl::compute_complex_path_homvec(graph_list[i], graph_list[i].number_of_vertices()+plus, -1);
                if (i % percent == percent-1) {
                    std::cout << "|";
                    std::cout.flush();
                }
            }

            std::cout << "\nDone with computation";
        }

        for (int i = 0; i < homvec_list.size(); i++) {
            classes["   "][wl::stringyfy_vector(homvec_list[i])] += 1;
            uint64_t n = homvec_list[i].size() -plus-1;
            for (int j = 0; j <= plus; ++j) {
                classes["N+"+std::to_string(j)][std::to_string(homvec_list[i][n+j].real()) + "|" + std::to_string(homvec_list[i][n+j].imag())] += 1;
            }
        }
    

        auto end = std::chrono::high_resolution_clock::now();
        std::cout << "\n\t " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms";
        for (auto &&[key, value] : classes) {
            std::cout << "\n\tcompute_complex_path_homvec(-) classes " << key << ": " << value.size() << "";
        }
        std::cout << "\n\n\n";
    }
    

    {
        auto classes = std::map<std::string, std::unordered_map<std::string, int>>{};
        // timer 
        auto start = std::chrono::high_resolution_clock::now();
        classes.clear();

        std::vector<std::vector<std::complex<long long>>> homvec_list{};

        {
            auto graph_list = get_graph_list();
            std::cout << "Read " << graph_list.size() << " graphs.\n";

            homvec_list.resize(graph_list.size());
            
            #pragma omp parallel for
            for (int i = 0; i < graph_list.size(); i++) {
                homvec_list[i] = wl::compute_complex_path_homvec_boost(graph_list[i], graph_list[i].number_of_vertices()+plus);
                if (i % percent == percent-1) {
                    std::cout << "|";
                    std::cout.flush();
                }
            }

            std::cout << "\nDone with computation";
        }

        for (int i = 0; i < homvec_list.size(); i++) {
            classes["   "][wl::stringyfy_vector(homvec_list[i])] += 1;
            uint64_t n = homvec_list[i].size() -plus-1;
            for (int j = 0; j <= plus; ++j) {
                classes["N+"+std::to_string(j)][std::to_string(homvec_list[i][n+j].real()) + "|" + std::to_string(homvec_list[i][n+j].imag())] += 1;
            }
        }
    

        auto end = std::chrono::high_resolution_clock::now();
        std::cout << "\n\t " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms";
        for (auto &&[key, value] : classes) {
            std::cout << "\n\tcompute_complex_path_homvec(B) classes " << key << ": " << value.size() << "";
        }
        std::cout << "\n\n\n";
    }
    




    {
        auto classes = std::map<std::string, std::unordered_map<std::string, int>>{};
        // timer 
        auto start = std::chrono::high_resolution_clock::now();
        classes.clear();

        std::vector<std::vector<long long>> homvec_list{};

        {
            auto graph_list = get_graph_list();
            std::cout << "Read " << graph_list.size() << " graphs.\n";

            homvec_list.resize(graph_list.size());
            
            #pragma omp parallel for
            for (int i = 0; i < graph_list.size(); i++) {
                homvec_list[i] = wl::compute_poly_path_homvec(graph_list[i], graph_list[i].number_of_vertices()+plus);
                if (i % percent == percent-1) {
                    std::cout << "|";
                    std::cout.flush();
                }
            }

            std::cout << "\nDone with computation";
        }

        for (int i = 0; i < homvec_list.size(); i++) {
            classes["   "][wl::stringyfy_vector(homvec_list[i])] += 1;
            uint64_t n = homvec_list[i].size() -plus-1;
            for (int j = 0; j <= plus; ++j) {
                classes["N+"+std::to_string(j)][std::to_string(homvec_list[i][n+j])] += 1;
            }
        }
    

        auto end = std::chrono::high_resolution_clock::now();
        std::cout << "\n\t " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms";
        for (auto &&[key, value] : classes) {
            std::cout << "\n\tcompute_poly_path_homvec classes " << key << ": " << value.size() << "";
        }
        std::cout << "\n\n\n";
    }
    


    // {
    //     auto rehash = std::unordered_map<std::string, unsigned long long>{};
    //     auto classes = std::unordered_map<std::string, std::vector<int>>{};
    //     auto start = std::chrono::high_resolution_clock::now();
    //     classes.clear();
    //     int i = 0;
    //     for (auto &&graph : graph_list) {
    //         auto vec = wl::colors_1(graph, rehash, graph.labels);
    //         classes[wl::stringyfy_vector(vec)].push_back(i);
    //         if (i % percent == percent-1) {
    //             std::cout << "|";
    //             std::cout.flush();
    //         }
    //         i += 1;
    //     }
    //     auto end = std::chrono::high_resolution_clock::now();
    //     std::cout << "\n\t " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms";
    //     std::cout << "\n\t #TW classes: " << classes.size() << "\n\n\n";
    // }

    // {
    //     auto rehash = wl::rehash_t{};
    //     auto start = std::chrono::high_resolution_clock::now();
    //     auto classes = std::unordered_map<std::string, std::vector<std::pair<int, wl::colvec_t>>>{};
    //     int i = 0;

    //     auto vecs = std::vector<wl::colvec_t>{};
    //     for (auto &&graph : graph_list) {
    //         auto vec = wl::colors_p1_05(graph, rehash, graph.labels);
    //         classes[wl::stringyfy_vector(vec)].emplace_back(i, vec);
    //         if (i % percent == percent-1) {
    //             std::cout << "|";
    //             std::cout.flush();
    //         }
    //         i += 1;
    //         vecs.push_back(vec);

    //         // std::cout << "\n\t";
    //         // std::cout << graph.labels.size() <<" -- "<< (int)graph.number_of_vertices() << "\t";
    //         // for (auto &&label : graph.labels) {
    //             // std::cout << (int)label << " ";
    //         // }

    //         // std::cout << "\n\t" << wl::stringyfy_vector(vec);
    //     }

    //     auto end = std::chrono::high_resolution_clock::now();
    //     std::cout << "\n\t " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms";
    //     std::cout << "\n\t #PW classes: " << classes.size() << "\n\n\n";

    //     // if (vecs.size() >= 2) {
    //     //     std::cout << "Explain:\n";
    //     //     std::cout << wl::explain(rehash, vecs[0], vecs[1]) << '\n';
    //     // }

    // }

    return 0;
}