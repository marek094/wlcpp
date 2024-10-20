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
    for (auto &&[key, value] : classes) {
        std::cout << "\n\t" << name << " classes " << key << ": " << value.size() << "";
    }
    std::cout << "\n\n\n";

}







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



    int plus = 9;
    

    std::vector<std::pair<std::string, std::function<void(wl::SmallGraph&, map_hashmap_t&)>>> runs;


    
    {
        auto func = [=]<typename Int>(auto&& graph, auto& classes) {
            auto homvec_out = wl::compute_path_homvec<Int>(graph, graph.number_of_vertices()+plus);
                
            classes["   "][wl::stringyfy_vector(homvec_out)] += 1;
            uint64_t n = homvec_out.size() -plus-1;
            for (int j = 0; j <= plus; ++j) {
                classes["N+"+to_string(j)][to_string(homvec_out[n+j])] += 1;
            }

            // int sum = 0;
            // for (int j = 0; j < homvec_out.size()-plus; ++j) {
            //     sum += homvec_out[j];
            // }
            // for (int j = 0; j <= plus; ++j) {
            //     sum += homvec_out[n+j];
            //     classes["xN+"+std::to_string(j)][std::to_string(sum)] += 1;
            // }

            Int sum_odd = 0;
            Int sum_even = 0;
            for (int j = 0; j < homvec_out.size()-plus; ++j) {
                if (j % 2 == 0) {
                    sum_even += homvec_out[j];
                } else {
                    sum_odd += homvec_out[j];
                }
            }
            for (int j = 0; j <= plus; ++j) {
                if ((n+j) % 2 == 0) {
                    sum_even += homvec_out[n+j];
                } else {
                    sum_odd += homvec_out[n+j];
                }
                classes["yN+"+to_string(j)][to_string(sum_even) + "|" + to_string(sum_odd)] += 1;
            }
            
        };



        runs.emplace_back("compute_path_homvec<4>", [=](auto&& g, auto& c) {return func.operator()<wl::crt::crtu64t<4>>(g, c);});
        runs.emplace_back("compute_path_homvec<1>", [=](auto&& g, auto& c) {return func.operator()<uint64_t>(g, c);});
    }


    runs.emplace_back("compute_complex_path_homvec", [=](auto&& graph, auto& classes) {
        auto homvec_out = wl::compute_complex_path_homvec(graph, graph.number_of_vertices()+plus);
                
        classes["   "][wl::stringyfy_vector(homvec_out)] += 1;
        uint64_t n = homvec_out.size() -plus-1;
        for (int j = 0; j <= plus; ++j) {
            classes["N+"+std::to_string(j)][std::to_string(homvec_out[n+j].real()) + "|" + std::to_string(homvec_out[n+j].imag())] += 1;
        }

        
    });


    // runs.emplace_back("compute_complex_path_homvec(-)", [=](auto&& graph, auto& classes) {
    //     auto homvec_out =  wl::compute_complex_path_homvec(graph, graph.number_of_vertices()+plus, -1);
                
    //     classes["   "][wl::stringyfy_vector(homvec_out)] += 1;
    //     uint64_t n = homvec_out.size() -plus-1;
    //     for (int j = 0; j <= plus; ++j) {
    //         classes["N+"+std::to_string(j)][std::to_string(homvec_out[n+j].real()) + "|" + std::to_string(homvec_out[n+j].imag())] += 1;
    //     }
    // });


    runs.emplace_back("compute_complex_path_homvec_boost(B)", [=](auto&& graph, auto& classes) {
        auto homvec_out =  wl::compute_complex_path_homvec_boost(graph, graph.number_of_vertices()+plus);
                
        classes["   "][wl::stringyfy_vector(homvec_out)] += 1;
        uint64_t n = homvec_out.size() -plus-1;
        for (int j = 0; j <= plus; ++j) {
            classes["N+"+std::to_string(j)][std::to_string(homvec_out[n+j].real()) + "|" + std::to_string(homvec_out[n+j].imag())] += 1;
        }
    });
    runs.pop_back();


    runs.emplace_back("compute_path_homvec_labeled", [=](auto&& graph, auto& classes) {
        auto homvec_out = wl::compute_path_homvec_labeled(graph, graph.number_of_vertices()+plus);
        classes["   "][wl::stringyfy_vector(homvec_out)] += 1;
    });


    {
        auto func = [=]<typename Int>(auto&& graph, auto& classes) {
            auto homvec_out = wl::compute_complex_pathwidth_one_homvec<Int>(graph, graph.number_of_vertices()+plus);

            for (int i = 0; i < homvec_out.size(); i++) {
                classes["   "][wl::stringyfy_vector(homvec_out)] += 1;
                uint64_t n = homvec_out.size() -plus-1;
                for (int j = 0; j <= plus; ++j) {
                    classes["N+"+to_string(j)][to_string(homvec_out[n+j].real()) + "|" + to_string(homvec_out[n+j].imag())] += 1;
                }
            }
        };
        
        runs.emplace_back("compute_complex_pathwidth_one_homvec<16>", [=](auto&& g, auto& c) {return func.operator()<wl::crt::crtu64t<16>>(g, c);});
        runs.emplace_back("compute_complex_pathwidth_one_homvec<1>", [=](auto&& g, auto& c) {return func.operator()<int64_t>(g, c);});
    }
    runs.pop_back();
    runs.pop_back();


    runs.emplace_back("compute_complex_pathwidth_one_homvec(D+Ax)<16>", [=](auto&& graph, auto& classes) {
        using Int = wl::crt::crtu64t<16>;
        auto homvec_out = wl::compute_complex_pathwidth_one_homvec<Int>(graph, graph.number_of_vertices()+plus, /*switch*/ true);

        for (int i = 0; i < homvec_out.size(); i++) {
            classes["   "][wl::stringyfy_vector(homvec_out)] += 1;
            uint64_t n = homvec_out.size() -plus-1;
            for (int j = 0; j <= plus; ++j) {
                classes["N+"+to_string(j)][to_string(homvec_out[n+j].real()) + "|" + to_string(homvec_out[n+j].imag())] += 1;
            }
        }
        
    });
    runs.pop_back();


    runs.emplace_back("compute_complex_pathwidth_one_homvec(B)<16>", [=](auto&& graph, auto& classes) {
        using Int = wl::crt::crtu64t<16>;
        auto homvec_out = wl::compute_complex_pathwidth_one_homvec<Int>(graph, graph.number_of_vertices()+plus, /*switch*/ false, true);

        for (int i = 0; i < homvec_out.size(); i++) {
            classes["   "][wl::stringyfy_vector(homvec_out)] += 1;
            uint64_t n = homvec_out.size() -plus-1;
            for (int j = 0; j <= plus; ++j) {
                classes["N+"+to_string(j)][to_string(homvec_out[n+j].real()) + "|" + to_string(homvec_out[n+j].imag())] += 1;
            }
        }
        
    });
    runs.pop_back();


    runs.emplace_back("compute_quatern_pathwidth_one_homvec(B)<16>", [=](auto&& graph, auto& classes) {
        using Int = wl::crt::crtu64t<16>;
        // using Int = int64_t;
        auto homvec_out = wl::compute_quatern_pathwidth_one_homvec<Int>(graph, graph.number_of_vertices()+plus);

        for (int i = 0; i < homvec_out.size(); i++) {
            classes["   "][wl::stringyfy_vector(homvec_out)] += 1;
            uint64_t n = homvec_out.size() -plus-1;
            for (int j = 0; j <= plus; ++j) {
                classes["N+"+to_string(j)][
                    (std::stringstream{} << homvec_out[n+j]).str()
                ] += 1;
            }
        }
        
    });





    // return 0;


    // runs.emplace_back("compute_pathwidth_one_homvec_bases", [=](auto&& graph, auto& classes) {
    //     long long square = graph.number_of_vertices() * graph.number_of_vertices() / 2;
    //     auto homvec_out = wl::compute_pathwidth_one_homvec_bases(graph, square, square);
            
    //     classes["   "][wl::stringyfy_vector(homvec_out)] += 1;
    //     uint64_t n = homvec_out.size() -plus-1;
    // });



    runs.emplace_back("compute_pathwidth_one_homvec_v2", [=](auto&& graph, auto& classes) {
        auto [homvec_out, homexpr] = wl::compute_pathwidth_one_homvec_v2(graph, graph.number_of_vertices()+5);
            
        classes["   "][wl::stringyfy_vector(homvec_out)] += 1;


        // std::map<std::string, std::string> sorted_exprs;
        // for (int i = 0; i < homvec_out.size(); ++i) {
        //     auto const& vec = homvec_out[i];
        //     auto const& expr = homexpr[i];

        //     std::string sorted_expr = expr;
        //     std::sort(sorted_expr.begin(), sorted_expr.end());
        //     sorted_exprs[sorted_expr] += std::to_string(vec) + " ";
        // }

        // for (auto &&[key, value] : sorted_exprs) {
        //     classes[key][value] += 1;
        // }

        
    });
    
    


    // list of all runs
    for (auto&& [name, func] : runs) {
        std::cout << name << std::endl;
    }

    for (auto&& [name, func] : runs) {
        test_run(name, get_graph_list, func);
    }
    





    return 0;
}