#include "read_g6.hpp"
#include "read_labeled_txt.hpp"
#include "read_tree_txt.hpp"
#include "homcounts.hpp"
#include "wl.hpp"
#include "explain.hpp"
#include "crtu64t.hpp"

#include <ska/unordered_map.hpp>
#include <omp.h>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/utility.hpp>


#include <iostream>
#include <vector>
#include <algorithm>
#include <ranges>
#include <chrono>
#include <filesystem>
#include <complex>
#include <map>
#include <functional>
#include <string>
#include <utility>
#include <sstream>
#include <ranges>
#include <unordered_set>
#include <regex>


using classes_update_t = std::vector<std::pair<std::string, std::string>>;
using hashmap_t = ska::unordered_map<std::string, std::vector<int>>;
using map_hashmap_t = std::map<std::string, hashmap_t>;
using vecmap_hashmap_t = std::vector<map_hashmap_t>;

using namespace std::string_literals;
using std::to_string;



template<typename InpFunc, typename ThreadFunc>
auto run_relation(std::string name, InpFunc&& get_graph_list, ThreadFunc&& thread_func, bool verbose = true, bool parallel = true) -> std::pair<std::string, hashmap_t> {
    std::ostringstream sout_quiet;
    auto& sout = (verbose ?  std::cout : sout_quiet);
    // timer 
    auto start = std::chrono::high_resolution_clock::now();

    auto par_classes = vecmap_hashmap_t{};
    unsigned num_pars = (parallel ? omp_get_max_threads() : 1u);
    par_classes.resize(num_pars);
    {
        auto graph_list = get_graph_list();
        sout << "Read " << graph_list.size() << " graphs.\n";
        int percent = std::max(graph_list.size() / 100, 1UL);

        #pragma omp parallel for if (parallel)
        for (int i = 0; i < graph_list.size(); i++) {
            auto classes_update = thread_func(graph_list[i]);

            for (auto&& [str1, str2] : classes_update) {
                unsigned par_idx = (parallel ? omp_get_thread_num() : 0);
                par_classes[par_idx][str1][str2].push_back(i);
            }

            if (i % percent == percent-1) {
                sout << "|" << std::flush;
            }
        }
        
        sout << "\nDone with computation" << std::endl;
    }
    
    // aggregate classes 
    auto& classes = par_classes.front();
    {
        for (int i = 1; i < par_classes.size(); ++i) {
            for (auto &&[key, value] : par_classes[i]) {
                for (auto &&[key2, value2] : value) {
                    auto &vec = classes[key][key2];
                    vec.insert(vec.end(), value2.begin(), value2.end());
                }
            }
            sout << "/" << std::flush;

            // make the space available in par_classes[i]
            par_classes[i] = map_hashmap_t{};
        }

        sout << "\nDone with aggregation" << std::endl;
    }

    auto end = std::chrono::high_resolution_clock::now();
    sout << "\n\t " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms";
    for (auto &&[key, value] : classes) {
        sout << "\n\t" << name << " classes " << key << ": " << value.size() << "";
    }
    sout << "\n\n\n";

    // return the main relation
    return {sout_quiet.str(), classes.begin()->second};
}




auto call_path_homvec(wl::SmallGraph& graph) -> classes_update_t {
    int plus = 9;
    using Int = int64_t;
    auto homvec_out = wl::compute_path_homvec(graph, graph.number_of_vertices()+ plus);

    auto classes_update = classes_update_t{};
    classes_update.emplace_back("   ", wl::stringyfy_vector(homvec_out));

    Int n = homvec_out.size() -plus-1;
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
        classes_update.emplace_back("yN+"+to_string(j), to_string(sum_even) + "|" + to_string(sum_odd));
    }

    return classes_update;
}

template<typename Int>
auto call_complex_path_homvec_log(wl::SmallGraph& graph) -> classes_update_t {
    int plus = 9;
    auto homvec_out = wl::compute_complex_path_homvec_log<Int>(graph, graph.number_of_vertices()+plus);
            
    auto classes_update = classes_update_t{};
    classes_update.emplace_back("   ", wl::stringyfy_vector(homvec_out));
    return classes_update;
}


template<typename Int>
auto call_complex_path_homvec(wl::SmallGraph& graph) -> classes_update_t {
    int plus = 9;
    auto homvec_out = wl::compute_complex_path_homvec<Int>(graph, graph.number_of_vertices()+14);

    auto classes_update = classes_update_t{};        
    // classes["   "][wl::stringyfy_vector(homvec_out)] += 1;
    classes_update.emplace_back("   ", wl::stringyfy_vector(homvec_out));
    uint64_t n = homvec_out.size() -plus-1;
    for (int j = 0; j <= plus; ++j) {
        // classes["N+"+std::to_string(j)][std::to_string(homvec_out[n+j].real()) + "|" + std::to_string(homvec_out[n+j].imag())] += 1;
        classes_update.emplace_back("N+"+to_string(j), to_string(homvec_out[n+j].real()) + "|" + to_string(homvec_out[n+j].imag()));
    }

    return classes_update;
}


template<typename InputFunc, typename SubFunc, typename SubSubFunc>
void check_relation(std::vector<std::vector<int>> const& relation_vec,
                    InputFunc&& get_graph_list,
                    std::string const& sub_name, SubFunc&& sub_func,
                    std::string const& subsub_name, SubSubFunc&& subsub_func,
                    bool do_outer) {
                        
    auto graph_list = get_graph_list();
    int percent = std::max(relation_vec.size() / 100, 1UL);
    #pragma omp parallel for if (do_outer)
    for (int i = 0; i < relation_vec.size(); ++i) {
        auto const& graph_idcs = relation_vec[i];

        if (graph_idcs.size() <= 1) {
            continue;
        }

        auto get_graph_sub_list = [&] {
            auto result = std::vector<wl::SmallGraph>{};
            result.reserve(graph_idcs.size());
            for (auto idx : graph_idcs) result.push_back(graph_list[idx]);
            return result;
        };

        auto [sub_output, sub_relation] = run_relation(sub_name, get_graph_sub_list, sub_func, false, !do_outer);

        if (sub_relation.size() != 1) {
            auto [subsub_output, subsub_relation] = run_relation(subsub_name, get_graph_sub_list, sub_func, false, false);
            if (subsub_relation.size() == 1) {
                std::cout << "\nWOW: a " << i << "-th non singleton class has " << sub_relation.size() << " subclasses:\n";
                std::cout << "BadGraphs\t"; for (auto idx : graph_idcs) std::cout << idx << "\t"; std::cout << "\n";
                std::cout << "rep: \n";
                std::cout << sub_output << "\n";
                std::cout << "subrep: \n";
                std::cout << subsub_output << "\n";
                std::cout << std::flush;
            } else {
                std::cout << "\nMOM: a right number of iteration fixed it!\n";
                // std::cout << "subrep: \n";
                // std::cout << subsub_output << "\n";
                std::cout << std::flush;
            }
        }

        if (i % percent == percent-1) {
            std::cout << "%" << std::flush;
        }
    }
    std::cout << std::endl;

}



int main(int argc, char* argv[]) {

    // using crt = wl::crt::crtu64t<2>;
    // crt crt1{5};
    // std::cout << crt1 << std::endl;

    // wl::crt::crtu64t<2> crt2{1ULL << 32};

    // std::cout << crt2 << std::endl;

    // std::cout << (crt2*crt2 - crt{2} * crt2 * crt1 + crt1*crt1)<< std::endl;
    // std::cout << (crt1 - crt2) * (- crt2 + crt1) << std::endl;
    

    // std::cout << ((crt1 - crt2) * (- crt2 + crt1) == (crt2*crt2 - crt{2} * crt2 * crt1 + crt1*crt1)) << std::endl;


    // Eigen::Matrix<crt, Eigen::Dynamic, Eigen::Dynamic> m(2,2);
    // std::cout << m << std::endl;

    // m(0,0) = crt{1};
    // m(0,1) = crt{2};
    // m(1,0) = crt{3};
    // m(1,1) = crt{4};

    // Eigen::Matrix<crt, Eigen::Dynamic, 1> ones = Eigen::Matrix<crt, Eigen::Dynamic, 1>::Ones(2);
    // std::cout << ones << std::endl;

    // std::cout << m * ones << std::endl;

    // for (int i = 0; i < 10; ++i) {
    //     m = m * m;
    // }
    // std::cout << m << std::endl;

    // std::complex<crt> c1{crt{1}, crt{2}};
    // std::complex<crt> unit{crt{0}, crt{1}};

    // std::cout << c1 << std::endl;
    // std::cout << unit << std::endl;
    // // 1 + 2i times i  == i + -2
    // auto r = c1 * unit;
    // std::cout << r << std::endl;

    // std::cout << r.real() + crt{2} << std::endl;



    // return 0;



    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <path_to_graph6_file>\n";
        return 1;
    }
    
    size_t threads = omp_get_max_threads();
    if (argc >= 3) {
        threads = std::atoi(argv[2]);
    }
    omp_set_num_threads(threads);
    
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
    std::string replaced_name = std::regex_replace(file_path_str, std::regex("/"), "__");

    auto get_graph_list = [=]() -> std::vector<wl::SmallGraph> {
        std::istringstream iss(file_path_str);
        auto graph_list = std::vector<wl::SmallGraph>{};
        for (std::string token; std::getline(iss, token, ' '); ) {
            std::filesystem::path file_path = token;
            
            if (file_path.extension() == ".txt" && file_path.stem().string().substr(0, 4) == "tree") {    
                auto r = wl::read_graph_from_tree_txt_file(file_path, limit, skip);
            } else
            if (file_path.extension() == ".txt") {
                auto r = wl::read_graph_from_labeled_txt_file(file_path, limit, skip);
            }
            assert(file_path.extension() == ".g6");
            auto r = wl::read_graph_from_graph6_file(file_path, limit, skip);
            if (graph_list.empty()) {
                graph_list = std::move(r);
            } else {
                graph_list.insert(graph_list.end(), r.begin(), r.end());
            }
        }
        return graph_list;
    };

    
    using lambda_func_t = std::function<classes_update_t(wl::SmallGraph&)>;

    auto name = "compute_complex_path_homvec_log<8>"s;
    lambda_func_t func = call_complex_path_homvec_log<wl::crt::crtu64t<8>>;

    auto sub_name = "compute_path_homvec"s;
    lambda_func_t sub_func = call_path_homvec;

    auto subsub_name = "compute_complex_path_homvec<8>"s;
    lambda_func_t subsub_func = call_complex_path_homvec<wl::crt::crtu64t<8>>;



    auto loaded_graph_list = [graphs = get_graph_list()]() mutable -> auto& { return graphs;};

    auto relation_vec = std::vector<std::pair<std::string, std::vector<int>>>{};
    {
        auto relation_cache_path = std::filesystem::path("cache") / std::format("cache___{}___{}.bin", replaced_name, name);
        if (! exists(relation_cache_path)) {
            auto [empty_output, relation] = run_relation(name, loaded_graph_list, func, true, true);
            auto relation_vec = std::vector(std::move_iterator(relation.begin()), std::move_iterator(relation.end()));
            std::ofstream ofs(relation_cache_path, std::ios::binary);
            std::cout << "\nFile " << relation_cache_path << " " << (ofs.is_open() ? "opened" : "not opened") << std::endl;
            boost::archive::binary_oarchive oarch(ofs);
            oarch << relation_vec;
            std::cout << "\nSaved " << relation_vec.size() << " classes." << std::endl;
        }
        
        assert(exists(relation_cache_path));

        {
            std::ifstream ifs(relation_cache_path, std::ios::binary);
            boost::archive::binary_iarchive iarch(ifs);
            iarch >> relation_vec;
            std::cout << "\nLoaded " << relation_vec.size() << " classes." << std::endl;
        }
    }

    auto big_small_th = omp_get_max_threads() * 2;
    auto relation_vec_bs = std::array<std::vector<std::vector<int>>, 2>{};
    auto relation_class_sizes = std::map<unsigned long long, unsigned long long>{};

    for (auto&& [class_string, graph_idcs] : relation_vec) {
        relation_class_sizes[graph_idcs.size()] += 1;
        if (graph_idcs.size() >= 2) {
            bool to_outer = (graph_idcs.size() < big_small_th);
            relation_vec_bs[(to_outer ? 0 : 1)].emplace_back(std::move(graph_idcs));
        }
    }

    relation_vec.clear();
    relation_vec.shrink_to_fit();



    std::cout << "These are the sizes of the classes of " << name << ":\n";
    for (auto&& [size, count] : relation_class_sizes) {
        std::cout << "\t" << size << " : " << count << std::endl;
    }

    check_relation(relation_vec_bs[0], loaded_graph_list, sub_name, sub_func, subsub_name, subsub_func, true);
    check_relation(relation_vec_bs[1], loaded_graph_list, sub_name, sub_func, subsub_name, subsub_func, false);


    std::cout << "Done." << std::endl;

    return 0;
}