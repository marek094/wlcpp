#include "read_g6.hpp"
#include "read_labeled_txt.hpp"
#include "read_tree_txt.hpp"
#include "logic.hpp"

#include <omp.h>

#include <iostream>
#include <filesystem>
#include <vector>
#include <iomanip>
#include <map>
#include <random>
#include <chrono>

int main(int argc, char* argv[]) {

    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <path_to_graph6_file>\n";
        return 1;
    }
    
    if (argc >= 3) {
        omp_set_num_threads(std::atoi(argv[2]));
    }

    int seed = 42;
    if (argc >= 4) {
        seed = std::atoi(argv[3]);
    }

    int rank = 0;
    // if (argc >= 4) {
    //     rank = std::atoi(argv[3]);
    // }
    
    size_t limit = -1;
    // if (argc >= 5) {
    //     limit = std::atoi(argv[3]);
    // }

    size_t skip = 0;
    // if (argc >= 6) {
    //     skip = std::atoi(argv[4]);
    // }


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


    {
        using namespace wl::logic;
        using namespace std;


        auto construct_logical_system = [seed] {
            auto f = std::map<std::string, shared_ptr<ElementBase>>{};

            f["true"] = make_shared<Fgt<1, 0>>(make_shared<True>());
            f["false"] = make_shared<Fgt<1, 0>>(make_shared<False>());
            auto adj_ = make_shared<Adj>();

            constexpr int max_c = 100;
            constexpr int max_l = 3;

            {
                auto childrenl1 = AtomicQuantifier::ChildrenArray{adj_, f["true"]};
                for (int c = 1; c <= max_c; ++c) {
                    f["phi_l1c" + to_string(c)] = make_shared<AtomicQuantifier>(childrenl1);
                    childrenl1.push_back(f["false"]);
                }
            }

            for (int l = 2; l <= max_l; ++l) {
                auto childrenl2 = AtomicQuantifier::ChildrenArray{adj_};
                for (int c = 1; c <= max_c; ++c) {
                    childrenl2.push_back(f["phi_l" + to_string(l-1) + "c" + to_string(c)]);
                    f["phi_l" + to_string(l) + "c" + to_string(c)] = make_shared<AtomicQuantifier>(childrenl2);
                }
            }

            
            // randomly select element from f
            auto rand_gen = mt19937(seed);
            for (int i = 0; i < 100; ++i) {
                auto children = AtomicQuantifier::ChildrenArray{adj_};
                for (int c = 0; c <= 200; ++c) {
                    auto it = f.begin();
                    std::advance(it, rand_gen() % f.size());
                    children.push_back(it->second);
                    // if (children.size() % 10 == 9) {
                        f["rand[" + to_string(i) + ", " + to_string(c) + "]"
                            ] = make_shared<AtomicQuantifier>(children);
                    // }
                }
            }

            return f;
        };


        // auto start_time = std::chrono::steady_clock::now();

        auto elements = construct_logical_system();
        
        std::cout << "Constructed logical system with " << elements.size() << " elements.\n";
        std::cout << std::string(30, '=') << "\n";

        auto results = std::vector< std::vector< std::vector<bool>>>(elements.size(), std::vector<std::vector<bool>>(graph_list.size()));

        #pragma omp parallel for
        for (int i = 0; i < graph_list.size(); ++i) {
            auto&& graph = graph_list[i]; 
            std::cout << "Graph " << i << std::endl;
            int j = 0;

            for (auto&& [name, phi] : construct_logical_system()) { // the cache is not used
                for (int node_i = 0; node_i < graph.number_of_vertices(); ++node_i) {
                    results[j][i].push_back(phi->evaluate(graph, node_i));
                }
                j += 1;
            }
            std::cout << "Graph " << i << " done\n";

        }



        // for (int j = 0; j < elements.size(); ++j) {
        int j = 0;
        for (auto &&[name, phi] : elements) {
            // std::ostringstream oss;
            auto counts = std::array<int, 2>{};
            for (int i = 0; i < graph_list.size(); ++i) {
                for (int k = 0; k < graph_list[i].number_of_vertices(); ++k) {
                    // if (graph) {
                        // oss << name << "\t";
                        // oss << "Graph " << std::setw(2) << i << "   Node " << std::setw(2) << k << "   ";
                        // oss << phi->to_string() << "\t" << results[j][i][k] << " ";
                        // oss << "\n";
                    // }
                    counts[i] += results[j][i][k];
                }
            }

            // std::cout << oss.str();
            if (counts[0] != counts[1]) {
                std::cout << "Counts: " << counts[0] << " " << counts[1] << "\n";
                std::cout << name << "\t";
                // phi->to_ostream(std::cout) << "\n";


            }

            j += 1;
        }

        std::cout << std::string(30, '=') << "\nBIG FAT DONE\n" << std::endl;



        // std::cout << std::endl;
        // for (int i = 0; i < graph_list.size(); ++i) {
        //     f["phi_l3c3"]->clear_cache();
        //     std::cout << "Graph" << i << ": " << f["phi_l3c3"]->evaluate(graph_list[i], 0) << std::endl;

        //     auto aq = dynamic_cast<AtomicQuantifier*>(f["phi_l3c3"].get());
        //     std::cout << aq->evaluation.size() << std::endl;
        //     // for (auto&& [va, result] : aq->evaluation) {
        //     //     std::cout << "EVAL " << (int)va[0] << " " << result << "\n";
        //     // }

        //     // for (auto&& child : aq->children) {
        //     //     // choose if dynamic_cast to Adj or AtomicQuantifier
        //     //     // in case of Adj, print the evaluation
        //     //     // in case of AtomicQuantifier, print the evaluation

        //     //     if (auto adj2 = dynamic_cast<Adj*>(child.get())) {
        //     //         for (auto&& [va, result] : adj2->evaluation) {
        //     //             std::cout << "->EVAL (A)" << (int)va[0] << " " << (int)va[1] << " " << result << "\n";
        //     //         }
        //     //         // std::cout << "EVAL " << (int)adj->evaluation[0] << " " << adj->evaluation[1] << "\n";
        //     //     } else if (auto aq = dynamic_cast<AtomicQuantifier*>(child.get())) {
        //     //         for (auto&& [va, result] : aq->evaluation) {
        //     //             std::cout << "->EVAL (P)" << (int)va[0] << " " << result << "\n";
        //     //         }
        //     //         // std::cout << "EVAL " << (int)aq->evaluation[0] << " " << aq->evaluation[1] << "\n";
        //     //     }
        //     // }

        //     std::cout << "\n\n";
        //     f["phi_l3c3"]->clear_cache();
        // }
        // return 0;

        // for (auto&& name : vector<string>{"phi_l2c1", "phi_l2c2", "phi_l2c3", "phi_l2c4"}) {
        //     for (int node_i = 0; node_i < graph_list[0].number_of_vertices(); ++node_i) {
        //         bool r1, r2;
        //         std::cout << name << " " << node_i << "\t" << (r1=f[name]->evaluate(graph_list[0], 0)) << "\n";
        //         f[name]->clear_cache();
        //         std::cout << name << " " << node_i << "\t" << (r2=f[name]->evaluate(graph_list[1], 0)) << "\n";
        //         f[name]->clear_cache();
        //         std::cout << "\n\n";
        //         if (r1 != r2) {
        //             std::cout << "ERROR\n";
        //             return 0;
        //         }
        //     }

            
        // }

    }






    return 0;
}