
#include "read_g6.hpp"
#include "read_labeled_txt.hpp"
#include "read_tree_txt.hpp"
#include "homcounts.hpp"
#include "unlabelled_graph.hpp"


#include <omp.h>

#include <Eigen/Dense>
#include <Eigen/Core>

#include <iostream>
#include <vector>
#include <filesystem>
#include <algorithm>
#include <unordered_map>
#include <map>
#include <compare>
#include <iomanip>


using EigenMatrixXull = Eigen::Matrix<uint64_t, Eigen::Dynamic, Eigen::Dynamic>;
using EigenVectorXull = Eigen::Matrix<uint64_t, Eigen::Dynamic, 1>;


// create cmp for vectors of EigenVectorXull
struct lex_cmp {
    auto operator()(EigenVectorXull const& a, EigenVectorXull const& b) const -> bool {
        return std::lexicographical_compare(a.data(), a.data()+a.size(), b.data(), b.data()+b.size());
    }
};


// returns integer 2D eigen matrix
auto compute_table(wl::SmallGraph const& G, int q_depth) -> std::pair<EigenMatrixXull, std::vector<std::vector<int>>> {
    // if depth is 0, return trivial
    if (q_depth == 0) {
        EigenMatrixXull table = EigenMatrixXull::Ones(1, 1);
        auto models = std::vector<std::vector<int>>{G.number_of_vertices(), std::vector<int>{0}};
        return {table, models};
    }

    auto [prev_table, prev_models] = compute_table(G, q_depth-1);
    assert (prev_table.rows() == q_depth);

    auto vectors_nodes = std::map<EigenVectorXull, std::vector<int>, lex_cmp>{};
    for (int node_i = 0; node_i < static_cast<int>(G.number_of_vertices()); ++node_i) {
        // maps counts indices to prev_table in neighbors of node_ik
        auto adj_models_counts = std::unordered_map<int, int>{};
        for (int node_j : G.all_adj(node_i)) {
            for (int prev_model : prev_models[node_j]) {
                adj_models_counts[prev_model] += 1;
            }
        }
        auto amc_ord = std::vector<std::pair<int, int>>{adj_models_counts.begin(), adj_models_counts.end()};
        std::sort(amc_ord.begin(), amc_ord.end(), [](auto&& a, auto&& b) { return a.second < b.second; });

        for (auto&& [prev_model, count] : amc_ord) {
            EigenVectorXull vector = EigenVectorXull::Zero(prev_table.rows()+1);
            // vector[0] = count;
            // vector.tail(prev_table.rows()) = prev_table.col(prev_model);
            // add the count to the end
            vector.head(prev_table.rows()) = prev_table.col(prev_model);
            vector[prev_table.rows()] = count;
            vectors_nodes[vector].push_back(node_i);
        }
    }

    EigenMatrixXull table = EigenMatrixXull::Zero(q_depth+1, vectors_nodes.size());
    auto models = std::vector<std::vector<int>>{G.number_of_vertices()};
    models.resize(G.number_of_vertices(), std::vector<int>{});

    int i = 0;
    for (auto&& [vector, nodes] : vectors_nodes) {
        table.col(i) = vector;
        for (int node : nodes) {
            models[node].push_back(i);
        }
        i += 1;
    }


    std::cout << "Table, v=0, q_depth = " << q_depth << "\n";
    for (auto m : models[0]) {
        std::cout << table.col(m).transpose().reverse() << "\n";
    }
    std::cout << "\n";

    return {table, models};
}
    

namespace wl {
void assert_(bool cond, std::string const& msg = "") {
    if (!cond) {
        std::cerr << msg << std::endl;
        throw std::runtime_error(msg);
    }
}
} // namespace wl

int main(int argc, char* argv[]) try {
    
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <path_to_graph6_file>\n";
        return 1;
    }
    
    if (argc >= 3) {
        omp_set_num_threads(std::atoi(argv[2]));
    }
    
    int q_depth = 1;
    if (argc >= 4) {
        q_depth = std::atoi(argv[3]);
    }


    std::filesystem::path file_path = argv[1];

    auto graph_list = [&]() {
        auto limit = ~0ULL;
        auto skip = 0;
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


    
    // for (int i = 0; auto&& G : graph_list) {
    //     std::cout << (i++) << "\n";
    //     auto homvec = wl::compute_path_homvec(G, q_depth);
    //     for (auto&& count : homvec) {
    //         std::cout << std::setw(5) << count << " ";
    //     }
    //     std::cout << "\n";
    // }
    // std::cout << "\n";



    int i = 0;
    EigenMatrixXull prev_table;

    // #pragma omp parallel for
    for (auto&& G : graph_list) {
        auto [table, models] = compute_table(G, q_depth);

        // std::cout << i << " ";
        // std::cout << "[" << table.cols() << ", " << table.rows() << "]\n";

        // for (int j = 0; auto&& model : models) {
        //     auto sum = 0ull;
        //     for (auto&& m : model) { sum += table.col(m).prod(); }
        //     std::cout << std::setw(5) << j << " " << std::setw(5) << sum << "\n"; 
        //     for (auto&& m : model) {
        //         std::cout << "\t" << table.col(m).reverse().transpose() << "\n";
        //     }
        //     j += 1;
        // }
        // std::cout << "\n";

        // std::cout << "Model Sums" << "\n";
        auto gsum = 0ULL;
        for (auto&& model : models) {
            auto model_sum = 0ULL;
            for (auto&& m : model) {
                model_sum += table.col(m).prod();
            }

            // std::cout << std::setw(5) << model_sum << " ";
            gsum += model_sum;
        }
        // std::cout << "\n";

        uint64_t table_max = table.maxCoeff();
        uint64_t max_degree = 0;
        for (int j = 0; j < G.number_of_vertices(); ++j) {
            max_degree = std::max(max_degree, G.adj_list[j].size());
        }


        wl::assert_(table_max <= max_degree, "Larger than max-degree");

        // unique values in the table
        auto vals = std::vector<bool>(max_degree, false);
        for (auto val : table.reshaped()) {
            vals[val-1] = true;
            // std::cout << val << std::endl;
        }

        std::cout << std::count(vals.begin(), vals.end(), false) << " ";
        

        // // wl::assert_(std::count(vals.begin(), vals.end(), false) == 0, "Not all values in table");
        // std::cout << "C" << std::count(vals.begin(), vals.end(), false) << "\n";
        
        // // // print table
        // std::cout << "Table\n" << std::setw(3) << table.transpose().rowwise().reverse() << "\n";
        // // if (i > 0) {
        // //     std::cout << "Table - prev_table\n" << (prev_table-table).cast<long long>().transpose().rowwise().reverse() << "\n";
        // // }


        // std::cout << "Number Of " << q_depth << "-Walks : " << gsum
        //     << "\tmax degree" << max_degree <<  "\ttable max"  << table_max << "\n";
        // std::cout << "\n";
        
        prev_table = table;
        i += 1;
    }



    auto end = std::chrono::steady_clock::now();
    std::cout << "\nElapsed time: " << std::chrono::duration<double>{end-start}.count() << "s\n";

    return 0;
} catch (std::exception& e) {
    std::cerr << "Error: " << e.what() << "\n";
    return 1;
} catch (...) {
    std::cerr << "Unknown error\n";
    return 1;
}