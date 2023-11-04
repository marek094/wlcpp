
#include "read_g6.hpp"
#include "read_labeled_txt.hpp"
#include "read_tree_txt.hpp"
#include "homcounts.hpp"
#include "unlabelled_graph.hpp"
#include "logic.hpp"

#include <omp.h>

#include <iostream>
#include <vector>
#include <filesystem>


using eq_classes_t = std::vector<std::vector<wl::SmallGraph>>;

auto compute_eq_classes(std::vector<wl::SmallGraph> graph_list) -> eq_classes_t {
    auto results = std::vector<std::vector<unsigned long long>>{};
    results.resize(graph_list.size());
    
    #pragma omp parallel for
    for (size_t i = 0; i < graph_list.size(); i++) {
        results[i] = wl::compute_path_homvec(graph_list[i], graph_list[i].number_of_vertices()+1);
    }

    auto eq_classes_vec = eq_classes_t{};
    auto eq_classes = std::unordered_map<std::vector<unsigned long long>, 
                        std::vector<wl::SmallGraph>, wl::hash_vec<unsigned long long>>{};
    auto nontriv_count = 0ULL;
    auto i = 0ULL;
    for (auto&& hom_counts : results) {
        auto& eq_class = eq_classes[hom_counts];
        eq_class.emplace_back(std::move(graph_list[i]));
        nontriv_count += (eq_class.size() == 2);
        i += 1;
    }

    std::cout << "Found " << eq_classes.size() << " equivalence classes.\n";

    eq_classes_vec.reserve(nontriv_count);
    for (auto &&[homvec, graphs] : eq_classes) {
        if (graphs.size() > 1) {
            eq_classes_vec.emplace_back(std::move(graphs));
        }
    }

    return eq_classes_vec;
}



template<wl::LogicQuantifier T>
class BigEvaluator {
public:
    BigEvaluator(eq_classes_t const& eq_classes_vec, int rank) : rank{rank} {
        evaluators_vec.reserve(eq_classes_vec.size());
        for (auto&& graphs : eq_classes_vec) {
            evaluators_vec.emplace_back();
            evaluators_vec.back().reserve(graphs.size());
            for (auto&& graph : graphs) {
                evaluators_vec.back().emplace_back(graph);
            }
        }
    } 


    auto evaluate(wl::LogicFormula<T> const& formula) -> bool {
        // auto result = std::vector<int>{};
        // result.reserve(evaluators_vec.size());
        
        bool all_true = true;
        bool all_false = true;

        for (auto&& evaluators : evaluators_vec) {
            
            for (auto&& evaluator : evaluators) {
                bool is_graph_model = evaluator.evaluate({~0ULL, formula});
                if (is_graph_model == true) {
                    all_false = false;
                } else {
                    all_true = false;
                }
                if (!all_true && !all_false) {
                    return false;
                }
            }

            // if (all_true) {
            //     result.push_back(1);
            // } else if (all_false) {
            //     result.push_back(0);
            // } else {
            //     result.push_back(-1);
            // }
        }

        return true;
    }


    int rank;
    std::vector< std::vector< wl::EvaluatorLogicPaths< T>>> evaluators_vec;
};



int main_freetypes(std::vector<wl::SmallGraph> graph_list, int rank) {

    auto eq_classes_vec = compute_eq_classes(std::move(graph_list));
    size_t m = eq_classes_vec.size();

    using QuantifierT = wl::LogicQuantifierBound;
    auto big_evaluator = BigEvaluator<QuantifierT>{eq_classes_vec, rank};

    auto formulas = wl::gen_cold_closed_formulas<QuantifierT>(rank);
    for (auto&& formula : formulas) {
        bool r = big_evaluator.evaluate(formula);
        std::cout << r << "\t" << formula << "\n";
    }


    return 0;
}



int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <path_to_graph6_file>\n";
        return 1;
    }
    
    if (argc >= 3) {
        omp_set_num_threads(std::atoi(argv[2]));
    }

    int rank = 0;
    if (argc >= 4) {
        rank = std::atoi(argv[3]);
    }
    
    size_t limit = -1;
    if (argc >= 5) {
        limit = std::atoi(argv[3]);
    }

    size_t skip = 0;
    if (argc >= 6) {
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


    return main_freetypes(std::move(graph_list), rank);
}