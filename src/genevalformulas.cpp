#include <iostream>
#include <vector>
#include <algorithm>
#include <ranges>
#include <string>
#include <chrono>
#include <filesystem>
#include <tuple>
#include <set>
#include <map>
#include <coroutine>
// #include <generator>
#include "ref_generator.hpp"

#include <omp.h>

#include "read_g6.hpp"
#include "read_labeled_txt.hpp"
#include "read_tree_txt.hpp"
#include "homcounts.hpp"
#include "wl.hpp"
#include "explain.hpp"
#include "unlabelled_graph.hpp"



struct logic_formula_step_t {
    int count;
    bool is_E;

    // ==
    auto operator==(logic_formula_step_t const& other) const -> bool {
        return count == other.count && is_E == other.is_E;
    }

    auto static forall_count() -> generator<int> {
        constexpr auto kMaxCount = 10;
        for (int i = 1; i < kMaxCount; ++i) {
            co_yield i;
        }
    }

    auto static forall_is_E() -> generator<bool> {
        co_yield true;
        // co_yield false;
    }

    auto static forall() -> generator<logic_formula_step_t> {
        for (auto &&count : forall_count()) {
            for (auto &&is_E : forall_is_E()) {
                co_yield logic_formula_step_t{.count = count, .is_E = is_E};
            }
        }
    }

    auto to_string(std::string var1 = "x", std::string var2 = "y") const -> std::string {
        std::stringstream os;
        if (this->count == 0) {
            os << "(E" << var1 << ")";
        } else {
            os << "(E<" << this->count << var1 << ")";
        }
        
        if (this->is_E) {
            os << "A" << var2 << var1;
        } else {
            os << "~A" << var2 << var1;
        }
        return std::move(os.str());
    }

    // print me as << "E()"
    friend auto operator<<(std::ostream& os, logic_formula_step_t const& step) -> std::ostream& {
        os << step.to_string();
        return os;
    }

    // sorting operator <=>
    auto operator<(logic_formula_step_t const& other) const -> bool {
        if (this->count < other.count) {
            return true;
        }
        if (this->count > other.count) {
            return false;
        }
        return this->is_E < other.is_E;
    }
};




using logic_formula_t = std::vector<logic_formula_step_t>;
auto operator<< (std::ostream& os, logic_formula_t const& formula) -> std::ostream& {
    os << "{ ";
    int i = 0;
    auto opts = std::array<std::string, 2>{"x", "y"};
    for (auto it = formula.rbegin(); it != formula.rend(); ++it) {
        auto str = it->to_string(opts[1-i%2], opts[i%2]);
        os << str << " ";
        i += 1;
    }
    os << "}";
    return os;
}



template<class T>
struct component {
    T const& value;
    component(T const& value) : value(value) {}
    auto operator,(uintmax_t n) const -> uintmax_t {
        n ^= std::hash<T>()(value);
        n ^= n << (sizeof(uintmax_t) * 4 - 1);
        return n ^ std::hash<uintmax_t>()(n);
    }
};


// specialization for vector<logic_formula_step_t>
template<>
struct component<logic_formula_t> {
    logic_formula_t const& value;
    component(logic_formula_t const& value) : value(value) {}
    auto operator,(uintmax_t n) const -> uintmax_t {
        for (auto &&step : value) {
            n ^= std::hash<int>()(step.count);
            n ^= std::hash<bool>()(step.is_E);
            n ^= n << (sizeof(uintmax_t) * 4 - 1);
        }
        return n ^ std::hash<uintmax_t>()(n);
    }
};
    

class hash_tuple {
public:
    template<class Tuple>
    auto operator()(Tuple const& tuple) const -> size_t {
        return std::hash<uintmax_t>()(
            std::apply([](auto const& ... xs) { return (component(xs), ..., 0); }, tuple));
    }
};


class EvaluatorLogicPaths {
public:
    using args_t = std::tuple<unsigned long long, std::vector<logic_formula_step_t>>;
    using result_t = std::pair<args_t, bool>;
    using evaluated_t = std::unordered_map<args_t, bool, hash_tuple>;

    EvaluatorLogicPaths(wl::SmallGraph const& graph) : graph(graph) {}

    auto evaluate(args_t args) -> bool {
        if (auto it = results.find(args); it != results.end()) {
            // std::cout << "\t" << std::get<1>(args) << " HIT \n" << std::flush;
            return it->second;
        }

        auto [node_i, formula] = args;

        // std::cout << "\t" << formula << "Miss \t" << std::flush;
        if (formula.empty()) {
            results[args] = true;
            return true;
        }

        auto [given_count, given_is_E] = formula.back();

        formula.pop_back();
        auto actual_count = 0;

        auto atp_func = [&](auto&& x) { 
            return given_is_E ? graph.all_adj(x) : graph.all_nonadj(x); };

        for (auto &&adj : atp_func(node_i)) {
            if (this->evaluate({adj, formula})) {
                ++actual_count;
                if (given_count == 0) {
                    // \exists
                    results[args] = true;
                    return true;
                } else if (actual_count >= given_count) {
                    // \exists <count
                    results[args] = false;
                    return false;
                }
            }
        }

        if (given_count == 0) {
            assert (actual_count == 0);
            // \exists 
            results[args] = false;
            return false;
        } else {
            assert (actual_count < given_count);
            results[args] = true;
            return true;
        }
    }

    wl::SmallGraph const& graph;
    std::vector<args_t> call_stack;
    std::unordered_map<args_t, bool, hash_tuple> results;
};


auto gen_cold_formulas(int quantifier_depth) -> generator<logic_formula_t> {
    if (quantifier_depth == 0) {
        co_yield logic_formula_t{};
        co_return;
    }

    for (logic_formula_t &&formula : gen_cold_formulas(quantifier_depth - 1)) {
        for (auto &&step : logic_formula_step_t::forall()) {
            auto formula_copy = formula;
            formula_copy.push_back(step);
            co_yield formula_copy;
        }
    }
}


auto check_formulas_on_graph(wl::SmallGraph const& graph, int max_q=3) -> EvaluatorLogicPaths::evaluated_t {
    auto evaluator = EvaluatorLogicPaths(graph);

    // std::cout << (int)graph.number_of_vertices() << "-- " << "\n";
    int i = 0;
    for (auto &&formula : gen_cold_formulas(max_q)) {
        // make here a progress bar where i is the progress and it is updated every 2 seconds (so that the previous number is overwritten)
        if (i % 500 == 0) {
            std::cout << "\r" << i << " " << std::flush;
        }

        for (int node_i = 0; node_i < graph.number_of_vertices(); ++node_i) {
            evaluator.evaluate({node_i, formula});
        }
        i += 1;
    }
    std::cout << "\n" << std::flush;

    return std::move(evaluator.results);
}




int main(int argc, char* argv[]) {
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
    
    assert (graph_list.size() == 2);

    auto b1 = check_formulas_on_graph(graph_list[0], q_depth);
    std::cout << "Graph 1 done.\n";

    // for (auto &&[args, result] : b1) {
    //     std::cout << std::get<0>(args) << " " << std::get<1>(args) << " " << result << "\n";
    // }

    auto b2 = check_formulas_on_graph(graph_list[1], q_depth);
    std::cout << "Graph 2 done.\n";


    // types
    auto logic_type1 = std::vector<std::vector<logic_formula_t>>{}; // vector of formulas
    logic_type1.resize(graph_list[0].number_of_vertices(), std::vector<logic_formula_t>{});
    
    for (int i = 0; i < graph_list[0].number_of_vertices(); ++i) {
        for (auto &&[args, result] : b1) {
            if (! result) {
                continue;
            }
            auto [node_i, formula] = args;
            if (node_i == i) {
                logic_type1[i].push_back(formula);
            }
        }
    }

    auto logic_type2 = std::vector<std::vector<logic_formula_t>>{}; // vector of formulas
    logic_type2.resize(graph_list[1].number_of_vertices(), std::vector<logic_formula_t>{});

    for (int i = 0; i < graph_list[1].number_of_vertices(); ++i) {
        for (auto &&[args, result] : b2) {
            if (! result) {
                continue;
            }
            auto [node_i, formula] = args;
            if (node_i == i) {
                logic_type2[i].push_back(formula);
            }
        }
    }

    auto sort_function = [](auto&& a, auto&& b) {
        if (a.size() < b.size()) {
            return true;
        }
        if (a.size() > b.size()) {
            return false;
        }

        for (int i = 0; i < a.size(); ++i) {
            if (a[i] < b[i]) {
                return true;
            }
            if (b[i] < a[i]) {
                return false;
            }
        }

        return false;
    };

    int i = 0;
    for (auto&& formulas1 : logic_type1) {
        std::sort(formulas1.begin(), formulas1.end(), sort_function);
        int n_found = 0;
        for (auto&& formulas2 : logic_type2) {
            std::sort(formulas2.begin(), formulas2.end(), sort_function);
            if (formulas1 == formulas2) {
                n_found += 1;
            }
        }
        std::cout << "#" << i << " " << n_found << "\n";
        i += 1;
    }



    return 0;

    // model1 -> model2
    const auto model_map1 = std::unordered_map<int, int>{
        {0, 0}, {1, 8}, {2, 1}, {3, 3}, {4, 2}, {5, 5}, {6, 4}, {7, 7}, {8, 6}, {9, 9}, {10, 10}
    };

    // model2 -> model1
    const auto model_map2 = std::unordered_map<int, int>{
        {0, 0}, {1, 2}, {2, 4}, {3, 3}, {4, 6}, {5, 5}, {6, 8}, {7, 7}, {8, 1}, {9, 9}, {10, 10}
    };

    for (auto &&[args1, result1] : b1) {
        auto [model1, formula] = args1;
        
        auto it = model_map1.find(model1);
        assert(it != model_map1.end());
        auto model2 = it->second;
        auto args2 = std::make_tuple(model2, formula);
        assert (b2.find(args2) != b2.end());
        auto&& result2 = b2[args2];
        if (result1 != result2) {
            std::cout << "DIFF: " << model1 << " " << model2 << " [] " << formula << " " << result1 << " " << result2 << "\n";
        }
    }

    std::cout << "\ndebug eval\n";

    auto formulalist1 = std::vector<EvaluatorLogicPaths::args_t>{
        {1, {logic_formula_step_t{.count = 3, .is_E = true}, logic_formula_step_t{.count = 2, .is_E = true} }},
        {3, {logic_formula_step_t{.count = 3, .is_E = true}, logic_formula_step_t{.count = 2, .is_E = true} }},
        {4, {logic_formula_step_t{.count = 3, .is_E = true}, logic_formula_step_t{.count = 2, .is_E = true} }},
        
        {6, {logic_formula_step_t{.count = 3, .is_E = true},}},
        {10, {logic_formula_step_t{.count = 3, .is_E = true},}},
        {9, {logic_formula_step_t{.count = 3, .is_E = true},}},

        {7, {logic_formula_step_t{.count = 3, .is_E = true},}},

        {8, {logic_formula_step_t{.count = 3, .is_E = true},}},

        // === === ===

    };

    for (auto &&args1 : formulalist1) { 
        auto model1 = std::get<0>(args1);
        auto formula = std::get<1>(args1);
        auto it = model_map1.find(model1);
        assert(it != model_map1.end());
        auto model2 = it->second;
        auto args2 = std::make_tuple(model2, formula);
        
        std::cout << model1 << " " << model2 << " [] " << formula << " " << b1[args1] << " " << b2[args2] << "\n";
    }

    

    std::cout << "\n";

    auto formulalist2 = std::vector<EvaluatorLogicPaths::args_t>{
        {2, {logic_formula_step_t{.count = 3, .is_E = true}, logic_formula_step_t{.count = 2, .is_E = true} }},
        {3, {logic_formula_step_t{.count = 3, .is_E = true}, logic_formula_step_t{.count = 2, .is_E = true} }},
        {6, {logic_formula_step_t{.count = 3, .is_E = true}, logic_formula_step_t{.count = 2, .is_E = true} }},

        {1, {logic_formula_step_t{.count = 3, .is_E = true},}},
        {2, {logic_formula_step_t{.count = 3, .is_E = true},}},
        {9, {logic_formula_step_t{.count = 3, .is_E = true},}},

        {7, {logic_formula_step_t{.count = 3, .is_E = true},}},

        {6, {logic_formula_step_t{.count = 3, .is_E = true},}},
    };

    for (auto &&args2 : formulalist2) {
        auto model2 = std::get<0>(args2);
        auto formula = std::get<1>(args2);
        auto it = model_map2.find(model2);
        assert(it != model_map2.end());
        auto model1 = it->second;
        auto args1 = std::make_tuple(model1, formula);
        
        std::cout << model1 << " " << model2 << " [] " << formula << " " << b1[args1] << " " << b2[args2] << "\n";
    }




    auto end = std::chrono::steady_clock::now();
    std::cout << "Time elapsed: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms\n";

    return 0;
}