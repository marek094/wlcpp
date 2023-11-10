#pragma once

#include "ref_generator.hpp"
#include "unlabelled_graph.hpp"


#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <utility>
#include <algorithm>
#include <ranges>
#include <set>
#include <memory>




// // enumerate relations
// enum struct Relation {
//     Equality = 0ULL,
//     Adjacency = 1ULL,
//     // Label etc
// };

// // template<int k>
// // struct Atp {
    
// // };

// // template<>
// // struct Atp<1> {
// //     auto static genall() -> generator<Atp<1>> {
// //         co_yield Atp<1>{};
// //     };
// // };



// // struct Atp1 {
// //     auto static genall() -> generator<Atp1> {
// //         co_yield Atp1{.values = {}};
// //     };

// //     std::array<bool, 0> values;
// // };

// // struct Atp2 {

// //     constexpr auto static listall() -> std::array<Atp2, 4> {
// //         return {
// //             Atp2{.values = {false, false}},
// //             Atp2{.values = {false, true}},
// //             Atp2{.values = {true, false}},
// //             Atp2{.values = {true, true}},
// //         };
// //     };

// //     std::array<bool, 2> values;
// // };

namespace wl {


struct hash_array {
    template<typename T, size_t N>
    auto operator()(std::array<T, N> const& arr) const -> size_t {
        size_t seed = arr.size();
        for (auto &&i : arr) {
            seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};



namespace logic {



struct ElementBase {
    virtual auto to_string() const -> std::string = 0;

    using sigma_struct_t = wl::SmallGraph;

    template<uint64_t NumFreeVariables>
    using variable_assignment_t = std::array<wl::SmallGraph::vertex_type, NumFreeVariables>;


    ElementBase() = default;

    explicit ElementBase(std::vector<std::shared_ptr<ElementBase>> children) : children(std::move(children)) {}


    auto evaluate(sigma_struct_t const& s) -> bool {
        return this->evaluate(s, ElementBase::variable_assignment_t<0>{});
    }

    virtual auto evaluate(sigma_struct_t const& s, variable_assignment_t<0> const& variable_assignment) -> bool {
        throw std::runtime_error("Not implemented (evaluate<0>) for `" + this->to_string() + "`");
        assert(false);
        return false;
    }


    virtual auto evaluate(sigma_struct_t const& s, variable_assignment_t<1> const& variable_assignment) -> bool {
        throw std::runtime_error("Not implemented (evaluate<1>) for `" + this->to_string() + "`");
        assert(false);
        return false;
    }


    virtual auto evaluate(sigma_struct_t const& s, variable_assignment_t<2> const& variable_assignment) -> bool {
        throw std::runtime_error("Not implemented (evaluate<2>) for `" + this->to_string() + "`");
        assert(false);
        return false;
    }

    std::vector<std::shared_ptr<ElementBase>> children;
};

template<uint64_t NumFreeVariables>
struct FormulaElementBase : public ElementBase {
    using sigma_struct_t = wl::SmallGraph;
    using variable_assignment_t = std::array<wl::SmallGraph::vertex_type, NumFreeVariables>;
    
    using ElementBase::ElementBase;
    
    static constexpr uint64_t kNumFreeVariables = NumFreeVariables;

    FormulaElementBase() = default;

    virtual auto evaluate_impl(sigma_struct_t const&, variable_assignment_t const&) -> bool {
        throw std::runtime_error("Not implemented (evaluate_impl)");
        assert(false);
        return false;
    }

    auto evaluate(sigma_struct_t const& s, variable_assignment_t const& variable_assignment) -> bool override {
        // std::cout << "Evaluate <" << NumFreeVariables << "> " << this->to_string() << std::endl;
        if (auto it = evaluation.find(variable_assignment); it != evaluation.end()) {
            return it->second;
        }

        auto result = this->evaluate_impl(s, variable_assignment);
        evaluation[variable_assignment] = result;
        return result;
    }

    auto num_free_variables() const -> uint64_t {
        return NumFreeVariables;
    }

    std::unordered_map<variable_assignment_t, bool, hash_array> evaluation;
};


template<uint64_t NumFreeVariables, typename Derived>
struct FormulaElement : public FormulaElementBase<NumFreeVariables> {
    using parent_t = FormulaElementBase<NumFreeVariables>;
    using variable_assignment_t = typename parent_t::variable_assignment_t;
    using sigma_struct_t = typename parent_t::sigma_struct_t;
    // using FormulaElementBase<NumFreeVariables>;
    static constexpr auto kNumFreeVariables = NumFreeVariables;

    using parent_t::parent_t;

    virtual auto evaluate_impl(sigma_struct_t const& s, variable_assignment_t const& variable_assignment) -> bool override {
        // std::cout << "Evaluate Impl <" << NumFreeVariables << "> " << this->to_string() << std::endl;
        return static_cast<Derived*>(this)->evaluate_(s, variable_assignment);
    }
};


struct True : public FormulaElement<0, True> {
    auto to_string() const -> std::string override {
        return "T";
    }

    auto evaluate_(sigma_struct_t const&, variable_assignment_t const&) -> bool {
        return true;
    }
};

struct False : public FormulaElement<0, False> {
    auto to_string() const -> std::string override {
        return "F";
    }

    auto evaluate_(sigma_struct_t const&, variable_assignment_t const&) -> bool {
        return false;
    }
};


struct Adj : public FormulaElement<2, Adj> {
    auto to_string() const -> std::string override {
        return "A";
    }

    auto evaluate_(sigma_struct_t const& s, variable_assignment_t const& variable_assignment) -> bool {
        auto [v1, v2] = variable_assignment;
        auto const& adjs = s.adj_list[v1];
        return std::find(adjs.begin(), adjs.end(), v2) != adjs.end(); 
    }
};

struct NonAdj : public FormulaElement<2, NonAdj> {

    auto to_string() const -> std::string override {
        return "~A";
    }

    auto evaluate_(sigma_struct_t const& s, variable_assignment_t const& variable_assignment) -> bool {
        auto [v1, v2] = variable_assignment;
        auto const& adjs = s.adj_list[v1];
        return std::find(adjs.begin(), adjs.end(), v2) == adjs.end(); 
    }

};

struct Eq : public FormulaElement<2, Eq> {

    auto to_string() const -> std::string override {
        return "=";
    }

    auto evaluate_(sigma_struct_t const&, variable_assignment_t const& variable_assignment) -> bool {
        auto [v1, v2] = variable_assignment;
        return v1 == v2;
    }

};


struct NonEq : public FormulaElement<2, NonEq> {
    
    auto to_string() const -> std::string override {
        return "~=";
    }

    auto evaluate_(sigma_struct_t const&, variable_assignment_t const& variable_assignment) -> bool {
        auto [v1, v2] = variable_assignment;
        return v1 != v2;
    }
    
};

struct And : public FormulaElement<2, And> {

    using FormulaElement<2, And>::FormulaElement;

    auto to_string() const -> std::string override {
        return this->children[0]->to_string() + "&" + this->children[1]->to_string();
    }

    auto evaluate_(sigma_struct_t const& s, variable_assignment_t const& variable_assignment) -> bool {
        return this->children[0]->evaluate(s, variable_assignment) && this->children[1]->evaluate(s, variable_assignment);
    }
    
};

struct Or : public FormulaElement<2, Or> {

    using FormulaElement<2, Or>::FormulaElement;

    auto to_string() const -> std::string override {
        return this->children[0]->to_string() + "|" + this->children[1]->to_string();
    }

    auto evaluate_(sigma_struct_t const& s, variable_assignment_t const& variable_assignment) -> bool {
        return this->children[0]->evaluate(s, variable_assignment) || this->children[1]->evaluate(s, variable_assignment);
    }
    
};

template<uint64_t NumFreeVariables, uint64_t... Variables>
struct Fgt : public FormulaElement<NumFreeVariables, Fgt<NumFreeVariables, Variables...>> {
    using parent_t = FormulaElement<NumFreeVariables, Fgt<NumFreeVariables, Variables...>>;
    using variable_assignment_t = typename parent_t::variable_assignment_t;
    using sigma_struct_t = typename parent_t::sigma_struct_t;

    using fgt_variable_assignment_t = FormulaElement<NumFreeVariables - sizeof...(Variables), void>::variable_assignment_t;
    
    static_assert(((Variables < NumFreeVariables) && ...) , "Variables must be smaller than NumFreeVariables");
    static_assert(sizeof...(Variables) <= NumFreeVariables, "Too many variables");


    using parent_t::parent_t;
    
    auto to_string() const -> std::string override {
        return this->children[0]->to_string();
    }

    auto evaluate_(sigma_struct_t const& s, variable_assignment_t const& variable_assignment) -> bool {
        return this->children[0]->evaluate(s, forget_variables(variable_assignment));
    }

private:
    constexpr static auto forget_variables(variable_assignment_t const& variable_assignment) -> fgt_variable_assignment_t {
        auto new_assignment = fgt_variable_assignment_t{};
        size_t newIndex = 0;
        for (size_t i = 0; i < NumFreeVariables; ++i) {
            if (((i != Variables) && ...)) {
                new_assignment[newIndex++] = variable_assignment[i];
            }
        }
        return new_assignment;
    }
};



struct Exists : public FormulaElement<1, Exists> {
    using parent_t = FormulaElement<1, Exists>;
    using variable_assignment_t = typename parent_t::variable_assignment_t;
    using sigma_struct_t = typename parent_t::sigma_struct_t;

    using parent_t::parent_t;

    auto to_string() const -> std::string override {
        return "Ex " + this->children[0]->to_string();
    }

    auto evaluate_(sigma_struct_t const& s, variable_assignment_t const& variable_assignment) -> bool {
        auto [v1] = variable_assignment;
        for (auto &&vtx : s.all({})) {
            if (this->children[0]->evaluate(s, {vtx, v1})) {
                return true;
            }
        }
        return false;
    }
};



struct ExistCount : public FormulaElement<1, ExistCount> {
    using parent_t = FormulaElement<1, ExistCount>;
    using variable_assignment_t = typename parent_t::variable_assignment_t;
    using sigma_struct_t = typename parent_t::sigma_struct_t;

    int count = 0;

    ExistCount(int count, std::vector<std::shared_ptr<ElementBase>> children) : parent_t(std::move(children)), count(count) {}

    auto to_string() const -> std::string override {
        return "E" + std::to_string(count) + "x " + this->children[0]->to_string();
    }

    auto evaluate_(sigma_struct_t const& s, variable_assignment_t const& variable_assignment) -> bool {
        auto [v1] = variable_assignment;
        
        int actual_count = 0;
        for (auto &&vtx : s.all({})) {
            actual_count += this->children[0]->evaluate(s, {vtx, v1});
            if (actual_count > count) {
                return false;
            }
        }
        return (actual_count == count);
    }
};



// template<uint64_t NumFreeVariables>
// struct Atp : public FormulaElement<NumFreeVariables, Atp<NumFreeVariables>> {
//     using variable_assignment_t = typename FormulaElementBase<NumFreeVariables>::variable_assignment_t;
//     using sigma_struct_t = typename FormulaElementBase<NumFreeVariables>::sigma_struct_t;

//     auto to_string() const -> std::string override {
//         return "atp_"s + NumFreeVariables;
//     }

//     auto evaluate_(sigma_struct_t const& s, variable_assignment_t const& variable_assignment_t) -> bool {
        
//     }
// };



} // namespace logic






template<typename Derived>
struct LogicQuantifierBase {
    // trivial version
    bool is_E = false;
    int k = 0;
    
    using self_type = LogicQuantifierBase<Derived>;

    LogicQuantifierBase() = default;
    LogicQuantifierBase(self_type const&) = default;
    LogicQuantifierBase(self_type&&) = default;
    LogicQuantifierBase& operator=(self_type const&) = default;
    LogicQuantifierBase& operator=(self_type&&) = default;


    LogicQuantifierBase(bool is_E, int k) : is_E{is_E}, k{k} {}

    auto static forall_is_E() -> generator<bool> {
        co_yield true;
        // co_yield false;
    
    }

    auto static forall() -> generator<Derived> {
        for (auto &&val : Derived::forall_impl()) {
            for (auto &&is_E : forall_is_E()) {
                auto cores = Derived{};
                cores.k = 2;
                cores.is_E = is_E;
                cores.val() = val;
                co_yield cores;
            }
        }
    }

    auto static forall_closing() -> generator<Derived> {
        for (auto &&val : Derived::forall_impl()) {
            auto cores = Derived{};
            cores.k = 1;
            cores.val() = val;
            co_yield cores;
        }
    }

    auto to_string(std::string var1 = "x", std::string var2 = "y") const -> std::string {
        std::stringstream os;
        os << "(E";
        static_cast<Derived const*>(this)->to_string_impl(os);
        os << var1 << ")";
        if (k == 2) {
            if (this->is_E) {
                os << "A" << var2 << var1;
            } else {
                os << "~A" << var2 << var1;
            }
        }
        return std::move(os.str());
    }

    friend auto operator<<(std::ostream& os, Derived const& step) -> std::ostream& {
        os << step.to_string();
        return os;
    }

    auto hash_impl() const -> size_t {
        auto val = static_cast<Derived const*>(this)->val();
        using hash_type = std::remove_cvref_t<decltype(val)>;
        return std::hash<hash_type>()(val);
    }

    friend auto hash_struct(self_type const& tuple) -> size_t {
        auto res = static_cast<Derived const&>(tuple).hash_impl();
        if (tuple.k == 1) {
            return res;
        }

        // assert( tuple.k == 2);
        res ^= res << (sizeof(uintmax_t) * 4 - 1);
        res ^= std::hash<bool>()(tuple.is_E);
        return res;
    }

    

    friend auto operator<=>(Derived const& lhs, Derived const& rhs) -> std::strong_ordering {
        if (auto cmp = lhs.val() <=> rhs.val(); cmp != std::strong_ordering::equal) return cmp;
        if (auto cmp = lhs.k <=> rhs.k; cmp != std::strong_ordering::equal) return cmp;
        return lhs.is_E <=> rhs.is_E;
    }

    friend auto operator==(Derived const& lhs, Derived const& rhs) -> bool {
        return (lhs <=> rhs) == std::strong_ordering::equal;
    }
};


template<typename T>
concept LogicQuantifier = requires(T a, const T b, std::ostream& os) {
    { a.val() } -> std::same_as<auto&>;
    { b.val() } -> std::same_as<const auto&>;
    { T::forall_impl() } -> std::same_as<generator<std::remove_reference_t<decltype(a.val())>>>;
    { a.to_string_impl(os) } -> std::same_as<void>;
};


template<int MinCount, int MaxCount>
struct LogicQuantifierBound : public LogicQuantifierBase<LogicQuantifierBound<MinCount, MaxCount>> {

    // expose constants
    static constexpr auto kMinCount = MinCount;
    static constexpr auto kMaxCount = MaxCount;

    int count = 0;

    LogicQuantifierBound() = default;

    auto val() const -> int const& {
        return this->count;
    }

    auto val() -> int& {
        return this->count;
    }

    auto static forall_impl() -> generator<int> {
        for (int i = MinCount; i <= MaxCount; ++i) {
            co_yield i;
        }
    }

    auto to_string_impl(std::ostream& os) const -> void {
        if (this->count == 0) {
            os << "";
        } else {
            os << "<" << this->count;
        }
    }

};



struct LogicQuantifierModulo : public LogicQuantifierBase<LogicQuantifierModulo> {
    int mod;
    
    LogicQuantifierModulo() = default;

    auto val() const -> int const& {
        return this->mod;
    }

    auto val() -> int& {
        return this->mod;
    }

    auto static forall_impl() -> generator<int> {
        constexpr auto kMaxMod = 4ULL;
        for (auto i = 1ULL; i <= kMaxMod; ++i) {
            co_yield static_cast<int>(1ULL << i);
        }
        // co_yield 3;
    }

    auto to_string_impl(std::ostream& os) const -> void {
        os << "%" << this->mod;
    }

};


struct LogicQuantifierBinMod : public LogicQuantifierBase<LogicQuantifierBinMod> {
    bool modulo_is_one;
    
    LogicQuantifierBinMod() = default;

    auto val() const -> bool const& {
        return this->modulo_is_one;
    }

    auto val() -> bool& {
        return this->modulo_is_one;
    }

    auto static forall_impl() -> generator<bool> {
        co_yield true;
        co_yield false;
    }

    auto to_string_impl(std::ostream& os) const -> void {
        os << (this->modulo_is_one ? "&" : "^");
    }

};



struct LogicQuantifierBinCount : public LogicQuantifierBase<LogicQuantifierBinCount> {
    bool is_digit_not_carry;
    
    LogicQuantifierBinCount() = default;
    
    auto val() const -> bool const& {
        return this->is_digit_not_carry;
    }

    auto val() -> bool& {
        return this->is_digit_not_carry;
    }

    auto static forall_impl() -> generator<bool> {
        co_yield true;
        co_yield false;
    }

    auto to_string_impl(std::ostream& os) const -> void {
        os << (this->is_digit_not_carry ? "&" : "^");
    }

};


template<int Length>
struct LogicQuantifierPartition : public LogicQuantifierBase<LogicQuantifierPartition<Length>> {
    int length;

    LogicQuantifierPartition() = default;

    auto val() const -> int const& {
        return this->length;
    }

    auto val() -> int& {
        return this->length;
    }

    auto static forall_impl() -> generator<int> {
        for (int i = 1; i <= Length; ++i) {
            co_yield i;
        }
    }

    auto to_string_impl(std::ostream& os) const -> void {
        os << "P" << this->length;
    }

};








template<LogicQuantifier T>
struct LogicFormula {

    using LogicQuantifierT = T;


    friend auto operator<< (std::ostream& os, LogicFormula const& formula) -> std::ostream& {
        os << "{ ";
        int i = 0;
        auto opts = std::array<std::string, 2>{"x", "y"};
        for (auto it = formula.quantifiers.rbegin(); it != formula.quantifiers.rend(); ++it) {
            auto str = it->to_string(opts[1-i%2], opts[i%2]);
            os << str << " ";
            i += 1;
        }
        os << "}";
        return os;
    }


    auto operator==(LogicFormula const& other) const -> bool {
        return this->quantifiers == other.quantifiers;
    }

    auto operator<=>(LogicFormula const& other) const -> std::strong_ordering {
        if (auto cmp = this->quantifiers.size() <=> other.quantifiers.size(); 
            cmp != std::strong_ordering::equal) return cmp;

        for (int i = 0; i < this->quantifiers.size(); ++i) {
            if (auto cmp = this->quantifiers[i] <=> other.quantifiers[i]; 
                cmp != std::strong_ordering::equal) return cmp;
        }

        return std::strong_ordering::equal;
    }


    std::vector<T> quantifiers;
};


template<class T>
struct hashed_component {
    T const& value;
    hashed_component(T const& value) : value(value) {}
    auto operator,(uintmax_t n) const -> uintmax_t {
        n ^= std::hash<T>()(value);
        n ^= n << (sizeof(uintmax_t) * 4 - 1);
        return n ^ std::hash<uintmax_t>()(n);
    }
};


// specialization for vector<logic_quantifier_bound_t>
template<LogicQuantifier T>
struct hashed_component<LogicFormula<T>> {
    LogicFormula<T> const& value;
    hashed_component(LogicFormula<T> const& value) : value(value) {}
    auto operator,(uintmax_t n) const -> uintmax_t {
        for (auto &&step : value.quantifiers) {
            n ^= hash_struct(step);
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
            std::apply([](auto const& ... xs) { return (hashed_component(xs), ..., 0); }, tuple));
    }
};

template<std::integral T>
struct hash_vec {
    size_t operator()(std::vector<T> const& vec) const {
        size_t seed = vec.size();
        for (auto &&i : vec) {
            seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};



} // namespace wl



namespace std {
template <wl::LogicQuantifier T>
struct hash<wl::LogicFormula<T>> {
    size_t operator()(wl::LogicFormula<T> const& f) const {
        return wl::hashed_component<wl::LogicFormula<T>>{f}.operator,(0);
    }
};
} // namespace std

namespace wl {


template<LogicQuantifier T>
class EvaluatorLogicPaths {
public:
    using QuantifierT = T;

    using args_t = std::tuple<unsigned long long, LogicFormula<T>>;
    using result_t = std::pair<args_t, bool>;
    using evaluated_t = std::unordered_map<args_t, bool, hash_tuple>;

    EvaluatorLogicPaths(wl::SmallGraph const& graph) : graph(graph) {}



    static auto atp_scope(LogicQuantifier auto const& quantifier) -> std::function<generator<wl::SmallGraph::vertex_type>(wl::SmallGraph const&, unsigned long long)> {
        assert(quantifier.k == 2 or quantifier.k == 1);

        if (quantifier.k == 2) {
            if (quantifier.is_E) {
                return &wl::SmallGraph::all_adj;
            } else {
                return &wl::SmallGraph::all_nonadj;
            }
        }

        assert(quantifier.k == 1);
        return &wl::SmallGraph::all;
    }
    
    template<int MinCount, int MaxCount>
    auto quantifier_semantics(
        LogicQuantifierBound<MinCount, MaxCount> const& quantifier, 
        unsigned long long node_i, LogicFormula<LogicQuantifierBound<MinCount, MaxCount>> const& formula) -> bool {
        auto actual_count = 0;

        for (auto &&adj : atp_scope(quantifier)(graph, node_i)) {
            if (this->evaluate({adj, formula})) {
                ++actual_count;
                if (quantifier.count == 0) {
                    // \exists
                    return true;
                } else if (actual_count >= quantifier.count) {
                    // \exists <count
                    return false;
                }
            }
        }

        if (quantifier.count == 0) {
            assert (actual_count == 0);
            // \exists 
            return false;
        } else {
            assert (actual_count < quantifier.count);
            return true;
        }
    }

    auto quantifier_semantics(LogicQuantifierModulo const& quantifier, unsigned long long node_i, LogicFormula<LogicQuantifierModulo> const& formula) -> bool {
        auto actual_count = 0;
        for (auto &&adj : atp_scope(quantifier)(graph, node_i)) {
            if (this->evaluate({adj, formula})) {
                ++actual_count;
            }
        }

        if (actual_count % quantifier.mod == 0) {
            return true;
        } else {
            return false;
        }
    }

    auto quantifier_semantics(LogicQuantifierBinMod const& quantifier, unsigned long long node_i, LogicFormula<LogicQuantifierBinMod> const& formula) -> bool {
        auto actual_count = 0;
        for (auto &&adj : atp_scope(quantifier)(graph, node_i)) {
            if (this->evaluate({adj, formula})) {
                ++actual_count;
            }
        }

        if (quantifier.modulo_is_one) {
            if (actual_count % 2 == 1) {
                return true;
            } else {
                return false;
            }
        } else {
            if (actual_count % 2 == 0) {
                return true;
            } else {
                return false;
            }
        }
    }

    auto quantifier_semantics(LogicQuantifierBinCount const& quantifier, unsigned long long node_i, LogicFormula<LogicQuantifierBinCount> const& formula) -> bool {
        auto actual_count = 0;
        for (auto &&adj : atp_scope(quantifier)(graph, node_i)) {
            if (this->evaluate({adj, formula})) {
                ++actual_count;
            }
        }

        if (quantifier.is_digit_not_carry) {
            if (actual_count % 2 > 0) {
                return true;
            } else {
                return false;
            }
        } else {
            if (actual_count / 2 > 0) {
                return true;
            } else {
                return false;
            }
        }
    }

    



    auto evaluate(args_t args) -> bool {
        if (auto it = results.find(args); it != results.end()) {
            // std::cout << "\t" << std::get<1>(args) << " HIT \n" << std::flush;
            return it->second;
        }

        auto [node_i, formula] = args;

        // std::cout << "\t" << formula << "Miss \t" << std::flush;
        if (formula.quantifiers.empty()) {
            results[args] = true;
            return true;
        }

        auto quantifier = formula.quantifiers.back();
        formula.quantifiers.pop_back();

        auto result = this->quantifier_semantics(quantifier, node_i, formula);
        results[args] = result;
        return result;
    }

    wl::SmallGraph const& graph;
    // std::vector<args_t> call_stack;
    std::unordered_map<args_t, bool, hash_tuple> results;
};


template<LogicQuantifier T>
auto gen_cold_formulas(int quantifier_depth) -> generator<LogicFormula<T>> {
    if (quantifier_depth == 0) {
        co_yield LogicFormula<T>{};
        co_return;
    }

    for (LogicFormula<T> &&formula : gen_cold_formulas<T>(quantifier_depth - 1)) {
        for (auto &&step : T::forall()) {
            auto formula_copy = formula;
            formula_copy.quantifiers.push_back(step);
            co_yield formula_copy;
        }
    }
}

template<LogicQuantifier T>
auto gen_cold_closed_formulas(int max_quantifier_depth) -> generator<LogicFormula<T>> {
    for (int quantifier_depth = 0; quantifier_depth < max_quantifier_depth; ++quantifier_depth) {
        for (LogicFormula<T> &&formula : gen_cold_formulas<T>(quantifier_depth)) {
            for (auto &&step : T::forall_closing()) {
                auto formula_copy = formula;
                formula_copy.quantifiers.push_back(step);
                co_yield formula_copy;
            }
        }
    }
}


template<LogicQuantifier T>
auto check_formulas_on_graph(wl::SmallGraph const& graph, int max_q=3) -> EvaluatorLogicPaths<T>::evaluated_t {
    auto evaluator = EvaluatorLogicPaths<T>(graph);

    // std::cout << (int)graph.number_of_vertices() << "-- " << "\n";
    int i = 0;
    for (auto &&formula : gen_cold_closed_formulas<T>(max_q)) {
        // make here a progress bar where i is the progress and it is updated every 2 seconds (so that the previous number is overwritten)
        // if (i % 500 == 0) {
        //     std::cout << "\r" << i << " " << std::flush;
        // }
        
        evaluator.evaluate({~0ULL, formula});
        i += 1;
    }
    // std::cout << "\n" << std::flush;

    return std::move(evaluator.results);
}

template<LogicQuantifier T>
auto compute_type(wl::SmallGraph const& graph, int rank) -> std::vector<LogicFormula<T>> {
    auto modeled_formulas = std::vector<LogicFormula<T>>{};

    auto evaluation = check_formulas_on_graph<T>(graph, rank);
    for (auto &&[args, result] : evaluation) {
        if (result && std::get<0>(args) == ~0ULL) {
            modeled_formulas.emplace_back(std::move(std::get<1>(args)));
        }
    }

    std::sort(modeled_formulas.begin(), modeled_formulas.end());
    return modeled_formulas;
}

template<LogicQuantifier T>
auto compute_type_set(wl::SmallGraph const& graph, int rank) -> std::unordered_set<LogicFormula<T>> {
    auto modeled_formulas = std::unordered_set<LogicFormula<T>>{};

    auto evaluation = check_formulas_on_graph<T>(graph, rank);
    for (auto &&[args, result] : evaluation) {
        if (result && std::get<0>(args) == ~0ULL) {
            modeled_formulas.emplace(std::move(std::get<1>(args)));
        }
    }

    return modeled_formulas;
}

} // namespace wl