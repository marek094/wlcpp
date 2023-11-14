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
#include <concepts>
#include <iomanip>

namespace wl {


struct hash_array {
    template<typename T, size_t N>
    auto operator()(std::array<T, N> const& arr) const -> size_t {
        size_t seed = arr.size();
        for (auto&& i : arr) {
            seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};


// template<uint64_t St, uint64_t Nd, uint64_t... Rest>
// struct is_nondecreasing {
//     static constexpr bool value = (St <= Nd) && is_nondecreasing<Nd, Rest...>::value;
// };


// template<uint64_t St, uint64_t Nd>
// struct is_nondecreasing<St, Nd> {
//     static constexpr bool value = (St <= Nd);
// };

// template<uint64_t St>
// struct is_nondecreasing<St> {
//     static constexpr bool value = true;
// };

// template<>
// struct is_nondecreasing<> {
//     static constexpr bool value = true;
// };


template<uint64_t... Rest>
struct is_nondecreasing;

template<uint64_t St, uint64_t Nd, uint64_t... Rest>
struct is_nondecreasing<St, Nd, Rest...> {
    static constexpr bool value = (St <= Nd) && is_nondecreasing<Nd, Rest...>::value;
};

template<uint64_t St>
struct is_nondecreasing<St> {
    static constexpr bool value = true;
};

template<>
struct is_nondecreasing<> {
    static constexpr bool value = true;
};





namespace logic {

constexpr uint64_t Dymanic = ~0ULL;

struct ElementBase {
    virtual auto to_string() const -> std::string = 0;

    virtual auto to_ostream(std::ostream& os) const -> std::ostream& {
        os << this->to_string();
        return os;
    }

    using sigma_struct_t = wl::SmallGraph;

    template<uint64_t NumFreeVariables>
    using variable_assignment_t = std::array<wl::SmallGraph::vertex_type, NumFreeVariables>;


    ElementBase() = default;


    template<typename ...VariableAssigment>
    auto evaluate(sigma_struct_t const& s, VariableAssigment... va) -> bool {
        return this->evaluate_impl(s, variable_assignment_t<sizeof...(va)>{static_cast<wl::SmallGraph::vertex_type>(va)...});
    }

    virtual auto evaluate_impl(sigma_struct_t const& s, variable_assignment_t<0> const& variable_assignment) -> bool {
        throw std::runtime_error("Not implemented (evaluate<0>) for `" + this->to_string() + "`");
        assert(false);
        return false;
    }


    virtual auto evaluate_impl(sigma_struct_t const& s, variable_assignment_t<1> const& variable_assignment) -> bool {
        throw std::runtime_error("Not implemented (evaluate<1>) for `" + this->to_string() + "`");
        assert(false);
        return false;
    }


    virtual auto evaluate_impl(sigma_struct_t const& s, variable_assignment_t<2> const& variable_assignment) -> bool {
        throw std::runtime_error("Not implemented (evaluate<2>) for `" + this->to_string() + "`");
        assert(false);
        return false;
    }

    virtual auto clear_cache() -> void {
        throw std::runtime_error("Not implemented (clear_cache) for `" + this->to_string() + "`");
        assert(false);
    }

};


template<uint64_t NumFreeVariables, uint64_t NumChildren, typename Derived>
struct FormulaElementBase : public ElementBase {
    using sigma_struct_t = wl::SmallGraph;
    using variable_assignment_t = std::array<wl::SmallGraph::vertex_type, NumFreeVariables>;

    static constexpr uint64_t kNumFreeVariables = NumFreeVariables;

    using ChildrenArray = std::conditional_t<NumChildren == Dymanic, 
        std::vector<std::shared_ptr<ElementBase>>, 
        std::array<std::shared_ptr<ElementBase>, NumChildren>>;

    FormulaElementBase() {
        children.fill(nullptr);
    }

    FormulaElementBase(ChildrenArray children) : children(std::move(children)) {
    }

    auto evaluate_impl(sigma_struct_t const& s, variable_assignment_t const& variable_assignment) -> bool override {
        if (auto it = evaluation.find(variable_assignment); it != evaluation.end()) {
            return it->second;
        }

        auto result = static_cast<Derived*>(this)->evaluate_(s, variable_assignment);
        evaluation[variable_assignment] = result;
        return result;
    }

    auto num_free_variables() const -> uint64_t {
        return NumFreeVariables;
    }

    auto clear_cache() -> void override {
        evaluation.clear();
        for (auto &&child : children) {
            child->clear_cache();
        }
    }

    std::unordered_map<variable_assignment_t, bool, hash_array> evaluation;
    ChildrenArray children;
};


struct True : public FormulaElementBase<0, 0, True> {
    auto to_string() const -> std::string override {
        return "T";
    }

    auto evaluate_(sigma_struct_t const&, variable_assignment_t const&) -> bool {
        return true;
    }
};

struct False : public FormulaElementBase<0, 0, False> {
    auto to_string() const -> std::string override {
        return "F";
    }

    auto evaluate_(sigma_struct_t const&, variable_assignment_t const&) -> bool {
        return false;
    }
};


struct Adj : public FormulaElementBase<2, 0, Adj> {
    auto to_string() const -> std::string override {
        return "A";
    }

    auto evaluate_(sigma_struct_t const& s, variable_assignment_t const& variable_assignment) -> bool {
        auto [v1, v2] = variable_assignment;
        auto const& adjs = s.adj_list[v1];
        return std::find(adjs.begin(), adjs.end(), v2) != adjs.end(); 
    }
};

struct NonAdj : public FormulaElementBase<2, 0, NonAdj> {

    auto to_string() const -> std::string override {
        return "~A";
    }

    auto evaluate_(sigma_struct_t const& s, variable_assignment_t const& variable_assignment) -> bool {
        auto [v1, v2] = variable_assignment;
        auto const& adjs = s.adj_list[v1];
        return std::find(adjs.begin(), adjs.end(), v2) == adjs.end(); 
    }

};

struct Eq : public FormulaElementBase<2, 0, Eq> {

    auto to_string() const -> std::string override {
        return "=";
    }

    auto evaluate_(sigma_struct_t const&, variable_assignment_t const& variable_assignment) -> bool {
        auto [v1, v2] = variable_assignment;
        return v1 == v2;
    }

};


struct NonEq : public FormulaElementBase<2, 0, NonEq> {
    
    auto to_string() const -> std::string override {
        return "~=";
    }

    auto evaluate_(sigma_struct_t const&, variable_assignment_t const& variable_assignment) -> bool {
        auto [v1, v2] = variable_assignment;
        return v1 != v2;
    }
    
};

struct And : public FormulaElementBase<2, 2, And> {

    using FormulaElementBase<2, 2, And>::FormulaElementBase;

    And(std::shared_ptr<ElementBase> child1, std::shared_ptr<ElementBase> child2) : FormulaElementBase<2, 2, And>({std::move(child1), std::move(child2)}) {}

    auto to_string() const -> std::string override {
        return this->children[0]->to_string() + "&" + this->children[1]->to_string();
    }

    auto evaluate_(sigma_struct_t const& s, variable_assignment_t const& variable_assignment) -> bool {
        return this->children[0]->evaluate_impl(s, variable_assignment) && this->children[1]->evaluate_impl(s, variable_assignment);
    }
    
};

struct Or : public FormulaElementBase<2, 2, Or> {

    using FormulaElementBase<2, 2, Or>::FormulaElementBase;

    Or(std::shared_ptr<ElementBase> child1, std::shared_ptr<ElementBase> child2) : FormulaElementBase<2, 2, Or>({std::move(child1), std::move(child2)}) {}

    auto to_string() const -> std::string override {
        return this->children[0]->to_string() + "|" + this->children[1]->to_string();
    }

    auto evaluate_(sigma_struct_t const& s, variable_assignment_t const& variable_assignment) -> bool {
        return this->children[0]->evaluate_impl(s, variable_assignment) || this->children[1]->evaluate_impl(s, variable_assignment);
    }
    
};

template<uint64_t NumFreeVariables, uint64_t... Variables>
struct Fgt : public FormulaElementBase<NumFreeVariables, 1, Fgt<NumFreeVariables, Variables...>> {
    using parent_t = FormulaElementBase<NumFreeVariables, 1, Fgt<NumFreeVariables, Variables...>>;
    using variable_assignment_t = typename parent_t::variable_assignment_t;
    using sigma_struct_t = typename parent_t::sigma_struct_t;

    using fgt_variable_assignment_t = FormulaElementBase<NumFreeVariables - sizeof...(Variables), 1, void>::variable_assignment_t;
    
    static_assert(((Variables < NumFreeVariables) && ...) , "Variables must be smaller than NumFreeVariables");
    static_assert(sizeof...(Variables) <= NumFreeVariables, "Too many variables");
    static_assert(is_nondecreasing<Variables...>::value, "Variables must be in nondecreasing order");

    using parent_t::parent_t;

    explicit Fgt(std::shared_ptr<ElementBase> child) : parent_t({std::move(child)}) {}

    auto to_string() const -> std::string override {
        return this->children[0]->to_string();
    }

    auto evaluate_(sigma_struct_t const& s, variable_assignment_t const& variable_assignment) -> bool {
        return this->children[0]->evaluate_impl(s, forget_variables(variable_assignment));
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



struct Exists : public FormulaElementBase<1, 1, Exists> {
    using parent_t = FormulaElementBase<1, 1, Exists>;
    using variable_assignment_t = typename parent_t::variable_assignment_t;
    using sigma_struct_t = typename parent_t::sigma_struct_t;

    using parent_t::parent_t;

    explicit Exists(std::shared_ptr<ElementBase> child) : parent_t({std::move(child)}) {}

    auto to_string() const -> std::string override {
        return "Ex " + this->children[0]->to_string();
    }

    auto to_ostream(std::ostream& os) const -> std::ostream& override {
        os << "Ex ";
        return this->children[0]->to_ostream(os);
    }

    auto evaluate_(sigma_struct_t const& s, variable_assignment_t const& variable_assignment) -> bool {
        auto [v1] = variable_assignment;
        for (auto &&vtx : s.all({})) {
            if (this->children[0]->evaluate(s, vtx, v1)) {
                return true;
            }
        }
        return false;
    }

    static auto forall(parent_t::ChildrenArray children) -> generator<std::shared_ptr<ElementBase>> {
        auto quantifier = std::make_shared<Exists>(std::move(children));
        co_yield quantifier;
    }
};


template<uint64_t Lo, uint64_t Hi>
struct ExistCount_ : public FormulaElementBase<1, 1, ExistCount_<Lo, Hi>> {
    using parent_t = FormulaElementBase<1, 1, ExistCount_<Lo,Hi>>;
    using variable_assignment_t = typename parent_t::variable_assignment_t;
    using sigma_struct_t = typename parent_t::sigma_struct_t;

    int count = 0;

    ExistCount_(int count, std::shared_ptr<ElementBase> child) : count(count), parent_t({std::move(child)}) {}

    ExistCount_(int count, parent_t::ChildrenArray children) : count(count), parent_t(std::move(children)) {}

    auto to_string() const -> std::string override {
        return "E" + std::to_string(count) + "x " + this->children[0]->to_string();
    }

    auto to_ostream(std::ostream& os) const -> std::ostream& override {
        os << "E" << count << "x ";
        return this->children[0]->to_ostream(os);
    }

    auto evaluate_(sigma_struct_t const& s, variable_assignment_t const& variable_assignment) -> bool {
        auto [v1] = variable_assignment;
        
        int actual_count = 0;
        for (auto &&vtx : s.all({})) {
            actual_count += this->children[0]->evaluate(s, vtx, v1);
            if (actual_count > count) {
                return false;
            }
        }
        return (actual_count == count);
    }

    static auto forall(parent_t::ChildrenArray children) -> generator<std::shared_ptr<ElementBase>> {
        for (uint64_t i = Lo; i <= Hi; ++i) {
            auto quantifier = std::make_shared<ExistCount_>(i, children);
            co_yield quantifier;
        }
    }
};

using ExistCount = ExistCount_<0, 10>;


struct AtomicQuantifier : public FormulaElementBase<1, Dymanic, AtomicQuantifier> {
    using parent_t = FormulaElementBase<1, Dymanic, AtomicQuantifier>;
    using variable_assignment_t = typename parent_t::variable_assignment_t;
    using sigma_struct_t = typename parent_t::sigma_struct_t;

    explicit AtomicQuantifier(parent_t::ChildrenArray children) : parent_t(std::move(children)) {
        assert (this->children.size() >= 1);
    }

    auto to_string() const -> std::string override {
        std::ostringstream oss;
        oss << "P" << (this->children.size()-1) << "x";
        if (this->children.size() >= 1) {
            oss << " " << this->children[0]->to_string() << "->";
            oss << "(";
            for (int i = 1; i < this->children.size(); ++i) {
                oss << this->children[i]->to_string();
                if (i < this->children.size()-1) {
                    oss << ",";
                }
            }
            oss << ")";
        }
        return std::move(oss.str());
    }

    auto to_ostream(std::ostream& os) const -> std::ostream& override {
        os << "P" << (this->children.size()-1) << "x";
        if (this->children.size() >= 1) {
            os << " ";
            this->children[0]->to_ostream(os);
            os << "->";
            os << "(";
            for (int i = 1; i < this->children.size(); ++i) {
                this->children[i]->to_ostream(os);
                if (i < this->children.size()-1) {
                    os << ",";
                }
            }
            os << ")";
        }
        return os;
    }


    auto evaluate_(sigma_struct_t const& s, variable_assignment_t const& variable_assignment) -> bool {
        auto [v1] = variable_assignment;
        assert (this->children.size() >= 1);
        auto&& atomic = this->children.front();

        int64_t wanted_c = this->children.size()-1;

        int64_t c_sum = 0;
        for (auto &&vtx : s.all({})) {
            if (atomic->evaluate(s, vtx, v1)) {
                bool is_disjunction = false;
                for (int ci = 1; ci < this->children.size(); ++ci) {
                    if (this->children[ci]->evaluate(s, vtx)) {
                        is_disjunction = true;
                        c_sum += ci;
                        if (c_sum > wanted_c) {
                            // std::cout << "false" << "\t" <<std::setw(2)<< (int)v1 << " (c_sum  > wanted_c) " << this->to_string() << " " << c_sum << " > " << wanted_c << std::endl;
                            return false;
                        }
                    }
                }
                if (!is_disjunction) {
                    // std::cout << "false" << "\t" <<std::setw(2)<< (int)v1 << " (!is_disjunction) " << this->to_string() << std::endl;
                    return false;
                }
            }
        }

        // std::cout << std::boolalpha << (c_sum == wanted_c) << "\t" << std::setw(2) << (int)v1 << " (c_sum == wanted_c) " << this->to_string() << " " << c_sum << " == " << wanted_c << std::endl;
        return (c_sum == wanted_c);
    }

    static auto forall(parent_t::ChildrenArray children) -> generator<std::shared_ptr<ElementBase>> {
        // auto&& atomic = children.front();
        // for (uint64_t subset = 1; subset < (1ULL<<(children.size())); ++subset) {
        //     std::cout << "\nf\t" << children.size() << "\t" << subset  << std::endl;   
        //     auto children2 = parent_t::ChildrenArray{};
        //     children2.push_back(atomic);
        //     for (uint64_t i = 1; i < children.size(); ++i) {
        //         if (subset & (1ULL << (i-1))) {
        //             children2.push_back(children[i]);
        //         }
        //     }
        //     std::shared_ptr<ElementBase> result = std::make_shared<AtomicQuantifier>(children2);
        //     std::cout << "\nff\t" << children2.size() << "\t" << result->to_string() <<  std::endl;
        //     co_yield result;
        // }

        for (uint64_t subset_size = 2; subset_size <= children.size(); ++subset_size) {
            co_yield std::make_shared<AtomicQuantifier>(parent_t::ChildrenArray{children.begin(), children.begin()+subset_size});
        }
    }
};



template<uint64_t FreeVarIndex, std::derived_from<ElementBase> Quantifier>
auto generate_formulas_helper(int rank) -> generator<std::shared_ptr<ElementBase>> {
    if (rank == 0) {
        co_yield std::make_shared<Fgt<1, 0>>(std::make_shared<True>());
        co_yield std::make_shared<Fgt<1, 0>>(std::make_shared<False>());
        co_return;
    }

    
    auto all_formulas = std::vector<std::shared_ptr<ElementBase>>{};
    all_formulas.push_back(std::make_shared<Adj>());

    for (auto&& formula : generate_formulas_helper<(FreeVarIndex+1)%2, Quantifier>(rank - 1)) {
        all_formulas.push_back(formula);
    }

    for (auto&& formula : Quantifier::forall(all_formulas)) {
        co_yield formula;
    }

    // if (rank == 3) {
    //     std::cout << "Z" << all_formulas.size() << std::endl;
    //     all_formulas.erase(all_formulas.begin()+10, all_formulas.end());
    //     for (auto&& formula : Quantifier::forall(all_formulas)) {
    //         std::cout << "R\t" << formula->to_string() << std::endl;
    //     }
    // }


    // for (auto&& formula : all_formulas) {
        
    //     co_yield formula;

    //     // edge
    //     auto fadj = std::make_shared<Adj>();
    //     auto fformula = std::make_shared<Fgt<2, FreeVarIndex>>(formula);
    //     auto fand = std::make_shared<And>(fadj, fformula);

    //     auto children = typename Quantifier::ChildrenArray{fand};
    //     assert( children.size() == 1);
    //     assert( children.front() != nullptr);
    //     for (auto quantifier : Quantifier::forall(children)) {
    //         co_yield quantifier;
    //     }
    // }

}

template<std::derived_from<ElementBase> Quantifier>
auto generate_formulas(int rank) -> generator<std::shared_ptr<ElementBase>> {
    return generate_formulas_helper<0, Quantifier>(rank);
}


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