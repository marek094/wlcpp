#include "logic.hpp"
#include "unlabelled_graph.hpp"


#include <iostream>
#include <string>
#include <exception>
#include <memory>
#include <source_location>
#include <map>
#include <string_view>


// Helper macro for correct expansion of __VA_ARGS__
#define EXPAND(x) x

// Macro to convert line number to string
#define STRINGIZE(x) STRINGIZE2(x)
#define STRINGIZE2(x) #x

// Updated test_true and test_false macros to include line number
#define test_val_base(cond, cond_s, spec, ...) EXPAND(test_true_(cond, spec, STRINGIZE(cond_s), std::source_location::current(), ##__VA_ARGS__))
#define test_true(cond, ...) EXPAND(test_val_base(cond, cond, "true", ##__VA_ARGS__))
#define test_false(cond, ...) EXPAND(test_val_base(!cond, cond, "false", ##__VA_ARGS__))


auto test_true_(bool cond, std::string spec, std::string cond_str, std::source_location location, std::string msg = "") -> void {
    if (! cond) {
        std::ostringstream oss;
        if (! msg.empty()) {
            oss << " [" << msg << "]: ";
        }
        oss << "`test_" << spec << "(" << cond_str << ")`";
        oss << " in function `" << location.function_name() << "`";
        oss << " at " << location.file_name() << ":" << location.line() << ":" << location.column();
        oss << std::endl;
        throw std::runtime_error(std::move(oss.str()));
    }
}


using std::make_shared;
using std::shared_ptr;
using std::vector;
using std::array;
using namespace wl::logic;


void test_basics() {

    wl::SmallGraph g;
    g.add_edge(0, 1);

    
    shared_ptr<ElementBase> ftrue = make_shared<True>();

    test_true( ftrue->to_string() == "T" );
    test_true( ftrue->evaluate(g) );

    shared_ptr<ElementBase> ffalse = make_shared<False>();

    test_true( ffalse->to_string() == "F" );
    test_false( ffalse->evaluate(g) );



    shared_ptr<ElementBase> fadj = make_shared<Adj>();
    
    test_true( fadj->to_string() == "A" );
    test_true( fadj->evaluate(g, 0, 1) );
    test_false( fadj->evaluate(g, 0, 2) );

    shared_ptr<ElementBase> fnonadj = make_shared<NonAdj>();

    test_true( fnonadj->to_string() == "~A" );
    test_false( fnonadj->evaluate(g, 0, 1) );
    test_true( fnonadj->evaluate(g, 0, 2) );

    
    shared_ptr<ElementBase> feq = make_shared<Eq>();

    test_true( feq->to_string() == "=" );
    test_true( feq->evaluate(g, 0, 0) );
    test_false( feq->evaluate(g, 0, 1) );


    shared_ptr<ElementBase> fneq = make_shared<NonEq>();

    test_true( fneq->to_string() == "~=" );
    test_false( fneq->evaluate(g, 0, 0) );
    test_true( fneq->evaluate(g, 0, 1) );



    shared_ptr<ElementBase> fatp = make_shared<And>(fadj, fneq);

    test_true( fatp->to_string() == "A&~=" );
    test_true( fatp->evaluate(g, 0, 1) );
    test_false( fatp->evaluate(g, 0, 0) );

    

    shared_ptr<ElementBase> ffgt_true = make_shared<Fgt<2, 0, 1>>(ftrue);

    test_true( (Fgt<2, 0, 1>::kNumFreeVariables == 2) , "x");
    test_true( ffgt_true->to_string() == "T" );
    test_true( ffgt_true->evaluate(g, 1, 1) );
    test_true( ffgt_true->evaluate(g, 0, 0) );
    test_true( ffgt_true->evaluate(g, 0, 1) );
    test_true( ffgt_true->evaluate(g, 1, 0) );
    

    shared_ptr<ElementBase> ftrue_and_true = make_shared<And>(ffgt_true, ffgt_true);
    test_true( ftrue_and_true->to_string() == "T&T" );
    test_true( ftrue_and_true->evaluate(g, 1, 1) );

}


void test_quantifier() {

    wl::SmallGraph g;
    g.add_edge(0, 1);
    g.add_edge(1, 2);

    wl::SmallGraph g2;
    g2.add_edge(0, 1);

    wl::SmallGraph g3;

    shared_ptr<ElementBase> ftrue = make_shared<True>();
    shared_ptr<ElementBase> ffgt_true = make_shared<Fgt<2, 0, 1>>(ftrue);

    {
        shared_ptr<ElementBase> fexists = make_shared<Exists>(ffgt_true);
        test_true( fexists->to_string() == "Ex T" );
        test_true( fexists->evaluate(g, 0) );

        fexists->clear_cache();
        test_true( fexists->evaluate(g2, 0) );

        fexists->clear_cache();
        test_false( fexists->evaluate(g3, 0) );
    }

    {
        shared_ptr<ElementBase> fexistcount = make_shared<ExistCount>(3, ffgt_true);
        test_true( fexistcount->to_string() == "E3x T" );
        test_true( fexistcount->evaluate(g, 0) );

        fexistcount->clear_cache();
        test_false( fexistcount->evaluate(g2, 0) );

        fexistcount->clear_cache();
        test_false( fexistcount->evaluate(g3, 0) );
    }

    {
        shared_ptr<ElementBase> fexistcount = make_shared<ExistCount>(2, ffgt_true);
        test_true( fexistcount->to_string() == "E2x T" );
        test_false( fexistcount->evaluate(g, 0) );

        fexistcount->clear_cache();
        test_true( fexistcount->evaluate(g2, 0) );

        fexistcount->clear_cache();
        test_false( fexistcount->evaluate(g3, 0) );
    }

}


void test_generators() {

    auto elements = std::vector<shared_ptr<ElementBase>>{};
    for (auto formula : wl::logic::generate_formulas(1)) {
        elements.push_back(formula);
    }

    test_true( elements.size() == 1 );
    test_true( elements[0]->to_string() == "Ex T" );


    // K_4
    wl::SmallGraph g;
    g.add_edge(0, 1);
    g.add_edge(0, 2);
    g.add_edge(0, 3);
    g.add_edge(1, 2);
    g.add_edge(1, 3);
    g.add_edge(2, 3);


    test_true( elements[0]->evaluate(g, 0) );

}





int main(int argc, char** argv) try {

    constexpr auto test_list = std::array{
        std::pair("basic", test_basics),
        std::pair("quantifier", test_quantifier),
        std::pair("generators", test_generators),
    };

    auto args = std::vector<std::string_view>{argv + 1, argv + argc};
    
    if (args.empty()) {
        for (auto&& [name, test] : test_list) {
            args.push_back(name);
        }
    } 
    
    for (auto&& arg : args) {
        auto it = std::find_if(test_list.begin(), test_list.end(), [&](auto&& p) {
            return p.first == arg;
        });
        if (it == test_list.end()) {
            std::cerr << "Unknown test: " << arg << "\n";
            return 1;
        }
        auto&& [name, test] = *it;
        std::cout << "Running test: " << name << "\n";
        test();
    }


    std::cout << "All tests passed.\n";
    return 0;

} catch (std::exception& e) {
    std::cerr << "Exception caught: " << e.what() << "\n";
    return 1;
} catch (...) {
    std::cerr << "Unknown exception caught.\n";
    return 1;
}
