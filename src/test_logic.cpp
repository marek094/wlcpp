#include "logic.hpp"
#include "unlabelled_graph.hpp"


#include <iostream>
#include <string>
#include <exception>
#include <memory>


// Helper macro for correct expansion of __VA_ARGS__
#define EXPAND(x) x

// Macro to convert line number to string
#define STRINGIZE(x) STRINGIZE2(x)
#define STRINGIZE2(x) #x

// Updated test_true and test_false macros to include line number
#define test_true(cond, ...) EXPAND(test_true_(cond, "test_true(" STRINGIZE(cond) ") at line " STRINGIZE(__LINE__), ##__VA_ARGS__))
#define test_false(cond, ...) EXPAND(test_true_(!cond, "test_false(" STRINGIZE(cond) ") at line " STRINGIZE(__LINE__), ##__VA_ARGS__))

auto test_true_(bool cond, std::string code_str, std::string msg = "") {
    if (! cond) {
        std::ostringstream oss;
        oss << code_str << " failed." << msg << "\n" << std::flush;
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



int main() try {
    test_basics();
    test_quantifier();

    std::cout << "All tests passed.\n";

} catch (std::exception& e) {
    std::cerr << "Exception caught: " << e.what() << "\n";
    return 1;
} catch (...) {
    std::cerr << "Unknown exception caught.\n";
    return 1;
}
