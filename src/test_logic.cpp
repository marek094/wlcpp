#include "logic.hpp"
#include "unlabelled_graph.hpp"


#include <iostream>
#include <string>
#include <exception>
#include <memory>


// unit tests
auto test_true(bool cond, std::string msg = "") {
    if (! cond) {
        std::cout << "Test failed: " << msg << "\n" << std::flush;
        throw std::runtime_error("Test failed.");
    }
}


void test_basics() {

    wl::SmallGraph g;
    g.add_edge(0, 1);

    std::shared_ptr<wl::logic::FormulaElementBase<0>> f_true = std::make_shared<wl::logic::True>();
    test_true( f_true->to_string() == "T" );
    test_true( f_true->evaluate(g, {}) );

    std::shared_ptr<wl::logic::FormulaElementBase<0>> f_false = std::make_shared<wl::logic::False>();
    test_true( f_false->to_string() == "F" );
    test_true( ! f_false->evaluate(g, {}) );

    std::shared_ptr<wl::logic::FormulaElementBase<2>> f_adj = std::make_shared<wl::logic::Adj>();
    test_true( f_adj->to_string() == "A" );
    test_true( f_adj->evaluate(g, {0, 1}) );
    test_true( ! f_adj->evaluate(g, {0, 2}) );


}



int main() try {
    test_basics();

} catch (std::exception& e) {
    std::cerr << "Exception caught: " << e.what() << "\n";
    return 1;
} catch (...) {
    std::cerr << "Unknown exception caught.\n";
    return 1;
}
