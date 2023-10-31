#include <utility>

#include "unlabelled_graph.hpp"


namespace wl::ops {
    
class BilabeledGraph {
public:
    // constants 




    BilabeledGraph(SmallGraph const& graph, colvec_t const& in_labels, colvec_t const& out_labels)
        : graph(graph), in_labels(in_labels), out_labels(out_labels) {}

    auto get_rank() -> std::pair<size_t, size_t> {
        return {in_labels.size(), out_labels.size()};
    }

    SmallGraph graph;
    colvec_t in_labels;
    colvec_t out_labels;
};


} // namespace wl
