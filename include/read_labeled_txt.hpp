
#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "unlabelled_graph.hpp"

namespace wl {

auto read_graph_from_labeled_txt_file(std::string const& filename, size_t limit = (~0), size_t skip = 0) -> std::vector<SmallGraph> {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open the file: " << filename << '\n';
        return {};
    }

    std::vector<SmallGraph> graphs;
    std::string line;

    unsigned long long n, m;
    int has_labels;

    while (file >> n >> m >> has_labels) {
        // std::cout << n << " " << m << " " << has_labels << '\n';
        auto graph = SmallGraph{};
        graph.num_vertices = n;
        
        if (has_labels) {
            graph.labels.reserve(n);
            for (unsigned long long i = 0; i < n; ++i) {
                unsigned long long label;
                file >> label;
                graph.labels.push_back(-label);
            }
        }

        for (unsigned long long i = 0; i < m; ++i) {
            unsigned long long u, v;
            file >> u >> v;
            graph.add_edge(u, v);
        }

        graphs.push_back(graph);

        if (graphs.size() >= limit) {
            break;
        }
    }

    return graphs;
}

} // namespace wl