#pragma once

#include <vector>
#include <utility>
#include <algorithm>
#include <ranges>

#include <Eigen/Dense>

namespace wl {


template <typename vertex_t>
class Graph {
public:
    using EigenMatrixXu64 = Eigen::Matrix<uint64_t, Eigen::Dynamic, Eigen::Dynamic>;

    auto to_adjacency_matrix() const -> EigenMatrixXu64 {
        auto n = this->number_of_vertices();
        EigenMatrixXu64 A = EigenMatrixXu64::Zero(n, n);
        for (auto &&[u, v] : edges) {
            A(u, v) = 1;
            A(v, u) = 1;
        }
        return A;
    }

    auto to_degree_matrix() const -> EigenMatrixXu64 {
        auto n = this->number_of_vertices();
        EigenMatrixXu64 D = EigenMatrixXu64::Zero(n, n);
        for (auto &&[u, v] : edges) {
            D(u, u) += 1;
            D(v, v) += 1;
        }
        return D;
    }

    auto number_of_vertices() const -> vertex_t {
        return num_vertices;
        // if (edges.empty()) {
        //     return 0;
        // }

        // // if (num_vertices == -1) {
        // auto projection = [](auto &&edge) { return std::max<vertex_t>(edge.first, edge.second); };
        // auto largest_label = std::ranges::max(edges, {}, projection);
        // auto num_vertices = projection(largest_label) + 1;
        // // }

        // return num_vertices;
    }


    using type = vertex_t;
    std::vector<std::pair<vertex_t, vertex_t>> edges;
    std::vector<unsigned long long> labels;
    unsigned long long num_vertices;
};


// graph with less than 256 vertices
using SmallGraph = Graph<unsigned char>;


} // namespace wl