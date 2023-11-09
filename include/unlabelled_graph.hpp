#pragma once

#include "ref_generator.hpp"

#include <Eigen/Dense>

#include <vector>
#include <utility>
#include <algorithm>
#include <ranges>
#include <set>
#include <iostream>


namespace wl {


template <typename vertex_t>
class Graph {
public:
    using EigenMatrixXu64 = Eigen::Matrix<uint64_t, Eigen::Dynamic, Eigen::Dynamic>;
    using vertex_type = vertex_t;

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

    auto number_of_vertices() const -> size_t {
        assert (num_vertices == adj_list.size());
        return static_cast<size_t>(num_vertices);
        // return adj_list.size();
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

    auto add_edge(vertex_t u, vertex_t v) -> void {
        adj_list.resize(num_vertices);

        edges.emplace_back(u, v);
        if (u >= adj_list.size()) {
            adj_list.resize(u+1, {});
        }
        adj_list[u].push_back(v);
        
        if (v >= adj_list.size()) {
            adj_list.resize(v+1, {});
        }
        
        adj_list[v].push_back(u);
    }

    auto all(vertex_t) const -> generator<vertex_t> {
        for (vertex_t node = 0; node < adj_list.size(); ++node) {
            co_yield node;
        }
    }

    auto all_adj(vertex_t node) const -> generator<vertex_t> {
        assert (node < adj_list.size());
        for (auto &&adj : adj_list[node]) {
            if (adj != node) {
                co_yield adj;
            }
        }
    }

    auto all_nonadj(vertex_t node) const -> generator<vertex_t> {
        assert (node < adj_list.size());
        auto set = std::set<vertex_t>{};
        for (auto &&adj : adj_list[node]) {
            set.insert(adj);
        }
        for (vertex_t nonadj = 0; nonadj < adj_list.size(); ++nonadj) {
            if ((nonadj != node) && (set.find(nonadj) == set.end())) {
                co_yield nonadj;
            }
        }
    }


    auto all_edges() const -> generator<std::pair<vertex_t, vertex_t>> {
        for (auto &&[u,v] : edges) {
            co_yield std::pair{u, v};
            co_yield std::pair{v, u};
        }
    }

    auto print_python() const -> void {
        std::cout << "nx.from_edgelist([";
        for (auto &&[u, v] : edges) {
            std::cout << "(" 
                << static_cast<int>(u) << ", " 
                << static_cast<int>(v) << "), ";
        }
        std::cout << "])\n";
    }

    auto print_as_txt_line(std::ostream& os) const -> void {
        assert(labels.empty());
        int n = this->number_of_vertices();
        int m = edges.size();
        int has_labels = false;
        os << n << " " << m << " " << has_labels << " ";
        for (auto &&[u, v] : edges) {
            os << static_cast<int>(u) << " " << static_cast<int>(v) << " ";
        }
        os << "\n";
    }

    auto degree_vector() const -> std::vector<unsigned long long> {
        auto degrees = std::vector<unsigned long long>(this->number_of_vertices(), 0);
        for (auto &&[u, v] : edges) {
            degrees[u] += 1;
            degrees[v] += 1;
        }
        return degrees;
    }

    using type = vertex_t;
    std::vector<unsigned long long> labels;
    unsigned long long num_vertices;

    std::vector<std::vector<vertex_t>> adj_list;
private:
    std::vector<std::pair<vertex_t, vertex_t>> edges;
};


// graph with less than 256 vertices
using SmallGraph = Graph<unsigned char>;


} // namespace wl