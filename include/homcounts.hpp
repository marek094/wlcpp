
// def compute_pathwidth_one_homvec(graph : nx.Graph, num_vertices : int):

//     A = nx.adjacency_matrix(graph).todense().astype(np.int64)
//     assert isinstance(A, np.ndarray)
//     D = np.diag(np.sum(A, axis=1)).astype(np.int64)

//     homvector = []
//     for n in range(num_vertices):
//         for bitvec in range(2**n+1):
//             # 0 corresponds to A and 1 corresponds to D
//             hom_count = np.eye(A.shape[0], dtype=np.int64)
//             for j in range(n):
//                 if (bitvec >> j) & 1:
//                     hom_count = hom_count @ D
//                 else:
//                     hom_count = hom_count @ A

//             homvector.append(np.sum(hom_count))
    
//     return np.array(homvector)


#pragma once

#include <vector>
#include <utility>

#include <Eigen/Dense>

#include "unlabelled_graph.hpp"

namespace wl {

auto compute_pathwidth_one_homvec(SmallGraph const& graph, SmallGraph::type num_vertices) -> std::vector<unsigned long long> {
    auto A = graph.to_adjacency_matrix();
    auto D = graph.to_degree_matrix();

    using eigen_matrix_t = decltype(A);

    std::vector<unsigned long long> homvector;
    for (int n = 0; n < num_vertices; ++n) {
        for (int bitvec = 0; bitvec < (1 << n); ++bitvec) {
            // 0 corresponds to A and 1 corresponds to D
            eigen_matrix_t hom_count = eigen_matrix_t::Identity(A.rows(), A.cols());
            for (int j = 0; j < n; ++j) {
                if ((bitvec >> j) & 1) {
                    hom_count = hom_count * D;
                } else {
                    hom_count = hom_count * A;
                }
            }
            homvector.push_back(hom_count.sum());
        }
    }
    
    return homvector;
}

auto compute_pathwidth_one_homvec_v2(SmallGraph graph, SmallGraph::type num_vertices) -> std::vector<unsigned long long> {
    auto A = graph.to_adjacency_matrix();
    auto D = graph.to_degree_matrix();

    using eigen_matrix_t = decltype(A);

    std::vector<eigen_matrix_t> matvector = {eigen_matrix_t::Identity(A.rows(), A.cols())};
    std::vector<eigen_matrix_t> matvector2;

    std::vector<unsigned long long> homvector;
    for (int n = 0; n < num_vertices; ++n) {
        for (auto &&mat : matvector) {
            auto optA = mat * A;
            auto optD = mat * D;
            if (n < num_vertices-1) {
                matvector2.push_back(optA);
                matvector2.push_back(optD);
            }
            homvector.push_back(optA.sum());
            homvector.push_back(optD.sum());
        }

        std::swap(matvector, matvector2);
        matvector2.clear();
    }
    
    return homvector;
}



} // namespace wl
