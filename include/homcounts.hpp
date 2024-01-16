
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

#include "unlabelled_graph.hpp"

#include <boost/math/quaternion.hpp>
#include <Eigen/Dense>

#include <vector>
#include <utility>
#include <complex>
#include <map>
#include <string>
#include <tuple>
#include <array>
#include <cmath>

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



auto compute_pathwidth_one_homvec_v2(SmallGraph graph, SmallGraph::type num_vertices) {
    auto A = graph.to_adjacency_matrix();
    auto D = graph.to_degree_matrix();

    using eigen_matrix_t = decltype(A);

    std::vector<std::tuple<eigen_matrix_t, std::string>> matvector = {std::make_tuple(eigen_matrix_t::Identity(A.rows(), A.cols()), std::string{})};
    std::vector<std::tuple<eigen_matrix_t, std::string>> matvector2;

    std::vector<unsigned long long> homvector;
    std::vector<std::string> homexpr;
    for (int n = 0; n < num_vertices; ++n) {
        for (auto &&[mat, expr] : matvector) {
            auto exprA = expr + "A";
            auto optA = mat * A;
            auto exprD = expr + "D";
            
            if (n < num_vertices-1) {
                auto optD = mat * D;
                matvector2.push_back(std::tuple{optA, exprA});
                matvector2.push_back(std::tuple{optD, exprD});

                homexpr.push_back(exprD);
                homvector.push_back(optD.sum());
            } else {
                homexpr.push_back(exprD);
                // A1 = D1
                homvector.push_back(optA.sum());
            }

            homexpr.push_back(exprA);
            homvector.push_back(optA.sum());
            
        }

        std::swap(matvector, matvector2);
        matvector2.clear();
    }
    
    return std::tuple{homvector, homexpr};
}


auto compute_pathwidth_one_homvec_bases(SmallGraph graph, SmallGraph::type plus,  SmallGraph::type max_Ds) -> std::vector<unsigned long long> {

    auto A = graph.to_adjacency_matrix();
    auto D = graph.to_degree_matrix();

    using eigen_matrix_t = decltype(A);

    std::vector<unsigned long long> homvector;

    std::vector<eigen_matrix_t> matvector = {eigen_matrix_t::Identity(A.rows(), A.cols())};
    {
        for (int n = 0; n < (graph.number_of_vertices() / 2 + graph.number_of_vertices() % 2 + plus/2); ++n) {
            matvector.push_back(matvector.back() * A);
        }
    }

    eigen_matrix_t optD = eigen_matrix_t::Identity(A.rows(), A.cols());
    for (int k = 0; k < max_Ds; ++k) {
        for (int i0 = 0; i0 < (int)matvector.size(); ++i0) {
            for (int i1 = i0; i1 < (int)matvector.size(); ++i1) {
                auto to_hom = matvector[i0] * optD * matvector[i1];
                homvector.push_back(to_hom.sum());
            }
        }
        optD = optD * D;
    }
    
    return homvector;
}






template<typename Int = uint64_t>
auto compute_path_homvec(SmallGraph graph, SmallGraph::type num_vertices) -> std::vector<Int> {
    auto A_was = graph.to_adjacency_matrix();

    using eigen_matrix_t = Eigen::Matrix<Int, Eigen::Dynamic, Eigen::Dynamic>;
    eigen_matrix_t A = eigen_matrix_t::Zero(A_was.rows(), A_was.cols());
    for (int i = 0; i < A_was.rows(); ++i) {
        for (int j = 0; j < A_was.cols(); ++j) {
            A(i, j) = A_was(i, j);
        }
    }

    std::vector<Int> homvector;
    homvector.reserve(num_vertices+1);
    
    eigen_matrix_t product = eigen_matrix_t::Identity(A.rows(), A.cols());
    homvector.push_back(product.sum());
    for (int n = 0; n < num_vertices; ++n) {
        product = product * A;
        homvector.push_back(product.sum());
    }

    return homvector;
}

auto compute_path_homvec_labeled(SmallGraph graph, SmallGraph::type num_vertices) -> std::vector<unsigned long long> {
    auto A = graph.to_adjacency_matrix();

    using eigen_matrix_t = decltype(A);

    std::vector<unsigned long long> homvector;
    homvector.reserve(num_vertices+1);
    
    eigen_matrix_t product = eigen_matrix_t::Identity(A.rows(), A.cols());
    homvector.push_back(product.sum());
    for (int n = 0; n < num_vertices; ++n) {
        product = product * A;
        for (auto val : product.colwise().sum()) homvector.push_back(val);
        std::sort(homvector.end() - product.rows(), homvector.end());
    }

    return homvector;
}



template<typename Int = int64_t>
auto compute_complex_path_homvec(SmallGraph graph, SmallGraph::type num_vertices, long long sign = 1) -> std::vector<std::complex<Int>> {
    auto A = graph.to_adjacency_matrix();

    using complex_matrix_t = Eigen::Matrix<std::complex<Int>, Eigen::Dynamic, Eigen::Dynamic>;

    // A to complex<uint64_t>
    complex_matrix_t Ai = A.cast<std::complex<Int>>();
    assert (Ai.rows() == Ai.cols());
    for (int i = 0; i < Ai.rows(); ++i) {
        Ai(i, i) = std::complex<Int>(0, sign*1);
    }

    std::vector<std::complex<Int>> homvector;
    homvector.reserve(num_vertices+1);
    
    complex_matrix_t product = complex_matrix_t::Identity(Ai.rows(), Ai.cols());
    homvector.push_back(product.sum());
    for (int n = 0; n < num_vertices; ++n) {
        product = product * Ai;
        homvector.push_back(product.sum());
    }

    return homvector;
}


template<typename Matrix>
Matrix fast_power(Matrix base, unsigned long long exponent) {
    Matrix result = Matrix::Identity(base.rows(), base.cols());

    while (exponent > 0) {
        if (exponent & 1) {
            result = result * base;
        }
        base = base * base;
        exponent >>= 1;
    }

    return result;
}



template<typename Int = int64_t>
auto compute_complex_path_homvec_log(SmallGraph graph, SmallGraph::type num_vertices) -> std::vector<std::complex<Int>> {
    auto A = graph.to_adjacency_matrix();
    auto n = A.rows();

    using complex_matrix_t = Eigen::Matrix<std::complex<Int>, Eigen::Dynamic, Eigen::Dynamic>;
    complex_matrix_t Ai = A.cast<std::complex<Int>>() * std::complex<Int>(0, 1);
    
    // compute ceil(nlog_2n)
    n = 2*n;
    unsigned long long nlog2n =  std::ceil((n) * std::log2(n));
    complex_matrix_t powAp0 = fast_power<complex_matrix_t>(Ai, nlog2n);
    complex_matrix_t powAp1 = powAp0 * Ai;
    complex_matrix_t powAp2 = powAp1 * Ai;
    complex_matrix_t powAp3 = powAp2 * Ai;

    // return {powAp0.sum()};
    // return {powAp0.sum(), powAp1.sum()};

    return {powAp0.sum(), powAp1.sum(), powAp2.sum(), powAp3.sum()};
}



auto compute_complex_path_homvec_boost(SmallGraph graph, SmallGraph::type num_vertices) -> std::vector<std::complex<long long>> {
    auto A = graph.to_adjacency_matrix();

    using complex_matrix_t = Eigen::Matrix<std::complex<int64_t>, Eigen::Dynamic, Eigen::Dynamic>;

    // A to complex<uint64_t>
    complex_matrix_t Ai = A.cast<std::complex<int64_t>>();
    assert (Ai.rows() == Ai.cols());
    for (int i = 0; i < Ai.rows(); ++i) {
        for (int j = 0; j < Ai.cols(); ++j) {
            if (A(i, j) != 0) {
                Ai(i, j) = std::complex<int64_t>(2, 0);
            } else if (i == j) {
                Ai(i, j) = std::complex<int64_t>(-1, 1);
            } else {
                Ai(i, j) = std::complex<int64_t>(-1, -1);
            }

        }
    }

    std::vector<std::complex<long long>> homvector;
    homvector.reserve(num_vertices+1);
    
    complex_matrix_t product = complex_matrix_t::Identity(Ai.rows(), Ai.cols());
    homvector.push_back(product.sum());
    for (int n = 0; n < num_vertices; ++n) {
        product = product * Ai;
        homvector.push_back(product.sum());
    }

    return homvector;
}


auto compute_poly_path_homvec(SmallGraph graph, SmallGraph::type num_vertices) -> std::vector<long long> {
    auto A = graph.to_adjacency_matrix();

    using double_matrix_t = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;

    double_matrix_t Ai = A.cast<double>();
    assert (Ai.rows() == Ai.cols());
    for (int i = 0; i < Ai.rows(); ++i) {
        Ai(i, i) = std::exp(static_cast<double>(1.0));
    }

    std::vector<long long> homvector;
    homvector.reserve(num_vertices+1);
    
    double_matrix_t product = double_matrix_t::Identity(Ai.rows(), Ai.cols());
    homvector.push_back(
        static_cast<long long>(std::round(product.sum()))
    );
    for (int n = 0; n < num_vertices; ++n) {
        product = product * Ai;
        homvector.push_back(
            static_cast<long long>(std::round(product.sum()*100000))
        );
    }

    return homvector;
}


template<typename Int = int64_t>
auto compute_complex_pathwidth_one_homvec(SmallGraph graph, SmallGraph::type num_vertices, bool x_with_A=false, bool boost = false) -> std::vector<std::complex<Int>> {
    auto A = graph.to_adjacency_matrix();
    auto D = graph.to_degree_matrix();

    using complex_matrix_t = Eigen::Matrix<std::complex<Int>, Eigen::Dynamic, Eigen::Dynamic>;

    // A to complex<uint64_t>
    complex_matrix_t Ai = A.cast<std::complex<Int>>();
    assert (Ai.rows() == Ai.cols());
    for (int i = 0; i < Ai.rows(); ++i) {
        for (int j = 0; j < Ai.cols(); ++j) {
            if (boost) {
                if (A(i, j) != 0) {
                    Ai(i, j) = std::complex<Int>(Int(2)*D(i, j), 0);
                } else if (i == j) {
                    Ai(i, j) = std::complex<Int>(-1, D(i, j));
                } else {
                    Ai(i, j) = std::complex<Int>(-1, -1);
                }
            } else {
                if (x_with_A) {
                    if (A(i, j) != 0) {
                        Ai(i, j) = std::complex<Int>(0, 1);
                    } else if (i == j) {
                        Ai(i, j) = std::complex<Int>(D(i, i), 0);
                    }
                } else {
                    if (A(i, j) != 0) {
                        Ai(i, j) = std::complex<Int>(1, 0);
                    } else if (i == j) {
                        Ai(i, j) = std::complex<Int>(0, D(i, i));
                    }
                }
            }
        }
    }

    std::vector<std::complex<Int>> homvector;
    homvector.reserve(num_vertices+1);
    
    complex_matrix_t product = complex_matrix_t::Identity(Ai.rows(), Ai.cols());
    homvector.push_back(product.sum());
    for (int n = 0; n < num_vertices; ++n) {
        product = product * Ai;
        homvector.push_back(product.sum());
    }

    return homvector;
}

template<typename Int = int64_t>
auto compute_quatern_pathwidth_one_homvec(SmallGraph graph, SmallGraph::type num_vertices) {
    auto A = graph.to_adjacency_matrix();
    auto D = graph.to_degree_matrix();

    using H = boost::math::quaternion<Int>;

    using complex_matrix_t = Eigen::Matrix<H, Eigen::Dynamic, Eigen::Dynamic>;

    // A to complex<uint64_t>
    complex_matrix_t Ai = A.cast<H>();
    assert (Ai.rows() == Ai.cols());
    for (int i = 0; i < Ai.rows(); ++i) {
        for (int j = 0; j < Ai.cols(); ++j) {

            if (A(i, j) != 0) {
                Ai(i, j) = H(1, 0, 0, 0);
            } else if (i == j) {
                Ai(i, j) = H(0, D(i, i), 0, 0);
            } else {
                Ai(i, j) = H(0, 0, 1, 0);
            }
        }
    }

    std::vector<H> homvector;
    homvector.reserve(num_vertices+1);
    
    complex_matrix_t product = complex_matrix_t::Identity(Ai.rows(), Ai.cols());
    homvector.push_back(product.sum());
    for (int n = 0; n < num_vertices; ++n) {
        product = product * Ai;
        homvector.push_back(product.sum());
    }

    return homvector;
}


} // namespace wl
