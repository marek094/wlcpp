
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

#include <Eigen/Dense>

#include <vector>
#include <utility>
#include <complex>
#include <map>
#include <string>
#include <tuple>
#include <array>

namespace wl {


namespace crt {
    using ull = unsigned long long;

    // prime source:) https://www.math.utah.edu/~pa/MDS/primes.html

    constexpr auto primes = std::array{
        2147483647ull,	
        2147483629ull,	
        2147483587ull,	
        2147483579ull, 
        2147483563ull,	
        2147483549ull,	
        2147483543ull, 
        2147483497ull, 
        2147483489ull,	
        2147483477ull,	
        2147483423ull, 
        2147483399ull,
        2147483353ull,
        2147483323ull, 
        2147483269ull, 
        2147483249ull,
        2147483237ull,
        2147483179ull,
        2147483171ull,
        2147483137ull,
    };

    template<ull R>
    struct crtu64t {
        static_assert(R < primes.size(), "R must be less than primes.size()");

        crtu64t() {
            for (auto& mod : mods) mod = 0ull;
        }

        crtu64t(ull val) {
            for (int r = 0; r < R; ++r) {
                mods[r] = val % primes[r];
            }
        }

        static auto min() -> crtu64t {
            return crtu64t(0ull);
        }

        static auto max() -> crtu64t {
            crtu64t res;
            for (int r = 0; r < R; ++r) {
                res.mods[r] = primes[r] - 1;
            }
            return res;
        }

        auto operator+=(crtu64t const& other) -> crtu64t& {
            for (int r = 0; r < R; ++r) {
                mods[r] = (mods[r] + other.mods[r]) % primes[r];
            }
            return *this;
        }

        friend auto operator+(crtu64t lhs, crtu64t const& rhs) -> crtu64t {
            lhs += rhs;
            return lhs;
        }

        auto operator*=(crtu64t const& other) -> crtu64t& {
            for (int r = 0; r < R; ++r) {
                mods[r] = (mods[r] * other.mods[r]) % primes[r];
            }
            return *this;
        }

        friend auto operator*(crtu64t lhs, crtu64t const& rhs) -> crtu64t {
            lhs *= rhs;
            return lhs;
        }

        auto operator-=(crtu64t const& other) -> crtu64t& {
            for (int r = 0; r < R; ++r) {
                mods[r] = (mods[r] + primes[r] - other.mods[r]) % primes[r];
            }
            return *this;
        }

        friend auto operator-(crtu64t lhs, crtu64t const& rhs) -> crtu64t {
            lhs -= rhs;
            return lhs;
        }

        auto operator-() const -> crtu64t {
            return crtu64t(0ULL) - *this;
        }

        auto operator<=>(crtu64t const& other) const -> std::strong_ordering {
            for (int r = 0; r < R; ++r) {
                if (mods[r] != other.mods[r]) {
                    return mods[r] <=> other.mods[r];
                }
            }
            return std::strong_ordering::equal;
        }

        auto operator==(crtu64t const& other) const -> bool {
            return (*this <=> other) == std::strong_ordering::equal;
        }

        friend auto to_string(crtu64t const& val) -> std::string {
            std::string res;
            for (int r = 0; r < R; ++r) {
                res += std::to_string(val.mods[r]);
                if (r != R-1) res += "x_";
            }
            return res;
        }

        friend auto operator<<(std::ostream& os, crtu64t const& val) -> std::ostream& {
            for (int r = 0; r < R; ++r) {
                os << val.mods[r];
                if (r != R-1) os << "x_";
            }
            return os;
        }

        friend auto log10(crtu64t const& val) -> double {
            double res = 0.0;
            for (int r = 0; r < R; ++r) {
                res += std::log10(static_cast<double>(val.mods[r]));
            }
            return res;
        }


        std::array<ull, R> mods;
    };


}

} // namespace wl


namespace Eigen {
    template<std::size_t R>
    struct NumTraits<wl::crt::crtu64t<R>> : GenericNumTraits<wl::crt::crtu64t<R>> {
        typedef wl::crt::crtu64t<R> Real;

        enum {
            IsComplex = 0,
            IsInteger = 1,
            IsSigned = 0,
            RequireInitialization = 1,
            ReadCost = 10 * R,
            AddCost = 10 * R,
            MulCost = 20 * R
        };

        // // Since it's an integral type, epsilon and precision are not applicable
        static inline Real epsilon() {
            return Real(0);
        }

        static inline Real dummy_precision() {
            return Real(0);
        }

        // For an integer type, digits10 can represent the maximum representable value
        static inline int digits10() {
            double log10_max = log10(Real::max());
            return static_cast<int>(std::ceil(log10_max));
        }

        static inline Real highest() {
            return Real::max();
        }

        static inline Real lowest() {
            return Real::min();
        }
    };
}




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






template<typename T = uint64_t>
auto compute_path_homvec(SmallGraph graph, SmallGraph::type num_vertices) -> std::vector<T> {
    auto A_was = graph.to_adjacency_matrix();

    using eigen_matrix_t = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    eigen_matrix_t A = eigen_matrix_t::Zero(A_was.rows(), A_was.cols());
    for (int i = 0; i < A_was.rows(); ++i) {
        for (int j = 0; j < A_was.cols(); ++j) {
            A(i, j) = A_was(i, j);
        }
    }

    std::vector<T> homvector;
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



template<typename T = int64_t>
auto compute_complex_path_homvec(SmallGraph graph, SmallGraph::type num_vertices, long long sign = 1) -> std::vector<std::complex<T>> {
    auto A = graph.to_adjacency_matrix();

    using complex_matrix_t = Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic>;

    // A to complex<uint64_t>
    complex_matrix_t Ai = A.cast<std::complex<T>>();
    assert (Ai.rows() == Ai.cols());
    for (int i = 0; i < Ai.rows(); ++i) {
        Ai(i, i) = std::complex<T>(0, sign*1);
    }

    std::vector<std::complex<T>> homvector;
    homvector.reserve(num_vertices+1);
    
    complex_matrix_t product = complex_matrix_t::Identity(Ai.rows(), Ai.cols());
    homvector.push_back(product.sum());
    for (int n = 0; n < num_vertices; ++n) {
        product = product * Ai;
        homvector.push_back(product.sum());
    }

    return homvector;
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


template<typename T = int64_t>
auto compute_complex_pathwidth_one_homvec(SmallGraph graph, SmallGraph::type num_vertices, bool x_with_A=false) -> std::vector<std::complex<T>> {
    auto A = graph.to_adjacency_matrix();
    auto D = graph.to_degree_matrix();

    using complex_matrix_t = Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic>;

    // A to complex<uint64_t>
    complex_matrix_t Ai = A.cast<std::complex<T>>();
    assert (Ai.rows() == Ai.cols());
    for (int i = 0; i < Ai.rows(); ++i) {
        for (int j = 0; j < Ai.cols(); ++j) {

            if (x_with_A) {
                if (A(i, j) != 0) {
                    Ai(i, j) = std::complex<T>(0, 1);
                } else if (i == j) {
                    Ai(i, j) = std::complex<T>(D(i, i), 0);
                }
            } else {
                if (A(i, j) != 0) {
                    Ai(i, j) = std::complex<T>(1, 0);
                } else if (i == j) {
                    Ai(i, j) = std::complex<T>(0, D(i, i));
                }
            }
        }
    }

    std::vector<std::complex<T>> homvector;
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
