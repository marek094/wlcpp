
#include "homcounts.hpp"
#include "unlabelled_graph.hpp"

#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>


auto get_path(size_t n) -> wl::SmallGraph {
    auto graph = wl::SmallGraph{};
    graph.num_vertices = n;
    for (size_t i = 0; i < n-1; ++i) {
        graph.add_edge(i, i+1);
    }
    return graph;
}


auto choose(int64_t n, int64_t k) -> int64_t {
    if (k < 0) {
        return 0;
    }

    if (n < 0) {
        return 0;
    }

    if (k > n) {
        return 0;
    }
    if (k == 0 || k == n) {
        return 1;
    }
    if (k > n/2) {
        return choose(n, n-k);
    }
    return n * choose(n-1, k-1) / k;
}


auto compute_formula(int64_t n, int64_t m) -> int64_t {
    // # P_m -> P_n

    int64_t sum_j = 0;

    for (int64_t j = 0; j < n; ++j) {
        int64_t l = std::max<int64_t>(0ll, (m - j - 1)/2 + (m - j - 1)%2 );
        int64_t u = std::min<int64_t>(m - 1, (m + n - j - 2)/2);

        for (int64_t i = l; i <= u; ++i) {

            int64_t const T = (m + n)/n;
            for (int64_t t = -T; t <= T; ++t) {
                sum_j += choose(m-1, i-t*(n+1)) - choose(m-1, i+j-t*(n+1)+1);
            }
        }
    }

    return sum_j;
}


int main(int argc, char const *argv[]) {
    int m = 10;
    int n = 10;

    if (argc > 1) {
        m = std::atoi(argv[1]);
    }

    if (argc > 2) {
        n = std::atoi(argv[2]);
    }


    // # P_m -> P_n
    auto path = get_path(n);

    std::cout << path.number_of_vertices() << "\n";

    auto homvec = wl::compute_path_homvec(path, m-1);


    // # P_m -> P_n
    for (int mi = 1; mi <= m; mi++) {
        std::cout 
            << "P" << std::setw(3) << std::setfill('_') << mi 
            << " -> P" << std::setw(3) << std::setfill('_') << n << ": "
            << std::setfill(' ')
            << std::setw(8) << compute_formula(n, mi) << " " 
            << std::setw(8) << homvec[mi-1] << " " << std::endl;
    }

    return 0;
}