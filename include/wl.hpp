// def __relabel(col, rehash):
//     colstr = str(col)
//     if colstr not in rehash:
//         rehash[colstr] = len(rehash)
//     return rehash[colstr]


// def __is_refined(vec1, vec2):
//     # return (vec1-vec1.min() != vec2-vec2.min()).any()
//     # return (vec1 != vec2).any()

//     rho1 = np.array(sorted([(vec1 == v).sum() for v in np.unique(vec1)]))
//     rho2 = np.array(sorted([(vec2 == v).sum() for v in np.unique(vec2)]))
//     return (rho1.shape[0] != rho2.shape[0]) or (rho1 != rho2).any()



// def colors_1(graph, rehash={}, labels=None):
//     """
//     latex_name: $1$-WL
//     """
//     A = __check_and_get_A(graph)
//     n = A.shape[0]
//     colvec = __check_and_get_labels(A, labels)
//     while True:
//         colvec2 = np.array([
//             __relabel((colvec[i], sorted(colvec[A[i] != 0])), rehash)
//             for i in range(n)
//         ])
//         if not __is_refined(colvec, colvec2):
//             colvec2.sort()
//             return colvec2, rehash
//         colvec = colvec2
//     return None



#pragma once

#include <sstream>
#include <string>
#include <unordered_map>

namespace wl {

using rehash_t = std::unordered_map<std::string, unsigned long long>;
using colvec_t = std::vector<unsigned long long>;
// auto format = [](auto&&... elems) { return ((std::ostringstream{} << ... << elems)).str(); };

template <typename T>
auto stringyfy_vector(std::vector<T> const& vec) -> std::string {
    std::ostringstream ss;
    for (auto &&value : vec) {
        ss << value << ",";
    }
    return std::move(ss).str();
}

auto relabel(std::string const& col, rehash_t& rehash) -> int {
    if (rehash.find(col) == rehash.end()) {
        rehash[col] = rehash.size();
    }
    return rehash[col];
}

auto relabel(colvec_t const& col, rehash_t& rehash) -> int {
    return relabel(stringyfy_vector(col), rehash);
}

auto is_refined(colvec_t const& vec1, colvec_t const& vec2) -> bool {
    std::unordered_map<unsigned long long, int> rho1, rho2;
    colvec_t sizes1, sizes2;
    
    for (auto &&value : vec1) {
        rho1[value] += 1;
    }

    for (auto &&value : vec2) {
        rho2[value] += 1;
    }

    for (auto &&[key, value] : rho1) {
        sizes1.push_back(value);
    }

    for (auto &&[key, value] : rho2) {
        sizes2.push_back(value);
    }

    if (sizes1.size() != sizes2.size()) {
        return true;
    }

    std::sort(sizes1.begin(), sizes1.end());
    std::sort(sizes2.begin(), sizes2.end());

    if (sizes1 != sizes2) {
        return true;
    }

    return false;
}

auto is_refined(std::vector<colvec_t> const& vec2d1, std::vector<colvec_t> const& vec2d2) -> bool {
    colvec_t vec1, vec2;
    auto rehash = rehash_t{};

    vec1.reserve(vec2d1.size());
    for (auto &&vec : vec2d1) {
        vec1.emplace_back(relabel(stringyfy_vector(vec), rehash));
    }
    vec2.reserve(vec2d2.size());
    for (auto &&vec : vec2d2) {
        vec2.emplace_back(relabel(stringyfy_vector(vec), rehash));
    }

    return is_refined(vec1, vec2);
}




auto colors_1(SmallGraph const& graph, rehash_t& rehash, colvec_t const& labels = {}) -> colvec_t {
    auto A = graph.to_adjacency_matrix();
    auto n = A.rows();
    
    auto colvec = (labels.empty()) ? colvec_t(n, -1) : labels;

    for (int t = 0;; ++t) {
        auto colvec2 = colvec_t{};
        colvec2.reserve(n);

        // std::cout << "t=" << t << ": ";
        for (int i = 0; i < n; ++i) {
            std::vector<std::string> col;
            col.reserve(n);
            for (int j = 0; j < n; ++j) {
                std::ostringstream ss;
                ss << (i==j) << A(i, j) << "&";
                ss << colvec[j] << "&" << colvec[i];
                col.emplace_back(std::move(ss).str());
            }
            std::sort(col.begin(), col.end());
            colvec2.push_back(relabel(stringyfy_vector(col), rehash));
        }
        // std::cout << "" << stringyfy_vector(colvec2) << " \n";

        if (!is_refined(colvec, colvec2)) {
            std::sort(colvec2.begin(), colvec2.end());
            // std::cout << "t=" << t;
            return colvec2;
        }
        colvec = std::move(colvec2);
    }

    return {};
}



auto colors_p1_03(SmallGraph const& graph, rehash_t& rehash, colvec_t const& labels = {}) -> colvec_t {
    auto A = graph.to_adjacency_matrix();
    auto n = A.rows();

    auto labels_ = (labels.empty() ? colvec_t(n, -1) : labels);
    auto colvec = labels_;

    while (true) {
        auto colvec2 = colvec_t{};
        colvec2.reserve(n);
        for (int i = 0; i < n; ++i) {
            std::vector<std::string> col;
            col.reserve(n);
            for (int j = 0; j < n; ++j) {
                std::ostringstream ss;
                ss << (i==j) << A(i, j) << " & ";
                ss << colvec[j];
                col.emplace_back(std::move(ss).str());
            }
            std::sort(col.begin(), col.end());
            colvec2.push_back(relabel(stringyfy_vector(col), rehash));
        }

        if (!is_refined(colvec, colvec2)) {
            std::sort(colvec2.begin(), colvec2.end());
            return colvec2;
        }

        colvec = std::move(colvec2);
    }
}


auto colors_p1_04(SmallGraph const& graph, rehash_t& rehash, colvec_t const& labels = {}) -> colvec_t {
    auto A = graph.to_adjacency_matrix();
    auto n = A.rows();

    auto labels_ = (labels.empty() ? colvec_t(n, -1) : labels);


    // every vertex has set of formulas that models 
    auto colvec2d = std::vector<colvec_t>{};
    colvec2d.reserve(n);
    for (auto label : labels_) {
        colvec2d.emplace_back(colvec_t{label});
    }

    for (int t = 0;; ++t) {
        auto colvec2d2 = std::vector<colvec_t>{};
        colvec2d2.reserve(n);
        
        // std::cout << "t=" << t << ": ";

        for (int i = 0; i < n; ++i) {
            auto sort_counter = std::unordered_map<std::string, unsigned long long>{};
            for (int j = 0; j < n; ++j) {
                std::ostringstream atp;
                atp << (i==j) << A(i, j) << -labels_[i] << -labels_[j] << " &";
                assert (colvec2d.size() > j);
                for (auto &&elem : colvec2d[j]) {
                    auto atp_color = atp.str() + std::to_string(elem);
                    sort_counter[atp_color] += 1;
                }
            }

            auto col = colvec_t{};
            col.reserve(n);
            for (auto &&[key, value] : sort_counter) {
                auto grouped_key = key + "X" + std::to_string(value);
                col.emplace_back(relabel(grouped_key, rehash));
            }

            std::sort(col.begin(), col.end());
            // std::cout << "(" << stringyfy_vector(col) << ") ";
            colvec2d2.emplace_back(std::move(col));
        }

        // std::cout << "\n";

        if (!is_refined(colvec2d, colvec2d2)) {
            auto colvec2 = colvec_t{};
            colvec2.reserve(n);
            for (auto &&vec : colvec2d2) {
                colvec2.emplace_back(relabel(stringyfy_vector(vec), rehash));
            }
            std::sort(colvec2.begin(), colvec2.end());
            // std::cout << "t=" << t;
            return colvec2;
        }

        colvec2d = std::move(colvec2d2);
    }
}


auto colors_p1_05(SmallGraph const& graph, rehash_t& rehash, colvec_t const& labels = {}) -> colvec_t {
    auto A = graph.to_adjacency_matrix();
    auto n = A.rows();

    auto labels_ = (labels.empty() ? colvec_t(n, -1) : labels);

    // every vertex has set of formulas that models 
    auto colvec2d = std::vector<colvec_t>{};
    colvec2d.reserve(n);
    for (auto label : colvec_t(n, -1)) {
        colvec2d.emplace_back(colvec_t{label});
    }

    for (int t = 0;; ++t) {
        auto colvec2d2 = std::vector<colvec_t>{};
        colvec2d2.reserve(n);
        
        // std::cout << "t=" << t << ": ";

        


        if (!is_refined(colvec2d, colvec2d2)) {
            auto colvec2 = colvec_t{};
            colvec2.reserve(n);
            for (auto &&vec : colvec2d2) {
                colvec2.emplace_back(relabel(stringyfy_vector(vec), rehash));
            }
            std::sort(colvec2.begin(), colvec2.end());
            // std::cout << "t=" << t;
            return colvec2;
        }

        colvec2d = std::move(colvec2d2);
    }
}


} // namespace wl
