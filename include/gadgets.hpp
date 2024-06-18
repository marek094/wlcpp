#pragma once

#include "unlabelled_graph.hpp"

#include <optional>


namespace wl {
    



auto gadgetize_small(SmallGraph const& graph, bool do_twist = false) -> std::optional<SmallGraph> {
    auto n = graph.number_of_vertices();
    
    ull nx = 0;

    auto pointervec = std::vector<ull>{};
    pointervec.reserve(n+1);
    pointervec.push_back(0);
    for (auto&& nbs : graph.adj_list) {
        ull deg = nbs.size();
        // 2^(k-1) + k + k
        nx += (1ULL<<(deg-1)) + 2*deg;

        // get max value of SmallGraph::vertex_t
        if (nx >= std::numeric_limits<SmallGraph::type>::max()-2) {
            return std::nullopt;
        }

        pointervec.push_back(nx);
    }
    // std::cout << "nx: " << nx << std::endl;
    // std::cout << "pointervec: " << pointervec.size() << std::endl;
    
    auto Xgraph = SmallGraph{};
    Xgraph.num_vertices = nx;
    Xgraph.adj_list.resize(nx);

    bool was_twisted = false;
    for (ull vtx = 0; vtx < graph.adj_list.size(); ++vtx) {
        auto vd = graph.adj_list[vtx].size();
        // std::cout << "vd: " << vd << std::endl;
        for (ull vidx = 0; vidx < graph.adj_list[vtx].size(); ++vidx) {
            auto wtx = graph.adj_list[vtx][vidx];
            if (wtx < vtx) continue;
            // std::cout << "wtx: " << (int)wtx << std::endl;
            auto wd = graph.adj_list[wtx].size();

            auto widx = 0; 
            for (auto ww : graph.adj_list[wtx]) {
                if (ww == vtx) {
                    break;
                }
                widx += 1;
            }
            // std::cout << "widx: " << widx << std::endl;

            if (do_twist && !was_twisted && vd >= 2 && wd >= 2) { // twist rule
                // a to b
                Xgraph.add_edge(pointervec[vtx]+2*vidx, pointervec[wtx]+2*widx+1);
                // b to a
                Xgraph.add_edge(pointervec[vtx]+2*vidx+1, pointervec[wtx]+2*widx);
                was_twisted = true;
            } else {
                // std::cout << "vtx: " << vtx << std::endl;
                // std::cout << "vidx: " << vidx << std::endl;
                // std::cout << "wtx: " << wtx << std::endl;
                // std::cout << "widx: " << widx << std::endl;
                // std::cout << "pointervec[vtx]+2*vidx: " << pointervec[vtx]+2*vidx << std::endl;
                // std::cout << "pointervec[wtx]+2*widx: " << pointervec[wtx]+2*widx << std::endl;
                // std::cout << "pointervec[vtx]+2*vidx+1: " << pointervec[vtx]+2*vidx+1 << std::endl;
                // std::cout << "pointervec[wtx]+2*widx+1: " << pointervec[wtx]+2*widx+1 << std::endl;

                // a to a
                Xgraph.add_edge(pointervec[vtx]+2*vidx, pointervec[wtx]+2*widx);
                // b to b
                Xgraph.add_edge(pointervec[vtx]+2*vidx+1, pointervec[wtx]+2*widx+1);
            }
        }

        auto middle_ptr = pointervec[vtx]+2*vd;
        for (ull S = 0; S < (1ULL<<vd); ++S) {
            // if sum of bits in S is odd, then a, else b
            if (__builtin_popcount(S) & 1) continue;
            for (ull part = 0; part < vd; ++part) {
                // std::cout << "\tmiddle_ptr: " << middle_ptr << std::endl;
                // std::cout << "\tpointervec[vtx]+2*part: " << pointervec[vtx]+2*part << std::endl;
                // std::cout << "\tpointervec[vtx]+2*part+1: " << pointervec[vtx]+2*part+1 << std::endl;

                if (S & (1ULL<<part)) {
                    // m_S to a
                    Xgraph.add_edge(middle_ptr, pointervec[vtx]+2*part);
                } else {
                    // m_S to b
                    Xgraph.add_edge(middle_ptr, pointervec[vtx]+2*part+1);
                }
            }
            middle_ptr += 1;
        }
        
    }

    // std::cout << "Xgraph: " << std::endl;
    // std::cout << "Xgraph.adj_list.size: " << Xgraph.adj_list.size() << std::endl;
    // std::cout << "num_vertices: " << nx << std::endl;
    // std::cout << "Xgraph.num_vertices: " << Xgraph.number_of_vertices() << std::endl;

    return Xgraph;
}


} // namespace wl
