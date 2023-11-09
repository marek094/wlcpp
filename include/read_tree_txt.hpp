#pragma once 



// There is one tree per line. The trees are given as an obvious list of edges, with vertices numbered from 0.

namespace wl {


auto read_graph_from_tree_txt_file(std::string const& filename, size_t limit = (~0), size_t skip = 0) -> std::vector<SmallGraph> {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open the file: " << filename << '\n';
        return {};
    }

    std::vector<SmallGraph> graphs;
    std::string line;

    unsigned long long n, m;
    int has_labels;

    while (getline(file, line)) {
        std::istringstream ss(line);

        auto graph = SmallGraph{};
        graph.num_vertices = 0;
        
        unsigned long long u, v; 
        while (ss >> u >> v) {
            graph.num_vertices = std::max(graph.num_vertices, u+1);
            graph.num_vertices = std::max(graph.num_vertices, v+1);
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