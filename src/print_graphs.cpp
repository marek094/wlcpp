

#include "read_g6.hpp"
#include "read_labeled_txt.hpp"
#include "read_tree_txt.hpp"
#include "unlabelled_graph.hpp"


#include <iostream>
#include <filesystem>


int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <path_to_graph6_file>\n";
        return 1;
    }


    std::filesystem::path file_path = argv[1];

    size_t limit = -1;
    size_t skip = 0;
    auto graph_list = [&]() {
        if (file_path.extension() == ".txt" && file_path.stem().string().substr(0, 4) == "tree") {    
            return wl::read_graph_from_tree_txt_file(file_path, limit, skip);
        } else
        if (file_path.extension() == ".txt") {
            return wl::read_graph_from_labeled_txt_file(file_path, limit, skip);
        }
        assert(file_path.extension() == ".g6");
        return wl::read_graph_from_graph6_file(file_path, limit, skip);
    }();

    

    for (int i = 2; i < argc; ++i) {
        
        int index = std::atoi(argv[i]);
        if (index < 0 || index >= graph_list.size()) {
            std::cerr << "Index out of bounds: " << index << "\n";
            return 1;
        }

        auto const& graph = graph_list[index];
        std::cout << "#Graph-" << index << "\t";
        graph.print_as_txt_line(std::cout);
        std::cout << "\n";
    }


    return 0;
}