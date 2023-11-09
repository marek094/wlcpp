#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <bitset>
#include <utility>

#include "unlabelled_graph.hpp"

namespace wl {

// General principles:

//   All numbers in this description are in decimal unless obviously 
//   in binary.

//   Apart from the header, there is one object per line. Apart from
//   the header, end-of-line characters, and the characters ":", ";"
//   and "&" which might start a line, all bytes have a value in the
//   range 63-126 (which are all printable ASCII characters). A file of
//   objects is a text file, so whatever end-of-line convention is
//   locally used is fine; however the C library input routines must
//   show the standard single-LF end of line to programs).

// Bit vectors:

//   A bit vector x of length k can be represented as follows.  
//       Example:  1000101100011100

//   (1) Pad on the right with 0 to make the length a multiple of 6.
//       Example:  100010110001110000

//   (2) Split into groups of 6 bits each.
//       Example:  100010 110001 110000

//   (3) Add 63 to each group, considering them as bigendian binary numbers.
//       Example:  97 112 111

//   These values are then stored one per byte.  
//   So, the number of bytes is ceiling(k/6).

//   Let R(x) denote this representation of x as a string of bytes.
      
// Small nonnegative integers:
 
//   Let n be an integer in the range 0-68719476735 (2^36-1).

//   If 0 <= n <= 62, define N(n) to be the single byte n+63.
//   If 63 <= n <= 258047, define N(n) to be the four bytes
//       126 R(x), where x is the bigendian 18-bit binary form of n.
//   If 258048 <= n <= 68719476735, define N(n) to be the eight bytes
//       126 126 R(x), where x is the bigendian 36-bit binary form of n.

//   Examples:  N(30) = 93
//              N(12345) = N(000011 000000 111001) = 126 66 63 120
//              N(460175067) = N(000000 011011 011011 011011 011011 011011)
//                           = 126 126 63 90 90 90 90 90


// Description of graph6 format.
// ----------------------------

// Data type:  
//    simple undirected graphs of order 0 to 68719476735.

// Optional Header: 
//    >>graph6<<     (without end of line!)

// File name extension:
//    .g6

// One graph:
//    Suppose G has n vertices.  Write the upper triangle of the adjacency
//    matrix of G as a bit vector x of length n(n-1)/2, using the ordering
//    (0,1),(0,2),(1,2),(0,3),(1,3),(2,3),...,(n-2,n-1).

//    Then the graph is represented as  N(n) R(x).

// Example:
//    Suppose n=5 and G has edges 0-2, 0-4, 1-3 and 3-4.
//        01, 02,12  03,13,23  04,14,24  34
//    x = 0   10     010       1001
    
//    Then N(n) = 68 and R(x) = R(010010 100100) = 81 99.
//    So, the graph is  68 81 99.



// Decoding a byte into a 6-bit binary
template<typename T>
auto decode_byte(T ch) -> int {
    return static_cast<int>(ch) - 63;
}

// Decode N(n)
auto decode_n(std::string const&input, size_t &index) -> unsigned long long {
    if (input[index] < 126) {
        return decode_byte(input[index++]);
    }

    assert(false); // this big is not supported yet
    
    int count = 0;
    while (input[index] == 126) {
        count += 3;
        index += 1;
    }
    std::cout << "count: " << count << '\n';
    unsigned long long value = 0;
    for (int i = 0; i < count; ++i) {
        // display 8 bites:
        std::cout << std::bitset<6>(decode_byte(input[index])) << '\t' << decode_byte(input[index]) << '\n';

        value = (value << 6) | decode_byte(input[index++]);
    }

    return value;
}

auto graph6_to_graph(std::string const&input) -> SmallGraph {
    // size of the graph <= 62
    SmallGraph graph;

    size_t index = 0;
    graph.num_vertices = decode_n(input, index);
    graph.adj_list.resize(graph.num_vertices);
    auto n = graph.num_vertices;

    //    Suppose G has n vertices.  Write the upper triangle of the adjacency
    //    matrix of G as a bit vector x of length n(n-1)/2, using the ordering
    //    (0,1),(0,2),(1,2),(0,3),(1,3),(2,3),...,(n-2,n-1).
    // Example:
    //    Suppose n=5 and G has edges 0-2, 0-4, 1-3 and 3-4.
    //    x = 0 10 010 1001
    
    // implementation
    auto pointer = 1ULL << 5;
    int col = 1;
    int row = 0;
    while(index < input.size()) {
        if (pointer == 0) {
            pointer = 1ULL << 5;
            index++;
        }

        // std::cout << "XX " << std::bitset<6>(decode_byte(input[index])) << "; " << std::bitset<6>(pointer) << "; ";

        if (decode_byte(input[index]) & pointer) {
            graph.add_edge(row, col);
            // std:: cout << 1 << '\t' << row << " " << col <<'\n';
        } else {
            // std:: cout << 0 << '\t' << row << " " << col <<'\n';
        }

        pointer >>= 1;
        
        row += 1;
        if (row == col) {
            col += 1;
            row = 0;
        }
    }

    return graph;
}

auto read_graph_from_graph6_file(std::string const& filename, size_t limit = (~0), size_t skip = 0) -> std::vector<SmallGraph> {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open the file: " << filename << '\n';
        return {};
    }

    std::vector<SmallGraph> graphs;
    std::string line;

    int skipper_i = 0;
    while (std::getline(file, line)) {
        if (skipper_i++ < skip) {
            continue;
        } else {
            skipper_i = 0;
        }

        if (line.substr(0, 10) == ">>graph6<<") {
            graphs.push_back(graph6_to_graph(line.substr(10)));
        } else {
            graphs.push_back(graph6_to_graph(line));
        }
        if (graphs.size() >= limit) {
            break;
        }

    }

    return graphs;
}



} // namespace wl
