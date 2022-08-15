#pragma once

#include "graph.hpp"
#include <boost/algorithm/string/predicate.hpp>
#include <regex>

template <typename EdgeData>
void Graph::load(std::string filename) {
    std::cout << "loading " << filename << '\n';

    if( boost::algorithm::ends_with(filename, "bel") ) 
        assert(loadBEL(filename)==true);
    else if( boost::algorithm::ends_with(filename, "csr") )
        assert(loadCSR(filename)==true);
    else if( boost::algorithm::ends_with(filename, "mix") )
        assert(loadMIX(filename)==true);   
    else
        std::cout << "unsupported graph format" << '\n';

    
    std::cout << "num vertices: " << num_vertices << " intV: " << sizeof(VerxId) << " Byte\n";
    std::cout << "num edges: " << num_edges << " intE: " << sizeof(VerxId) << " Byte\n";
    if(!std::is_same<EdgeData, Empty>::value) std::cout << "weighted graph" << '\n';
    
    if(params::output) 
        storeMIX(filename);
}


template <typename EdgeData>
bool Graph::loadBEL(std::string filename) {
    std::ifstream input_file (filename, std::ios::binary);
    if(!input_file.is_open()) {
        std::cout << "cannot open the input bel file!" << '\n';
        return false;
    }
    input_file.read(reinterpret_cast<char*>(&num_vertices), sizeof(VerxId));
    input_file.read(reinterpret_cast<char*>(&num_edges), sizeof(EdgeId));

    size_t edge_unit = std::is_same<EdgeData, Empty>::value ? 0 : sizeof(EdgeData);
    edge_unit += 2 * sizeof(VerxId);
    #ifdef WEIGHTED
    edge_unit += sizeof(float);
    #endif
    EdgeUnit<EdgeData> * edge_buffer = new EdgeUnit<EdgeData> [num_edges];    
 //   input_file.read(reinterpret_cast<char*>(edge_buffer), num_edges * edge_unit);

    out_degree = new VerxId[num_vertices];
    in_degree = new VerxId[num_vertices];
    csr_offset = new EdgeId[num_vertices+1];
    csc_offset = new EdgeId[num_vertices+1];
    csr_index = new VerxId[num_edges];
    csc_index = new VerxId[num_edges];

    for(int i = 0; i < num_edges; i++) {
        out_degree[edge_buffer[i].src]++;
        in_degree[edge_buffer[i].dst]++;
    }

    for(int i = 0; i < num_vertices; i++) {
        csr_offset[i+1] = csr_offset[i] + out_degree[i];
        csc_offset[i+1] = csc_offset[i] + in_degree[i];
    }


    std::vector<VerxId> csr_count(num_vertices, 0);
    std::vector<VerxId> csc_count(num_vertices, 0);

    for(EdgeId i = 0; i < num_edges; i++) {
        auto dst = edge_buffer[i].dst;
        auto src = edge_buffer[i].src;
        csr_index[csr_offset[src] + csr_count[src]++]  = dst ;
        csc_index[csc_offset[dst] + csc_count[dst]++]  = src ;
    }

    // sort the edges of each vertex
    for(VerxId i = 0; i < num_vertices; i++) {
        __gnu_parallel::sort(csr_index + csr_offset[i], csr_index + csr_offset[i+1]);
        __gnu_parallel::sort(csc_index + csc_offset[i], csc_index + csc_offset[i+1]);
    }

    input_file.close();
    delete[] edge_buffer;
    return true;
}

template <typename EdgeData>
bool Graph::loadCSR(std::string filename) {
    std::ifstream input_file (filename, std::ios::binary);
    if(!input_file.is_open()) {
        std::cout << "cannot open the input csr file!" << '\n';
        return false;
    }

    input_file.read(reinterpret_cast<char*>(&num_vertices), sizeof(VerxId));
    input_file.read(reinterpret_cast<char*>(&num_edges), sizeof(EdgeId));
   
    csr_offset = new EdgeId[num_vertices+1];
    csr_index = new VerxId[num_edges];

    input_file.read(reinterpret_cast<char*>(csr_offset), num_vertices * sizeof(EdgeId));
    input_file.read(reinterpret_cast<char*>(csr_index), num_edges * sizeof(VerxId));

    csr_offset[num_vertices] = num_edges;

#ifdef WEIGHTED
    std::vector<VerxId> local_wei(num_edges);
    input_file.read(reinterpret_cast<char*>(local_wei), num_edges * sizeof(VerxId));
    edge_weight = move(local_wei);
#endif
    input_file.close();

    out_degree = new VerxId[num_vertices];
    in_degree = new VerxId[num_vertices];
    csc_offset = new EdgeId[num_vertices+1];
    csc_index = new VerxId[num_edges];

    int count = 0;
    #pragma omp parallel for
    for(VerxId i = 0; i < num_vertices; i++) {
        out_degree[i] = (VerxId) (csr_offset[i + 1] - csr_offset[i]);
        for(VerxId j = csr_offset[i]; j < csr_offset[i+1]; j++) {
            #pragma omp atomic
            in_degree[csr_index[j]]++;
        }
    }
    
    __gnu_parallel::partial_sum(in_degree, in_degree+num_vertices, csc_offset + 1);

    std::vector<VerxId> csc_count(num_vertices, 0);
    #pragma omp parallel for
    for(VerxId i = 0; i < num_vertices; i++) {
        for(VerxId j = csr_offset[i]; j < csr_offset[i+1]; j++) {
            auto dst = csr_index[j];
            csc_index[csc_offset[dst] + csc_count[dst]++] = i;
        }
    }
    #pragma omp parallel for
    for(VerxId i = 0; i < num_vertices; i++) {
    __gnu_parallel::sort(csc_index + csc_offset[i], csc_index + csc_offset[i+1]);
    }
    return true;
}

template <typename EdgeData>
bool Graph::loadMIX(std::string filename) {
    std::ifstream input_file (filename, std::ios::binary);
    if(!input_file.is_open()) {
        std::cout << "cannot open the input mix file!" << '\n';
        return false;
    }

    input_file.read(reinterpret_cast<char*>(&num_vertices), sizeof(VerxId));
    input_file.read(reinterpret_cast<char*>(&num_edges), sizeof(EdgeId));

    out_degree = new VerxId[num_vertices];
    in_degree = new VerxId[num_vertices];
    csr_offset = new EdgeId[num_vertices+1];
    csc_offset = new EdgeId[num_vertices+1];
    csr_index = new VerxId[num_edges];
    csc_index = new VerxId[num_edges];
    
    input_file.read(reinterpret_cast<char*>(csr_offset), num_vertices * sizeof(EdgeId));
    input_file.read(reinterpret_cast<char*>(csr_index), num_edges * sizeof(VerxId));
    input_file.read(reinterpret_cast<char*>(csc_offset), num_vertices * sizeof(EdgeId));
    input_file.read(reinterpret_cast<char*>(csc_index), num_edges * sizeof(VerxId));
    input_file.close();

    csr_offset[num_vertices] = num_edges;
    csc_offset[num_vertices] = num_edges;

    #pragma omp for
    for(VerxId i = 0; i < num_vertices; i++) {
        out_degree[i] = (VerxId) (csr_offset[i + 1] - csr_offset[i]);
        in_degree[i] = (VerxId) (csc_offset[i + 1] - csc_offset[i]);
    }
    return true;
}

template <typename EdgeData>
bool Graph::storeMIX(std::string filename) {
    filename = std::regex_replace(filename, std::regex("csr|bel"), "mix"); 
    std::cout << "storing " << filename << '\n';
    std::ofstream output_file (filename, std::ios::binary);
    if(!output_file.is_open()) {
        std::cout << "cannot open the output mix file!" << '\n';
        return false;
    }

    output_file.write(reinterpret_cast<char*>(&num_vertices), sizeof(EdgeId));
    output_file.write(reinterpret_cast<char*>(&num_edges), sizeof(VerxId));
 
    output_file.write(reinterpret_cast<char*>(csr_offset), num_vertices * sizeof(EdgeId));
    output_file.write(reinterpret_cast<char*>(csr_index),num_edges * sizeof(VerxId));

    output_file.write(reinterpret_cast<char*>(csc_offset), num_vertices * sizeof(EdgeId));
    output_file.write(reinterpret_cast<char*>(csc_index), num_edges * sizeof(VerxId));

    output_file.close();   
    return true;    
}
