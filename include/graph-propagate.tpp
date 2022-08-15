#pragma once
/////////////////////////
//! scatter and gather
////////////////////////
#include <fstream>
#include <boost/timer/timer.hpp>

template<class Graph_Func>
void Graph::run(Graph_Func& func) {

    computeBlock(func);

    computePull(func);
   
    func.swapFunc();
}

template<class Graph_Func>
void Graph::computeBlock(Graph_Func& func) {
    #pragma omp parallel for schedule(dynamic,1) num_threads(params::threads) 
    for(int i = 0; i < grid_col; i++) {    
        for(int j = 0; j < grid_row; j++) 
          scatterBlock(func, buffer[i][j], inters[i][j]);           
    }

    #pragma omp parallel for schedule(dynamic,1) num_threads(params::threads) 
    for(int i = 0; i < grid_row; i++) {
        for(auto v = row_block[i].start; v < row_block[i].stop; v++)
            func.resetFunc(v);  
        for(int j = 0; j < grid_col; j++) 
            gatherBlock(func, buffer[j][i], dstver[j][i]);
        for(int v = row_block[i].start; v < row_block[i].stop; v++) 
            func.applyFunc(v);
    }
}

/////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!! Dense mode !!!!!!!!!!!!!!!!!!!!!
///////////////////
template<class Graph_Func>
void Graph::scatterBlock(Graph_Func& func, std::vector<Attr_t>& buf_bin, std::vector<VerxId>& inter_bin) {
    unsigned index = 0;
    for(const auto src_ver: inter_bin) // auto const src_ver: inter_bin) 
        buf_bin[index++] = func.scatterFunc(src_ver);
}

template<class Graph_Func>
void Graph::gatherBlock(Graph_Func& func, std::vector<Attr_t>& buf_bin, std::vector<VerxId>& dst_bin) {
    unsigned index = MASK::MAX_UINT;
    for(const auto curr_ver: dst_bin) { 
        index += (curr_ver >> MASK::MSB);
        func.gatherFunc(curr_ver & MASK::MAX_POS, buf_bin[index]); // pagerank[vertex] += update;
    }
}

