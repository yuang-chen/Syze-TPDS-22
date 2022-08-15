#pragma once

#include "graph.hpp"
//#include "env.hpp"

void Graph::partition(bool dynamic) {
    //////////////////////////////////////////
    dynamicSplit(dynamic);
    //////////////////////////////////////////
    setBlock();
}
///////////////
// dynamic caching of hot vertices;
// the time cost can be ignored
///////////////
void Graph::dynamicSplit(bool dynamic) {
    std::cout << "--------dynamic split--------" << '\n';
    #ifdef INFO
    boost::timer::cpu_timer timer;
    #endif

    volume = params::subgraph_size * 1024 / sizeof(VerxId);
    num_subs = (num_vertices - 1)/ volume + 1; 
    offset = (unsigned) log2 (volume); 
    grid_row = num_hot ==0? num_subs: (num_hot - 1)/volume + 1;
    grid_col = num_cold==0? num_subs: (num_hot + num_cold - 1)/volume + 1; 
    row_verx = num_hot ==0? num_vertices: num_hot;
    col_verx = num_cold==0? num_vertices: num_hot + num_cold;
    #ifdef INFO
    std::cout << "verx bounds at Row Verx * Col Verx: " << row_verx << " X " << col_verx
                                                        <<  " (" << (float) row_verx/num_vertices 
                                                        <<  " X "<< (float) col_verx/num_vertices << ")\n";
    std::cout << "static grid: Row X Col: " << grid_row << " X " << grid_col << " | " << num_subs << "\n";
    #endif
   // auto& degree = out_degree; //params::dtype == DType::in? in_degree:
    auto averdeg = csc_offset[row_verx] / grid_row; //params::dtype == DType::in? csc_offset[row_verx+1] / grid_row:
    auto stc_grid_row = grid_row; // static grid_col
    
    std::vector<unsigned> degsum(stc_grid_row);
    std::vector<unsigned> div(stc_grid_row, 1);
    std::vector<unsigned> dvolume(stc_grid_row, volume);

    /////////////////////////////
    // perform the dynamic partitioning based on the out-degrees 
    // when a subgraph is overflowed, i.e., degsum > 2 * averdeg,
    // such subgraph is further subdivided until the degsumre <= 2 * averdeg
    /////////////////////////////
    #pragma omp parallel for reduction(+:grid_row)
    for(unsigned n = 0; n < stc_grid_row; n++) {
        auto start = n * volume;
        auto end = (n + 1) * volume > row_verx? 
                    row_verx: (n + 1) * volume; 
        degsum[n] = csc_offset[end] - csc_offset[start];
        if(degsum[n] >= 2 * averdeg && dynamic) {
            div[n] = pow(2, round(log2(degsum[n]/averdeg)));
            dvolume[n] = dvolume[n]/ div[n] ;
            grid_row += div[n] -1;
        } 
    }
    std::vector<VerxId> size(grid_row, volume);
    row_block.resize(grid_row);

   // dyn_offset = std::vector<unsigned>(stc_grid_row, offset);
   // dyn_map = std::vector<unsigned>(stc_grid_row);
    dyn_offset.resize(stc_grid_row);
    dyn_map.resize(stc_grid_row);

   // dyn_map = new unsigned [stc_grid_row];
    #pragma omp parallel for 
        for(VerxId n= 0; n < stc_grid_row; n++) {
        unsigned n_map = std::accumulate(div.begin(), div.begin() + n, 0);
        dyn_map[n] = n_map;
        dyn_offset[n] = offset - log2(div[n]);
        for(VerxId i = 0; i < div[n]; i++) {
            size[i+n_map] = size[i+n_map] / div[n];
            row_block[i+n_map].start = n *volume + i*size[n_map];
            row_block[i+n_map].stop = std::min(n*volume + (i+1)*size[n_map], col_verx);
        }
    }
    /////////////////////////////////
    // row-wise dynamic spliting delivers better result than col-wise one
    /////////////////////////////////
    col_block.resize(grid_col);
    #pragma omp parallel for 
        for(unsigned n = 0; n < grid_col; n++) {
        col_block[n].id = n;
        col_block[n].start = n * volume;
        col_block[n].stop = (n + 1) * volume > col_verx? 
                             col_verx: (n + 1) * volume; 
    }
    #ifdef INFO
    std::cout << "dynamic grid: Row X Col: " << grid_row << " X " << grid_col << '\n';
    #endif
}


////////////////////////////////////////
// construct the blocks
// encode the destination vertices
////////////////////////////////////////

void Graph::setBlock() {
    std::cout << "--------block building--------" << '\n';
    #ifdef INFO
    boost::timer::cpu_timer timer;
    #endif
    dstver = alloc2AoV<VerxId>(grid_col, grid_row, volume);
    inters = alloc2AoV<VerxId>(grid_col, grid_row, volume);
    buffer = alloc2AoV<Attr_t>(grid_col, grid_row, volume);


    // * dstver: outgoing inter-out_edges
    // * encode inter-out_edges: curr_ver |= MASK::MAX_NEG;
        /////////////////////////////////////
    #pragma omp parallel for schedule(dynamic,1) //num_threads(params::threads)
    for(int n = 0; n < grid_col; n++) { 
        VerxId prev_sub = 0;
        std::vector<unsigned> count_dstver(grid_row, 0);
        unsigned vertex_count = 0;
        /////////////////////////////////////////
        // compute the sub this->csr_index
        ///////////////////////////////////////        
        std::vector<unsigned> count_sub(grid_row,0);
        for(auto i = col_block[n].start; i < col_block[n].stop; i++) {
            prev_sub = grid_row;
            for(auto j = csr_offset[i]; j < csr_offset[i+1]; j++) {
                if(csr_index[j]>=row_verx)
                    continue;
                auto dst_ver = csr_index[j];
                auto curr_sub = locateSub(csr_index[j]); //csr_index[j] >> offset;
                if(prev_sub != curr_sub) {                  // subgraph[n].inters[s][x] = i
                    inters[n][curr_sub].push_back(i); //
                    dst_ver |= MASK::MAX_NEG; 
                    prev_sub = curr_sub;
                }
                dstver[n][curr_sub].push_back(dst_ver);
            }
        }
    }
    #ifdef INFO
    std::cout << "+" << timer.elapsed().wall/(1e9) << "s: dst vertices is encoded for edge compression"  << '\n';
    #endif

    #pragma omp parallel for
    for(unsigned i = 0; i < grid_col; i++) {
        for(unsigned j = 0; j < grid_row; j++) {
            buffer[i][j].resize(inters[i][j].size());
        }
    }
}