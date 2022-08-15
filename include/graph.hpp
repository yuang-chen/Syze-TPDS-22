#pragma once

#include "global.hpp"
#include "utils.hpp"
//#include "env.hpp"

class Graph{
public:
    /////////////////////
    // basic graph components
    /////////////////////
    VerxId num_vertices = 0;
    EdgeId num_edges = 0;

    unsigned num_subs = 0;
    unsigned volume = 0;
    unsigned offset = 0;
    ////////////////
    // dynamic split
    ////////////////
 //   unsigned* dyn_offset = nullptr;
  //  unsigned* dyn_map = nullptr;
    std::vector<unsigned> dyn_offset;
    std::vector<unsigned> dyn_map;
    ///////////////////
    // hot, cold vertices
    //////////////////
    unsigned num_hot = 0;
    unsigned num_cold = 0;
    // a grid is established according to the hot and cold vertices.
    unsigned grid_row = 0;   // grid_row the number of row-wise blocks
    unsigned grid_col = 0;  // grid_col the number of column-wise blocks
    unsigned row_verx = 10000; // row_v the boundary vertex for defining the grid_row
    unsigned col_verx = 0;   //  col_v border vertex defining the grid_col
    std::vector<SubMat> row_block;  // [new num subs] new number = params::num_subgraphs, after dynamic partitioning
    std::vector<SubMat> col_block; 
    
    VerxId* new_id = nullptr; //std::unique_ptr<VerxId[]> new_id; smart pointers are good, but slower than C-style ptrs
    EdgeId* csr_offset = nullptr; // vertex index 
    VerxId* csr_index = nullptr;  // edge index
    VerxId* out_degree = nullptr;

    EdgeId* csc_offset = nullptr; // vertex index 
    VerxId* csc_index = nullptr;   // edge index
    VerxId* in_degree = nullptr;
   
    /////////////////////
    // components for subgraphs
    ////////////////////
    AoV<Attr_t> buffer = nullptr;
    AoV<VerxId> inters = nullptr;
    AoV<VerxId> dstver = nullptr;

    Graph(VerxId _vertices = 0, EdgeId _out_edges = 0):
            num_vertices(_vertices),
            num_edges(_out_edges) {
    }

    ~Graph() {
      delete[] csr_offset; 
      delete[] csr_index;
      delete[] out_degree;

      delete[] csc_offset;
      delete[] csc_index;
      delete[] in_degree;

      delete2AoV(dstver, grid_col, grid_row);
      delete2AoV(inters, grid_col, grid_row);
      delete2AoV(buffer, grid_col, grid_row);

    }
    //////////////////
    // load
    // functions for load the dataset and constructing the basic graph data structure
    ///////////////
    template <typename EdgeData = Empty>
    void load(std::string filename); 
    template <typename EdgeData = Empty>
    bool loadBEL(std::string filename);
    template <typename EdgeData = Empty>
    bool loadCSR(std::string filename);
    template <typename EdgeData = Empty>
    bool loadMIX(std::string filename);
    template <typename EdgeData = Empty>
    bool storeMIX(std::string filename);  
    /////////////////
    // reordering methods
    /////////////////
    void reorder(RAlgo ralgo, DType dtype);
    void BiHub(DType dtype);
    void PoP(DType dtype);
    void TriHub(DType dtype);
    void FBC(DType dtype);
    void DBG(DType dtype);
    void buildNewGraph();
    //////////////
    // subdivide
    // functions for partitioning and building the subgraphs
    //////////////
    void partition(bool dynamic);
    void dynamicSplit(bool dynamic);
    void setBlock();
    //////////////////////////
    // propagate.tpp
    // function for propagating the data
    // is it good to put every thing inside a class?
    //////////////////////////
    //////////////////////////
    template<class Graph_Func> void run(Graph_Func& func);
    template<class Graph_Func> void run(Graph_Func& func, unsigned rounds, unsigned iterations);
    //! block-based GAS
    template<class Graph_Func> void computeBlock(Graph_Func& func);
    template<class Graph_Func> void scatterBlock(Graph_Func& func, std::vector<Attr_t>& buf_bin, std::vector<VerxId>& inter_bin);
    template<class Graph_Func> void gatherBlock(Graph_Func& func, std::vector<Attr_t>& buf_bin, std::vector<VerxId>& dst_bin);

    /////////
    // helper function
    ////////

    unsigned inline locateSub(unsigned vertex) {
      //! if static caching
        auto const orig_sub = vertex >> offset;
        auto const curr_off = dyn_offset[orig_sub];
        auto const div_sub = (((1 << (offset-curr_off))) - 1) & (vertex >> curr_off);
        return dyn_map[orig_sub] + div_sub; 
      //! if dynamic caching
      // return vertex >> offset;    
    }

    inline int at(unsigned row, unsigned col) {
        return row * num_subs + col;
    }

};
///////////////////
// TODO: this thing is ugly
// maybe we should define new classes for these components
//////////////////
#include "graph-load.hpp"
#include "graph-subdivide.hpp"
#include "graph-propagate.tpp"    
