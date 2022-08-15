
#include <pthread.h>
#include <time.h>
#include <string>
#include <boost/timer/timer.hpp>
#include "option.hpp"
#include "graph.hpp"

using namespace boost::timer;
using namespace boost::program_options;


class PR_Func{
public:
    VerxId size = 0;
    Attr_t const damp = 0.85;
    Attr_t* pr_src = nullptr; 
    Attr_t* pr_dst = nullptr;
    VerxId* out_degree = nullptr;

    ~PR_Func() {
        delete[] pr_src;
        delete[] pr_dst;
    }

    PR_Func(VerxId* _out_degree, VerxId _size): out_degree(_out_degree), size(_size) {
        pr_src = new Attr_t[size]();
        pr_dst = new Attr_t[size]();
        std::cout << "PCPM-based pagerank" << "\n";
    }
    ////////////////////////
    //user-defined functions
    ////////////////////////
    inline float scatterFunc(VerxId vertex) // this one is used in both modes
    {
        return pr_src[vertex];
    }
    inline bool resetFunc(VerxId vertex) {
        pr_dst[vertex] = 0;
        return true;
    }
    inline bool gatherFunc(VerxId vertex, Attr_t update) {
        pr_dst[vertex] += update;
        return true;
    }
    inline bool applyFunc(VerxId vertex)
    {
        pr_dst[vertex] = 1 - damp + damp * pr_dst[vertex];
        if (out_degree[vertex] > 0)
            pr_dst[vertex] = pr_dst[vertex] / out_degree[vertex];        // pagerank[vertex] = 1 - damp + damp * pagerank[vertex];
        return true;
    }
    /////////////////
    // push and pull
    /////////////////
     inline void pushFunc(VerxId dst, VerxId src) {
        atomicAdd(&pr_dst[dst], pr_src[src]);
    }   
    inline void pullFunc(VerxId dst, VerxId src) {
        pr_dst[dst] += pr_src[src];
    }
    inline void swapFunc() {
        std::swap(pr_src, pr_dst);
    }
 
    ///////////////////
    // helper functions
    ///////////////////
    void inline init(){  
        #pragma omp for 
        for (VerxId i=0; i < size; i++) {
            pr_src[i] = out_degree[i] > 0? (Attr_t) 1 / out_degree[i]: 1.0;
        }  
    }
    /////////////////////
    // verify the pagerank value
    //////////////////////
    void verify(){
        #pragma omp parallel for schedule(static, 1024)     
        for(VerxId i = 0; i < size; i++) {
            if(out_degree[i] > 0) 
                pr_src[i] = pr_src[i] * out_degree[i];
        }

        std::vector<unsigned> idx(size, 0);
        std::iota(idx.begin(), idx.end(), 0);
        std::stable_sort(idx.begin(), idx.end(),
           [this](size_t i1, size_t i2) {return pr_src[i1] > pr_src[i2];});
        std::stable_sort(pr_src, pr_src + size, std::greater<Attr_t>());
        auto filename = params::input_file.substr(params::input_file.find_last_of("/") + 1);
        std::ofstream output("./logs/sorted_"+ filename +".txt");
        if(!output.is_open()) {
            std::cout << "cannot open txt file!" << std::endl;
            exit(1);
        }
        // only write 100 values
        for(int i = 0; i < 100; i++) 
            output << std::left << std::setw(8) << idx[i] <<" "
                   << std::right << std::setw(10) << pr_src[i] << '\n';
        output.close();    
        std::cout << "the sorted pagerank values are printed out in: ./logs/sorted_"
                  << filename +".txt\n";  
    }
};

//////////////////////////////////////////
//main function
//////////////////////////////////////////
int main(int argc, char** argv)
{
    options(argc, argv); 
    Graph graph;
    // Compute the preprocessing time
    cpu_timer timer;
    //////////////////////////////////////////
    // read csr file
    //////////////////////////////////////////
    graph.load(params::input_file);

    std::cout << std::fixed << std::setprecision(3);   
    std::cout << "@" << timer.elapsed().wall/(1e9) << "s: graph is loaded from " << params::input_file << '\n';
    std::cout << "restart the global timer\n";   
    timer.start();

 //   graph.reorder(params::ralgo, params::dtype);
    //////////////////////
    //! Graph Partitioning
    //////////////////////
    graph.num_hot = 0;
    graph.num_cold = 0;
    graph.partition(params::dynamic);
    std::cout << "@" << timer.elapsed().wall/(1e9) << "s: graph is cache-blocked" << '\n';

    // //std::vector<float> pagerank(graph.num_vertices);
    PR_Func pr(graph.out_degree, graph.num_vertices);
    
    std::cout << "--------program executing--------" << '\n';
    float total_time = 0, current_time = 0;
    #ifdef INFO
    float block_time = 0, pull_time = 0, total_bt = 0, total_pt = 0;
    #endif
    ////////////////////////////////////////
    // iterative execution
    ////////////////////////////////////////
   // graph.run(pr, params::rounds, params::iters);
    cpu_timer iter_timer;
    for(int r = 0; r < params::rounds; r++) {
        pr.init();
        timer.start();
        for(int i = 0; i < params::iters; i++) {
            graph.computeBlock(pr);
            pr.swapFunc();
        }
     
        current_time = timer.elapsed().wall/(1e9);
   //    if(params::is_filter) graph.postloop(pr);
        std::cout << "round " << r + 1 << ": " <<  params::input_file << ", processing time: " << current_time << '\n';
        total_time += current_time;
    }
    std::cout << "average time: " << total_time / params::rounds 
    << "\n----------------------\n";
    if(params::verify) pr.verify();
    return 0;

}
