//for asynchronous update propagation//
//converges faster//
#define ASYNCH

//#define TEST
#include "../include/pcp.h"
#include <numeric>


//#define THREAD
struct CC_F{
    intV* label;
    CC_F(intV* _label):label(_label){}

    inline intV scatterFunc (intV node)
    {
        return label[node];
    }

    inline bool initFunc(intV node)
    {
        return false;
    }

    inline bool gatherFunc (intV updateVal, intV destId)
    {
        bool cond = (updateVal < label[destId]);
        if (cond)
            label[destId] = updateVal;
        return cond;
    }  
    
    inline bool filterFunc(intV node)
    {
        return true;
    } 

};


int main(int argc, char** argv)
{
    graph<intV> G;
    initialize(&G, argc, argv);
    initBin<intV>(&G);
    intV n = G.numVertex;
    intV* label = new intV [n]();
    intV initFrontierSize = n;
    intV* initFrontier = new intV [initFrontierSize];
    for (intV i=0; i<initFrontierSize; i++)
        initFrontier[i] = i;

    struct timespec start, end, half;
    float time, aver_time=0;

    int ctr =0, numIter = 0;
    while(ctr < G.rounds){
        resetFrontier(&G);
        std::iota (label, label+n, 0);
        numIter = 0;

        loadFrontier(&G, initFrontier, initFrontierSize);    
    
        if( clock_gettime(CLOCK_REALTIME, &start) == -1) { perror("clock gettime");}

        while(G.frontierSize > 0)
        {
            numIter++;
             scatter_and_gather<intV>(&G, CC_F(label));
        }
        if( clock_gettime( CLOCK_REALTIME, &end) == -1 ) { perror("clock gettime");}
        time = (end.tv_sec - start.tv_sec)+ (int)(end.tv_nsec - start.tv_nsec)/1e9;
        if(ctr!=0) aver_time+=time;
        printf("cc, %d, %s, %lf\n", NUM_THREADS, input, time);
        ctr++;
    }
    printf("\n");
    std::cout << "The average running time of " << ctr << " rounds is: " << aver_time/(ctr-1) << std::endl;

    return 0;
}



