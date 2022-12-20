
#define THREAD
#undef THREAD
#define TEST
//unsigned int numIter = 0;
//#define DUMP
#include "../include/pcp.h"
#include "../include/sort.hpp"


float damping =0.15;
float epsilon = 0.01;

struct PR_F{
    float* pageRank;
    float* oldPageRank;
    intV* deg;
    PR_F(float* _pcurr, float* _old, intV* _outDeg):
        pageRank(_pcurr),  oldPageRank(_old), deg(_outDeg){}
    inline float scatterFunc (intV node)
    {
        return pageRank[node];
    }
#ifndef DENSE
    inline bool initFunc(intV node)
    {
        oldPageRank[node] = pageRank[node];
        pageRank[node]=0;
        return true;
    }
    inline bool gatherFunc (float updateVal, intV destId)
    {
        //sum[destId] +=  updateVal / (float) deg[destId];
        pageRank[destId] += updateVal; 
        return true;
    }
    inline bool filterFunc(intV node)
    {
        //delta[node] = damping * sum[node] + constant;
        pageRank[node] = ((damping) + (1-damping)*pageRank[node]);
        if (deg[node]>0)
            pageRank[node] = pageRank[node]/deg[node];

        return fabs((oldPageRank[node] - pageRank[node]) > epsilon/deg[node]);
    } 
#else
    inline void initFunc(intV node)
    {
        pageRank[node]=0;
    }
    inline void gatherFunc (float updateVal, intV destId)
    {
        pageRank[destId] += updateVal;
    }
    inline void filterFunc(intV node)
    {
 
        pageRank[node] = ((damping_const) + (1-damping)*pageRank[node]);
        if (deg[node]>0)
            pageRank[node] = pageRank[node]/deg[node];
    }
#endif
};




int main(int argc, char** argv)
{
    graph<float> G;
    initialize(&G, argc, argv);    
    initBin<float>(&G);    
    float* pcurr = new float [G.numVertex]();
    float* old_pr = new float [G.numVertex]();
    std::cout << "delta threshold: " << epsilon << std::endl;
#ifndef DENSE
    intV* initFrontier = new intV [G.numVertex];
    intV initFrontierSize = G.numVertex;
    #pragma omp parallel for
    for (intV i=0; i<G.numVertex; i++)
        initFrontier[i] = i;
    loadFrontier (&G, initFrontier, G.numVertex);
#endif

    intV n = G.numVertex;

    struct timespec start, end;
    float time, average_time=0;

    int ctr =0;
/*   
  */  while(ctr < G.rounds)
    {
        resetFrontier(&G);

        int numIter = 0;
        intV initFrontierSize = n;

        #pragma omp parallel for
        for(int i=0;i<G.numVertex;i++){
            if (G.outDeg[i] > 0)
                pcurr[i] = 1.0/G.outDeg[i];
            else
                pcurr[i] = 1.0;
            //old_pr[i] = pcurr[i];
        }
        loadFrontier (&G, initFrontier, initFrontierSize);

        if( clock_gettime(CLOCK_REALTIME, &start) == -1) { perror("clock gettime");}

        while(numIter < MAX_ITER && G.frontierSize > 0)
        {
            scatter_and_gather<float>(&G, PR_F(pcurr, old_pr, G.outDeg));
            numIter++;
        }

        if( clock_gettime( CLOCK_REALTIME, &end) == -1 ) { perror("clock gettime");}
        time = (end.tv_sec - start.tv_sec)+ (int)(end.tv_nsec - start.tv_nsec)/1e9;
        printf("prdelta, %d, %s, %lf\n", numIter, input, time);
        average_time += time;
        ctr++;
    }
    std::cout << "The average running time of " << ctr << " rounds is: " << average_time/ctr << std::endl;

    printf("\n");

    return 0;
}



