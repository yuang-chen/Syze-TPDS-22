#define DENSE
//#undef DENSE
//#define DUMP
#define THREAD
#undef THREAD
unsigned int numIter = 0;

#include "../include/pcp.h"
#include "../include/sort.hpp"

float damping =0.15;

struct PR_F{
    float* pageRank;
    intV* deg;
    PR_F(float* _pcurr, intV* _outDeg):pageRank(_pcurr), deg(_outDeg){}
    inline float scatterFunc (intV node)
    {
        return pageRank[node];
    }
#ifndef DENSE
    inline bool initFunc(intV node)
    {
        pageRank[node]=0;
        return true;
    }
    inline bool gatherFunc (float updateVal, intV destId)
    {
        pageRank[destId] += updateVal;
        return true;
    }
    inline bool filterFunc(intV node)
    {
        pageRank[node] = ((damping) + (1-damping)*pageRank[node]);
        if (deg[node]>0)
            pageRank[node] = pageRank[node]/deg[node];
        return true;
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
        pageRank[node] = ((damping) + (1-damping)*pageRank[node]);
        if (deg[node]>0)
            pageRank[node] = pageRank[node]/deg[node];
    }
#endif
};




int main(int argc, char** argv)
{
    graph<float> G;
    initialize(&G, argc, argv);    
    struct timespec pre_begin, pre_end;
    if( clock_gettime(CLOCK_REALTIME, &pre_begin) == -1) { perror("clock gettime");}
    initBin<float>(&G);    
    if( clock_gettime( CLOCK_REALTIME, &pre_end) == -1 ) { perror("clock gettime");}
    float pre_time = (pre_end.tv_sec - pre_begin.tv_sec)+ (int)(pre_end.tv_nsec - pre_begin.tv_nsec)/1e9;
    printf("preprocessing@ %lf\n", pre_time);
    float* pcurr = new float [G.numVertex]();
    #pragma omp parallel for
    for(int i=0;i<G.numVertex;i++){
        if (G.outDeg[i] > 0)
            pcurr[i] = 1.0/G.outDeg[i];
        else
            pcurr[i] = 1.0;
    }
#ifndef DENSE
    intV* initFrontier = new intV [G.numVertex];
#pragma omp parallel for
    for (intV i=0; i<G.numVertex; i++)
        initFrontier[i] = i;
    loadFrontier (&G, initFrontier, G.numVertex);
#endif

    struct timespec start, end;
    float time;
    float sum_time = 0;

    int ctr =0;
    while(ctr < G.rounds)
    {
        numIter = 0;
        for(int i=0;i<G.numVertex;i++){
            if (G.outDeg[i] > 0)
                pcurr[i] = 1.0/G.outDeg[i];
            else
                pcurr[i] = 1.0;
        }
        if( clock_gettime(CLOCK_REALTIME, &start) == -1) { perror("clock gettime");}

        while(numIter < MAX_ITER)
        {
            scatter_and_gather<float>(&G, PR_F(pcurr, G.outDeg));
            numIter++;
        }

        if( clock_gettime( CLOCK_REALTIME, &end) == -1 ) { perror("clock gettime");}
        time = (end.tv_sec - start.tv_sec)+ (int)(end.tv_nsec - start.tv_nsec)/1e9;
        printf("pr_dense, %d, %s, %lf\n",NUM_THREADS, input, time);
        if(ctr!=0) 
            sum_time += time;
        ctr++;

    }
    std::cout << "The average running time of " << ctr << " rounds is: " << sum_time/(ctr-1) << std::endl;
    
    printf("\n");
#ifdef DUMP
    for (unsigned int i=0; i<G.numVertex; i++)
    {
        if (G.outDeg[i]>0)
            pcurr[i] = pcurr[i]*G.outDeg[i];
    }
    mergeSortWOkey<float>(pcurr, 0, G.numVertex-1);
    FILE* fdump = fopen("dumpPR.txt", "w");
    if (fdump == NULL)
    {
        fputs("file error\n", stderr);
        exit(1);
    }
    int printVertices = (G.numVertex>1000) ? 1000 : G.numVertex;
    for (int i=0; i<printVertices; i++)
        fprintf(fdump, "%lf\n", pcurr[i]);
    fclose(fdump);
#endif
    return 0;
}



