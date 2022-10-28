/**
 * Author: Kartik Lakhotia
           Sourav Pati
 * Email id: klakhoti@usc.edu
             spati@usc.edu
 * Date: 27-Feb-2018
 *
 * Compute PageRank using Partition-centric graph processing
 *
 */

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
    #ifdef THREAD
    double smax = 0, smin = 0, saver=0, gmax=0, gmin=0, gaver=0;
    double s_t_max, s_t_min, s_t_avr, g_t_max, g_t_min, g_t_avr;
    double s_t_max_total = 0, s_t_min_total = 0, s_t_avr_total = 0, g_t_max_total = 0, g_t_min_total = 0, g_t_avr_total = 0;
    #endif
    int ctr =0;
    while(ctr < G.rounds)
    {
         #ifdef THREAD
            s_t_max_total = 0;
            s_t_min_total = 0;
            s_t_avr_total = 0; 
            g_t_max_total = 0; 
            g_t_min_total = 0; 
            g_t_avr_total = 0; 
        #endif
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
         #ifndef THREAD
            scatter_and_gather<float>(&G, PR_F(pcurr, G.outDeg));
            #else
           std::tie(s_t_max, s_t_min, s_t_avr, g_t_max, g_t_min, g_t_avr) = scatter_and_gather_thread<float>(&G, PR_F(pcurr, G.outDeg));

            s_t_max_total += s_t_max;
            s_t_min_total += s_t_min;
            s_t_avr_total += s_t_avr; 
            g_t_max_total += g_t_max; 
            g_t_min_total += g_t_min; 
            g_t_avr_total += g_t_avr;        
            #endif
            numIter++;
        }

        if( clock_gettime( CLOCK_REALTIME, &end) == -1 ) { perror("clock gettime");}
        time = (end.tv_sec - start.tv_sec)+ (int)(end.tv_nsec - start.tv_nsec)/1e9;
        printf("pr_dense, %d, %s, %lf\n",NUM_THREADS, input, time);
         #ifdef THREAD
                std::cout << "max, min, average time for scatter and gather threads (in ms): "
                  <<  s_t_max_total/MAX_ITER  << " " <<  s_t_min_total/MAX_ITER  << " " <<  s_t_avr_total/MAX_ITER  << " " << 
                   g_t_max_total/MAX_ITER  << " " <<  g_t_min_total/MAX_ITER  << " " <<  g_t_avr_total/MAX_ITER << std::endl;
        
        if(ctr!=0) {
            smax+=s_t_max_total/numIter;
            smin+=s_t_min_total/numIter;
            saver+=s_t_avr_total/numIter;
            gmax+=g_t_max_total/numIter;
            gmin+=g_t_min_total/numIter;
            gaver+=g_t_avr_total/numIter;
        }
        #endif
        if(ctr!=0) 
            sum_time += time;
        ctr++;

    }
    std::cout << "The average running time of " << ctr << " rounds is: " << sum_time/(ctr-1) << std::endl;
    #ifdef THREAD
    std::cout << "THREAD: " << smax/(ctr-1) << " " << smin/(ctr-1) << " " << saver/(ctr-1)
                << " " << gmax/(ctr-1) << " " << gmin/(ctr-1) << " " << gaver/(ctr-1) << '\n';
    #endif
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



