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

        #ifdef THREAD
        double s_t_max, s_t_min, s_t_avr, g_t_max, g_t_min, g_t_avr;
        double s_t_max_total = 0, s_t_min_total = 0, s_t_avr_total = 0, g_t_max_total = 0, g_t_min_total = 0, g_t_avr_total = 0;
        double activity_max = 0.0;
        double activity_min = 100000;
        double activity_avr = 0.0;
        #endif 
        if( clock_gettime(CLOCK_REALTIME, &start) == -1) { perror("clock gettime");}

        while(numIter < MAX_ITER && G.frontierSize > 0)
        {
         #ifndef THREAD
            scatter_and_gather<float>(&G, PR_F(pcurr, old_pr, G.outDeg));
            numIter++;
        #else
            std::tie(s_t_max, s_t_min, s_t_avr, g_t_max, g_t_min, g_t_avr) = scatter_and_gather_thread<float>(&G, PR_F(pcurr, sum, delta, G.outDeg, one_over_n));
            s_t_max_total += s_t_max;
            s_t_min_total += s_t_min;
            s_t_avr_total += s_t_avr; 
            g_t_max_total += g_t_max; 
            g_t_min_total += g_t_min; 
            g_t_avr_total += g_t_avr; 
 
            activity_max = fmax(activity_max, (double) G.frontierSize/n);
            activity_min = fmin(activity_min, (double) G.frontierSize/n);
            activity_avr += (double) G.frontierSize/n;        
        #endif
        }

        if( clock_gettime( CLOCK_REALTIME, &end) == -1 ) { perror("clock gettime");}
        time = (end.tv_sec - start.tv_sec)+ (int)(end.tv_nsec - start.tv_nsec)/1e9;
        printf("prdelta, %d, %s, %lf\n", numIter, input, time);
        average_time += time;

        #ifdef THREAD
        std::cout << "max min and average activity: " << activity_max << " " << activity_min << " " << activity_avr/numIter <<  std::endl;         std::cout << "max, min, average time for scatter and gather threads: "
                 <<  s_t_max_total/MAX_ITER  << " " <<  s_t_min_total/MAX_ITER  << " " <<  s_t_avr_total/MAX_ITER  << " " << 
                   g_t_max_total/MAX_ITER  << " " <<  g_t_min_total/MAX_ITER  << " " <<  g_t_avr_total/MAX_ITER << std::endl;
        #endif
        ctr++;
    }
    //if(ctr > 3)
    std::cout << "The average running time of " << ctr << " rounds is: " << average_time/ctr << std::endl;
#ifdef DUMP
/*     for (unsigned int i=0; i<G.numVertex; i++)
    {
        if (G.outDeg[i]>0)
            pcurr[i] = pcurr[i]*G.outDeg[i];
    } */
    mergeSortWOkey<float>(pcurr, 0, G.numVertex-1);
    FILE* fdump = fopen("app/PRD.txt", "w");
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
    printf("\n");

    return 0;
}



