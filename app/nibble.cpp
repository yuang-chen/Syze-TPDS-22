/**
 * Author: Kartik Lakhotia
           Sourav Pati
 * Email id: klakhoti@usc.edu
             spati@usc.edu
 * Date: 27-Feb-2018
 *
 * This code implements nibble algorithm
 * for computing probability distribution of 
 * a seeded random walk
 */


#define DUMP
#undef DUMP

#define TEST
unsigned int numIter = 0;


float threshold = 0.000000001;

#include "../include/pcp.h"



struct PR_F{
    float* pageRank;
    float* pageRankScat;
    intV* deg;
    PR_F(float* _pcurr, float* _pscat, intV* _outDeg):pageRank(_pcurr), pageRankScat(_pscat), deg(_outDeg){}
    inline float scatterFunc (intV node)
    {
        return pageRankScat[node];
    }
    inline bool initFunc(intV node)
    {
        pageRank[node]=pageRank[node]/2;
        pageRankScat[node] = 0;
        return (pageRank[node] >= threshold*deg[node]);
    }
    inline bool gatherFunc (float updateVal, intV destId)
    {
        pageRank[destId] += updateVal;
        return (updateVal > 0);
    }
    inline bool filterFunc(intV node)
    {
        bool cond = (pageRank[node] >= threshold*deg[node]);
        if (!cond)
            pageRank[node] = 0;
        if (cond && (deg[node]>0))
            pageRankScat[node] = pageRank[node]/(2*deg[node]);
        return cond; 
    } 
};




int main(int argc, char** argv)
{

    // for(int i = 0; i < 3; i++) { 
    // boost::timer::cpu_timer timer;
    // boost::timer::cpu_times times;
    // timer.start();



    graph<float> G;
    initialize(&G, argc, argv);    
    initBin<float>(&G);    
    float* pcurr = new float [G.numVertex]();
    float* pscat = new float [G.numVertex]();
    for (intV i=0; i<G.numVertex; i++)
    {
        pcurr[i] = 0;
        pscat[i] = 0;
    }
    intV initFrontierSize = 1;
    intV* initFrontier = new intV [initFrontierSize];
    intV actVertex;
    intV n = G.numVertex;

    for (intV i=0; i<initFrontierSize; i++)
        initFrontier[i] = G.start;

    struct timespec start, end;
    float time, aver_time=0;

    int ctr =0;
    while(ctr < G.rounds)
    {
    /* double activity_max = 0.0;
        double activity_min = 100000;
        double activity_avr = 0.0;
    */
        resetFrontier(&G);

        for (intV i=0; i<initFrontierSize; i++)
        {
            actVertex = initFrontier[i];
            pcurr[actVertex] = 1;
            if (G.outDeg[actVertex] > 0)
                pscat[actVertex] = pcurr[actVertex]/(2*G.outDeg[actVertex]);
        }
        loadFrontier (&G, initFrontier, initFrontierSize);

    // times = timer.elapsed();
    // double pre_time = (double) times.wall/(1e9);
    // std::cout << " preprocessing time (s): " << pre_time << std::endl; 
    // }
    
      //  double s_t_max, s_t_min, s_t_avr, g_t_max, g_t_min, g_t_avr;
      //  double s_t_max_total = 0, s_t_min_total = 0, s_t_avr_total = 0, g_t_max_total = 0, g_t_min_total = 0, g_t_avr_total = 0;
        numIter = 0;

        if( clock_gettime(CLOCK_REALTIME, &start) == -1) { perror("clock gettime");}

        while(numIter < MAX_ITER)//G.frontierSize > 0
        {
            scatter_and_gather<float>(&G, PR_F(pcurr, pscat, G.outDeg)); 
            numIter++;

/*            std::tie(s_t_max, s_t_min, s_t_avr, g_t_max, g_t_min, g_t_avr) = scatter_and_gather<float>(&G, PR_F(pcurr, pscat, G.outDeg));   
            numIter++;
            if(s_t_max == 0 || s_t_min == 0 || g_t_max == 0 || g_t_min == 0)
                break;
           // std::cout << (double) G.frontierSize << std::endl;
            s_t_max_total += s_t_max;
            s_t_min_total += s_t_min;
            s_t_avr_total += s_t_avr; 
            g_t_max_total += g_t_max; 
            g_t_min_total += g_t_min; 
            g_t_avr_total += g_t_avr; */
        }   

        getFrontier(&G);
        #pragma omp parallel for
        for (intV i=0; i<G.frontierSize; i++)
        {
            pcurr[G.frontier[i]]=0;
            pscat[G.frontier[i]]=0;
        }
 
        if( clock_gettime( CLOCK_REALTIME, &end) == -1 ) { perror("clock gettime");}
        time = (end.tv_sec - start.tv_sec)+ (int)(end.tv_nsec - start.tv_nsec)/1e9;
        printf("appr, %d, %s, %lf\n",NUM_THREADS, argv[1], time);
        aver_time+=time;
      /*  std::cout << "max, min, average time for scatter and gather threads: "
                  <<  s_t_max_total/MAX_ITER  << " " <<  s_t_min_total/MAX_ITER  << " " <<  s_t_avr_total/MAX_ITER  << " " << 
                   g_t_max_total/MAX_ITER  << " " <<  g_t_min_total/MAX_ITER  << " " <<  g_t_avr_total/MAX_ITER << std::endl;   
                   */
        ctr++;
    }
    std::cout << "The average running time of " << ctr << " rounds is: " << aver_time/ctr << std::endl;

    printf("\n");

    return 0;
}
