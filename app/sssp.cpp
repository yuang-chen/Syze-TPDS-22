/**
 * Author: Kartik Lakhotia
           Sourav Pati
 * Email id: klakhoti@usc.edu
             spati@usc.edu
 * Date: 27-Feb-2018
 *
 * Single source shortest path using Bellman-Ford 
 *
 */

unsigned int numIter = 0;

#define WEIGHTED
//#define TEST
//for asynchronous update propagation//
//converges faster//
#define ASYNCH
//#define THREAD
#include "../include/pcp.h"



struct SSSP_F{
    unsigned int* distance;
    SSSP_F(unsigned int* _distance):distance(_distance){}

    inline unsigned int scatterFunc (intV node)
    {
        return distance[node];
    }

    inline bool initFunc(intV node)
    {
        return false;
    }

    inline bool gatherFunc (unsigned int updateVal, intV destId)
    {
        if(updateVal < distance[destId])
        {
            distance[destId] = updateVal;
            return true;
        }
        else
            return false;
    }  
    
    inline bool filterFunc(intV node)
    {
        return true;
    } 

    inline unsigned int applyWeight (unsigned int updateVal, unsigned int weight)
    {
        return updateVal + weight;
    }

};



int main(int argc, char** argv)
{
    graph<unsigned int> G;
    initialize(&G, argc, argv);
    initBin<unsigned int>(&G);    
    intV n = G.numVertex;
    unsigned int* distance = new unsigned int [n]();
    intV initFrontierSize = 1;
    intV* initFrontier = new intV [initFrontierSize];
    for (intV i=0; i<initFrontierSize; i++)
        initFrontier[i] = G.start;  

    loadFrontier(&G, initFrontier, initFrontierSize);
   
    struct timespec start, end;
    float time,aver_time=0;
          #ifdef THREAD
    double smax = 0, smin = 0, saver=0, gmax=0, gmin=0, gaver=0;
    double s_t_max, s_t_min, s_t_avr, g_t_max, g_t_min, g_t_avr;
    double s_t_max_total = 0, s_t_min_total = 0, s_t_avr_total = 0, g_t_max_total = 0, g_t_min_total = 0, g_t_avr_total = 0;
    #endif  
    int ctr = 0;
    while(ctr < G.rounds){
        resetFrontier(&G);
         numIter=0;
         std::fill(distance, distance+n, 1<<31);
         distance[G.start] = 0;

         loadFrontier(&G, initFrontier, initFrontierSize);
        

         if( clock_gettime(CLOCK_REALTIME, &start) == -1) { perror("clock gettime");}

         while((G.frontierSize > 0) && (numIter < G.numVertex))
         {
            #ifndef THREAD
             scatter_and_gather<unsigned int>(&G, SSSP_F(distance));
            #else
            std::tie(s_t_max, s_t_min, s_t_avr, g_t_max, g_t_min, g_t_avr) = scatter_and_gather_thread<unsigned int>(&G, SSSP_F(distance));
         //   activity_max = fmax(activity_max, (double) G.frontierSize/n);
         //   activity_min = fmin(activity_min, (double) G.frontierSize/n);
         //   activity_avr += (double) G.frontierSize/n;
            if(s_t_max == 0 || s_t_min == 0 || g_t_max == 0 || g_t_min == 0)
                break;
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
        printf("sssp, %d, %s, %lf\n", NUM_THREADS, input, time);
        if(ctr!=0) aver_time+=time;
        #ifdef THREAD
     //   std::cout << "max min and average activity: " << activity_max << " " << activity_min << " " << activity_avr/numIter <<  std::endl;
        std::cout << "max, min, average time for scatter and gather threads: "
                  <<  s_t_max_total/numIter  << " " <<  s_t_min_total/numIter  << " " <<  s_t_avr_total/numIter  << " " << 
                   g_t_max_total/numIter  << " " <<  g_t_min_total/numIter  << " " <<  g_t_avr_total/numIter << std::endl;    
                  if(ctr!=0) {
            smax+=s_t_max_total/numIter;
            smin+=s_t_min_total/numIter;
            saver+=s_t_avr_total/numIter;
            gmax+=g_t_max_total/numIter;
            gmin+=g_t_min_total/numIter;
            gaver+=g_t_avr_total/numIter;
        }
                 
        #endif
        ctr++;
    }
    printf("\n");
    std::cout << "The average running time of " << ctr << " rounds is: " << aver_time/(ctr-1) << std::endl;
    #ifdef THREAD
    std::cout << "THREAD: " << smax/(ctr-1) << " " << smin/(ctr-1) << " " << saver/(ctr-1)
                << " " << gmax/(ctr-1) << " " << gmin/(ctr-1) << " " << gaver/(ctr-1) << '\n';
    #endif
   // FILE* fp = fopen("dump.bin", "wb");
   // fwrite(distance, sizeof(unsigned int), n, fp);
  //  fclose(fp);
    
    return 0;
}