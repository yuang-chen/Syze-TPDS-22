/**
 * Author: Kartik Lakhotia
           Sourav Pati
 * Email id: klakhoti@usc.edu
             spati@usc.edu
 * Date: 27-Feb-2018
 *
 * Weakly connected components
 * 
 */

unsigned int numIter = 0;

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
    //for(int i = 0; i < 3; i++) { 
    boost::timer::cpu_timer timer;
    boost::timer::cpu_times times;
    timer.start();

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
        #ifdef THREAD
    double smax = 0, smin = 0, saver=0, gmax=0, gmin=0, gaver=0;
    double s_t_max, s_t_min, s_t_avr, g_t_max, g_t_min, g_t_avr;
    double s_t_max_total = 0, s_t_min_total = 0, s_t_avr_total = 0, g_t_max_total = 0, g_t_min_total = 0, g_t_avr_total = 0;
    #endif
    while((G.frontierSize > 0))
    {
         scatter_and_gather<intV>(&G, CC_F(label));
         numIter++;
    }
    int ctr =0;
    while(ctr < G.rounds){
        resetFrontier(&G);
        std::iota (label, label+n, 0);
        numIter = 0;
 /*       double activity_max = 0.0;
        double activity_min = 100000;
        double activity_avr = 0.0;
   */ 
        loadFrontier(&G, initFrontier, initFrontierSize);

    // times = timer.elapsed();
    // double pre_time = (double) times.wall/(1e9);
    // std::cout << " preprocessing time (s): " << pre_time << std::endl; 
    
    
        if( clock_gettime(CLOCK_REALTIME, &start) == -1) { perror("clock gettime");}

        while(G.frontierSize > 0)
        {
            numIter++;
         #ifdef THREAD

            std::tie(s_t_max, s_t_min, s_t_avr, g_t_max, g_t_min, g_t_avr) = scatter_and_gather_thread<intV>(&G, CC_F(label));
        //    activity_max = fmax(activity_max, (double) G.frontierSize/n);
        //    activity_min = fmin(activity_min, (double) G.frontierSize/n);
        //    activity_avr += (double) G.frontierSize/n;
            numIter++;

            if(s_t_max == 0 || s_t_min == 0 || g_t_max == 0 || g_t_min == 0)
                break;
            s_t_max_total += s_t_max;
            s_t_min_total += s_t_min;
            s_t_avr_total += s_t_avr; 
            g_t_max_total += g_t_max; 
            g_t_min_total += g_t_min; 
            g_t_avr_total += g_t_avr;
            #else
             scatter_and_gather<intV>(&G, CC_F(label));

            #endif
        }
       // std::cout << "max min and average activity: " << activity_max << " " << activity_min << " " << activity_avr/numIter <<  std::endl;
        if( clock_gettime( CLOCK_REALTIME, &end) == -1 ) { perror("clock gettime");}
        time = (end.tv_sec - start.tv_sec)+ (int)(end.tv_nsec - start.tv_nsec)/1e9;
        if(ctr!=0) aver_time+=time;
        printf("cc, %d, %s, %lf\n", NUM_THREADS, input, time);
         #ifdef THREAD
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

    return 0;
}



