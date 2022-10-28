/**
 * Author: Kartik Lakhotia
           Sourav Pati
 * Email id: klakhoti@usc.edu
             spati@usc.edu
 * Date: 27-Feb-2018
 *
 * This code implements work efficient BFS
 * 
 */

unsigned int numIter = 0;

//#define TEST

#include "../include/pcp.h"

//#define THREAD

struct BFS_F{
    intV* parent;
    bool* visited;
    BFS_F(intV* _parent, bool* _visited):parent(_parent), visited(_visited){}
    inline intV scatterFunc (intV node)
    {
        return (((!visited[node])<<MSB_ROT) | node);
    }

    inline bool initFunc(intV node)
    {
        return false;
    }

    inline bool gatherFunc (intV updateVal, intV destId)
    {
//        if((!visited[destId]) && (!(updateVal>>MSB_ROT)))
//        {
//            parent[destId] = updateVal;
//            visited[destId] = true;
//            return true;
//        }
        if (!visited[destId]) //if destination vertex is not yet visited
        {
            parent[destId] = updateVal; //set its parent
            visited[destId] = (!(updateVal>>MSB_ROT)); //new visited status depends on parent's status
            return visited[destId]; //active if it is now visited
        }
        return false;
    }  
    
    inline bool filterFunc(intV node)
    {
        return true;
    } 

};

struct BFS_F2{
    intV* parent;
    BFS_F2(intV* _parent):parent(_parent){}
    inline intV scatterFunc (intV node)
    {
        return ((parent[node] & MAX_NEG)| node);
    }

    inline bool initFunc(intV node)
    {
        return false;
    }

    inline bool gatherFunc (intV updateVal, intV destId)
    {
        if((parent[destId] & MAX_NEG) && (!(updateVal>>MSB_ROT))) //if parent of destID is negative and received update from a visited node
        {
            parent[destId] = updateVal; //update parent
            return true;
        }
        return false;
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
    intV* parent = new intV [n]();
    bool* visited = new bool [n]();
    intV initFrontierSize = 1;
    intV* initFrontier = new intV [initFrontierSize];

    std::vector<float> vertex_time;


  //  for(auto start_ver: start_vertex) {
    float sum_time = 0;
    for (intV i=0; i<initFrontierSize; i++)
        initFrontier[i] = G.start;  

    struct timespec start, end;
    float time;

    #ifdef THREAD
    double smax = 0, smin = 0, saver=0, gmax=0, gmin=0, gaver=0;
    double s_t_max, s_t_min, s_t_avr, g_t_max, g_t_min, g_t_avr;
    double s_t_max_total = 0, s_t_min_total = 0, s_t_avr_total = 0, g_t_max_total = 0, g_t_min_total = 0, g_t_avr_total = 0;
    #endif

    int ctr = 0;
    while(ctr < G.rounds){
        resetFrontier(&G);

        std::fill(parent, parent+n, 1<<MSB_ROT);
        std::fill(visited, visited+n, false);

/*         double activity_max = 0.0;
        double activity_min = 100000;
        double activity_avr = 0.0; */
         #ifdef THREAD
            s_t_max_total = 0;
            s_t_min_total = 0;
            s_t_avr_total = 0; 
            g_t_max_total = 0; 
            g_t_min_total = 0; 
            g_t_avr_total = 0; 
        #endif

        for(int i=0;i<n;i++)
        {
            parent[i] = 1<<MSB_ROT;
            visited[i] = false;
        }
        visited[G.start] = true;
        parent[G.start] = G.start;

        if( clock_gettime(CLOCK_REALTIME, &start) == -1) { perror("clock gettime");}

//        loadFrontier(&G, initFrontier, initFrontierSize);
        loadFrontierPar(&G, initFrontier, initFrontierSize);

        numIter=0;
        while((G.frontierSize > 0))
        {    
            #ifdef THREAD
            std::tie(s_t_max, s_t_min, s_t_avr, g_t_max, g_t_min, g_t_avr) = scatter_and_gather_thread<intV>(&G, BFS_F(parent, visited));
        //    activity_max = fmax(activity_max, (double) G.frontierSize/n);
        //    activity_min = fmin(activity_min, (double) G.frontierSize/n);
          //  activity_avr += (double) G.frontierSize/n;
            s_t_max_total += s_t_max;
            s_t_min_total += s_t_min;
            s_t_avr_total += s_t_avr; 
            g_t_max_total += g_t_max; 
            g_t_min_total += g_t_min; 
            g_t_avr_total += g_t_avr; 
            #else
            scatter_and_gather<intV>(&G, BFS_F(parent, visited));
            #endif
            numIter++;
        }   

        if( clock_gettime( CLOCK_REALTIME, &end) == -1 ) { perror("clock gettime");}
        time = (end.tv_sec - start.tv_sec)+ (int)(end.tv_nsec - start.tv_nsec)/1e9;
        if(ctr!=0) sum_time += time;
        printf("bfs, %d, %s, %lf, %d\n", NUM_THREADS, input, time, numIter);
      
        #ifdef THREAD
      //  std::cout << "max min and average activity: " << activity_max << " " << activity_min << " " << activity_avr/numIter <<  std::endl;
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
    /*   printf("average time for vertex %d: %lf\n", start_ver, sum_time/G.rounds);

        vertex_time.push_back(sum_time/G.rounds);
   }
    float aver_time = 0;
    for(auto t: vertex_time)
        aver_time += t;
    printf("average time for all: %lf\n", aver_time/start_vertex.size());*/
    std::cout << "The average running time of " << ctr << " rounds is: " << sum_time/(ctr-1) <<   std::endl;
    #ifdef THREAD
    std::cout << "THREAD: " << smax/(ctr-1) << " " << smin/(ctr-1) << " " << saver/(ctr-1)
                << " " << gmax/(ctr-1) << " " << gmin/(ctr-1) << " " << gaver/(ctr-1) << '\n';
    #endif
    printf("\n");
    return 0;
}



