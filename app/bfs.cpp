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

    int ctr = 0;
    while(ctr < G.rounds){
        resetFrontier(&G);
        std::fill(parent, parent+n, 1<<MSB_ROT);
        std::fill(visited, visited+n, false);

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
            scatter_and_gather<intV>(&G, BFS_F(parent, visited));
            numIter++;
        }   

        if( clock_gettime( CLOCK_REALTIME, &end) == -1 ) { perror("clock gettime");}
        time = (end.tv_sec - start.tv_sec)+ (int)(end.tv_nsec - start.tv_nsec)/1e9;
        if(ctr!=0) sum_time += time;
        printf("bfs, %d, %s, %lf, %d\n", NUM_THREADS, input, time, numIter);
        ctr++; 
    }

    std::cout << "The average running time of " << ctr << " rounds is: " << sum_time/(ctr-1) <<   std::endl;

    printf("\n");
    return 0;
}



