/**
 * @brief this code is incorrect
 *
 */

#include <sys/_types/_dev_t.h>
#define DENSE
// #undef DENSE
#define DUMP
#define THREAD
#undef THREAD
unsigned int numIter = 0;

#include "../include/pcp.h"
#include "../include/sort.hpp"

float damping = 0.15;

struct HITS_F1
{
  float* auth{ nullptr };    // [num_regular]
  float* hub{ nullptr };     // [num_regular]

  HITS_F1(float* _auth, float* _hub) : auth(_auth), hub(_hub) {}
  // ~HITS_F1() { delete[] auth_cache;}s
  inline float scatterFunc(intV node)
  {
    return hub[node];
  }

  inline void initFunc(intV node)
  {
    hub[node] = 0.0;
  }

  inline void gatherFunc(float updateVal, intV vertex)
  {
    auth[vertex] += updateVal;
  }

  inline void filterFunc(intV node) {}

  void reset(intV num_verx)
  {
    std::fill(auth, auth + num_verx, 0.00000001);
    std::fill(hub, hub + num_verx, 0.00000001);
  }
};

struct HITS_F2
{
  float* auth{ nullptr };    // [num_regular]
  float* hub{ nullptr };     // [num_regular]

  HITS_F2(float* _auth, float* _hub) : auth(_auth), hub(_hub) {}

  inline float scatterFunc(intV node)
  {
    return auth[node];
  }

  inline void initFunc(intV node)
  {
    auth[node] = 0.0;
  }

  inline void gatherFunc(float updateVal, intV vertex)
  {
    hub[vertex] += updateVal;    // here vertex is the dst, but we want it to be src
  }

  inline void filterFunc(intV node) {}
};

int main(int argc, char** argv)
{
  graph<float> G;
  initialize(&G, argc, argv);
  struct timespec pre_begin, pre_end;
  if(clock_gettime(CLOCK_REALTIME, &pre_begin) == -1)
  {
    printf("clock gettime");
  }
  initBin<float>(&G);
  if(clock_gettime(CLOCK_REALTIME, &pre_end) == -1)
  {
    printf("clock gettime");
  }
  float pre_time =
      (pre_end.tv_sec - pre_begin.tv_sec) + ( int )(pre_end.tv_nsec - pre_begin.tv_nsec) / 1e9;
  printf("preprocessing@ %lf\n", pre_time);

  const int   dim    = 1;    // assume dim = 1
  const float lambda = 0.001;
  const float step   = 0.00000035;
  float*      auth   = new float[G.numVertex * dim]();
  float*      hub    = new float[G.numVertex * dim]();

  struct timespec start, end;
  float           time;
  float           sum_time = 0;

  HITS_F1 hits_f1(auth, hub);
  HITS_F2 hits_f2(auth, hub);

  int ctr = 0;
  while(ctr < G.rounds)
  {
    numIter = 0;

    hits_f1.reset(G.numVertex);

    if(clock_gettime(CLOCK_REALTIME, &start) == -1)
    {
      printf("clock gettime");
    }

    while(numIter < MAX_ITER)
    {
      scatter_and_gather<float>(&G, hits_f1);
      scatter_and_gather<float>(&G, hits_f2);
      numIter++;
      // todo: pull
    }
    if(clock_gettime(CLOCK_REALTIME, &end) == -1)
    {
      printf("clock gettime");
    }
    time = (end.tv_sec - start.tv_sec) + ( int )(end.tv_nsec - start.tv_nsec) / 1e9;
    printf("pr_dense, %d, %s, %lf\n", NUM_THREADS, input, time);

    if(ctr != 0)
      sum_time += time;
    ctr++;
  }
  std::cout << "The average running time of " << ctr << " rounds is: " << sum_time / (ctr - 1)
            << std::endl;

  printf("\n");

#ifdef DUMP

  mergeSortWOkey<float>(auth, 0, G.numVertex - 1);
  FILE* fdump = fopen("dumpHITS.txt", "w");
  if(fdump == nullptr)
  {
    fputs("file hub\n", stderr);
    exit(1);
  }
  int printVertices = (G.numVertex > 1000) ? 1000 : G.numVertex;
  for(int i = 0; i < printVertices; i++)
    fprintf(fdump, "%lf\n", auth[i]);
  fclose(fdump);
  std::cout << "result is dumped into dumpHTIS.txt" << '\n';

#endif

  delete[] auth;
  delete[] hub;
  return 0;
}
