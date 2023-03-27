/**
 * eigenvector centrality
 */

#include <numeric>
#define DENSE
// #undef DENSE
#define DUMP
#define THREAD
#undef THREAD
unsigned int numIter = 0;

#include "../include/pcp.h"
#include "../include/sort.hpp"

struct EC_F
{
  float* centrality;
  float  norm;
  intV   size;
  EC_F(float* _cen, intV _size) : centrality(_cen), size(_size) {}
  inline float scatterFunc(intV node)
  {
    return centrality[node];
  }

  inline void initFunc(intV node)
  {
    centrality[node] = 0;
  }

  inline void gatherFunc(float updateVal, intV destId)
  {
    centrality[destId] += updateVal;
  }

  inline void filterFunc(intV node)
  {
    //centrality[node] *= centrality[node];
  }

  inline void reset()
  {
    std::fill(centrality, centrality + size, 0.000000001 / size);
    norm = 1;
  }
  inline void normalize()
  {
    norm = std::inner_product(centrality, centrality + size, centrality, 0.0);
    norm = std::sqrt(norm);
  }
};

int main(int argc, char** argv)
{
  graph<float> G{};
  initialize(&G, argc, argv);
  struct timespec pre_begin;
  struct timespec pre_end;
  if(clock_gettime(CLOCK_REALTIME, &pre_begin) == -1)
  {
    perror("clock gettime");
  }
  initBin<float>(&G);
  if(clock_gettime(CLOCK_REALTIME, &pre_end) == -1)
  {
    perror("clock gettime");
  }
  float pre_time =
      (pre_end.tv_sec - pre_begin.tv_sec) + ( int )(pre_end.tv_nsec - pre_begin.tv_nsec) / 1e9;
  printf("preprocessing@ %lf\n", pre_time);
  auto* centrality = new float[G.numVertex]();
  EC_F  ec_f(centrality, G.numVertex);

  struct timespec start;
  struct timespec end;
  float           time;
  float           sum_time = 0;

  int ctr = 0;
  while(ctr < G.rounds)
  {
    numIter = 0;
    ec_f.reset();
    if(clock_gettime(CLOCK_REALTIME, &start) == -1)
    {
      perror("clock gettime");
    }

    while(numIter < MAX_ITER)
    {
      scatter_and_gather<float>(&G, ec_f);
      ec_f.normalize();
      numIter++;
    }

    if(clock_gettime(CLOCK_REALTIME, &end) == -1)
    {
      perror("clock gettime");
    }
    time = (end.tv_sec - start.tv_sec) + ( int )(end.tv_nsec - start.tv_nsec) / 1e9;
    printf("pr_dense, %d, %s, %lf\n", NUM_THREADS, input, time);
    if(ctr != 0)
    {
      sum_time += time;
    }
    ctr++;
  }
  std::cout << "The average running time of " << ctr << " rounds is: " << sum_time / (ctr - 1)
            << std::endl;

  printf("\n");

#ifdef DUMP

  mergeSortWOkey<float>(centrality, 0, G.numVertex - 1);
  FILE* fdump = fopen("dumpEC.txt", "w");
  if(fdump == nullptr)
  {
    fputs("file error\n", stderr);
    exit(1);
  }
  int printVertices = (G.numVertex > 1000) ? 1000 : G.numVertex;
  for(int i = 0; i < printVertices; i++)
    fprintf(fdump, "%lf\n", centrality[i]);
  fclose(fdump);
  std::cout << "result is dumped into dumpEC.txt" << '\n';
#endif

  delete[] centrality;
  return 0;
}