#define DENSE
// #undef DENSE
#define DUMP
#define THREAD
#undef THREAD
unsigned int numIter = 0;

#include "../include/pcp.h"
#include "../include/sort.hpp"

struct INDEG_F
{
  intV* sum;
  INDEG_F(intV* _sum) : sum(_sum) {}
  inline float scatterFunc(intV node)
  {
    return sum[node];
  }

  inline void initFunc(intV node)
  {
    sum[node] = 0;
  }
  inline void gatherFunc(intV updateVal, intV destId)
  {
    sum[destId] += updateVal;
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
  intV* sum = new intV[G.numVertex]();

  struct timespec start, end;
  float           time;
  float           sum_time = 0;

  INDEG_F indeg_f(sum);
  int     ctr = 0;
  while(ctr < G.rounds)
  {
    numIter = 0;
    std::fill(sum, sum + G.numVertex, 1);
    if(clock_gettime(CLOCK_REALTIME, &start) == -1)
    {
      perror("clock gettime");
    }

    while(numIter < MAX_ITER)
    {
      std::fill(sum, sum + G.numVertex, 1);
      scatter_and_gather<float>(&G, indeg_f);
      numIter++;
    }

    if(clock_gettime(CLOCK_REALTIME, &end) == -1)
    {
      perror("clock gettime");
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

  mergeSortWOkey<intV>(sum, 0, G.numVertex - 1);
  FILE* fdump = fopen("dumpIN.txt", "w");
  if(fdump == NULL)
  {
    fputs("file error\n", stderr);
    exit(1);
  }
  int printVertices = (G.numVertex > 1000) ? 1000 : G.numVertex;
  for(int i = 0; i < printVertices; i++)
    fprintf(fdump, "%d\n", sum[i]);
  fclose(fdump);
  std::cout << "result is dumped into dumpIN.txt" << '\n';

#endif
  delete[] sum;
  return 0;
}
