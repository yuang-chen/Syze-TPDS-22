/**
 * collaborative filtering (matrix factorization)
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

struct CF_F
{
  float* latent{ nullptr };    // [num_regular]
  float* error{ nullptr };     // [num_regular]
  int    dim{ 1 };             // dim = 1
  float  lambda{ 0 };
  float  step{ 0 };

  CF_F(float* _latent, float* _error, int _dim, float _lambda, float _step)
      : latent(_latent), error(_error), dim(_dim), lambda(_lambda), step(_step)
  {
  }
  // ~CF_F() { delete[] latent_cache;}s
  inline float scatterFunc(intV node)
  {
    return latent[node];
  }

  inline void initFunc(intV node)
  {
    error[node] = 0.0;
  }

  inline void gatherFunc(float updateVal, intV vertex)
  {
    const float estimate = latent[vertex] * updateVal;
    const float err      = 1.0 - estimate;
    error[vertex] += updateVal * err;
  }

  inline void filterFunc(intV node)
  {
    latent[node] += step * (-lambda * latent[node] + error[node]);
  }
  void reset(intV num_verx)
  {
    std::fill(latent, latent + num_verx, 0.5);
    std::fill(error, error + num_verx, 0);
  }
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

  const int   dim    = 1;    // assume dim = 1
  const float lambda = 0.001;
  const float step   = 0.00000035;
  float*      latent = new float[G.numVertex * dim]();
  float*      error  = new float[G.numVertex * dim]();

  struct timespec start, end;
  float           time;
  float           sum_time = 0;

  CF_F cf_func(latent, error, dim, lambda, step);

  int ctr = 0;
  while(ctr < G.rounds)
  {
    numIter = 0;

    cf_func.reset(G.numVertex);

    if(clock_gettime(CLOCK_REALTIME, &start) == -1)
    {
      perror("clock gettime");
    }

    while(numIter < MAX_ITER)
    {
      scatter_and_gather<float>(&G, cf_func);
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

  mergeSortWOkey<float>(latent, 0, G.numVertex - 1);
  FILE* fdump = fopen("dumpCF.txt", "w");
  if(fdump == nullptr)
  {
    fputs("file error\n", stderr);
    exit(1);
  }
  int printVertices = (G.numVertex > 1000) ? 1000 : G.numVertex;
  for(int i = 0; i < printVertices; i++)
    fprintf(fdump, "%lf\n", latent[i]);
  fclose(fdump);
#endif

  delete[] latent;
  delete[] error;
  return 0;
}
