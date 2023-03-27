/**
 * Author: Kartik Lakhotia
           Sourav Pati
 * Email id: klakhoti@usc.edu
             spati@usc.edu
 * Date: 27-Feb-2018
 *
 * This code implements work optimized propagation blocking with
 * transposed bin graph to reduce cache misses in scatter
 */

#include <boost/timer/timer.hpp>
#include <float.h>

#include <cmath>
#include <iostream>
#include <numeric>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <vector>

#include "../include/gas.h"
#include <algorithm>
#include <assert.h>
#include <omp.h>
#include <string.h>

intV         binWidth      = (256 * 1024) / sizeof(float);    // 512kB
unsigned int binOffsetBits = ( unsigned int )std::log2(( double )binWidth);
intV         NUM_BINS      = 10000000 / binWidth;

bool     countDense = false;
unsigned iteration  = 0;

// #undef DENSE

//////////////////////////////////////////
// performance monitoring via PCM
//////////////////////////////////////////
#define PERF_MON
#undef PERF_MON

//////////////////////////////////////////
// level 2 debugging - asserts enabled
//////////////////////////////////////////
#define DEBUGL2
#undef DEBUGL2

#define ITERTIME
#undef ITERTIME

int          NUM_THREADS = 20;
unsigned int MAX_ITER    = 20;
bool         dynamic     = false;
char*        input       = NULL;

template < class graph >
void initialize(graph* G, int argc, char** argv)
{
  int count = 0;
  G->start  = 1;
  G->rounds = 3;
  for(int i = 1; i < argc; i++)
  {
    if(i + 1 != argc)
    {
      if(strcmp(argv[i], "-d") == 0)    // This is your parameter name
      {
        input = argv[i + 1];    // The next value in the array is your value
        i++;                    // Move to the next flag
      }
      if(strcmp(argv[i], "-v") == 0)    // This is your parameter name
      {
        G->start = ( intV )atoi(argv[i + 1]);    // The next value in the array is your value
        i++;                                     // Move to the next flag
      }
      if(strcmp(argv[i], "-s") == 0)    // This is your parameter name
      {
        binWidth = (( unsigned int )atoi(argv[i + 1]) * 1024)
                   / sizeof(float);    // The next value in the array is your value
        count = i;
        i++;    // Move to the next flag
      }
      if(strcmp(argv[i], "-i") == 0)    // This is your parameter name
      {
        MAX_ITER =
            ( unsigned int )atoi(argv[i + 1]);    // The next value in the array is your value
        i++;                                      // Move to the next flag
      }
      if(strcmp(argv[i], "-r") == 0)    // This is your parameter name
      {
        G->rounds =
            ( unsigned int )atoi(argv[i + 1]);    // The next value in the array is your value
        i++;                                      // Move to the next flag
      }
      if(strcmp(argv[i], "-t") == 0)    // This is your parameter name
      {
        NUM_THREADS =
            ( unsigned int )atoi(argv[i + 1]);    // The next value in the array is your value
        i++;                                      // Move to the next flag
      }
      if(strcmp(argv[i], "-y") == 0)    // This is your parameter name
      {
        dynamic =
            ( unsigned int )atoi(argv[i + 1]) == 1;    // The next value in the array is your value
        i++;                                           // Move to the next flag
      }
    }
  }
  if(argc < 2)
  {
    printf("Usage : %s -d <filename> -v <start vertex>-s <part size> -i "
           "<#iterations> -t <threads> -r <rounds> -y <dynamic or not> \n",
           argv[0]);
    exit(1);
  }

  omp_set_num_threads(NUM_THREADS);

  //////////////////////////////////////////
  // read csr file
  //////////////////////////////////////////
  if(read_mix(input, G) == -1)
  {
    printf("couldn't read %s\n", argv[1]);
    exit(1);
  }

  // if (strcmp(argv[count], "-s") != 0) // This is your parameter name
  {
    intV numVerticesPerBin = (G->numVertex / (NUM_THREADS * 4));
    numVerticesPerBin      = (numVerticesPerBin < binWidth) ? numVerticesPerBin : binWidth;
    intV pow2              = 1;
    while(pow2 <= numVerticesPerBin)
      pow2 *= 2;
    pow2 /= 2;
    if(pow2 == 0)
      binWidth = 4;
    else
      binWidth = pow2;
  }

  NUM_BINS = (G->numVertex - 1) / binWidth + 1;
  // G->numBins = NUM_BINS;
  printf("root vertex %d, threads %d, max iters %d, rounds %d, is dynamic %d\n", G->start,
         NUM_THREADS, MAX_ITER, G->rounds, dynamic);

  printf("number of partitions %d, size of partitions %ld KB\n", NUM_BINS,
         ( unsigned )binWidth * sizeof(intV) / 1024);
  binOffsetBits = ( unsigned int )std::log2(( double )binWidth);
  //////////////////////////////////////////
  // initialize graph frontier, degree etc.//
  //////////////////////////////////////////
  initGraph(G);
}

template < class type, class graph >
void initBin(graph* G)
{
  dynamicSplit(G, dynamic);
  //////////////////////////////////////////
  // static work allocation to threads
  // equal no. of edges to all bins
  //////////////////////////////////////////
  partition(G->TD, G);

  printf("partitioning successful\n");

  //////////////////////////////////////////////////
  // compute storage space required for each bin and
  // offsets for storage in bins for a partition
  // 1 column -> 1 gather bin; 1 row -> 1 scatter bin
  // bin[i][j] -> stores what i sends to j
  //////////////////////////////////////////////////
  G->updateBinAddrSize = allocateBinMat< intE >(NUM_BINS, NUM_BINS);
  G->destIdBinAddrSize = allocateBinMat< intE >(NUM_BINS, NUM_BINS);
  G->binFlag           = allocateBinMat< bool >(NUM_BINS, NUM_BINS);
  G->activeBins        = allocateBinMat< intV >(NUM_BINS, NUM_BINS);

  //    struct timespec preStart, preEnd;
  //    float preTime;
  //    if( clock_gettime(CLOCK_REALTIME, &preStart) == -1) { perror("clock
  //    gettime");}

//////////////////////////////////////////
//// transpose and compute offsets ///////
//////////////////////////////////////////
#pragma omp parallel for schedule(dynamic, 1)
  for(intV i = 0; i < NUM_BINS; i++)
    transposePartition(G, &(G->TD[i]), G->updateBinAddrSize[i], G->destIdBinAddrSize[i]);

  printf("PNG construction successful\n");

  //    if( clock_gettime( CLOCK_REALTIME, &preEnd) == -1 ) { perror("clock
  //    gettime");} preTime = (preEnd.tv_sec - preStart.tv_sec)+
  //    (int)(preEnd.tv_nsec - preStart.tv_nsec)/1e9; printf("%s, preprocessing
  //    time - %lf\n", argv[1], preTime);

  //////////////////////////////////////////
  //////////////// BINNING ////////////////
  //////////////////////////////////////////

  //////////////////////////////////////////
  ////individual bins to->fro each partition //////
  //////////////////////////////////////////
  G->indUpdateBins    = allocateBinMatPtr< type >(NUM_BINS, NUM_BINS);
  G->indDestIdBins    = allocateBinMatPtr< intV >(NUM_BINS, NUM_BINS);
  G->sparseDestIdBins = allocateBinMatPtr< intV >(NUM_BINS, NUM_BINS);
#ifdef WEIGHTED
  G->indWeightBins = allocateBinMatPtr< unsigned int >(NUM_BINS, NUM_BINS);
#endif
#pragma omp parallel for num_threads(NUM_THREADS) schedule(dynamic, 1)
  for(intV i = 0; i < NUM_BINS; i++)
  {
    for(intV j = 0; j < NUM_BINS; j++)
    {
      G->indUpdateBins[i][j]    = new type[G->destIdBinAddrSize[i][j]];
      G->indDestIdBins[i][j]    = new intV[G->destIdBinAddrSize[i][j]];
      G->sparseDestIdBins[i][j] = new intV[G->destIdBinAddrSize[i][j]];
#ifdef WEIGHTED
      G->indWeightBins[i][j] = new unsigned int[G->destIdBinAddrSize[i][j]];
#endif
    }
  }

  // pointers for each (i,j) bin for later use //
  G->updateBinPointers = allocateBinMat< intE >(NUM_BINS, NUM_BINS);
  G->destIdBinPointers = allocateBinMat< intE >(NUM_BINS, NUM_BINS);
  struct timespec preStart, preEnd;
  float           preTime;
  if(clock_gettime(CLOCK_REALTIME, &preStart) == -1)
  {
    perror("clock gettime");
  }
#pragma omp parallel for num_threads(NUM_THREADS) schedule(dynamic, 1)
  for(intV i = 0; i < NUM_BINS; i++)
  {
#ifdef WEIGHTED
    writeDestIds(G, &G->TD[i], G->indDestIdBins[i], G->indWeightBins[i], G->destIdBinPointers[i]);
#else
    writeDestIds(G, &G->TD[i], G->indDestIdBins[i], G->destIdBinPointers[i]);
#endif
  }
  if(clock_gettime(CLOCK_REALTIME, &preEnd) == -1)
  {
    perror("clock gettime");
  }
  preTime = (preEnd.tv_sec - preStart.tv_sec) + ( int )(preEnd.tv_nsec - preStart.tv_nsec) / 1e9;
  printf("writing destid time - %lf\n", preTime);

#ifdef DEBUG
  printf("binning complete\n");
#endif

  //////////////////////////////////////////
  //////////// BINNING COMPLETE ////////////
  //////////////////////////////////////////
  unsigned total_edges = 0;
  for(intV i = 0; i < NUM_BINS; i++)
  {
    for(intV j = 0; j < NUM_BINS; j++)
    {
      total_edges += G->updateBinAddrSize[i][j];
    }
  }
  std::cout << "total_edges " << total_edges << '\n';
}

template < class graph >
void dynamicSplit(graph* G, bool is_dynamic)
{
  unsigned average_deg = G->numEdges / NUM_BINS * 2;
  unsigned old_num_sub = NUM_BINS;
  unsigned new_num_sub = NUM_BINS;

  // auto ver_deg = G->outDeg;
  std::vector< unsigned > sub_deg(old_num_sub);
  std::vector< unsigned > sub_div(old_num_sub, 1);
  std::vector< unsigned > sub_size(old_num_sub, binWidth);

  unsigned overflow_count  = 0;
  unsigned underflow_count = 0;
  std::cout << "dynamic splitting threshold: " << average_deg / (G->numEdges / NUM_BINS) << '\n';
#pragma omp parallel for schedule(static) num_threads(NUM_THREADS) reduction(+ : NUM_BINS)
  for(unsigned n = 0; n < old_num_sub; n++)
  {
    unsigned start = n * binWidth;
    unsigned end   = (n + 1) * binWidth > G->numVertex ? G->numVertex : (n + 1) * binWidth;

    unsigned curr_sub, prev_sub;
    sub_deg[n] = G->VI[end] - G->VI[start];
    if(sub_deg[n] > average_deg && is_dynamic)
    {
      sub_div[n]  = pow(2, round(log2(sub_deg[n] / average_deg)));
      sub_size[n] = sub_size[n] / sub_div[n];
      NUM_BINS += sub_div[n] - 1;
    }
  }
  G->numBins = NUM_BINS;
  // NUM_BINS = new_num_sub;
  auto                    x = sub_div.size();
  std::vector< unsigned > div_index;
  std::vector< unsigned > div_value;

  // for(int i = 0; i < x; i++) {
  //     if(sub_div[i]!=1) {
  //         div_index.push_back(i);
  //         div_value.push_back(sub_div[i]);
  //         std::cout <<  div_index.back() << " " << div_value.back() << "\n";

  //     }
  // }
  // std::cout << (float) div_value.size()/x << '\n' ;
  // std::cout << "the largest element: " << max_element(sub_div.begin(),
  // sub_div.end()) << '\n';
  G->TD = ( partitionData* )malloc(sizeof(partitionData) * NUM_BINS);

  std::vector< unsigned > size(NUM_BINS, binWidth);
  G->sub_offset = std::vector< unsigned >(old_num_sub, binOffsetBits);
  G->sub_map    = std::vector< unsigned >(old_num_sub);

#pragma omp parallel for schedule(static) num_threads(NUM_THREADS)
  for(auto n = 0; n < old_num_sub; n++)
  {
    auto n_map    = std::accumulate(sub_div.begin(), sub_div.begin() + n, 0);
    G->sub_map[n] = n_map;
    G->sub_offset[n] -= log2(sub_div[n]);
    for(unsigned i = 0; i < sub_div[n]; i++)
    {
      size[i + n_map]              = size[i + n_map] / sub_div[n];
      G->TD[i + n_map].startVertex = n * binWidth + i * size[n_map];
      G->TD[i + n_map].endVertex   = std::min(n * binWidth + (i + 1) * size[n_map], G->numVertex);
    }
  }

  std::cout << "after dynamic split: " << NUM_BINS << '\n';
}

template < class type, class graph, class userArg >
void scatter_and_gather(graph* G, userArg UA)
{
#ifdef ITERTIME
  float           time;
  struct timespec scatterStart, scatterEnd, gatherStart, gatherEnd, start, end;
  float           scatterTime = 0.0, gatherTime = 0.0;
  struct timespec iterStart, iterEnd;
  float           iterTime;
  if(clock_gettime(CLOCK_REALTIME, &start) == -1)
  {
    perror("clock gettime");
  }
#endif
  intV numActiveBins;
  ///////////////////////////////////////
  ////Set FLAG For Scatter and Gather////
  ///////////////////////////////////////

#ifndef DENSE
  G->frontierSize = 0;
#endif
//    printf("\n");
//	for(intV i=0;i<NUM_BINS; i++)
//    {
//        for (intV j=0; j<G->TD[i].frontierSize; j++)
//            printf("%d, ", G->TD[i].frontier[j]);
//    }
//    printf("\n");
#ifdef ITERTIME
  if(clock_gettime(CLOCK_REALTIME, &scatterStart) == -1)
  {
    perror("clock gettime");
  }
#endif

  numActiveBins      = G->partListPtr;
  unsigned NumTreads = NUM_THREADS;

#ifndef DENSE
  G->partListPtr = 0;
#ifndef ASYNCH

  /*     if(countDense) {
          unsigned count = 0;
          #pragma omp parallel for num_threads(NUM_THREADS) schedule (dynamic,
     1) reduction(+:count) for (intV i=0; i<numActiveBins; i++) {
              if(G->TD[G->activeScatter[i]].isDense)
                  count++;
          }
          if(count > numActiveBins/2) {
              NumTreads = 20;
          } else {
              NumTreads = 36;
          }
          if(count == 0 & firstTime)
      } else { */
#pragma omp parallel for num_threads(NUM_THREADS) schedule(dynamic, 1)
  for(intV i = 0; i < numActiveBins; i++)
  {
    densityCheck(&G->TD[G->activeScatter[i]]);
  }
  // }

#endif
#endif

#ifdef TEST
  // if cond(!dense & !asyn & test)
  unsigned count = 0;
#pragma omp parallel for num_threads(NUM_THREADS) schedule(dynamic, 1) reduction(+ : count)
  for(intV i = 0; i < numActiveBins; i++)
  {
    if(G->TD[G->activeScatter[i]].isDense)
      count++;
  }
  std::cout << numActiveBins << " " << count << "\n";
#endif

#pragma omp parallel for schedule(dynamic, 1) num_threads(NumTreads)
  for(intV ptr = 0; ptr < numActiveBins; ptr++)
  {
    intV i = G->activeScatter[ptr];
#ifdef ASYNCH
    sgMix(G, &G->TD[i], G->indUpdateBins, G->indDestIdBins, G->sparseDestIdBins, G->TD,
          G->destIdBinAddrSize, G->destIdBinPointers, G->updateBinPointers, G->scatterDone, UA);
#else
    scatter< type >(G, &G->TD[i], G->indUpdateBins[i], G->sparseDestIdBins[i],
                    G->updateBinPointers[i], G->destIdBinPointers[i], UA);
#endif
  }

#pragma omp parallel for
  for(intV ptr = 0; ptr < numActiveBins; ptr++)
    G->scatterDone[G->activeScatter[ptr]] = false;    // reset scatter done status of partitions

#ifdef ITERTIME

  if(clock_gettime(CLOCK_REALTIME, &scatterEnd) == -1)
  {
    perror("clock gettime");
  }

  scatterTime += (scatterEnd.tv_sec - scatterStart.tv_sec)
                 + ( float )(scatterEnd.tv_nsec - scatterStart.tv_nsec) / 1e9;

  if(clock_gettime(CLOCK_REALTIME, &gatherStart) == -1)
  {
    perror("clock gettime");
  }
#endif
#pragma omp parallel for schedule(dynamic, 1) num_threads(NumTreads)
  for(intV ptr = 0; ptr < G->partListPtr; ptr++)
  {
    intV i = G->activeGather[ptr];
    gather< type >(G, &G->TD[i], G->indUpdateBins, G->indDestIdBins, G->sparseDestIdBins, G->TD,
                   G->destIdBinAddrSize, G->destIdBinPointers, G->updateBinPointers, UA);
    G->activeScatter[ptr] = i;
  }
#ifdef ITERTIME
  if(clock_gettime(CLOCK_REALTIME, &gatherEnd) == -1)
  {
    perror("clock gettime");
  }
  gatherTime += (gatherEnd.tv_sec - gatherStart.tv_sec)
                + ( float )(gatherEnd.tv_nsec - gatherStart.tv_nsec) / 1e9;
#endif
  //  std::cout << "number of active subgraphs GATHER: " << G->partListPtr <<
  //  "\n";

#ifdef ITERTIME

  if(clock_gettime(CLOCK_REALTIME, &iterEnd) == -1)
  {
    perror("clock gettime");
  }
  iterTime =
      (iterEnd.tv_sec - iterStart.tv_sec) + ( int )(iterEnd.tv_nsec - iterStart.tv_nsec) / 1e9;
  printf("scatter time= %lf gather time = %lf, total time = %lf\n", scatterTime, gatherTime,
         scatterTime + gatherTime);
#endif
  // std::cout << G->frontierSize << '\n';

  // free allocated memory//
  //    freeMem(&G);
  //    free(TD);
  //    freeMat<intE>(updateBinAddrSize, NUM_BINS+1);
  //    freeMat<intE>(destIdBinAddrSize, NUM_BINS+1);
  //    freeMat<intE>(updateBinPointers, NUM_BINS);
  //    freeMat<intE>(destIdBinPointers, NUM_BINS);
  //    freeMatPtr<type>(indUpdateBins, NUM_BINS,NUM_BINS);
  //    freeMatPtr<intV>(indDestIdBins, NUM_BINS,NUM_BINS);
  // #ifdef WEIGHTED
  //    freeMatPtr<unsigned int>(indWeightBins, NUM_BINS, NUM_BINS);
  // #endif
}
