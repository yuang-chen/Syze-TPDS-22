#pragma once

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <immintrin.h>
#include <assert.h>
#include <unordered_map>
#include <algorithm>
#include <iomanip>
#include <numeric>
#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include <atomic>
#include <mutex>
#include <memory>
#include <deque>
#include <thread>
#include <pthread.h>
#include <fstream>
#include <boost/timer/timer.hpp>
#include <list>

//#define HUGE_E
//#define WEIGHTED
// sizeof(VerxID) <= sizeof(EdgeID)
#ifdef HUGE_V
using EdgeId = uint64_t;
#else
using VerxId = uint32_t;
#endif

#if defined (HUGE_E) || defined (HUGE_V)
using EdgeId = uint64_t;
#else
using EdgeId = uint32_t;
#endif

//typedef uint64_t EdgeId;
using Attr_t = float;
//typedef double Attr_t;
using bool_byte = uint8_t;


enum DType {
  in=0,
  out=1
};

enum RAlgo {
    original = 0,
    pop = 1,
    bihub = 2,
    trihub = 3,
    fbc = 4,
    dbg = 5
};

std::string DType_Name[2] = {"InDegree", "OutDegree"};
std::string RAlgo_Name[6] = {"Origin", "PoP", "BiHub", "TirHub", "FBC", "DBG"};
std::string Bool_Str[2] = {"No", "Yes"};


namespace MASK {
    const unsigned MAX_NEG = 0x80000000;
    const unsigned MAX_POS = 0x7fffffff;
    const unsigned MAX_UINT = 0xffffffff;
    const unsigned MSB = 31;
    const unsigned PAGESIZE = 1 << 12;

};


namespace params {
    unsigned subgraph_size = (1024 * 1024)/sizeof(unsigned);; 
    unsigned threads = 20;
    unsigned iters = 20;
    unsigned rounds = 5;
    unsigned root_vertex = MASK::MAX_UINT;
    bool dynamic = true;
    bool output = false;
    bool verify = false;
    std::string input_file;
};


class SubMat {
public:
    ///////////////
    //! basic submatrix info.
    ///////////////
    VerxId id = 0;
    VerxId start = 0;
    VerxId stop = 0;
    SubMat(){id=0; start=0; stop=0;}
};

struct Empty { };

template <typename EdgeData>
struct EdgeUnit {
  VerxId src;
  VerxId dst;
  EdgeData edge_data;
} __attribute__((packed));

template <>
struct EdgeUnit <Empty> {
  VerxId src;
  union {
    VerxId dst;
    Empty edge_data;
  };
} __attribute__((packed));