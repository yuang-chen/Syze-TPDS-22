#pragma once

#include <boost/multi_array.hpp>
#include <boost/smart_ptr.hpp>
#include <sys/mman.h>
#include <vector>
#include <random>
#include "global.hpp"

template<class T>
using AoV = std::vector<T>**;

template<class T>
T** allocArray2d(size_t row, size_t col) {
    T** array2d;
    array2d = new T* [row];
   // #pragma omp parallel for
    for (size_t i=0; i<row; i++)
        array2d[i] = new T [col]();
    return array2d;
}

template <class T>
void deleteArray2d (T** array2d, size_t row)
{
    
    for (size_t i=0; i<row; i++)
        delete[] array2d[i];
    delete[] array2d;
}

template<class T>
T*** allocArray2dPtr(size_t row, size_t col) {
    T*** array2dptr;
    array2dptr = new T** [row];
    for (size_t i=0; i<row; i++)
        array2dptr[i] = new T* [col];
    return array2dptr;
}

template <class T>
void deleteArray3d (T*** array3d, size_t row, size_t col)
{
    for (size_t i=0; i<row; i++)
      for (size_t j=0; j<col; j++)
        delete[] array3d[i][j];
    delete[] array3d;
}

// 2-d array of vector -> simulate 3d array
template<class T> 
std::vector<T>** alloc2AoV(size_t row, size_t col, size_t volume) {
    std::vector<T>** array2dvec;
    array2dvec = new std::vector<T>* [row];
  //  #pragma omp parallel for
    for (size_t i=0; i<row; i++) {
        array2dvec[i] = new std::vector<T>[col];
        for(size_t j = 0; j < col; j++) {
          array2dvec[i][j].reserve(volume);
      }
    }
    return array2dvec;
}

template<class T>
void delete2AoV(std::vector<T>** array2dvec, size_t row, size_t col) {
    if(array2dvec==nullptr) return;
    for (size_t i=0; i<row; i++) {
      delete[] array2dvec[i];
    }
    delete[] array2dvec;
}



template <class ET>
inline bool CAS(ET *ptr, ET oldv, ET newv) {
  if (sizeof(ET) == 1) {
    return __sync_bool_compare_and_swap((bool*)ptr, *((bool*)&oldv), *((bool*)&newv));
  } else if (sizeof(ET) == 4) {
    return __sync_bool_compare_and_swap((int*)ptr, *((int*)&oldv), *((int*)&newv));
  } else if (sizeof(ET) == 8) {
    return __sync_bool_compare_and_swap((long*)ptr, *((long*)&oldv), *((long*)&newv));
  }
  else {
    std::cout << "CAS bad length : " << sizeof(ET) << std::endl;
    abort();
  }
}

template <class ET>
inline void atomicAdd(ET *a, ET b) {
  volatile ET newV, oldV;
  do {oldV = *a; newV = oldV + b;}
  while (!CAS(a, oldV, newV));
}


bool CheckCommon(std::vector<long> const& inVectorA, std::vector<long> const& inVectorB)
{
    return std::find_first_of (inVectorA.begin(), inVectorA.end(),
                               inVectorB.begin(), inVectorB.end()) != inVectorA.end();
}



template<class Iterator>
void checkRepeatition(const Iterator begin, const Iterator end) {
  auto size = std::distance(begin, end); //auto  = typename std::iterator_traits<Iterator>::difference_type 
  using type = typename std::iterator_traits<Iterator>::value_type;
  std::vector<type> dup(size);
  std::copy(begin, end, dup.begin());
  std::sort(dup.begin(), dup.end() );
  dup.erase(std::unique( dup.begin(), dup.end() ), dup.end() );
  std::cout<<"repeated elements: " << size - dup.size() <<'\n';
}


template<class Iterator>
static void fill_random(Iterator begin, Iterator end){
    using type = typename std::iterator_traits<Iterator>::value_type;
    static std::uniform_real_distribution<type> distribution(0, 1); // value ranging from 0 to 1
    static std::default_random_engine generator;
    std::generate(begin, end, []() { return distribution(generator); });
}