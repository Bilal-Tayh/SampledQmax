#include <cstdio>
#include <cstdlib>
#include <list>
#include <fstream>
#include <iostream>

#include <stdio.h>
#include <stdlib.h>
#include <sys/timeb.h>
#include "../SampledQMax_LV.hpp"
#include "../SampledQMax_LV_V2.hpp"
// #include "../SampledQMax_MC.hpp"
// #include "../SampledQMax_MC_V2.hpp"
#include "../Qmax.hpp"
#include "../Heap.hpp"
#include "../Skiplist.hpp"
#include "Utils.hpp"

#define N 150000000
#define CLK_PER_SEC CLOCKS_PER_SEC

using namespace std;

void benchmark_skiplist(int q, int** data, ofstream &ostream) {
  int *elements = *data;
  struct timeb begintb, endtb;
  clock_t begint, endt;
  double time;
  Skiplist sl(q);
  begint = clock();
  ftime(&begintb);
  for (int i = 0; i < N; i++) {
         sl.add(elements[i]);
  }
  endt = clock();
  ftime(&endtb);
  time = ((double)(endt-begint))/CLK_PER_SEC;
  ostream << "random,SkipList," << N << "," << q << ",," << time << endl;
}

void benchmark_heap(int q, int** data, ofstream &ostream) {
  int *elements = *data;
  struct timeb begintb, endtb;
  clock_t begint, endt;
  double time;
  Heap heap(q);
  begint = clock();
  ftime(&begintb);
  for (int i = 0; i < N; i++) {
    heap.add(elements[i]);
  }
  endt = clock();
  ftime(&endtb);
  time = ((double)(endt-begint))/CLK_PER_SEC;
  ostream << "random,Heap," << N << "," << q << ",," << time << endl;
}

void benchmark_qmax(int q, double gamma, int** data, ofstream &ostream) {
  int* elements = *data;
  struct timeb begintb, endtb;
  clock_t begint, endt;
  double time;
  QMax qmax(q, gamma);
  begint = clock();
  ftime(&begintb);
  for (int i = 0; i < N; i++) {
    qmax.insert(elements[i]);
  }
  endt = clock();
  ftime(&endtb);
  time = ((double)(endt-begint))/CLK_PER_SEC;
  ostream << "random,AmortizedQMax," << N << "," << q << "," << gamma << "," << time << endl;
}


void benchmark_qmax_LV(int q, double gamma, int** data, ofstream& ostream) {
    int* elements = *data;
    struct timeb begintb, endtb;
    clock_t begint, endt;
    double time;
    SampledQMax_LV qmax = SampledQMax_LV(q, gamma);
    begint = clock();
    ftime(&begintb);
    for (int i = 0; i < N; i++) {
        qmax.insert(elements[i]);
    }
    endt = clock();
    ftime(&endtb);
    time = ((double)(endt - begint)) / CLK_PER_SEC;
    ostream << "random,AmortizedSampledQMax_LV," << N << "," << q << "," << gamma << "," << time << endl;
}

// void benchmark_qmax_MC(int q, double gamma, int** data, ofstream& ostream) {
//     int* elements = *data;
//     struct timeb begintb, endtb;
//     clock_t begint, endt;
//     double time;
//     SampledQMax_MC qmax = SampledQMax_MC (q, gamma);
//     begint = clock();
//     ftime(&begintb);
//     for (int i = 0; i < N; i++) {
//         qmax.insert(elements[i]);
//     }
//     endt = clock();
//     ftime(&endtb);
//     time = ((double)(endt - begint)) / CLK_PER_SEC;
//     ostream << "random,AmortizedSampledQMax_MC," << N << "," << q << "," << gamma << "," << time << endl;
// }

template<int param> void myfunction()
{
    switch (param)
    {
    case 1:
        cout << "asd";
        break;

    case 2:
        cout << "asd2";
        break;

    case 100:                 // 100 is taken as example
        break;
    }
}



template<int q, int actualSize>
void benchmark_qmax_LV_V2(int** data, ofstream& ostream) {
    int* elements = *data;
    struct timeb begintb, endtb;
    clock_t begint, endt;
    double time;
    SampledQMax_LV_V2<q,actualSize> qmax;
    begint = clock();
    ftime(&begintb);
    for (int i = 0; i < N; i++) {
        qmax.insert(elements[i]);
    }
    endt = clock();
    ftime(&endtb);
    time = ((double)(endt - begint)) / CLK_PER_SEC;
    ostream << "random,AmortizedSampledQMax_LV_V2U1," << N << "," << q << "," << (double)actualSize/q-1<< "," << time << endl;
}

/*
template<int q, int actualSize>
void benchmark_qmax_MC_V2(int** data, ofstream& ostream) {
    int* elements = *data;
    struct timeb begintb, endtb;
    clock_t begint, endt;
    double time;
    SampledQMax_MC_V2<q,actualSize> qmax;
    begint = clock();
    ftime(&begintb);
    for (int i = 0; i < N; i++) {
        qmax.insert(elements[i]);
    }
    endt = clock();
    ftime(&endtb);
    time = ((double)(endt - begint)) / CLK_PER_SEC;
    ostream << "random,AmortizedSampledQMax_MC_V2U1," << N << "," << q << "," << (double)actualSize/q-1<< "," << time << endl;
}
*/




int main() {
  ofstream ostream;
  setupOutputFile("../results/timing_random.raw_res", ostream, false);
  for (int run = 0; run < 3; run++) {
    int* data = (int*) malloc(sizeof(int) * N);
    for (int i = 0; i< N; ++i){
      data[i] = std::rand();
    }

      list<double> gammas = {/*0.005,0.01,0.05,0.1,0.25,0.5,1,2,*/4};
      list<unsigned int> qs = {/*10000, 100000, 1000000,*/ 10000000};
      for (unsigned q: qs) {
    
        for (double g : gammas) {
            benchmark_qmax(q, g, &data, ostream);
            benchmark_qmax_LV(q, g, &data, ostream);
        }
      }
        
        benchmark_qmax_LV_V2<10000000,(int) (10000000*(1+4))>(&data, ostream);

  }
  ostream.close();
  return 0;
}
