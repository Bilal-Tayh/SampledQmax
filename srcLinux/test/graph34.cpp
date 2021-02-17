#include <cstdio>
#include <cstdlib>
#include <list>
#include <fstream>
#include <iostream>

#include <stdio.h>
#include <stdlib.h>
#include <sys/timeb.h>
#include "../SampledQMax.hpp"
#include "../SampledQMaxV2.hpp"
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


void benchmark_sqmax(int q, double gamma, int** data, ofstream& ostream) {
    int* elements = *data;
    struct timeb begintb, endtb;
    clock_t begint, endt;
    double time;
    SampledQMax qmax(q, gamma);
    begint = clock();
    ftime(&begintb);
    for (int i = 0; i < N; i++) {
        qmax.insert(elements[i]);
    }
    endt = clock();
    ftime(&endtb);
    time = ((double)(endt - begint)) / CLK_PER_SEC;
    ostream << "random,AmortizedSampledQMax," << N << "," << q << "," << gamma << "," << time << endl;
}

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
void benchmark_sqmaxV2(int** data, ofstream& ostream) {
    int* elements = *data;
    struct timeb begintb, endtb;
    clock_t begint, endt;
    double time;
    SampledQMaxV2<q,actualSize> qmax;
    begint = clock();
    ftime(&begintb);
    for (int i = 0; i < N; i++) {
        qmax.insert(elements[i]);
    }
    endt = clock();
    ftime(&endtb);
    time = ((double)(endt - begint)) / CLK_PER_SEC;
    ostream << "random,AmortizedSampledQMaxV2U1," << N << "," << q << "," << (double)actualSize/q-1<< "," << time << endl;
}




int main() {
  ofstream ostream;
  setupOutputFile("../results/timing_random.raw_res", ostream, false);
  for (int run = 0; run < 3; run++) {
    int* data = (int*) malloc(sizeof(int) * N);
    for (int i = 0; i< N; ++i){
      data[i] = std::rand();
    }
//     list<int> qs = {10000, 100000, 1000000, 10000000};
//     int q = 10000;
//     int q = 100000;
//     int q = 1000000;
//     int q = 10000000;
//     for (int q : qs) {
      //benchmark_heap(q, &data, ostream);
      //benchmark_skiplist(q, &data, ostream);
//       list<double> gammas = {0.005,0.01,0.05,0.1,0.25,0.5,1,2,4};
      
//       list<double> gammas = { 0.005};
//       list<double> gammas = { 0.01};
//       list<double> gammas = { 0.05};
//       list<double> gammas = { 0.1};
//       list<double> gammas = { 0.25};
//       list<double> gammas = { 0.5};
//       list<double> gammas = { 1};
//       list<double> gammas = { 2};
//       list<double> gammas = { 4};
//       for (double g : gammas) {
//         benchmark_qmax(q, g, &data, ostream);
//         benchmark_sqmax(q, g, &data, ostream);
//         benchmark_sqmaxV2<10000000,(int) (10000000*(1+0.005))>(&data, ostream);
//         benchmark_sqmaxV2<10000000,(int) (10000000*(1+0.01))>(&data, ostream);
//         benchmark_sqmaxV2<10000000,(int) (10000000*(1+0.05))>(&data, ostream);
//         benchmark_sqmaxV2<10000000,(int) (10000000*(1+0.1))>(&data, ostream);
//         benchmark_sqmaxV2<10000000,(int) (10000000*(1+0.25))>(&data, ostream);
//         benchmark_sqmaxV2<10000000,(int) (10000000*(1+0.5))>(&data, ostream);
//         benchmark_sqmaxV2<10000000,(int) (10000000*(1+1))>(&data, ostream);
//         benchmark_sqmaxV2<10000000,(int) (10000000*(1+2))>(&data, ostream);
        benchmark_sqmaxV2<10000000,(int) (10000000*(1+4))>(&data, ostream);
//       }
//     }
  }
  ostream.close();
  return 0;
}












/*



#include <cstdio>
#include <cstdlib>
#include <list>
#include <fstream>
#include <iostream>

#include <stdio.h>
#include <stdlib.h>
#include <sys/timeb.h>

#include "../Qmax.hpp"
#include "../Heap.hpp"
#include "../Skiplist.hpp"
#include "../SampledQMax.hpp"
#include "../SampledQMaxV2.hpp"
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

double totalQmaxTime = 0;
double totalSampledQmaxTime = 0;
double totalSampledQmaxV2Time = 0;

void benchmark_qmax(int q, double gamma, int** data, ofstream& ostream) {
    int* elements = *data;
    struct timeb begintb, endtb;
    clock_t begint, endt;
    double time;
    QMax qmax(q, gamma);
    begint = clock();
    ftime(&begintb);
    for (int i = 0; i < N; i++) {
        qmax.insert((int)elements[i]);
    }
    endt = clock();
    ftime(&endtb);
    time = ((double)(endt - begint)) / CLK_PER_SEC;
    totalQmaxTime += time;
    ostream << "random,AmortizedQMax," << N << "," << q << "," << gamma << "," << time << endl;
}



void benchmark_sqmax(int q, double gamma, int** data, ofstream& ostream) {
    int* elements = *data;
    struct timeb begintb, endtb;
    clock_t begint, endt;
    double time;
    SampledQMax qmax(q, gamma);
    begint = clock();
    ftime(&begintb);
    for (int i = 0; i < N; i++) {
        qmax.insert(elements[i]);
    }
    endt = clock();
    ftime(&endtb);
    time = ((double)(endt - begint)) / CLK_PER_SEC;
    totalSampledQmaxTime += time;
    ostream << "random,AmortizedSampledQMax," << N << "," << q << "," << gamma << "," << time << endl;
}

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
void benchmark_sqmaxV2(int** data, ofstream& ostream) {
    int* elements = *data;
    struct timeb begintb, endtb;
    clock_t begint, endt;
    double time;
    SampledQMaxV2<q,actualSize> qmax;
    begint = clock();
    ftime(&begintb);
    for (int i = 0; i < N; i++) {
        qmax.insert(elements[i]);
    }
    endt = clock();
    ftime(&endtb);
    time = ((double)(endt - begint)) / CLK_PER_SEC;
    totalSampledQmaxV2Time += time;
    ostream << "random,AmortizedSampledQMaxV2," << N << "," << q << "," << (double)actualSize/q-1<< "," << time << endl;
}


int main() {
  ofstream ostream;
  setupOutputFile("timing_random.raw_res", ostream, false);
  int* data = (int*)malloc(sizeof(int) * N);
  if (false) {
      for (int i = 0; i < N; ++i) {
          data[i] = std::rand();
      }
      if (false) {
          auto myfile = std::fstream("file.binary", std::ios::out | std::ios::binary);
          myfile.write((char*)&data[0], N * sizeof(int));
          myfile.close();
      }
  }
  else {
      fstream file;
      file.open("file.binary", ios::in | ios::binary);
      file.read((char*)(data), N * sizeof(int));
      file.close();
  }
  for (int run = 0; run < 1; run++) {
    

    //list<int> qs = {10000, 100000, 1000000, 10000000};
    const int qs[] = { 100000};
    for (int q : qs) {
      //benchmark_heap(q, &data, ostream);
      //benchmark_skiplist(q, &data, ostream);
      //list<double> gammas = {/*0.005,0.01,0.05,0.1,0.25,0.5,0.999,1.999,3.999};
      const list<double> gammas = {0.5};
      for (double g : gammas) {
          benchmark_qmax(q, g, &data, ostream);
//           benchmark_sqmax(q, g, &data, ostream);
          //myfunction<2>();
//           benchmark_sqmaxV2<100000,(int) (100000*(1+0.5))>(&data, ostream);
          //benchmark_sqmaxV2<20, (int)(20 * (1 + 0.5))>(&data, ostream);
      }
    }
  }
  cout << "totalSampledQmaxV2Time:" << totalSampledQmaxV2Time << " " << "totalSampledQmaxTime:" << totalSampledQmaxTime << " " << "totalQmaxTime:" << totalQmaxTime << " " << endl;
  ostream.close();
  return 0;
}*/

