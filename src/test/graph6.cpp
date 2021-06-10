#include <cstdio>
#include <cstdlib>
#include <list>
#include <fstream>
#include <iostream>
#include <cstring>

#include <stdio.h>
#include <stdlib.h>
#include <sys/timeb.h>

#include "../QmaxKV.hpp"
#include "../HeapKV.hpp"
#include "../SkiplistKV.hpp"
#include "../xxhash.h"
#include "Utils.hpp"
#include "../SampledQMaxKV_LV_V2.hpp"

#define CLK_PER_SEC CLOCKS_PER_SEC
#define CAIDA16_SIZE 152197439
#define CAIDA18_SIZE 175880896
#define UNIV1_SIZE 17323447
using namespace std;

void benchmark_hhskiplist(int q, key** keys, ofstream &ostream, string dataset, int numKeys) {
  key *elements = *keys;
  struct timeb begintb, endtb;
  clock_t begint, endt;
  double time;
  SkiplistKV sl(q);
  begint = clock();
  ftime(&begintb);
  for (int i = 0; i < numKeys; i++) {
    val hash = XXH64((void *) &(elements[i]), sizeof(key), 1235);
    sl.add(pair<key, val>(elements[i], hash));
  }
  endt = clock();
  ftime(&endtb);
  time = ((double)(endt-begint))/CLK_PER_SEC;
  ostream << dataset << ",SkipList," << numKeys << "," << q << ",," << time << endl;
}

void benchmark_hhheap(int q, key** keys, ofstream &ostream, string dataset, int numKeys) {
  key *elements = *keys;
  struct timeb begintb, endtb;
  clock_t begint, endt;
  double time;
  HeapKV heap(q);
  begint = clock();
  ftime(&begintb);
  for (int i = 0; i < numKeys; i++) {
    val hash = XXH64((void *) &(elements[i]), sizeof(key), 1235);
    heap.add(pair<key,val>(elements[i], hash));
  }
  endt = clock();
  ftime(&endtb);
  time = ((double)(endt-begint))/CLK_PER_SEC;
  ostream << dataset << ",Heap," << numKeys << "," << q << ",," << time << endl;
}

void benchmark_hhqmax(int q, double gamma, key** keys, ofstream &ostream, string dataset, int numKeys) {
  key *elements = *keys;
  struct timeb begintb, endtb;
  clock_t begint, endt;
  double time;
  QMaxKV qmax = QMaxKV(q, gamma);
  begint = clock();
  ftime(&begintb);
  for (int i = 0; i < numKeys; i++) {
    val hash = XXH64((void *) &(elements[i]), sizeof(key), 1235);
    qmax.insert(elements[i], hash);
  }
  endt = clock();
  ftime(&endtb);
  time = ((double)(endt-begint))/CLK_PER_SEC;
  ostream << dataset << ",LVAmortizedQMax," << numKeys << "," << q << "," << gamma << "," << time << endl;
}


template<int q, int actualSize>
void benchmark_hhqmax_LV(key** keys, ofstream &ostream, string dataset, int numKeys) {
  key *elements = *keys;
  struct timeb begintb, endtb;
  clock_t begint, endt;
  double time;
  SampledQMaxKV_LV_V2<q,actualSize> lvqmax = SampledQMaxKV_LV_V2<q,actualSize>();
  begint = clock();
  ftime(&begintb);
  for (int i = 0; i < numKeys; i++) {
    val hash = XXH64((void *) &(elements[i]), sizeof(key), 1235);
    lvqmax.insert(elements[i], hash);
  }
  endt = clock();
  ftime(&endtb);
  time = ((double)(endt-begint))/CLK_PER_SEC;
  ostream << dataset << ",AmortizedQMax," << numKeys << "," << q << "," << (double)actualSize/q-1 << "," << time << endl;
}

void getKeysFromFile(string filename, vector<key*> &keys, int size) {
  ifstream stream;
  stream.open(filename, fstream::in | fstream::out | fstream::app);
  if (!stream) {
    throw invalid_argument("Could not open " + filename + " for reading.");
  }

  key* data = (key*) malloc(sizeof(key) * size);
  string line;
  for (int i = 0; i< size; ++i){
    getline(stream, line);
    try {
      data[i] = stoull(line);
    } catch (const invalid_argument& ia) {
      cerr << "Invalid argument: " << ia.what() << " at line " << i << endl;
      cerr << line << endl;
      --i;
    }
  }

  keys.push_back(data);

  stream.close();
}

int main() {
  vector<ofstream*> streams;
  vector<key*> keys;
  vector<int> sizes;
  vector<string> datasets;

  ofstream univ1stream;
  setupOutputFile("../results/hh_univ1.xxh.raw_res", univ1stream, false);
  streams.push_back(&univ1stream);
  getKeysFromFile("../datasets/UNIV1/mergedIPID_Srcip", keys, UNIV1_SIZE);
  sizes.push_back(UNIV1_SIZE);
  datasets.push_back("univ1");

  ofstream caida16stream;
  setupOutputFile("../results/hh_caida.xxh.raw_res", caida16stream, false);
  streams.push_back(&caida16stream);
  getKeysFromFile("../datasets/CAIDA16/mergedIPID_Srcip", keys, CAIDA16_SIZE);
  sizes.push_back(CAIDA16_SIZE);
  datasets.push_back("caida");

  ofstream caida18stream;
  setupOutputFile("../results/hh_caida18.xxh.raw_res", caida18stream, false);
  streams.push_back(&caida18stream);
  getKeysFromFile("../datasets/CAIDA18/mergedIPID_Srcip", keys, CAIDA18_SIZE);
  sizes.push_back(CAIDA18_SIZE);
  datasets.push_back("caida18");

  list<unsigned int> chis = {1000000,10000000};
  for (int run = 0; run < 2; run++) {
    vector<key*>::iterator k_it = keys.begin();
    vector<int>::iterator s_it = sizes.begin();
    vector<string>::iterator d_it = datasets.begin();
    
    for (auto& stream : streams) {
      
        key* k = *k_it;
        int size = *s_it;
        string dataset = *d_it;
      
    
      

        
       
        benchmark_hhqmax_LV<1000000,(int) (1000000*(1+0.25))>(&k, *stream, dataset, size);

        benchmark_hhqmax_LV<1000000,(int) (1000000*(1+0.1))>(&k, *stream, dataset, size);

        benchmark_hhqmax_LV<1000000,(int) (1000000*(1+0.05))>(&k, *stream, dataset, size);

        benchmark_hhqmax_LV<10000000,(int) (10000000*(1+0.25))>(&k, *stream, dataset, size);

        benchmark_hhqmax_LV<10000000,(int) (10000000*(1+0.1))>(&k, *stream, dataset, size);

        benchmark_hhqmax_LV<10000000,(int) (10000000*(1+0.05))>(&k, *stream, dataset, size);
        
        
        
        
         for (unsigned int chi : chis ) {
            list<double> gammas = { 0.25, 0.1, 0.05};
            for (double g : gammas) {
                benchmark_hhqmax(chi, g, &k, *stream, dataset, size);
                
            }
        }
    
      ++k_it;
      ++s_it;
      ++d_it;
    }
  }
  univ1stream.close();
  caida16stream.close();
  caida18stream.close();
  return 0;
}
