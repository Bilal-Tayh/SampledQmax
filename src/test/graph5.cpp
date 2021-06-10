#include <cstdio>
#include <cstdlib>
#include <list>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cstring>
#include <random>

#include <stdio.h>
#include <stdlib.h>
#include <sys/timeb.h>

#include "../QmaxKV.hpp"
#include "../HeapKV.hpp"
#include "../SkiplistKV.hpp"
#include "Utils.hpp"
#include "../SampledQMaxKV_LV_V2.hpp"

#define CLK_PER_SEC CLOCKS_PER_SEC
#define CAIDA16_SIZE 152197437
#define CAIDA18_SIZE 175880808
#define UNIV1_SIZE 17323447

using namespace std;

void benchmark_psskiplist(int q, key** keys, val** vals, ofstream &ostream, string dataset, int numKeys) {
  std::random_device _rd;
  std::mt19937 _e2(_rd());
  std::uniform_real_distribution<double> _dist(0,1);
  key *elements = *keys;
  val *weights = *vals;
  struct timeb begintb, endtb;
  clock_t begint, endt;
  double time;
  SkiplistKV sl(q);
  begint = clock();
  ftime(&begintb);
  for (int i = 0; i < numKeys; i++) {
    val priority = weights[i] / (1-_dist(_e2));
    sl.add(pair<key, val>(elements[i], priority));
  }
  endt = clock();
  ftime(&endtb);
  time = ((double)(endt-begint))/CLK_PER_SEC;
  ostream << dataset << ",SkipList," << numKeys << "," << q << ",," << time << endl;
}

void benchmark_psheap(int q, key** keys, val** vals, ofstream &ostream, string dataset, int numKeys) {
  std::random_device _rd;
  std::mt19937 _e2(_rd());
  std::uniform_real_distribution<double> _dist(0,1);
  key *elements = *keys;
  val *weights = *vals;
  struct timeb begintb, endtb;
  clock_t begint, endt;
  double time;
  HeapKV heap(q);
  begint = clock();
  ftime(&begintb);
  for (int i = 0; i < numKeys; i++) {
    val priority = weights[i] / (1-_dist(_e2));
    heap.add(pair<key,val>(elements[i], priority));
  }
  endt = clock();
  ftime(&endtb);
  time = ((double)(endt-begint))/CLK_PER_SEC;
  ostream << dataset << ",Heap," << numKeys << "," << q << "," << time << endl;
}

void benchmark_psqmax(int q, double gamma, key** keys, val** vals, ofstream &ostream, string dataset, int numKeys) {
  std::random_device _rd;
  std::mt19937 _e2(_rd());
  std::uniform_real_distribution<double> _dist(0,1);
  key *elements = *keys;
  val *weights = *vals;
  struct timeb begintb, endtb;
  clock_t begint, endt;
  double time;
  QMaxKV qmax = QMaxKV(q, gamma);
  begint = clock();
  ftime(&begintb);
  for (int i = 0; i < numKeys; i++) {
    val priority = weights[i] / (1-_dist(_e2));
    qmax.insert(elements[i], priority);
  }
  endt = clock();
  ftime(&endtb);
  time = ((double)(endt-begint))/CLK_PER_SEC;
  ostream << dataset << ",AmortizedQMax," << numKeys << "," << q << "," << gamma << "," << time << endl;
}


template<int q, int actualSize>
void benchmark_psqmax_LV(key** keys, val** vals, ofstream &ostream, string dataset, int numKeys) {
  std::random_device _rd;
  std::mt19937 _e2(_rd());
  std::uniform_real_distribution<double> _dist(0,1);
  key *elements = *keys;
  val *weights = *vals;
  struct timeb begintb, endtb;
  clock_t begint, endt;
  double time;
  SampledQMaxKV_LV_V2<q,actualSize> lvqmax = SampledQMaxKV_LV_V2<q,actualSize>();
  begint = clock();
  ftime(&begintb);
  for (int i = 0; i < numKeys; i++) {
    val priority = weights[i] / (1-_dist(_e2));
    lvqmax.insert(elements[i], priority);
  }
  endt = clock();
  ftime(&endtb);
  time = ((double)(endt-begint))/CLK_PER_SEC;
  ostream << dataset << ",LVAmortizedQMax," << numKeys << "," << q << "," << (double)actualSize/q-1 << "," << time << endl;
}

void getKeysAndValsFromFile(string filename, vector<key*> &keys, vector<val*> &vals, int size) {
  ifstream stream;
  stream.open(filename, fstream::in | fstream::out | fstream::app);
  if (!stream) {
    throw invalid_argument("Could not open " + filename + " for reading.");
  }

  key* file_keys = (key*) malloc(sizeof(key) * size);
  val* file_vals = (val*) malloc(sizeof(val) * size);

  string line;
  string len;
  string id;
  for (int i = 0; i < size; ++i){
    getline(stream, line);
    std::istringstream iss(line);
    iss >> len;
    iss >> id;
    try {
      file_keys[i] = stoull(id);
      
      file_vals[i] = stoull(len);
    } catch (const std::invalid_argument& ia) {
      cerr << "Invalid argument: " << ia.what() << " at line " << i << endl;
      cerr << len << " " << id << endl;;
      --i;
      exit(1);
    }
  }
  keys.push_back(file_keys);
  vals.push_back(file_vals);

  stream.close();
}

int main() {
  vector<ofstream*> streams;
  vector<key*> keys;
  vector<val*> vals;
  vector<int> sizes;
  vector<string> datasets;

  ofstream univ1stream;
  setupOutputFile("../results/ps_univ1.raw_res", univ1stream, false);
  streams.push_back(&univ1stream);
  getKeysAndValsFromFile("../datasets/UNIV1/mergedPktlen_Srcip", keys, vals, UNIV1_SIZE);
  sizes.push_back(UNIV1_SIZE);
  datasets.push_back("univ1");

  ofstream caida16stream;
  setupOutputFile("../results/ps_caida.raw_res", caida16stream, false);
  streams.push_back(&caida16stream);
  getKeysAndValsFromFile("../datasets/CAIDA16/mergedPktlen_Srcip", keys, vals, CAIDA16_SIZE);
  sizes.push_back(CAIDA16_SIZE);
  datasets.push_back("caida");

  ofstream caida18stream;
  setupOutputFile("../results/ps_caida18.raw_res", caida18stream, false);
  streams.push_back(&caida18stream);
  getKeysAndValsFromFile("../datasets/CAIDA18/mergedPktlen_Srcip", keys, vals, CAIDA18_SIZE);
  sizes.push_back(CAIDA18_SIZE);
  datasets.push_back("caida18");
  

  

  list<unsigned int> qs = {1000000, 10000000};
  for (int run = 0; run < 1; run++) {
      
      
      
      
      
      for (int run = 0; run < 1; run++) {
        vector<key*>::iterator k_it = keys.begin();
        vector<val*>::iterator v_it = vals.begin();
        vector<int>::iterator s_it = sizes.begin();
        vector<string>::iterator d_it = datasets.begin();
        
        for (auto& stream : streams) {
        
            key* k = *k_it;
            val* v = *v_it;
            int size = *s_it;
            string dataset = *d_it;
      
    
      

        
       
        benchmark_psqmax_LV<1000000,(int) (1000000*(1+0.25))>(&k, &v, *stream, dataset, size);

        benchmark_psqmax_LV<1000000,(int) (1000000*(1+0.1))>(&k, &v, *stream, dataset, size);

        benchmark_psqmax_LV<1000000,(int) (1000000*(1+0.05))>(&k, &v, *stream, dataset, size);

        benchmark_psqmax_LV<10000000,(int) (10000000*(1+0.25))>(&k, &v, *stream, dataset, size);

        benchmark_psqmax_LV<10000000,(int) (10000000*(1+0.1))>(&k, &v, *stream, dataset, size);

        benchmark_psqmax_LV<10000000,(int) (10000000*(1+0.05))>(&k, &v, *stream, dataset, size);
        
        
        
        
         for (unsigned int q : qs ) {
            list<double> gammas = { 0.25, 0.1, 0.05};
            for (double g : gammas) {
                benchmark_psqmax(q, g, &k, &v, *stream, dataset, size);
                
            }
        }
    
      ++k_it;
      ++v_it;
      ++s_it;
      ++d_it;
    }
  }  
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
/*      
      
  for (unsigned q: qs) {
    vector<key*>::iterator k_it = keys.begin();
    vector<val*>::iterator v_it = vals.begin();
    vector<int>::iterator s_it = sizes.begin();
    vector<string>::iterator d_it = datasets.begin();
    for (auto& stream : streams) {
      key* k = *k_it;
      val* v = *v_it;
      int size = *s_it;
      string dataset = *d_it;
      list<double> gammas = {0.25,0.1 , 0.05};
      for (double g : gammas) {
        benchmark_psqmax(q, g, &k, &v, *stream, dataset, size);
      }
      ++k_it;
      ++v_it;
      ++s_it;
      ++d_it;
    }
  }*/
  }
  
  
  
  
  univ1stream.close();
  caida16stream.close();
  caida18stream.close();
  return 0;
}

