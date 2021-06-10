#include <iostream>
#include <stdlib.h>
#include <random>
#include <chrono>
#include <math.h>
#include <fstream>
#include <assert.h>
#include <time.h>
#include <sstream>
#include <string>
#include <cstring>
#include <unordered_map>
#include <stdint.h>
#include <xmmintrin.h>
#include <smmintrin.h>
#include <immintrin.h>
#include <cstdio>
#include <cstdlib>
#include <list>
#include <sys/timeb.h>
#include "Utils.hpp"


#include "catch.hpp"


#include "frequent_items_sketch.hpp"


#include <memory>
#include <algorithm>
#include <iterator>
#include <cmath>

#include "MurmurHash3.h"


#define CLK_PER_SEC CLOCKS_PER_SEC
#define CAIDA16_SIZE 152197439
#define CAIDA18_SIZE 175880896
#define UNIV1_SIZE 17323447

typedef unsigned long long key;
typedef double val;

using namespace std;



void getKeysAndWeightsFromFile(string filename, vector<key*> &keys, vector<val*> &value, int size) {
    ifstream stream;
    stream.open(filename, fstream::in | fstream::out | fstream::app);
    if (!stream) {
        throw invalid_argument("Could not open " + filename + " for reading.");
    }

    key* file_keys = (key*) malloc(sizeof(key) * size);
    val* file_ws = (val*) malloc(sizeof(val) * size);

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
        file_ws[i] = stod(len);
        } catch (const std::invalid_argument& ia) {
        cerr << "Invalid argument: " << ia.what() << " at line " << i << endl;
        cerr << len << " " << id << endl;;
        --i;
        exit(1);
        }
    }
    keys.push_back(file_keys);
    value.push_back(file_ws);

    stream.close();
}


namespace datasketches {
    




    static int a()
    {
        
        
        vector<key*> keys;
        vector<val*> values;
        vector<int> sizes;
        vector<std::string> datasets;
        

        getKeysAndWeightsFromFile("../../datasets/UNIV1/mergedAggregatedPktlen_Srcip", keys, values, UNIV1_SIZE);
        sizes.push_back(UNIV1_SIZE);
        datasets.push_back("univ1");
        


        getKeysAndWeightsFromFile("../../datasets/CAIDA16/mergedAggregatedPktlen_Srcip", keys, values, CAIDA16_SIZE);
        sizes.push_back(CAIDA16_SIZE);
        datasets.push_back("caida");



        getKeysAndWeightsFromFile("../../datasets/CAIDA18/mergedAggregatedPktlen_Srcip", keys, values, CAIDA18_SIZE);
        sizes.push_back(CAIDA18_SIZE);
        datasets.push_back("caida18");

        
        
        
        int k = 3;
        int skSize = 23;
        ofstream ostream;
        setupOutputFile("time.raw_res", ostream, false);
        ofstream ostream1;
        setupOutputFile("nrmse.raw_res", ostream1, false);

        
        for (int run = 0; run < k; run++) {
            vector<key*>::iterator k_it = keys.begin();
            vector<val*>::iterator v_it = values.begin();
            vector<int>::iterator s_it = sizes.begin();
            vector<std::string>::iterator d_it = datasets.begin();
            
            for (int trc = 0; trc < 3; trc++) {
                key* kk = *k_it;
                val* vv = *v_it;
                int size = *s_it;
                std::string dataset = *d_it;
                
                struct timeb begintb, endtb;
                clock_t begint, endt;
                
                //construct
                frequent_items_sketch<key> sk(skSize);
                
                begint = clock();
                ftime(&begintb);
                for (int i = 0; i < size; ++i) {
                    //insert
                    sk.update(kk[i], vv[i]);
                }
                endt = clock();
                ftime(&endtb);
                double time = ((double)(endt-begint))/CLK_PER_SEC;

                
                    
//                 double c =0.0; 
//                 uint64_t vol =0;
//                 unordered_map<key, val> map;
//                 frequent_items_sketch<key> sk1(skSize);
//                 for (int i = 0; i < size; ++i) {
//                         map[kk[i]] +=vv[i];
//                         sk1.update(kk[i], vv[i]);
//                         uint32_t err = map[kk[1]] - ((double)sk1.get_estimate(kk[i]));
//                         c+= err*err;
//                         vol+=vv[i];
//                 }
//                 double mse = c/size;
//                 double rmse = sqrt(mse);
//                 double nrmse = rmse/vol;
//                 
                
                

                
                
                ostream << "Fmeg    " << dataset<<   "     " << pow(2,skSize) <<  "    " << time << endl;
//                 ostream1 << "Fmeg   " << dataset<<  "     " <<  pow(2,skSize)  << "     " << nrmse << endl;
                
                
                ++k_it;
                ++v_it;
                ++s_it;
                ++d_it;
                
            }
            
        }
        ostream.close();
        ostream1.close();
        return 0;   
        
    }

}


int main(int argc, char* argv[])
{
    datasketches::a();
    
	return 0;
}
