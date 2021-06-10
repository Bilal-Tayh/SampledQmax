

#include "catch.hpp"

#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>
#include <random>
#include <sys/timeb.h>
#include "frequent_items_sketch.hpp"
#include <iostream>
#include <stdlib.h>
#include <chrono>
#include <math.h>
#include <assert.h>
#include <time.h>
#include <cstring>
#include <unordered_map>
#include <stdint.h>
#include <xmmintrin.h>
#include <smmintrin.h>
#include <immintrin.h>
#include <cstdio>
#include <cstdlib>
#include <list>
#include "Utils.hpp"

#define CLK_PER_SEC CLOCKS_PER_SEC

using namespace std;

namespace datasketches {
    static frequent_items_sketch<int> create_unweighted_sketch(uint64_t N) {
        struct timeb begintb, endtb;
        
        double time = 0;
            
        frequent_items_sketch<int> sk(14);
        
        
        for (int run = 0; run < 5; run++) {
            int* data = (int*) malloc(sizeof(int) * N);
            for (int i = 0; i< N; ++i){
            data[i] = std::rand();
            }
            
            
            clock_t begint, endt;
            begint = clock();
            ftime(&begintb);
            for (uint64_t i = 0; i < N; ++i) {
                sk.update(i,data[i]);
            }
            endt = clock();
            ftime(&endtb);
            time += ((double)(endt-begint))/CLK_PER_SEC;
            
            
        }
        cout << time/5.0 << endl;
        return sk;
    }

    static void a()
    {
            struct timeb begintb, endtb;
            clock_t begint, endt;
            double time;
            frequent_items_sketch<int> a = create_unweighted_sketch( 150000000);
            
//             cout << a.get_num_samples() ;
    }

}

int main(int argc, char* argv[])
{
    datasketches::a();
    
	return 0;
}
