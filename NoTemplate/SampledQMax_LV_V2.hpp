#pragma once
#ifndef S_QMAXLV_H_V2
#define S_QMAXLV_H_V2
#include <fstream>
#include <cstdio>
#include <stdio.h>
#include <assert.h>
#include <iostream>
#include <time.h>
#include <immintrin.h>
#include <iomanip> 
#include "RngFast.hpp"
#include <queue>

__m256i m1;                                            // multiplier used in fast division
__m256i d1;
uint32_t shf1;
uint32_t shf2;

union U256 {
    int32_t i[8];
    __m256i ymm;
};


void initFastMod(uint32_t d) {                                 // Set or change divisor, calculate parameters
    uint32_t L, L2, sh1, sh2, m;
    switch (d) {
    case 0:
        m = sh1 = sh2 = 1 / d;                         // provoke error for d = 0
        break;
    case 1:
        m = 1; sh1 = sh2 = 0;                          // parameters for d = 1
        break;
    case 2:
        m = 1; sh1 = 1; sh2 = 0;                       // parameters for d = 2
        break;
    default:                                           // general case for d > 2
        L = ceil(log2(d));//bit_scan_reverse(d - 1) + 1;                  // ceil(log2(d))
        L2 = L < 32 ? 1 << L : 0;                      // 2^L, overflow to 0 if L = 32
        m = 1 + uint32_t((uint64_t(L2 - d) << 32) / d); // multiplier
        sh1 = 1;  sh2 = L - 1;                         // shift counts
    }
    m1 = _mm256_set1_epi32(m);
    d1 = _mm256_set1_epi32(d);
    shf1 = sh1;
    shf2 = sh2;
}


class SampledQMax_LV_V2
{
    int _actualsize;
    int q;
	int* _A;
	int _curIdx;
	int _actualsizeMinusOne;
	int _qMinusOne;
	float _gamma;
	int _nminusq;
	int _phi;
	std::default_random_engine _generator;
	void maintenance();
    int PartitionAroundPivotValue(int left, int right, int pivot_val, int* nums);
    int PartitionAroundPivot(int left, int right, int pivot_idx, int* nums);
	double _delta;
	double _alpha;
	double _psi;
	int _k;
	int _Z;
	rng::rng128 gen_arr;
	int counter;
    U256 eight_samples;
	std::priority_queue <int, std::vector<int>, std::greater<int>> testp;
	int testSize;
    int _correctness_threshold;
    int _min_maintenance_free;
public:
	void reset();
	int findKthLargestAndPivot();
	SampledQMax_LV_V2(int a, int b);
	void insert(int id);
	int* largestQ();
	void print();
	int checkPivot(int value);
	int findValueIndex(int value);
	int test();
	int findReplaceValueIndex(int value);
	uint32_t mm256_extract_epi32_var_indx(int i);
    void get8samples();
};

#include "SampledQMax_LV_V2.hpp"
#include <iostream>
#include <random>
#include <math.h> 

#include <vector>

#include <bitset>
#include <cstring>
#include <algorithm>
#include <memory.h>
#include <stdlib.h>
#include <chrono>
#include <fstream>
#include <assert.h>
#include <time.h>
#include <immintrin.h>
#include <iomanip> 


void SampledQMax_LV_V2::print() {
    for (int i = 0; i < _actualsize; ++i)
        std::cout << _A[i] << " ";
    std::cout << std::endl;
    std::cout << "_phi  = " << _phi << std::endl;
}



SampledQMax_LV_V2::SampledQMax_LV_V2(int q1, int actualsize) {
    int q=q1;
    int _actualsize = actualsize;
    //_actualsize = q * (1 + gamma);
    _actualsizeMinusOne = _actualsize - 1;
    _curIdx = _actualsize;
    _A = (int*)malloc(sizeof(int) * _actualsize);
    if (!_A)
        exit(1);
    _gamma = (double)_actualsize / q - 1;
    //_q = q;
    _qMinusOne = q - 1;
    _nminusq = _actualsize - q;
    _phi = -1;
    _delta =0.6;
    _alpha = 0.83;
    _psi = 1.0 / 5.0;
    _k = ceil(((_alpha * _gamma * (2 + _gamma - _alpha * _gamma)) / (pow(_gamma - _alpha * _gamma, 2))) * log(1 / _delta));
    _Z = (int)((_k * (1 + _gamma)) / (_alpha * _gamma));
    if (_Z & 0b111) // integer multiple of 8
        _Z = _Z - (_Z & 0b111) + 8;
    _k = _Z * (_alpha * _gamma) / (1 + _gamma)+0.5;
    
    
    if(q <= 100000){
        _delta =0.99;
        _k = ceil(((_alpha * _gamma * (2 + _gamma - _alpha * _gamma)) / (pow(_gamma - _alpha * _gamma, 2))) * log(1 / _delta));
        _Z = (int)((_k * (1 + _gamma)) / (_alpha * _gamma));
        if (_Z & 0b111) // integer multiple of 8
            _Z = _Z - (_Z & 0b111) + 8;
        _k = _Z * (_alpha * _gamma) / (1 + _gamma)+0.5;
    }
    
//     printf("%d\n",_k);
    
  
    
    gen_arr();
    //rand_bits = _mm256_set_epi64x(gen_arr(), gen_arr(), gen_arr(), gen_arr());
    counter = 0;
    initFastMod(_actualsize);
    _correctness_threshold = int(_gamma * q);
    _min_maintenance_free = int(_gamma * q * _psi);
}


 
void SampledQMax_LV_V2::insert(int v) {
    if (v < _phi) {
        return;
    }
    else {
        _A[--_curIdx] = v;
        if (_curIdx) {
            return;
        }
        else {
            maintenance();
        }
    }
}
 
void SampledQMax_LV_V2::maintenance() {
    _phi = findKthLargestAndPivot();
}

 
int* SampledQMax_LV_V2::largestQ() {
    maintenance();
    return _A + _nminusq;
}

inline void swap(int& x, int& y) {
    int z = x;
    x = y;
    y = z;
}

 
int SampledQMax_LV_V2::PartitionAroundPivotValue(int left, int right, int pivot_val, int* nums) {
    //int pivot_value = nums[pivot_idx];
    //int new_pivot_idx = right;
    
    //bool valFound = false;
    int valLoc = -1;
    int curLeft = left;
    int curRight = right;
    while (curRight >= curLeft) {
        if (nums[curRight] > pivot_val) {
            --curRight;
        }
        else {
            swap(nums[curRight], nums[curLeft++]);
            if (nums[curLeft-1] == pivot_val) {
                valLoc = curLeft - 1;
            }
        }
    }
    //swap(nums[pivot_idx], nums[right]);

    /*
    while (valFound == false){
        if (nums[curRight] > pivot_val) {
            --curRight;
        }
        else if (nums[curRight] < pivot_val) {
            swap(nums[curRight], nums[curLeft++]);
        }
        else {
            valFound = true;
            if (curRight < right){
                swap(nums[curRight--], nums[right]); //put pivot at [right]
            }
            break;
        }
    }
    while (curRight >= curLeft) {
        if (nums[curRight] > pivot_val) {
            --curRight;
        }
        else {
            swap(nums[curRight], nums[curLeft++]);
        }
    }*/
    swap(nums[valLoc], nums[curRight]);
    return curRight;
}

 
int SampledQMax_LV_V2::PartitionAroundPivot(int left, int right, int pivot_idx, int* nums) {
    int pivot_value = nums[pivot_idx];
    int new_pivot_idx = right;
    swap(nums[pivot_idx], nums[right]);
    for (int i = right - 1; i >= left; --i) {
        if (nums[i] > pivot_value) {
            swap(nums[i], nums[--new_pivot_idx]);
        }
    }
    swap(nums[right], nums[new_pivot_idx]);
    return new_pivot_idx;
}

/*
 
uint32_t SampledQMax_LV_V2::mm256_extract_epi32_var_indx(int i) {
    __m128i indx = _mm_cvtsi32_si128(i);
    __m256i val = _mm256_permute2f128_si256(rand_bits, _mm256_castsi128_si256(indx), 0);
    return _mm_cvtsi128_si32(_mm256_castsi256_si128(val));
}
 
int SampledQMax_LV_V2::GenerateRandom() {
    //     return std::rand() % max; 
    int indx = 0;
    if (counter >= 8) {
        rand_bits = _mm256_set_epi64x(gen_arr(), gen_arr(), gen_arr(), gen_arr());
        counter = 0;
    }
    indx = mm256_extract_epi32_var_indx(counter) % _actualsize;
    counter++;
    return indx;
}
*/
 
void SampledQMax_LV_V2::get8samples() {
    int indx = 0;
    eight_samples.ymm = _mm256_set_epi64x(gen_arr(), gen_arr(), gen_arr(), gen_arr());
    __m256i t1 = _mm256_mul_epu32(eight_samples.ymm, m1);                           // 32x32->64 bit unsigned multiplication of even elements of a
    __m256i t2 = _mm256_srli_epi64(t1, 32);                         // high dword of even numbered results
    __m256i t3 = _mm256_srli_epi64(eight_samples.ymm, 32);                          // get odd elements of a into position for multiplication
    __m256i t4 = _mm256_mul_epu32(t3, m1);                          // 32x32->64 bit unsigned multiplication of odd elements
    __m256i t5 = _mm256_set_epi32(-1, 0, -1, 0, -1, 0, -1, 0);      // mask for odd elements
    __m256i t7 = _mm256_blendv_epi8(t2, t4, t5);                    // blend two results
    __m256i t8 = _mm256_sub_epi32(eight_samples.ymm, t7);                           // subtract
    __m256i t9 = _mm256_srli_epi32(t8, shf1);                       // shift right logical
    __m256i t10 = _mm256_add_epi32(t7, t9);                         // add
    __m256i t11 = _mm256_srli_epi32(t10, shf2);                     // shift right logical 
    __m256i t12 = _mm256_mullo_epi32(t11, d1);                      // multiply quotient with divisor
    eight_samples.ymm = _mm256_sub_epi32(eight_samples.ymm, t12);                   // subtract
}



 
// if value dont exist in _A return -1, otherwise return the index of the value in _A
int SampledQMax_LV_V2::findValueIndex(int value) {
    for (int i = 0; i < _actualsize; i++) {
        if (_A[i] == value) {
            return i;
        }
    }
    return -1;
}


// check if the conditions holds for the possible pivot "value"
// return the pivot index in _A if it hold otherwise return -1
 
int SampledQMax_LV_V2::checkPivot(int value) {
    int left = 0, right = _actualsizeMinusOne;
    int new_pivot_idx;

#if true
       int pivot_idx = findValueIndex(value);
    
       if (pivot_idx == -1) {
           return -1;
       }
       
       new_pivot_idx = PartitionAroundPivot(left, right, pivot_idx, _A);
#else
        new_pivot_idx = PartitionAroundPivotValue(left, right, value, _A);
#endif

    

    if (new_pivot_idx <= _correctness_threshold) {
//         if (new_pivot_idx >= _min_maintenance_free) {
            return new_pivot_idx;
//         }
//         else{
// //                 printf("2:  %d\n",new_pivot_idx);
//                 return -1;
//         }
    }
//     printf("1:  %d\n",new_pivot_idx);
    return -1*new_pivot_idx;
}



 
int SampledQMax_LV_V2::findKthLargestAndPivot() {
    int tries = 2;
    while (tries != 0) {

        // p should contain the smallest _k values from the _Z sampled random values from _A
        std::priority_queue <int> p;
        int left_to_fill = _k;
        while (left_to_fill > 0) {
            get8samples();
            for (int i = 0; i < 8; ++i) {
                p.push(_A[eight_samples.i[i]]);
            }
            left_to_fill -= 8;
        }

        int left_to_sample = _Z-_k+ left_to_fill;
        while (left_to_fill++ < 0) {
            p.pop();
        }
        
        int top = p.top();

        while (left_to_sample > 0) {
            get8samples();
            for (int i = 0; i < 8; ++i) {
                if (top > _A[eight_samples.i[i]]) {
                    p.pop();
                    p.push(_A[eight_samples.i[i]]);
                    top = p.top();
                }
            }
            left_to_sample -= 8;
        }

        // check if the conditions holds for the possible pivot B[Kth_minimal_idx]
        int idx = checkPivot(top);
        if (idx > 0) {
            _curIdx = idx;
            return _A[idx];
        }
        else if(q>100000 and idx < -1){
            int l = -idx - _correctness_threshold;
//             printf("l = %d\n",l);
            

            int pop = ( (_k-1) - _k * ( (q*_gamma*_alpha)/(l+(q*_gamma)) ) );
            
//             printf("per = %f\n", (_k-1) -_k *( (q*_gamma*_alpha)/(l+(q*_gamma)) ));
//             printf("pop = %d\n",pop);
            for(int i=0;i<pop;i++){
                p.pop();
            }
            top = p.top();
            idx = checkPivot(top);
//             printf("newidx = %d\n",idx);
            if(idx < 0){
//                 int s;
//                 printf("wrong\n");
//                 scanf("wrong %d\n",s);
            }
            else{
                _curIdx = idx;
                return top;
            }
        }
            
        tries--;
    }

    // sampling didn't work twice so we are doing the old SampledQMax_LV_V2 maintanance
    int left = 0, right = _actualsizeMinusOne;
    while (left <= right) {
        int pivot_idx = left;
        int new_pivot_idx = PartitionAroundPivot(left, right, pivot_idx, _A);
        if (new_pivot_idx == _nminusq) {
            _curIdx = new_pivot_idx;
            return _A[new_pivot_idx];
        }
        else if (new_pivot_idx > _nminusq) {
            right = new_pivot_idx - 1;
        }
        else {  // new_pivot_idx < _q - 1.
            left = new_pivot_idx + 1;
        }
    }

}
 
void SampledQMax_LV_V2::reset() {
    _phi = -1;
    _curIdx = _actualsize;
}

#endif
