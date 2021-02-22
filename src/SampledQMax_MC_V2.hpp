#pragma once
#ifndef S_QMAXM_H_V2
#define S_QMAXM_H_V2
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

template<int q, int _actualsize>
class SampledQMax_MC_V2
{
	float* _A;
	int _curIdx;
	int _actualsizeMinusOne;
	int _qMinusOne;
	float _gamma;
	int _nminusq;
    __m256 phi_x8;
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
public:
	void reset();
	int findKthLargestAndPivot();
	SampledQMax_MC_V2();
	void insert(int id);
	float* largestQ();
	void print();
	int checkPivot(int value);
	int findValueIndex(int value);
	int test();
	uint32_t mm256_extract_epi32_var_indx(int i);
    void get8samples();
    void replaceOneOf8(int v);
    void updateZK();
};

#include "SampledQMax_MC_V2.hpp"
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

template<int q, int _actualsize>
void SampledQMax_MC_V2<q, _actualsize>::print() {
    for (int i = 0; i < _actualsize; ++i)
        std::cout << _A[i] << " ";
    std::cout << std::endl;
    std::cout << "_phi  = " << _phi << std::endl;
}


template<int q, int _actualsize>
SampledQMax_MC_V2<q, _actualsize>::SampledQMax_MC_V2() {
    _actualsizeMinusOne = _actualsize - 1;
    _curIdx = 0;
    _A = (float*)malloc(sizeof(float) * _actualsize);
    if (!_A)
        exit(1);
    _gamma = (double)_actualsize / q - 1;
    //_q = q;
    _qMinusOne = q - 1;
    _nminusq = _actualsize - q;
    _phi = 0;
    _delta = 1-0.001;
    _alpha = 0.83;
    _psi = 2.0 / 3.0;
    _k = ceil(((_alpha * _gamma * (2 + _gamma - _alpha * _gamma)) / (pow(_gamma - _alpha * _gamma, 2))) * log(1 / _delta));
    _Z = (int)((_k * (1 + _gamma)) / (_alpha * _gamma));
    if (_Z & 0b111) // integer multiple of 8
        _Z = _Z - (_Z & 0b111) + 8;
    _k = _Z * (_alpha * _gamma) / (1 + _gamma)+0.5;
    
//     printf("%d\n",_k);
    
    gen_arr();
    //rand_bits = _mm256_set_epi64x(gen_arr(), gen_arr(), gen_arr(), gen_arr());
    counter = 0;
    initFastMod(_actualsize);
    _correctness_threshold = int(_gamma * q);
}



template<int q, int _actualsize>
void SampledQMax_MC_V2<q, _actualsize>::insert(int v) {
    if (v <= _phi){
		return;
	}
	else{
        
        while(true){
            
            // if there is less than 8 elemnts that we need to check we run on them one by one
            if(_curIdx+7>=_actualsize){
                for(int j=_curIdx;j<_actualsize;j++){
                    if(_A[j]<=_phi){
                        _A[j] = float(v);
                        _curIdx = j+1;
                        return;
                    }
                }
                maintenance();   
            }
            
            
            const __m256 item2 = _mm256_loadu_ps(_A+_curIdx);
            
            
            const __m256 match1 = _mm256_cmp_ps(item2,phi_x8, _CMP_LE_OS);
            
            const int mask1 = _mm256_movemask_ps(match1);
            
            // if one of the 8 elemnts we checked are lower than phi we replace it
            if(mask1!=0){
                replaceOneOf8(v);
                return;
            }
            else{
                _curIdx+=8;
            }
        }
	}

}


// finding which one of the 8 elemnts (starting from _curIdx) is lower than phi
// this function to be called only in case we know for sure that one of the 8 elemnts from _curIdx are lower than phi.
template<int q, int _actualsize>
void SampledQMax_MC_V2<q, _actualsize>::replaceOneOf8(int v){
    for(int j=_curIdx;j<_curIdx+8;j++){
        if(_A[j]<=_phi){
            _A[j] = float(v);
            _curIdx=j+1;
            break;
        }
    }
    
}




template<int q, int _actualsize>
void SampledQMax_MC_V2<q, _actualsize>::maintenance() {
    _phi = findKthLargestAndPivot();
    phi_x8 = _mm256_set1_ps((float)_phi);
    _curIdx = 0;
}

template<int q, int _actualsize>
float* SampledQMax_MC_V2<q, _actualsize>::largestQ() {
    maintenance();
    return _A + _nminusq;
}


// inline void swap(int& x, int& y) {
//     int z = x;
//     x = y;
//     y = z;
// }

template<int q, int _actualsize>
int SampledQMax_MC_V2<q, _actualsize>::PartitionAroundPivotValue(int left, int right, int pivot_val, int* nums) {
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
            
//             swap(nums[curRight], nums[curLeft++]);
            {
                curLeft++;
                int k = nums[curRight];
                nums[curRight] = nums[curLeft];
                nums[curLeft] = k;
            }
            
            
            if (nums[curLeft-1] == pivot_val) {
                valLoc = curLeft - 1;
            }
        }
    }
    
//     swap(nums[valLoc], nums[curRight]);
    {
        int k = nums[valLoc];
        nums[curRight] = nums[valLoc];
        nums[valLoc] = k;
    }
    
    return curRight;
}

template<int q, int _actualsize>
int SampledQMax_MC_V2<q, _actualsize>::PartitionAroundPivot(int left, int right, int pivot_idx, int* nums) {
    int pivot_value = nums[pivot_idx];
    int new_pivot_idx = right;
    
//     swap(nums[pivot_idx], nums[right]);
    {
        int k = nums[pivot_idx];
        nums[pivot_idx] = nums[right];
        nums[right] = k;
    }
    
    for (int i = right - 1; i >= left; --i) {
        if (nums[i] > pivot_value) {
            
//             swap(nums[i], nums[--new_pivot_idx]);
            {
                new_pivot_idx --;
                int k = nums[new_pivot_idx];
                nums[new_pivot_idx] = nums[i];
                nums[i] = k;
            }
        }
    }
    
//     swap(nums[right], nums[new_pivot_idx]);
    {
        int k = nums[new_pivot_idx];
        nums[new_pivot_idx] = nums[right];
        nums[right] = k;
    }
    
    return new_pivot_idx;
}


template<int q, int _actualsize>
void SampledQMax_MC_V2<q, _actualsize>::get8samples() {
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
template<int q, int _actualsize>
int SampledQMax_MC_V2<q, _actualsize>::findValueIndex(int value) {
    for (int i = 0; i < _actualsize; i++) {
        if (_A[i] == value) {
            return i;
        }
    }
    return -1;
}


// check if the conditions holds for the possible pivot "value"
// return the pivot index in _A if it hold otherwise return -1
template<int q, int _actualsize>
int SampledQMax_MC_V2<q, _actualsize>::checkPivot(int value) {
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
        //if (new_pivot_idx >= int(_gamma * q * _psi)) {
            return new_pivot_idx;
        //}
    }
    return -1*new_pivot_idx;
}



template<int q, int _actualsize>
int SampledQMax_MC_V2<q, _actualsize>::findKthLargestAndPivot() {

        // p should contain the smallest _k values from the _Z sampled random values from _A
        std::priority_queue <float> p;
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
        /*
        for (int i = 0; i < _k; i++) {
            int j = GenerateRandom();
            p.push(_A[j]);
        }*/
        float top = p.top();

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
    updateZK();
    return top;
        
	
}


template<int q, int _actualsize>
void SampledQMax_MC_V2<q, _actualsize>::updateZK(){
    _delta=_delta/2.0;
    _k=ceil( ((_alpha*_gamma*(2+_gamma - _alpha*_gamma)) / (pow(_gamma -_alpha*_gamma,2))) * log(1/_delta));
    _Z = (int)  ( (_k*(1+_gamma)) / (_alpha*_gamma) );
}

template<int q, int _actualsize>
void SampledQMax_MC_V2<q, _actualsize>::reset() {
    _phi = 0;
    _curIdx = 0;
}

#endif
