#pragma once
#ifndef S_QMAXLVKV_H_V2
#define S_QMAXLVKV_H_V2
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
#include <bits/stdc++.h>

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


class SampledQMaxKV_LV_V2
{
    int q;
    int _actualsize;
	key* _K;
	val* _V;
	int _curIdx;
	int _actualsizeMinusOne;
	int _qMinusOne;
	float _gamma;
	int _nminusq;
	val _phi;
    key _k_phi;
	std::default_random_engine _generator;
	void maintenance();
    int PartitionAroundPivot(int left, int right, int pivot_idx, val* nums);
    void swap(int a, int b);
	double _delta;
	double _alpha;
	double _psi;
	int _k;
	int _Z;
	rng::rng128 gen_arr;
	int counter;
    U256 eight_samples;
// 	std::priority_queue <int, std::vector<int>, std::greater<int>> testp;
	int testSize;
    int _correctness_threshold;
    int _min_maintenance_free;
public:
	void reset();
	val findKthLargestAndPivot();
	SampledQMaxKV_LV_V2(int a, int b);
	void insert(key k, val v);
	outputkv largestQ();
	void print();
	int checkPivot(val value);
	int findValueIndex(val value);
    int findValueDup(val value,int idx, int check);
	int test();
// 	int findReplaceValueIndex(int value);
	uint32_t mm256_extract_epi32_var_indx(int i);
    void get8samples();
};

#include "SampledQMaxKV_LV_V2.hpp"
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


void SampledQMaxKV_LV_V2::print() {
    for (int i = 0; i < _actualsize; ++i)
        std::cout << _K[i] << ":"<< _V[i]<< " ";
    std::cout << std::endl;
    std::cout << "_phi  = " << _phi << std::endl;
}


SampledQMaxKV_LV_V2::SampledQMaxKV_LV_V2(int q1, int actualsize) {
    q=q1;
    _actualsize = actualsize;
    //_actualsize = q * (1 + gamma);
    _actualsizeMinusOne = _actualsize - 1;
    _curIdx = _actualsize;
    _K = (key*) malloc(sizeof(key) * _actualsize);
	_V = (val*) malloc(sizeof(val) * _actualsize);
    if (!_V || !_K) {
		exit(1);
	}
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



void SampledQMaxKV_LV_V2::insert(key k, val v) {
//     printf("in2\n");
    if (v < _phi) {
        return;
    }
    else {
        _V[--_curIdx] = v;
        _K[_curIdx] = k;
        if (_curIdx) {
            return;
        }
        else {
            maintenance();
        }
    }
}



void SampledQMaxKV_LV_V2::maintenance() {
    _phi = findKthLargestAndPivot();
}


outputkv SampledQMaxKV_LV_V2::largestQ() {
	outputkv out;
	maintenance();
	out.keyArr = _K + _nminusq;
	out.valArr = _V + _nminusq;
	return out;
}


void SampledQMaxKV_LV_V2::swap(int a, int b) {
	if (a==b) return;
	key k = _K[a];
	_K[a] = _K[b];
	_K[b] = k;
	val v = _V[a];
	_V[a] = _V[b];
	_V[b] = v;
}




int SampledQMaxKV_LV_V2::PartitionAroundPivot(int left, int right, int pivot_idx, val* nums) {
    val pivot_value = nums[pivot_idx];
	int new_pivot_idx = right;
	swap(pivot_idx, right);
	for (int i = right-1; i >= left; --i) {
		if (nums[i] > pivot_value) { // DIRECT COMPARSION: compare if nums[i] is bigger than pivot_value
			swap(i, --new_pivot_idx);
		}
	}
	swap(right, new_pivot_idx);
	return new_pivot_idx;
}




 
void SampledQMaxKV_LV_V2::get8samples() {
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



 
// if value dont exist in _K return -1, otherwise return the index of the value in _K
int SampledQMaxKV_LV_V2::findValueIndex(val value) {
    for (int i = 0; i < _actualsize; i++) {
        if (_V[i] == value) {
            return i;
        }
    }
    return -1;
}



 
// if value dont exist in _K return -1, otherwise return the index of the value in _K
int SampledQMaxKV_LV_V2::findValueDup(val value,int idx, int check) {
    int sum=check;
    for (int i = idx; i < _actualsize; i++) {
        if (_V[i] == value) {
            sum--;
            if(sum==0){
                return 1;
            }
        }
    }
    return 0;
}


// check if the conditions holds for the possible pivot "value"
// return the pivot index in _A if it hold otherwise return -1
 
int SampledQMaxKV_LV_V2::checkPivot(val value) {
    int left = 0, right = _actualsizeMinusOne;
    int new_pivot_idx;


    int pivot_idx = findValueIndex(value);
//     int dup = findValueDup(value ,pivot_idx);
    if (pivot_idx == -1) {
        return -1;
    }
    
    new_pivot_idx = PartitionAroundPivot(left, right, pivot_idx, _V);

    
    
    

    if (new_pivot_idx <= _correctness_threshold) {
//         if (new_pivot_idx >= _min_maintenance_free) {
            return new_pivot_idx;
//         }
//         else{
// //                 printf("2:  %d\n",new_pivot_idx);
//                 return -1;
//         }
    }
    else if ( findValueDup(value ,0, new_pivot_idx - _correctness_threshold)==1  ) {
//         printf("dup %d\n", dup);
        return _correctness_threshold;
    }
    
    
//     printf("1:  %d\n",new_pivot_idx);
    return -1*new_pivot_idx;
}



 
val SampledQMaxKV_LV_V2::findKthLargestAndPivot() {
    int tries = 2;
    while (tries != 0) {

        // p should contain the smallest _k values from the _Z sampled random values from _A
        std::priority_queue <val> p;
        int left_to_fill = _k;
        while (left_to_fill > 0) {
            get8samples();
            for (int i = 0; i < 8; ++i) {
                p.push(_V[eight_samples.i[i]]);
            }
            left_to_fill -= 8;
        }

        int left_to_sample = _Z-_k+ left_to_fill;
        while (left_to_fill++ < 0) {
            p.pop();
        }
        
        val top = p.top();

        while (left_to_sample > 0) {
            get8samples();
            for (int i = 0; i < 8; ++i) {
                if (top > _V[eight_samples.i[i]]) {
                    p.pop();
                    p.push(_V[eight_samples.i[i]]);
                    top = p.top();
                }
            }
            left_to_sample -= 8;
        }
        
        
//         printf("full\n");

        // check if the conditions holds for the possible pivot B[Kth_minimal_idx]
        int idx = checkPivot(top);
//         printf("ind %d\n", idx);
        
        if (idx > 0) {
            _curIdx = idx;
            _k_phi = _K[idx];
            return _V[idx];
        }
        else if(q>100000 and idx < -1){
//             printf("pop in\n");
            int l = -idx - _correctness_threshold;
//             printf("l = %d\n",l);
            

            int pop = ( (_k-1) - _k * ( (q*_gamma*_alpha)/(l+(q*_gamma)) ) );
//              printf("poping :  %d\n",pop);
            
//             printf("per = %f\n", (_k-1) -_k *( (q*_gamma*_alpha)/(l+(q*_gamma)) ));
//             printf("pop = %d\n",pop);
            for(int i=0;i<pop;i++){
                p.pop();
            }
            top = p.top();
//             printf("pvt in\n");
            idx = checkPivot(top);
//             printf("pvt out\n");
//             printf("newidx = %d\n",idx);
            if(idx < 0){
//                 int s;
//                 printf("wrong\n");
//                 scanf("wrong %d\n",s);
            }
            else{
                _curIdx = idx;
                _k_phi = _K[idx];
                return _V[idx];
            }
//             printf("pop out\n");
        }
        
        tries--;
//         printf("tries : %d\n", tries);
    }
    
//     printf("fail\n");

    // sampling didn't work twice so we are doing the old SampledQMaxKV_LV_V2 maintanance
    int left = 0, right = _actualsizeMinusOne;
	while (left <= right) {
	
//         printf("left %d, right %d\n", left, right);
        
		int pivot_idx = q;
		int new_pivot_idx = PartitionAroundPivot(left, right, pivot_idx, _V);
		if (new_pivot_idx == _nminusq) {
			_k_phi = _K[new_pivot_idx];
			return _V[new_pivot_idx];
		} else if (new_pivot_idx > _nminusq) {
			right = new_pivot_idx - 1;
		} else {  // new_pivot_idx < _q - 1.
			left = new_pivot_idx + 1;
		}
	}
// 	printf("org out\n");
}



 
void SampledQMaxKV_LV_V2::reset() {
    _k_phi = -1;
    _phi = -1;
    _curIdx = _actualsize;
}

#endif
