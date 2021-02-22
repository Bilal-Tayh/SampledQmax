#include "SampledQMax_LV.hpp"
#include <iostream>
#include <random>
#include <math.h> 

#include <vector>

#include <bitset>
//#include "unistd.h"
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
#include "RngFast.hpp"

void SampledQMax_LV::print(){
	for (int i = 0; i < _actualsize; ++i)
		std::cout << _A[i] << " ";
	std::cout << std::endl;
	std::cout << "_phi  = " << _phi << std::endl;
}



SampledQMax_LV::SampledQMax_LV(int q, float gamma){
	_actualsize = q * (1+gamma);
	_actualsizeMinusOne = _actualsize - 1;
	_curIdx = _actualsize;
	_A = (int*) malloc(sizeof(int) * _actualsize);
	if (!_A)
		exit(1);
	_gamma = gamma;
	_q = q;
	_qMinusOne = q - 1;
	_nminusq = _actualsize - q;
	_phi = -1;
    _delta = 0.001;
    _alpha=0.83;
    _psi = 2.0/3.0;
    _k=ceil( ((_alpha*_gamma*(2+_gamma - _alpha*_gamma)) / (pow(_gamma -_alpha*_gamma,2))) * log(1/_delta));
    _Z = (int)  ( (_k*(1+_gamma)) / (_alpha*_gamma) );
    gen_arr();
    rand_bits = _mm256_set_epi64x(gen_arr(), gen_arr(), gen_arr(), gen_arr());
    counter=0;
}



void SampledQMax_LV::insert(int v){
	if (v < _phi){
		return;
	}
	else{
		_A[--_curIdx] = v;
		if (_curIdx){
			return;
		}
		else {
			maintenance();
		}
	}
}

void SampledQMax_LV::maintenance(){
	_phi = findKthLargestAndPivot();
}

int* SampledQMax_LV::largestQ(){
	maintenance();
	return _A + _nminusq;
}

inline void swap(int &x, int &y){
	int z = x;
	x = y;
	y = z;
}


int SampledQMax_LV::PartitionAroundPivot(int left, int right, int pivot_idx, int* nums) {
	int pivot_value = nums[pivot_idx];
	int new_pivot_idx = right;
	swap(nums[pivot_idx], nums[right]);
	for (int i = right-1; i >= left; --i) {
		if (nums[i] > pivot_value) {
			swap(nums[i], nums[--new_pivot_idx]);
		}
	}
	swap(nums[right], nums[new_pivot_idx]);
	return new_pivot_idx;
}



uint32_t SampledQMax_LV::mm256_extract_epi32_var_indx(int i){   
    __m128i indx = _mm_cvtsi32_si128(i);
    __m256i val  = _mm256_permute2f128_si256(rand_bits, _mm256_castsi128_si256(indx),0);
    return _mm_cvtsi128_si32(_mm256_castsi256_si128(val));
}





int SampledQMax_LV::GenerateRandom(int max){
//     return std::rand() % max; 
    int indx = 0;
    if(counter >= 8){
        rand_bits =_mm256_set_epi64x(gen_arr(), gen_arr(), gen_arr(), gen_arr());
        counter=0;
    }
    indx = mm256_extract_epi32_var_indx(counter)%max;
    counter++;
    return indx;
}




// if value dont exist in _A return -1, otherwise return the index of the value in _A
int SampledQMax_LV::findValueIndex(int value){
    for(int i=0;i<_actualsize;i++){
        if(_A[i]==value){
            return i;
        }
    }
    return -1;
}


// check if the conditions holds for the possible pivot "value"
// return the pivot index in _A if it hold otherwise return -1
int SampledQMax_LV::checkPivot(int value){
    int left = 0, right = _actualsizeMinusOne;
    int pivot_idx = findValueIndex(value);
    
    
    
     if(pivot_idx==-1){
        return -1;
    }
    int new_pivot_idx = PartitionAroundPivot(left, right, pivot_idx, _A);
    

    
    if (new_pivot_idx <= int(_gamma*_q)) {
        if(new_pivot_idx >= int(_gamma*_q*_psi)) {
            return new_pivot_idx;
        }
    }
    return -1;
}




int SampledQMax_LV::findKthLargestAndPivot(){
    int tries=2;
    while(tries!=0){
        
        // p should contain the smallest _k values from the _Z sampled random values from _A
       std::priority_queue <int> p;

       for(int i=0;i<_k;i++){
            int j=GenerateRandom(_actualsize);
            p.push(_A[j]);
        }
        int top = p.top();
        
        for(int i=_k;i<_Z;i++){
            int j=GenerateRandom(_actualsize);
            if(top>_A[j]){
                p.pop();
                p.push(_A[j]); 
                top = p.top();
            }
        }
    
        
        
        
        // check if the conditions holds for the possible pivot B[Kth_minimal_idx]
        int idx = checkPivot(top);
        if(idx!=-1){
            _curIdx = idx;
            return _A[idx];
        }

        // if the conditions dont hold try sample Z elemnts from _A again...
        tries--;
    }
    
    // sampling didn't work twice so we are doing the old SampledQMax maintanance
    int left = 0, right = _actualsizeMinusOne;
	while (left <= right) {
		int pivot_idx = left;
		int new_pivot_idx = PartitionAroundPivot(left, right, pivot_idx, _A);
		if (new_pivot_idx == _nminusq) {
            _curIdx = new_pivot_idx;
			return _A[new_pivot_idx];
		} else if (new_pivot_idx > _nminusq) {
			right = new_pivot_idx - 1;
		} else {  // new_pivot_idx < _q - 1.
			left = new_pivot_idx + 1;
		}
	}
	
}

void SampledQMax_LV::reset(){
	_phi = -1;
	_curIdx = _actualsize;
}
