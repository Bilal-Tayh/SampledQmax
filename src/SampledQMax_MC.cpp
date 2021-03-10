#include "SampledQMax_MC.hpp"
#include <iostream>
#include <random>
#include <math.h> 

#include <vector>

#include <bitset>
#include "unistd.h"
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

void SampledQMax_MC::print(){
	for (int i = 0; i < _actualsize; ++i)
		std::cout << _A[i] << " ";
	std::cout << std::endl;
	std::cout << "_phi  = " << _phi << std::endl;
}



SampledQMax_MC::SampledQMax_MC(int q, float gamma){    
    _actualsize = q * (1+gamma);
	_actualsizeMinusOne = _actualsize - 1;
	_curIdx = 0;
	_A = (float*) malloc(sizeof(float) * _actualsize);
	if (!_A)
		exit(1);
	_gamma = gamma;
	_q = q;
	_qMinusOne = q - 1;
	_nminusq = _actualsize - q;
	_phi = 0;
    phi_x8 = _mm256_set1_ps((float)_phi);
    _delta = 1.0-0.001;
    _alpha=0.83;
    _k=ceil( ((_alpha*_gamma*(2+_gamma - _alpha*_gamma)) / (pow(_gamma -_alpha*_gamma,2))) * log(1/_delta));
    _Z = (int)  ( (_k*(1+_gamma)) / (_alpha*_gamma));
    gen_arr();
    rand_bits = _mm256_set_epi64x(gen_arr(), gen_arr(), gen_arr(), gen_arr());
    counter=0;
    Z_bound = q/10;
    switch_To_LV_Flag = 0;
    mask=0;
}




void SampledQMax_MC::insert(int v){
    if (v <= _phi){
		return;
	}
	else{
        
        if(switch_To_LV_Flag == 1){
            return  insert_LV(v);  
            
        }
        

        if(mask !=0){
            unsigned int emptySlot = _tzcnt_u32(mask);
            _A[_curIdx + emptySlot] = float(v);
            mask = mask & (~ int(pow(2,emptySlot)));
            if(mask == 0){
                _curIdx+=8;
            }
        }
                
        
        
        
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
                // to make transition from MC to LV we need to Partition around the pivot and update the current index
                if(switch_To_LV_Flag == 1){
                    int left = 0, right = _actualsizeMinusOne;
                    int pivot_idx = findValueIndex(_phi);
                    int new_pivot_idx = PartitionAroundPivot(left, right, pivot_idx, _A);
                    _curIdx = new_pivot_idx;
                    return  insert_LV(v);  
                    
                }
            }
            
            
            const __m256 item2 = _mm256_loadu_ps(_A+_curIdx);
            
            
            const __m256 match1 = _mm256_cmp_ps(item2,phi_x8, _CMP_LE_OS);
            
            int mask1 = _mm256_movemask_ps(match1);
            
            // if one of the 8 elemnts we checked are lower than phi we replace it
            if(mask1!=0){
                    unsigned int emptySlot = _tzcnt_u32(mask1);
                    _A[_curIdx + emptySlot] = float(v);
                    
                    mask = mask1 & (~ (int(pow(2,emptySlot))));
                    if(mask == 0){
                        _curIdx+=+8;
                    }
                return;
            }
            else{
                _curIdx+=8;
            }
        }
	}

}


void SampledQMax_MC::insert_LV(int v){
    _A[--_curIdx] = v;
    if (_curIdx){
        return;
    }
    else {
        maintenance_LV();
    }
}



void SampledQMax_MC::maintenance(){
	_phi = findKthLargestAndPivot();
    phi_x8 = _mm256_set1_ps((float)_phi);
    _curIdx = 0;
}


void SampledQMax_MC::maintenance_LV(){
	_phi = findKthLargestAndPivot_LV();
}

float* SampledQMax_MC::largestQ(){
	maintenance();
	return _A + _nminusq;
}

inline void swap(float &x, float &y){
	float z = x;
	x = y;
	y = z;
}


int SampledQMax_MC::PartitionAroundPivot(int left, int right, int pivot_idx, float* nums) {
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



uint32_t SampledQMax_MC::mm256_extract_epi32_var_indx(int i){   
    __m128i indx = _mm_cvtsi32_si128(i);
    __m256i val  = _mm256_permute2f128_si256(rand_bits, _mm256_castsi128_si256(indx),0);
    return _mm_cvtsi128_si32(_mm256_castsi256_si128(val));
}





int SampledQMax_MC::GenerateRandom(int max){
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
int SampledQMax_MC::findValueIndex(int value){
    for(int i=0;i<_actualsize;i++){
        if(_A[i]==value){
            return i;
        }
    }
    return -1;
}


// if value dont exist return -1
int SampledQMax_MC::findReplaceValueIndex(int value){
    int INT_MIN = -2147483648;
    for(int i=0;i<_actualsize;i++){
        if(_A[i]==value){
            _A[i]=float(INT_MIN);
            return i;
        }
    }
    return -1;
}


// check if the conditions holds for the possible pivot "value"
// return the pivot index in _A if it hold otherwise return -1
int SampledQMax_MC::checkPivot(int value){
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












int SampledQMax_MC::findKthLargestAndPivot_LV(){
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









int SampledQMax_MC::findKthLargestAndPivot(){
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
    
    updateZK();
    if(_Z >= Z_bound){
        switch_To_LV_Flag = 1;
    }

    return top;
        
	
}


void SampledQMax_MC::updateZK(){
    _delta=_delta/2.0;
    _k=ceil( ((_alpha*_gamma*(2+_gamma - _alpha*_gamma)) / (pow(_gamma -_alpha*_gamma,2))) * log(1/_delta));
    _Z = (int)  ( (_k*(1+_gamma)) / (_alpha*_gamma) );
}

void SampledQMax_MC::reset(){
	_phi = 0;
    phi_x8 = _mm256_set1_ps((float)_phi);
	_curIdx = 0;
}

