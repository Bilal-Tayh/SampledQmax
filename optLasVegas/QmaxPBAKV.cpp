#include "QmaxPBAKV.hpp"
#include <iostream>
#include <unordered_map>

using namespace std;

QMaxPBAKV::QMaxPBAKV(int q, float gamma) {
	_actualsize = q * (1+gamma);
	_actualsizeMinusOne = _actualsize - 1;
	_curIdx = _actualsize;
	_K = (key*) malloc(sizeof(key) * _actualsize);
	_V = (val*) malloc(sizeof(val) * _actualsize);
	if (!_V || !_K) {
		exit(1);
	}
	_gamma = gamma;
	_q = q;
	_qMinusOne = q - 1;
	_nminusq = _actualsize - q;
	_phi = 0;
	z_star = 0;
}

QMaxPBAKV::~QMaxPBAKV() {
	free(_K);
	free(_V);
}

void QMaxPBAKV::insert(key k, weight x, double u) {
	auto it = hash_table.find(k);
	if (it != hash_table.end()) {
	//outputkv out = largestQ();
	//for(int i=0; i <_q; ++i) {
	//	if(out.keyArr[i] == k) {
			unsigned int i = it->second;
			_V[i].a *= _V[i].q;
			if(z_star != 0) {
				double w_by_z_star = _V[i].w / z_star;
				_V[i].q = _V[i].q < (w_by_z_star) ? _V[i].q : w_by_z_star;
			}
			if(_V[i].q != 0) {
				_V[i].a /= _V[i].q;
			}
			_V[i].a += x;
			_V[i].w += x;
			_V[i].w_by_u = _V[i].w / _V[i].u;
			return;
	//	}
	}
	pba_val v = {x, x, 1, u, x/u};
	if (v.w_by_u < _phi) {
		z_star = (z_star > v.w_by_u) ? z_star : v.w_by_u;
		return;
	} else {
		key key_star = _K[_curIdx];
		weight w_by_u_star = _V[_curIdx].w_by_u;
		z_star = (z_star > w_by_u_star) ? z_star : w_by_u_star;
		_K[--_curIdx] = k;
		_V[_curIdx] = v;
		hash_table[k] = _curIdx;
		if (_curIdx){
			return;
		} else {
			outputkv out = largestQ(); // performs maintenance implicitly
			hash_table.clear();
			for(int i=0; i<_q; ++i) {
				hash_table[out.keyArr[i]] = _nminusq + i;
			}
			//maintenance();
		}
	}
}

void QMaxPBAKV::maintenance() {
	_phi = findKthLargestAndPivot();
	_curIdx = _nminusq;
}

outputkv QMaxPBAKV::largestQ() {
	outputkv out;
	maintenance();
	out.keyArr = _K + _nminusq;
	out.valArr = _V + _nminusq;
	return out;
}

inline void QMaxPBAKV::swap(int a, int b) {
	key k = _K[a];
	_K[a] = _K[b];
	_K[b] = k;
	val v = _V[a];
	_V[a] = _V[b];
	_V[b] = v;
}


inline void QMaxPBAKV::swap1(int a, int b, val* num) {
	if (a==b) return;
	val v = num[a];
	num[a] = num[b];
	num[b] = v;
}



int QMaxPBAKV::PartitionAroundPivot(int left, int right, int pivot_idx, val* nums) {
	weight pivot_value = nums[pivot_idx].w_by_u;
	int new_pivot_idx = right;
	swap(pivot_idx, right);
	for (int i = right-1; i >= left; --i) {
		if (nums[i].w_by_u > pivot_value) {
			swap(i, --new_pivot_idx);
		}
	}
	swap(right, new_pivot_idx);
	return new_pivot_idx;
}



int QMaxPBAKV::PartitionAroundPivot1(int left, int right, int pivot_idx, val* nums) {
	weight pivot_value = nums[pivot_idx].w_by_u;
	int new_pivot_idx = right;
	swap1(pivot_idx, right,nums);
	for (int i = right-1; i >= left; --i) {
		if (nums[i].w_by_u > pivot_value) {
			swap1(i, --new_pivot_idx,nums);
		}
	}
	swap1(right, new_pivot_idx,nums);
	return new_pivot_idx;
}



int QMaxPBAKV::GenerateRandom(int min,int max){
  // uniform_real_distribution documentation
  // http://www.cplusplus.com/reference/random/uniform_real_distribution/
  std::uniform_real_distribution<double> distribution(min,max+1);
  double u = distribution(_generator);
  return int(u);
}





int QMaxPBAKV::findValueIndex(val value){
    for(int i=0;i<_actualsize;i++){
        if(_V[i]==value){
            return i;
        }
    }
    return -1;
}





// check if the conditions holds for the possible pivot "value"
// return the pivot index in _V if it hold otherwise return -1
int QMaxPBAKV::checkPivot(val value, double psi){
//     int index=0;
//     int bigger=0;
//     int smaller=0;
//     for(int i=0;i<_actualsize;i++){
//         if(_V[i]==value){
//             index=i;
//         }
//         else if(_V[i]>value){
//             bigger++;
//         }
//         else smaller++;
//     } 
//     
    
    
    
    
    int left = 0, right = _actualsizeMinusOne;
    int pivot_idx = findValueIndex(value);
    
     if(pivot_idx==-1){
        return -1;
    }
    
    int new_pivot_idx = PartitionAroundPivot(left, right, pivot_idx, _V);
    
    
//         std::cout<<"*****" <<new_pivot_idx <<"*******"<<std::endl;
//         std::cout<<"bigger= "<< bigger << ",   n-q= "<< _nminusq<<std::endl;
//         std::cout<<"should= "<< _q << ",   should = "<< _gamma*_q*psi<<std::endl;
    
    
    if (new_pivot_idx <= _nminusq) {
        if(new_pivot_idx >= _gamma*_q*psi) {  // new_pivot_idx < _q - 1.
            return _V[new_pivot_idx];
        }
    }
    return -1;
}




val QMaxPBKV::findKthLargestAndPivot() {
    double delta = 1.0-0.999;
    double alpha=0.83;
    double psi = 2.0/3.0;
    int k=ceil( ((alpha*_gamma*(2+_gamma - alpha*_gamma)) / (pow(_gamma -alpha*_gamma,2))) * log(1/delta));
    int Z = (int)  ( (k*(1+_gamma)) / (alpha*_gamma) );
    
    int tries=2;
    while(tries!=0){
        
        // B should contain Z random values from _V
        val *B = (val*) malloc(sizeof(val) * Z);
        for(int i=0;i<Z;i++){
            int j=GenerateRandom(0,_actualsize);
            B[i]=_V[j];
        }
        
        
        int left1 = 0, right1 = Z-1;
        int Kth_minimal_idx=0;
        //find kth minimal value in B
        while (left1 <= right1) {
            int pivot_idx = left1;
            int new_pivot_idx = PartitionAroundPivot1(left1, right1, pivot_idx, B);
            if (new_pivot_idx == k) {
                Kth_minimal_idx = new_pivot_idx;
                break;
            } else if (new_pivot_idx > k) {
                right1 = new_pivot_idx - 1;
            } else {  // new_pivot_idx < k - 1.
                left1 = new_pivot_idx + 1;
            }
        }


        
        // check if the conditions holds for the possible pivot B[Kth_minimal_idx]
        int idx = checkPivot(B[Kth_minimal_idx],psi);
        
        free(B);
        if(idx!=-1){
            return _V[idx].w_by_u;   
        }
        // if the conditions dont hold try sample Z elemnts from _V again...
        tries--;
    }
    int left = 0, right = _actualsizeMinusOne;
	while (left <= right) {
		int pivot_idx = left;
		int new_pivot_idx = PartitionAroundPivot(left, right, pivot_idx, _V);
		if (new_pivot_idx == _nminusq) {
			return _V[new_pivot_idx].w_by_u;
		} else if (new_pivot_idx > _nminusq) {
			right = new_pivot_idx - 1;
		} else {  // new_pivot_idx < _q - 1.
			left = new_pivot_idx + 1;
		}
	}
	
}

