#ifndef S_QMAXLV_H
#define S_QMAXLV_H
#include <fstream>
#include <assert.h>
#include <time.h>
#include <immintrin.h>
#include <iomanip> 
#include "RngFast.hpp"
#include <queue>

class SampledQMax_LV
{
	int *_A;
	int _curIdx;
	int _q;
	int _actualsize;
	int _actualsizeMinusOne;
	int _qMinusOne;
	float _gamma;
	int _nminusq;
	int _phi;
    std::default_random_engine _generator;
	void maintenance();
	int PartitionAroundPivot(int left, int right, int pivot_idx, int* nums);
    double _delta;
    double _alpha;
    double _psi;
    int _k;
    int _Z;
    rng::rng128 gen_arr;
    int counter;
    __m256i rand_bits;
    std::priority_queue <int, std::vector<int> ,std::greater<int>> testp;
    int testSize;
public:
	void reset();
	int findKthLargestAndPivot();
	SampledQMax_LV(int q, float gamma);
	void insert(int id);
	int* largestQ();
	void print();
    int checkPivot(int value);
    int findValueIndex(int value);
    int GenerateRandom(int max);
    int test();
    int findReplaceValueIndex(int value);
    uint32_t mm256_extract_epi32_var_indx(int i);
};
#endif

