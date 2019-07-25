#include "utils.cpp" 
#include "parse.cpp"
#include <iostream>

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

using namespace std;

__device__ long gcd(long a, long b) {
	long temp;
	while (b != 0) {
		temp = a % b;
		a = b;
		b = temp;
	}
	return a;
}

__device__ __managed__ Lookup_Table LT;

__global__ void run(long a_start, long a_end, long bs_bound, long search_bound) {
	long index = blockIdx.x * blockDim.x + threadIdx.x;
        long stride = (blockDim.x * gridDim.x);
	for (long a = a_start + index; a < a_end; a += stride) {
		if (a % 100000000 == 3 || a % 100000000 == -3)  {
			// Print progress report
			printf("progress %ld\n", a);
		}
		for (long b_s_ = 1; b_s_ <= bs_bound; b_s_++) {
			long b_s = b_s_*2;
			long b = (b_s * b_s) * (b_s * b_s);
			// Check analytic bound
			float fa = __ll2float_rd(a);
			float fb = __ll2float_rd(b);
			if (fa / fb > -0.913942) continue;
			// Check gcd
			long g = gcd(a, b);
			if (g != 1 && g != -1) continue;
			// Perform actual computation
			bool possibles = possible_periods(a, b, search_bound, LT);
			if (possibles) {
				printf("manually check %ld / %ld \n", a, b);
			}
		}
	}
}

int main(void) {
	// Get the number of available GPU's
	int dev_ct;
	gpuErrchk( cudaGetDeviceCount(&dev_ct) );
	
	long a_bd = 100000000000, bs_bd = 281;
	long N = a_bd / dev_ct;
	cout << "Starting search" << endl;
	Lookup_Table table = raw_lt();
	
	//LT = table;
	
	long blockSize = 512;
	long numBlocks = (N + blockSize - 1) / blockSize;
	
	for (int i = 0; i < dev_ct; i++) {
		cout << "Initializing gpu " << i << endl;
		cudaSetDevice(i);
		LT = table; // Re-assign the lookup table on the new device
		cout << i << ": " << (-a_bd + i*((a_bd) / dev_ct)) << ", " << -a_bd + (i + 1)*((a_bd) / dev_ct) << endl;
		run<<<numBlocks, blockSize>>>(-a_bd + i*((a_bd) / dev_ct), -a_bd + (i + 1)*((a_bd) / dev_ct), bs_bd, 3);
	}
	
	cout << "Waiting for gpu's to finish computations" << endl;
		
	gpuErrchk( cudaDeviceSynchronize() );
	
	cout << "Finished search" << endl;
	return 0;
}

