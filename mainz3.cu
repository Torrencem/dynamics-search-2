#include "utils.cpp" 
#include "parse.cpp"
#include <iostream>

#include "eisred.cpp"

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

__global__ void run(long a_bound, long bs_bound, long search_bound, long a_min, long a_max) {
	long index = blockIdx.x * blockDim.x + threadIdx.x;
        long stride = (blockDim.x * gridDim.x);
	//printf("%ld\n", stride);
	for (long a = a_min + index; a < a_max; a += stride) {
		if (a % 100 == 3 || a % 100 == -3) {
			printf("progress: %ld\n", a);
		}
		for (long a2 = -a_bound; a2 < a_bound; a2++) {
			//printf("%ld , %ld", a, a2);
			for (long b = -bs_bound; b < bs_bound; b++) {
				for (long b2 = -bs_bound; b2 < bs_bound; b2++) {
					eis_int bt_(b, b2);
					if (b2 != 0 && bt_.phase_angle() >= (3.141593 / 9.0)) {
						continue;
					}
					eis_int bt = bt_*bt_*bt_;
					eis_int at(a, a2);
					
					if (bt.is_zero()) {
						continue;
					}
					
					// Check gcd
					if (!at.gcd(bt).is_unit()) {
						continue;
					}

					qw_elem r(at, bt);
					
					bool possibles = possible_periods_qw(r, search_bound, LT);
					if (possibles) {
						printf("manually check (%ld + w*%ld)/(%ld + w*%ld) (btph = %f)\n", at.a, at.b, bt.a, bt.b, bt_.phase_angle(true));
					}
				}
			}
		}
	}
}

int main(void) {
	// Get the number of available GPU's
	int dev_ct;
	gpuErrchk( cudaGetDeviceCount(&dev_ct) );
	
	long a_bd = 20000, bs_bd = 27;
	long N = 2 * a_bd / dev_ct;
	cout << "Starting search" << endl;
	Lookup_Table table = raw_lt_z3();
	
	//LT = table;
	
	long blockSize = 256;
	long numBlocks = (N + blockSize - 1) / blockSize;
	
	for (int i = 0; i < dev_ct; i++) {
		cout << "Initializing gpu " << i << endl;
		cudaSetDevice(i);
		LT = table; // Re-assign the lookup table on the new device
		run<<<numBlocks, blockSize>>>(a_bd, bs_bd, 3, -a_bd + i*(2 * a_bd / dev_ct), -a_bd + (i + 1)*(2 * a_bd / dev_ct));
	}
	
	cout << "Waiting for gpu's to finish computations" << endl;
	
	cudaDeviceSynchronize();
	
	cout << "Finished search" << endl;
	return 0;
}
