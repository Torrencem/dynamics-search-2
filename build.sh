nvcc -std=c++11 -arch=compute_30 -Xptxas -dlcm=ca main.cu -O3 -o main
