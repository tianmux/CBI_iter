tmake: ./src/main.cpp 
	g++ ./src/main.cpp -o CBI -fopenmp -std=c++14 -O3  -g #-I$(HOME)/parallelstl/include -tbb 

