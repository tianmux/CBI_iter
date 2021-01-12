tmake: ./src/main.cpp 
	icpc ./src/main.cpp -o CBI -qopenmp -qopt-report5 -qopt-report-phase=openmp,vec  -qopt-zmm-usage=high -std=c++14 -O3  -g -xCORE-AVX512 `gsl-config --cflags --libs` -I /opt/intel/advisor/include #-I$(HOME)/parallelstl/include -tbb 

