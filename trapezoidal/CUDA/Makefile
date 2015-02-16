# Matheus da Silva Serpa - Ciência da Computação (2015)
# Universidade Federal do Pampa - Campus Alegrete
# matheusserpa@gmail.com
# https://github.com/matheusserpa/applications

# Flags
FLAGS := -DUSE_DOUBLE -m64 -arch=sm_35 --use_fast_math --ptxas-options=-v

# Target rules
all: build

trapezoidal_cuda.o: trapezoidal_cuda.cu
	nvcc $(FLAGS) -o $@ -c $<

trapezoidal_cuda: ../../boltzmann/CUDA/helper.o trapezoidal_cuda.o
	nvcc $(FLAGS) -o $@ $+

build: trapezoidal_cuda

run: build
	./trapezoidal_cuda debug check n=1000000002

clean:
	rm -f trapezoidal_cuda *.o 