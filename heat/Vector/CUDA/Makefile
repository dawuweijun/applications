# Matheus da Silva Serpa - Ciência da Computação (2015)
# Universidade Federal do Pampa - Campus Alegrete
# matheusserpa@gmail.com
# https://github.com/matheusserpa/applications

# Flags
FLAGS := -DUSE_DOUBLE -m64 -arch=sm_35 --use_fast_math --ptxas-options=-v

# Target rules
all: build

heat_cuda.o: heat_cuda.cu
	nvcc $(FLAGS) -o $@ -c $<

heat_cuda: ../../../boltzmann/CUDA/helper.o heat_cuda.o
	nvcc $(FLAGS) -o $@ $+

build: heat_cuda

run: build
	./heat_cuda size=10000

clean:
	rm -f heat_cuda *.o 