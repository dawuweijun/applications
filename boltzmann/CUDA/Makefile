# Matheus da Silva Serpa - Ciência da Computação (2015)
# Universidade Federal do Pampa - Campus Alegrete
# matheusserpa@gmail.com
# https://github.com/matheusserpa/applications

# Flags
FLAGS := -DUSE_DOUBLE -m64 -arch=sm_35 --use_fast_math --ptxas-options=-v

# Target rules
all: build

lb_cuda.o: lb_cuda.cu
	nvcc $(FLAGS) -o $@ -c $<

lb_3D.o: lb_3D.cu
	nvcc $(FLAGS) -o $@ -c $<

helper.o: helper.cu
	nvcc $(FLAGS) -o $@ -c $<

lb_cuda: helper.o lb_3D.o lb_cuda.o
	nvcc $(FLAGS) -o $@ $+

build: lb_cuda

run: build
	./lb_cuda obstacle="../in/anb-128.obs" debug out output=saida.out

clean:
	rm -f lb_cuda *.o