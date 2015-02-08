# Matheus da Silva Serpa - Ciência da Computação (2015)
# Universidade Federal do Pampa - Campus Alegrete
# matheusserpa@gmail.com
# https://github.com/matheusserpa/applications

# Flags
FLAGS := -DUSE_DOUBLE -fopenmp -Wall -Wextra

# Target rules
all: build

lb_mp.o: lb_mp.c
	gcc $(FLAGS) -o $@ -c $<

lb_mp: ../D3Q19/helper.o ../D3Q19/lb_3D.o lb_mp.o
	gcc $(FLAGS) -o $@ $+

build: lb_mp

run: build
	./lb_mp obstacle="../../in/anb-128.obs" debug out output=saida.out threads=4

clean:
	rm -f lb_mp *.o 