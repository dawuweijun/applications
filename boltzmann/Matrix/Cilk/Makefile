# Matheus S. Serpa
# matheusserpa@gmail.com
# https://github.com/matheusserpa/boltzmann

# Flags
FLAGS := -DUSE_DOUBLE -fcilkplus -lcilkrts -Wall -Wextra

# Target rules
all: build

lb_cilk.o: lb_cilk.c
	~/cilkplus-install/bin/gcc $(FLAGS) -o $@ -c $<

lb_cilk: ../D3Q19/helper.o ../D3Q19/lb_3D.o lb_cilk.o
	~/cilkplus-install/bin/gcc $(FLAGS) -o $@ $+

build: lb_cilk

run: build
	./lb_cilk obstacle="../../in/anb-128.obs" debug out output=saida.out threads=4

clean:
	rm -f lb_cilk *.o 