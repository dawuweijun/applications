# Matheus da Silva Serpa - Ciência da Computação (2015)
# Universidade Federal do Pampa - Campus Alegrete
# matheusserpa@gmail.com
# https://github.com/matheusserpa/applications

# Flags
FLAGS := -DUSE_DOUBLE -pthread -Wall -Wextra

# Target rules
all: build

lb_pthread.o: lb_pthread.c
	gcc $(FLAGS) -o $@ -c $<

lb_pthread: ../../../helper.o ../D3Q19/lb_3D.o lb_pthread.o
	gcc $(FLAGS) -o $@ $+

build: lb_pthread

run: build
	./lb_pthread obstacle="../../in/anb-128.obs" debug out output=saida.out threads=3

clean:
	rm -f lb_pthread *.o 