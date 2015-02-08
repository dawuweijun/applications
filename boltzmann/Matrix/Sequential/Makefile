# Matheus da Silva Serpa - Ciência da Computação (2015)
# Universidade Federal do Pampa - Campus Alegrete
# matheusserpa@gmail.com
# https://github.com/matheusserpa/applications

# Flags
FLAGS := -DUSE_DOUBLE -Wall -Wextra

# Target rules
all: build

lb_seq.o: lb_seq.c
	gcc $(FLAGS) -o $@ -c $<

lb_seq: ../D3Q19/helper.o ../D3Q19/lb_3D.o lb_seq.o
	gcc $(FLAGS) -o $@ $+

build: lb_seq

run: build
	./lb_seq obstacle="../../in/anb-128.obs" debug out output=saida.out

clean:
	rm -f lb_seq *.o 