# Matheus da Silva Serpa - Ciência da Computação (2015)
# Universidade Federal do Pampa - Campus Alegrete
# matheusserpa@gmail.com
# https://github.com/matheusserpa/applications

# Flags
FLAGS := -DUSE_DOUBLE -fopenmp -Wall -Wextra

# Target rules
all: build

heat_mp.o: heat_mp.c
	gcc $(FLAGS) -o $@ -c $<

heat_mp: ../../../helper.o heat_mp.o
	gcc $(FLAGS) -o $@ $+

build: heat_mp

run: build
	./heat_mp size=15000 threads=4
	#./heat_mp threads=4 check debug

clean:
	rm -f heat_mp *.o 