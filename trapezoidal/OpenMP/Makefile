# Matheus da Silva Serpa - Ciência da Computação (2015)
# Universidade Federal do Pampa - Campus Alegrete
# matheusserpa@gmail.com
# https://github.com/matheusserpa/applications

# Flags
FLAGS := -DUSE_DOUBLE -fopenmp -Wall -Wextra

# Target rules
all: build

trapezoidal_mp.o: trapezoidal_mp.c
	gcc $(FLAGS) -o $@ -c $<

trapezoidal_mp: ../../helper.o trapezoidal_mp.o
	gcc $(FLAGS) -o $@ $+

build: trapezoidal_mp

run: build
	./trapezoidal_mp debug check n=1000000002 threads=4

clean:
	rm -f trapezoidal_mp *.o 