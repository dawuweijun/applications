# Matheus da Silva Serpa - Ciência da Computação (2015)
# Universidade Federal do Pampa - Campus Alegrete
# matheusserpa@gmail.com
# https://github.com/matheusserpa/applications

# Flags
FLAGS := -DUSE_DOUBLE -fcilkplus -lcilkrts -Wall -Wextra

# Target rules
all: build

trapezoidal_cilk.o: trapezoidal_cilk.c
	gcc $(FLAGS) -o $@ -c $<

trapezoidal_cilk: ../../helper.o trapezoidal_cilk.o
	gcc $(FLAGS) -o $@ $+

build: trapezoidal_cilk

run: build
	./trapezoidal_cilk debug check n=1000000002 threads=4

clean:
	rm -f trapezoidal_cilk *.o 