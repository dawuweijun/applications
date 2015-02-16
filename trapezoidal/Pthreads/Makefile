# Matheus da Silva Serpa - Ciência da Computação (2015)
# Universidade Federal do Pampa - Campus Alegrete
# matheusserpa@gmail.com
# https://github.com/matheusserpa/applications

# Flags
FLAGS := -DUSE_DOUBLE -pthread -Wall -Wextra

# Target rules
all: build

trapezoidal_pthread.o: trapezoidal_pthread.c
	gcc $(FLAGS) -o $@ -c $<

trapezoidal_pthread: ../../helper.o trapezoidal_pthread.o
	gcc $(FLAGS) -o $@ $+

build: trapezoidal_pthread

run: build
	./trapezoidal_pthread debug check n=1000000002 threads=4

clean:
	rm -f trapezoidal_pthread *.o 