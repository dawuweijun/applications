# Matheus da Silva Serpa - Ciência da Computação (2015)
# Universidade Federal do Pampa - Campus Alegrete
# matheusserpa@gmail.com
# https://github.com/matheusserpa/applications

# Flags
FLAGS := -DUSE_DOUBLE -Wall -Wextra

# Target rules
all: build

trapezoidal_seq.o: trapezoidal_seq.c
	gcc $(FLAGS) -o $@ -c $<

trapezoidal_seq: ../../helper.o trapezoidal_seq.o
	gcc $(FLAGS) -o $@ $+

build: trapezoidal_seq

run: build
	./trapezoidal_seq debug check n=1000000002

clean:
	rm -f trapezoidal_seq *.o 