# Matheus da Silva Serpa - Ciência da Computação (2015)
# Universidade Federal do Pampa - Campus Alegrete
# matheusserpa@gmail.com
# https://github.com/matheusserpa/applications

# Flags
FLAGS := -DUSE_DOUBLE -Wall -Wextra

# Target rules
all: build

seidel_seq.o: seidel_seq.c
	gcc $(FLAGS) -o $@ -c $<

seidel_seq: ../../helper.o seidel_seq.o
	gcc $(FLAGS) -o $@ $+

build: seidel_seq

run: build
	./seidel_seq debug check n=1000000002

clean:
	rm -f seidel_seq *.o 