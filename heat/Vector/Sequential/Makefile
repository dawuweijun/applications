# Matheus da Silva Serpa - Ciência da Computação (2015)
# Universidade Federal do Pampa - Campus Alegrete
# matheusserpa@gmail.com
# https://github.com/matheusserpa/applications

# Flags
FLAGS := -DUSE_DOUBLE -Wall -Wextra

# Target rules
all: build

heat_seq.o: heat_seq.c
	gcc $(FLAGS) -o $@ -c $<

heat_seq: ../../../helper.o heat_seq.o
	gcc $(FLAGS) -o $@ $+

build: heat_seq

run: build
	./heat_seq size=10000

clean:
	rm -f heat_seq *.o 