# Matheus da Silva Serpa - Ciência da Computação (2015)
# Universidade Federal do Pampa - Campus Alegrete
# matheusserpa@gmail.com
# https://github.com/matheusserpa/applications

# Flags
FLAGS := -DUSE_DOUBLE -pthread -Wall -Wextra

# Target rules
all: build

heat_pthread.o: heat_pthread.c
	gcc $(FLAGS) -o $@ -c $<

heat_pthread: ../../../helper.o heat_pthread.o
	gcc $(FLAGS) -o $@ $+

build: heat_pthread

run: build
	./heat_pthread size=15000 threads=4

clean:
	rm -f heat_pthread *.o 