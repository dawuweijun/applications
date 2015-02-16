# Matheus da Silva Serpa - Ciência da Computação (2015)
# Universidade Federal do Pampa - Campus Alegrete
# matheusserpa@gmail.com
# https://github.com/matheusserpa/applications

# Flags
FLAGS := -DUSE_DOUBLE -fcilkplus -lcilkrts -Wall -Wextra

# Target rules
all: build

heat_cilk.o: heat_cilk.c
	~/cilkplus-install/bin/gcc $(FLAGS) -o $@ -c $<

heat_cilk: ../../../helper.o heat_cilk.o
	~/cilkplus-install/bin/gcc $(FLAGS) -o $@ $+

build: heat_cilk

run: build
	./heat_cilk size=15000 threads=4
	#./heat_cilk threads=4 check debug

clean:
	rm -f heat_cilk *.o 