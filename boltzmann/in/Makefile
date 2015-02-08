# Matheus S. Serpa
# matheusserpa@gmail.com
# https://github.com/matheusserpa/boltzmann

# Flags
FLAGS := -Ofast -Wall -Wextra

# Target rules
all: build

generator.o: generator.c
	gcc $(FLAGS) -o $@ -c $<

generator: generator.o
	gcc $(FLAGS) -o $@ $+

build: generator

clean:
	rm -f generator *.o 