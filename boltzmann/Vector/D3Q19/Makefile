# Matheus da Silva Serpa - Ciência da Computação (2015)
# Universidade Federal do Pampa - Campus Alegrete
# matheusserpa@gmail.com
# https://github.com/matheusserpa/applications

# Flags
FLAGS := -DUSE_DOUBLE -Wall -Wextra

# Target rules
all: build

helper.o: helper.c
	gcc $(FLAGS) -o $@ -c $<

lb_3D.o: lb_3D.c
	gcc $(FLAGS) -o $@ -c $<

build: ../../../helper.o lb_3D.o

clean:
	rm -f *.o