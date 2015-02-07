# Matheus S. Serpa
# matheusserpa@gmail.com
# https://github.com/matheusserpa/boltzmann

# Flags
FLAGS := -DUSE_DOUBLE -Wall -Wextra

# Target rules
all: build

helper.o: helper.c
	gcc $(FLAGS) -o $@ -c $<

lb_3D.o: lb_3D.c
	gcc $(FLAGS) -o $@ -c $<

build: helper.o lb_3D.o

clean:
	rm -f *.o