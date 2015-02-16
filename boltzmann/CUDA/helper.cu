/* Matheus da Silva Serpa - Ciência da Computação (2015)
 * Universidade Federal do Pampa - Campus Alegrete
 * matheusserpa@gmail.com
 * https://github.com/matheusserpa/applications */
 
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include "helper.cuh"

/* Used by HANDLE_ERROR for show file and line where error occurred */
void HandleError(const char *file, const unsigned int line){
	fprintf(stderr, "Error in %s at line %u\n", file, line);
    exit(EXIT_FAILURE);
}

/* Return time in a specific moment */
double crono(){
	struct timeval tv;
	gettimeofday(&tv, 0);
	return (double)tv.tv_sec + (double)tv.tv_usec / 1e6;
}

/* Calculate the length of a string. */
unsigned int mystrlen(const char *s){
	unsigned int size;

	for(size = 0; s[size] != '\0'; ++size);

	return size;
}

/* Compare two strings. Delimiter is the final caracter. */
unsigned short int mystrcmp(const char *s1, const char *s2, const char delimiter){
	unsigned int i;

	for(i = 0; s1[i] != delimiter; ++i)
		if(s1[i] != s2[i])
			return false;

	return true;
}

/* Check if the flag exists. */
unsigned short int checkFlag(char **argv, int argc, char *flag){
	int i;

	for(i = 1; i < argc; ++i)
		if(mystrcmp(argv[i], flag, '\0'))
			return true;
	
	return false;
}

/* Return the flag in a integer value. */
int flagValueInt(char **argv, int argc, char *flag){
	int i;

	for(i = 1; i < argc; ++i)
		if(mystrcmp(argv[i], flag, '='))
			return atoi(&argv[i][mystrlen(flag)+1]);

	return -1;
}

/* Return the flag in a long long value. */
long long int flagValueLong(char **argv, int argc, char *flag){
	long long int i;

	for(i = 1; i < argc; ++i)
		if(mystrcmp(argv[i], flag, '='))
			return atoll(&argv[i][mystrlen(flag)+1]);

	return -1;
}

/* Return the flag in a text value. */
char *flagValueText(char **argv, int argc, char *flag){
	int i;

	for(i = 1; i < argc; ++i)
		if(mystrcmp(argv[i], flag, '='))
			return &argv[i][mystrlen(flag)+1];

	return '\0';
}

/* Return the flag in a real value. */
#ifdef USE_DOUBLE
double flagValueReal(char **argv, int argc, char *flag)
#else
float flagValueReal(char **argv, int argc, char *flag)
#endif
{
	int i;

	for(i = 1; i < argc; ++i)
		if(mystrcmp(argv[i], flag, '='))
			return atof(&argv[i][mystrlen(flag)+1]);

	return -1.0;
}

/* Show flags and their respective values. */
void flagStatus(char **argv, int argc){
	int i;

	printf("Flags\n");

	for(i = 1; i < argc; ++i)
		printf("%s\n", argv[i]);

	printf("\n");
}

