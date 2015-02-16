/* Matheus da Silva Serpa - Ciência da Computação (2015)
 * Universidade Federal do Pampa - Campus Alegrete
 * matheusserpa@gmail.com
 * https://github.com/matheusserpa/applications */
 
#ifndef HELPER
#define HELPER

#ifndef true
#define true 1
#endif

#ifndef false
#define false 0
#endif

enum Boltzmann{REDISTRIBUTE = 0, PROPAGATE, BOUNCEBACK, RELAXATION, EXECUTION, AUX};

/* Call this, if one error probably occurred. */
#define HANDLE_ERROR() (HandleError(__FILE__, __LINE__))

/* Used by HANDLE_ERROR for show file and line where error occurred. */
void HandleError(const char *file, const unsigned int line);

/* Return time in a specific moment. */
double crono();

/* Calculate the length of a string. */
unsigned int mystrlen(const char *s);

/* Compare two strings. Delimiter is the final caracter. */
unsigned short int mystrcmp(const char *s1, const char *s2, const char delimiter);

/* Check if the flag exists. */
unsigned short int checkFlag(char **argv, int argc, char *flag);

/* Return the flag in a integer value. */
int flagValueInt(char **argv, int argc, char *flag);

/* Return the flag in a long long value. */
long long int flagValueLong(char **argv, int argc, char *flag);

/* Return the flag in a text value. */
char *flagValueText(char **argv, int argc, char *flag);

/* Return the flag in a real value. */
#ifdef USE_DOUBLE
double flagValueReal(char **argv, int argc, char *flag);
#else
float flagValueReal(char **argv, int argc, char *flag);
#endif


/* Show flags and their respective values. */
void flagStatus(char **argv, int argc);

#endif