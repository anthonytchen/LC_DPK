/*
 parse the arguments to the main program
 and others
*/

#ifndef _H_ARG_
#define _H_ARG_

#include <stdio.h>
#include <stdlib.h>


typedef enum { NOTYPE, SEPERATOR,BYTE, SHORT, INT, LONG,U_BYTE, U_SHORT, U_INT, U_LONG,FLOAT, DOUBLE,BOOL, CHAR, STRING
} arg_t;

typedef char* caddr_t;
typedef unsigned char u_char;
typedef unsigned short u_short;
typedef unsigned int u_int;
typedef unsigned long u_long;


typedef union _ptr {
   char *CharP;
   char *ByteP;
   u_char  *UByteP;
   short   *ShortP;
   u_short *UShortP;
   int  *IntP;
   u_int   *UIntP;
   long *LongP;
   u_long  *ULongP;
   float   *FloatP;
   double  *DoubleP;
   int      *BoolP;
   char **StringP;
} ptr_t;

typedef struct _config {
   char  *LongName;     /* Long Name */
   char  *Doc;       /* Documentation string */
   char  *swtch;        /* Switch Name */
   arg_t arg_type;      /* Argument Type */
   caddr_t  var;        /* Pointer to the variable */
} config_t;

typedef struct _InternalConfig {
   char  *LongName;     /* Long Name */
   char  *Doc;       /* Documentation string */
   char  *swtch;        /* Switch Name */
   arg_t arg_type;      /* Argument Type */
   ptr_t var;        /* Pointer to the variable */
} Config_t;

int required(char **variable, Config_t *vparam, int n);
void pusage (char *prog, Config_t *cp);
int pconf (
             int    argc,    // arg count
             char  *argv[], // arguments
             Config_t *config_p   // argument configuratiuon
             );
int ppconf (
              int    argc,    // arg count
              char  *argv[], // arguments
              Config_t *config_p,  // argument configuratiuon
              char   last    // flag indicating that all remain arguments should be parsed
              );

// count the number of rows & columns of a data file
// if not consistent exit
void countf (char *fn, int *nr, int *nc);

// parse the range of mcmc samples for averaging from the string <str>
void parse_range(char *str, int *low, int *high, int *modulus);
#endif


