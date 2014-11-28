/* 
  parse the arguments to the main program
 */

#include <string.h>
#include <ctype.h>
#include "arg.h"
#include "gsl/gsl_errno.h"

static int SetVal (Config_t *cp, char *str);
static void SPrintVal (Config_t *cp, char *str);

int pconf
/*------------------------------------------------------------*
 * DESCRIPTION
 *    Configure arguments according to configuration structure.
 * Ignore arguments that don't match the configuration list.
 *
 * RETURN
 * non zero if an argument value can't be parsed or if
 * the special argument -? or -help appears in the argument
 * list.
 *
 * Return 0 on success.
 */
(
 int   argc,    // arg count
 char    *argv[], // arguments
 Config_t *config_p  // argument configuratiuon
 )
{
   int i, parsed;
   int bad_usage = false;
   Config_t *cp;

   for (i = 1; i < argc; i++) {
      for (parsed = false, cp = config_p; cp->arg_type != NOTYPE; cp++) {
         if ( strcmp (argv[i], cp->swtch) == 0) {
            parsed = true;
            if (++i < argc)
               bad_usage |= SetVal (cp, argv[i]);
            else
               bad_usage = true;

            break;
            }
         }

      if (parsed)
         continue;

      if ( (strcmp ("-?", argv[i]) == 0) || (strcmp ("-help", argv[i]) == 0))
         bad_usage = true;
   }
   return (bad_usage);
}

int ppconf
/*------------------------------------------------------------*
 * DESCRIPTION
 *    Configure arguments according to configuration structure.
 * - Arguments that are parsed are deleted.
 * - Ignore arguments that don't match the configuration list unless
 * last is true in which case bad_usage is set and and we print an
 * error.
 * - The effect is to allow incremental argument parsing.
 *
 * RETURN
 * non zero if an argument value can't be parsed or if
 * the special argument -? or -help appears in the argument
 * list.
 *
 * Return 0 on success.
 */
(
 int   argc,    // arg count
 char    *argv[], // arguments
 Config_t *config_p, // argument configuratiuon
 char  last    // flag indicating that all remain arguments should be parsed
 )
{
   int i, parsed;
   int bad_usage = false;
   Config_t *cp;

   for (i = 1; i < argc; i++) {
      /* argument has been processed already */
      if (argv[i][0] == '\0')
         continue;

      if ((strcmp ("-?", argv[i]) == 0) || (strcmp ("-help", argv[i]) == 0))
         bad_usage = true;

      for (parsed = false, cp = config_p; cp->arg_type != NOTYPE; cp++) {
		 if (cp->arg_type == SEPERATOR)	continue;
         if (strcmp (argv[i], cp->swtch) == 0) {
            parsed = true;
            /* remove this switch from consideration */
            argv[i][0] = '\0';
            if (++i < argc) {
               bad_usage |= SetVal (cp, argv[i]);
               argv[i][0] = '\0';
               }
            else
               bad_usage = true;

            break;  // From inner for.
         }
      }

      if (parsed || !last)
         continue;

      printf("%s: Unrecognized argument, %s\n", argv[0], argv[i]);
      bad_usage = true;
   }
   return (bad_usage);
}

void pusage
/*------------------------------------------------------------*
 * DESCRIPTION
 * Print the usage of a configuration.
 */
(
 char *prog,
 Config_t *cp
 )
{
   char valstr[2560];

   printf("Usage: %s\n", prog);
   for (; cp->arg_type != NOTYPE && cp->arg_type != SEPERATOR; cp++) {
      SPrintVal (cp, valstr);
      //      printf(" [%s %s] %s - %s\n", cp->swtch, valstr, cp->LongName, cp->Doc);
      printf(" [%s %s] %s\n", cp->swtch, valstr, cp->Doc);
   }
}

static char *env_scan(char *str)
/*------------------------------------------------------------*
 * DESCRIPTION: substitutes environment values into a string
 */
{
   char buf[1024];             /* buffer for temp use */
   register char *p = buf;     /* holds place in the buffer */
   char var[50];               /* holds the name of the env variable */
   char *val;

   while (*str)
      if (*str == '$') {
         if (*++str == '$') {
            *p++ = *str++;
            continue;
            }
         val = var;
         while (isalnum(*str) || *str == '_')
            *val++ = *str++;
         *p = '\0';
	 if ( val == var ) val = "$";
	 else val = getenv(var);
	 //val = (val == var) ? "$" : getenv(var);
         if (val) {
            strcat(p, val);
            p += strlen(val);
            }
         } else
   *p++ = *str++;
   *p = '\0';

   if (strcmp(buf, "NULL") == 0)
      return NULL;

   //return CM_strdup(buf);
   //return _strdup(buf);
   char *ret = new char [strlen(buf)+1];
   strcpy ( ret, buf );
   return ret;
}

static int SetVal (Config_t *cp, char *str)
{
   switch (cp->arg_type) {
   case CHAR:
      *(cp->var.CharP) = (char) str[0];
      break;

   case BYTE:
      *(cp->var.ByteP) = (char) atoi (str);
      break;

   case U_BYTE:
      *(cp->var.UByteP) = (u_char) atoi (str);
      break;

   case SHORT:
      *(cp->var.ShortP) = (int) atoi (str);
      break;

   case U_SHORT:
      *(cp->var.UShortP) = (unsigned int) atoi (str);
      break;

   case INT:
      *(cp->var.IntP) = (int) atoi (str);
      break;

   case U_INT:
      *(cp->var.UIntP) = (u_int) atoi (str);
      break;

   case FLOAT:
      *(cp->var.FloatP) = (float) atof (str);
      break;

   case DOUBLE:
      *(cp->var.DoubleP) = (double) atof (str);
      break;

   case BOOL:
      if ((strcmp ("on", str) == 0) || (strcmp ("1", str) == 0) ||
            (strcmp ("t", str) == 0) ||  (strcmp ("true", str) == 0) ||
            (strcmp ("y", str) == 0) ||  (strcmp ("yes", str) == 0))
         *(cp->var.BoolP) = true;
      else if ((strcmp ("off", str) == 0) || (strcmp ("0", str) == 0) ||
                 (strcmp ("f", str) == 0) || (strcmp ("false", str) == 0) ||
                 (strcmp ("n", str) == 0) || (strcmp ("no", str) == 0))
         *(cp->var.BoolP) = false;
      else
         return (-1);
      break;

   case STRING:
      //sprintf (*(cp->var.StringP), "%s", str);
      //char buf[1024];
	  //sprintf(buf, "%s", str);
	   *(cp->var.StringP) = env_scan(str);
      break;

   default:
      printf("bad param type %d\n", cp->arg_type);
      return (-1);
   }
   return (0);
}

static void SPrintVal (Config_t *cp, char *str)
{
   switch (cp->arg_type) {
   case CHAR:
      sprintf (str, "%c", *(cp->var.CharP));
      break;

   case BYTE:
      sprintf (str, "%d", *(cp->var.ByteP));
      break;

   case U_BYTE:
      sprintf (str, "%u", *(cp->var.UByteP));
      break;

   case SHORT:
      sprintf (str, "%d", *(cp->var.ShortP));
      break;

   case U_SHORT:
      sprintf (str, "%u", *(cp->var.UShortP));
      break;

   case INT:
      sprintf (str, "%d", *(cp->var.IntP));
      break;

   case U_INT:
      sprintf (str, "%u", *(cp->var.UIntP));
      break;

   case FLOAT:
      sprintf (str, "%f", *(cp->var.FloatP));
      break;

   case DOUBLE:
      sprintf (str, "%g", *(cp->var.DoubleP));
      break;

   case BOOL:
      sprintf (str, "%s", (*(cp->var.BoolP) ? "TRUE" : "FALSE"));
      break;

   case STRING:
      sprintf (str, "%s", *(cp->var.StringP));
      break;

   default:
      sprintf (str, "bad param type %d\n", cp->arg_type);
   }
}


char *find_switch (char **variable, Config_t *vparam, int n)
{
   int  i;

   for (i=0; i < n; i++) {
      if (vparam[i].var.StringP == variable)
         return(vparam[i].swtch);
   }
   return NULL;
}

int required(char **variable, Config_t *vparam, int n)
{
   if ((*variable) == NULL) {
      printf("need to specify %s\n", find_switch(variable, vparam, n));
      return(true);
   }
   return (false);
}

void countf(char *fn, int *nr, int *nc)
{
  int i, j, c, n, obs, gsl_errno;
  double l;
  char tmp[1024], item[101];

  FILE *fp = NULL;
  if ( ( fp = fopen(fn, "rt") ) == NULL )
    gsl_error("Cannot open file to read", __FILE__, __LINE__, gsl_errno);

  i = 0; c = getc(fp);
  while (c!=EOF) { /* count number of lines */
    i++;
    while ( c!=EOF && c!='\n' ) c = getc(fp);
    c = getc(fp);
  }
  *nr = i; // number of observations

  rewind(fp);
  
  n = 0; 
  *nc = 0;

  for ( i = 0; i < *nr;  ) {

    do {
      c = getc(fp); 
      //      if ( c == '\n' ) i ++;
      if ( c=='\n' ) {
	if ( *nc == 0 ) {
	  *nc = n;
	} else {
	  if ( *nc != n ) {
	    gsl_error("Number of columns is not consistent in data file", __FILE__, __LINE__, gsl_errno);
	    //	  exit (-1);
	  }
	}

	i ++ ;
	n = 0;
      }

    } while (c==' ' || c=='\t' || c==',' || c==';' || c=='\n' || c=='\r');

    if (c==EOF) break;
    j = 0;
    while (c!=EOF && c!='\n' && c!=' ' && c!='\t' && c!=',' && c!=';') {
      if (j<100) item[j++] = c;
      c = getc(fp);
    }
    item[j] = 0;
    if ( sscanf(item, "%lf", &l) != 1 ){
      gsl_error("Bad numeric item on datafile", __FILE__, __LINE__, gsl_errno);
      exit(1);
    }
    n ++;

    if ( c=='\n' ) {

      if ( *nc == 0 ) {
	*nc = n;
      } else {
	if ( *nc != n ) 
	  gsl_error("Number of columns is not consistent in data file", __FILE__, __LINE__, gsl_errno);
      }

      i ++ ;
      n = 0;
    }

  }
  fclose(fp);
#if 0
  printf("nros = %d ncol = %d\n", *nr, *nc);
#endif

}

void parse_range(char *str, int *low, int *high, int *modulus)
/* Parses a string of the following form: [@][low][:[high]][+modulus], where
   'low' and 'high' must be a non-negative integers; 'modulus' must be a
   positive integer.  If 'low' is omitted, it defaults to one (not zero, even
   though zero is a legal value).  If 'modulus' is omitted, it defaults to
   zero.  If 'high' is omitted but the colon or a modulus is present, 'high'
   is set to -1 as a signal to the caller.  If just a single number appears,
   it is interpreted as a value for 'low', and 'high' is set to -2. If an
   initial @ appears, '-low' is returned instead of 'low'. Something must
   appear (i.e. the empty string is not legal).

   The caller provides pointers to places to store 'low', 'high', and
   'modulus'.  The pointer for 'modulus' may be null, in which case a modulus
   specification is illegal.

   Incorrect specifications lead to an error message being displayed, and the
   program being terminated.  The checks include verifying that low<=high
   when 'high' is specified, but the caller may need to check this also after
   in cases where 'high' is left to default. */
{
  char *s;

  s = str;
  if (*s==0) goto error;
  if (*s=='@') s += 1;                                /* skip an initial '@' */

  *low = 1;                   /* look for value for 'low', or let it default */
  if (*s>='0' && *s<='9') {
    *low = 0;
    while (*s>='0' && *s<='9') { *low = 10*(*low) + (*s-'0'); s += 1; }
  }
  if (*str=='@') *low = -*low;                 /* negate low iff leading '@' */

  *high = *s==0 ? -2 : -1;                     /* look for value for 'high', */
  if (*s==':') {                                    /* or signal its absence */
    s += 1;
    if (*s!=0 && *s>='0' && *s<='9') {
      *high = 0;
      while (*s>='0' && *s<='9') { *high = 10*(*high) + (*s-'0'); s += 1; }
      if (*high<*low){
        fprintf(stderr,"High end of range is less than low end: %s\n",str);
        exit(1); }
    }
  }

  if (modulus!=0) *modulus = 0;              /* look for value for 'modulus' */
  if (*s=='+') {                         /* if it's legal, or let it default */
    if (modulus==0) goto error;
    s += 1;
    *modulus = 0;
    while (*s>='0' && *s<='9') { *modulus = 10*(*modulus) + (*s-'0'); s += 1; }
    if (*modulus==0) goto error;
  }

  if (*s!=0) goto error;                         /* check for garbage at end */
  return;

error:                                                       /* report error */
  fprintf(stderr,"Bad range specification: %s\n",str); exit(1);
}
