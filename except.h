// error message facility
// $Log: /Speaker_Adaptation/lib/L_ui.h $
// 
// 3     02-04-10 17:26 Tchen

#include <stdio.h>
#include <stdlib.h>

#ifndef EXCEPT_H
#define EXCEPT_H


#define SayBye(msg) { fprintf (stderr, "Error - %s (%d): %s\n", __FILE__, __LINE__, msg); exit(-1); }

#define LookOut(msg) { fprintf (stderr,"Warning - %s (%d): %s\n", __LINE__, __FILE__, msg ); }

#endif