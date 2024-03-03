#ifndef MYGEMMH_
#define MYGEMMH_

#define MR 8
#define NR 6
#define MC 480
#define NC 1440
#define KC 48
#include <time.h>

/*
  Any #define or function declaration must be provided in this header.
*/

void fiveloops( int, int, int, double *, int, int, double *, int, int,  double *, int, int );

#endif
