#ifndef MYGEMMH_
#define MYGEMMH_

#define MR 8
#define NR 6
#define MC 144
#define NC 216
#define KC 216

/*
  Any #define or function declaration must be provided in this header.
*/

void fiveloops( int, int, int, double *, int, int, double *, int, int,  double *, int, int );

#endif
