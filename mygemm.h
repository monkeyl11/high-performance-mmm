#ifndef MYGEMMH_
#define MYGEMMH_

#define MR 8
#define NR 6
#define MC 240
#define NC 960
#define KC 240

/*
  Any #define or function declaration must be provided in this header.
*/

void fiveloops( int, int, int, double *, int, int, double *, int, int,  double *, int, int );

#endif
