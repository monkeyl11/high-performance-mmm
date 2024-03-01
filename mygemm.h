#ifndef MYGEMMH_
#define MYGEMMH_

#define MR 8
#define NR 6
#define MC 16
#define NC 24
#define KC 12

/*
  Any #define or function declaration must be provided in this header.
*/

void fiveloops( int, int, int, double *, int, int, double *, int, int,  double *, int, int );

#endif
