#ifndef FEM_H
#define FEM_H

#include "types.h"

double phi(int a, double xi);          // a = 0,1
double dphi_dxi(int a, double xi);
double globalNodeCoord(const Data& system, int id);

double sigma(double u);
double dsigma_du(double u);
double rhs_f(double x, double t);      // не зависит от u
double exactSolution(double x, double t); // для тестов

void buildMatrixSimpleIter(Data& system, double dt, double t);
void buildMatrixNewton(Data& system, double dt, double t);
void applyBC(Data& system, int type0, double val0, int type1, double val1);

#endif