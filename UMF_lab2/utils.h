#ifndef UTILS_H
#define UTILS_H

#include "types.h"

void fillElems(Data& system);
void fillTime(Data& system);
void input(Data& system);
void output(const Data& system);
void nodalError(const Data& system);
void elemError(const Data& system);

#endif