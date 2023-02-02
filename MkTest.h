//---------------------------------------------------------------------------
#ifndef MkTestH
#define MkTestH
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <GL/glut.h>
#include "rvgs.h"
#include "MkMatrix.hpp"
#include "MkSmoothFunc.h"
#include "MkParticle.hpp"
#include "MkGrid.h"
//---------------------------------------------------------------------------
#endif

int particles();
void drawBox(void);
void display(void);
void init(void);
float Wijab(MkParticle &mp1, MkParticle &mp2);
