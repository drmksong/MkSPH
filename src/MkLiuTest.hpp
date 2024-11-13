//---------------------------------------------------------------------------
#ifndef MkLiuTestH
#define MkLiuTestH
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <GL/glut.h>
#include <boost/smart_ptr.hpp>
#include <iostream>
#include <vector>
#include <string>
#include "rvgs.h"
#include "MkMatrix.hpp"
#include "MkLiuKernel.hpp"
#include "MkLiuParticle.hpp"
#include "MkLiuGrid.hpp"
#include "MkLiuSPH_mk2.hpp"
#include "MkArc.hpp"
#include "MkCube.hpp"
#include "MkCircle.hpp"
#include "MkCylinder.hpp"
#include "MkLine.hpp"
#include "MkRect.hpp"
//---------------------------------------------------------------------------
void initGL(void);
void InitLiuSPH(void);
void display(void);
void update(void);

#endif
