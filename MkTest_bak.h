//---------------------------------------------------------------------------

#ifndef MkTestH
#define MkTestH
#include <math.h>
#include "MkFloat.h"
#include "MkInt.h"
#include "MkMatrix.h"
#include "MkPolygon.h"
#include "MkNurbs.h"
#include "MkEntity.h"
#include "MkLayer.h"
#include "MkPile.h"
#include "MkSection.h"
//---------------------------------------------------------------------------
#endif

void nurbtest();
void exceptiontest();
void polytest();
void simpletest();
void setlayer(MkLayers &layer,MkLoads &load,MkWall &w1, MkWall &w2);
void setpile(MkPiles &pile);
void setcut(MkCuts &cut,MkPiles &pile);
void setbc(MkBndConds &bc);
void setbcrange(MkRangeTree &rt);
void setload(MkLoads &load,MkWall &w1,MkWall &w2);
void setlloadrange(MkRangeTree &rt);
void setrloadrange(MkRangeTree &rt);
void setrankine(MkLoads &load,MkLayers &lay, MkCuts &cut, MkFills &fill, MkWall &w0);
void setsubreact(MkSubreacts &sub, MkLayers &lay, MkCuts &cut, MkFills &fill, MkPile &pile0);

