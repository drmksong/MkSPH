//---------------------------------------------------------------------------
#include "MkDimUnit.h"
//---------------------------------------------------------------------------

MkDimUnit::MkDimUnit()
{
  LenType = utMeter;
  TimeType = utSec;
  PressType = utTPM2;
  ForceType = utTONF;
}

float MkDimUnit::len()
{
  static float length;
  switch(LenType){
    case utMeter: length = 1.0; break;
    case utCentimeter: length = 100.0; break;
    case utMilimeter: length = 1000.0;break;
    default:length=1.0;
  }
  return length;
}

float MkDimUnit::area()
{
  return len()*len();
}

float MkDimUnit::volume()
{
  return len()*len()*len();
}

float MkDimUnit::time()
{
  static float t;
  switch(TimeType){
    case utSec: t = 1.0; break;
    case utMinute: t=60.0; break;
    case utHour: t = 3600.0; break;
    case utDay:t = 86400.0; break; 
    case utMonth: t = 2592000.0;break;
    case utYear: t = 946080000.0;break;
    default:t=1.0;
  }
  return t;

}

float MkDimUnit::press()
{
  return 1.0;
}

float MkDimUnit::stress()
{
  return press();
}

float MkDimUnit::force()
{
  return 1.0;
}

float MkDimUnit::weight()
{
  return force();
}

