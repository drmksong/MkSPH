//---------------------------------------------------------------------------
// This module is general purposed simple graphic class to store, draw,
// manipulate object. It is well suited to VCL component, but not restricted.
// It forms the base for the higher level class, such as tunnel component.
//
// Copyright (c) 1999 Myung Kyu Song, ESCO Consultant Co., Ltd.#include <vcl.h>
#include "MkShape.hpp"

#include "MkLine.hpp"
#include "MkCircle.hpp"
#include "MkArc.hpp"
#include "MkCube.hpp"
#include "MkCylinder.hpp"
#include "MkRect.hpp"
#include "MkSphere.hpp"
#include "MkTriangle.hpp"

bool MkShape::operator==(MkShape &ms)
{
  if (isTriangle() && ms.isTriangle())
  {
    MkTriangle *ta, *tb;
    ta = (MkTriangle *)this;
    tb = (MkTriangle *)&ms;
    return *ta == *tb;
  }
  else if (isRect() && ms.isRect())
  {
    MkRect *ra, *rb;
    ra = (MkRect *)this;
    rb = (MkRect *)&ms;
    return *ra == *rb;
  }
  else if (isCube() && ms.isCube())
  {
    MkCube *ca, *cb;
    ca = (MkCube *)this;
    cb = (MkCube *)&ms;
    return *ca == *cb;
  }
  else if (isCircle() && ms.isCircle())
  {
    MkCircle *ca, *cb;
    ca = (MkCircle *)this;
    cb = (MkCircle *)&ms;
    return *ca == *cb;
  }
  else if (isSphere() && ms.isSphere())
  {
    MkSphere *sa, *sb;
    sa = (MkSphere *)this;
    sb = (MkSphere *)&ms;
    return *sa == *sb;
  }
  else if (isCylinder() && ms.isCylinder())
  {
    MkCylinder *ca, *cb;
    ca = (MkCylinder *)this;
    cb = (MkCylinder *)&ms;
    return *ca == *cb;
  }

  MkDebug("MkShape::operator==() is called.\n");
  return false;
}

bool MkShape::operator!=(MkShape &ms)
{
  if (isTriangle() && ms.isTriangle())
  {
    MkTriangle *ta, *tb;
    ta = (MkTriangle *)this;
    tb = (MkTriangle *)&ms;
    return *ta != *tb;
  }
  else if (isRect() && ms.isRect())
  {
    MkRect *ra, *rb;
    ra = (MkRect *)this;
    rb = (MkRect *)&ms;
    return *ra != *rb;
  }
  else if (isCube() && ms.isCube())
  {
    MkCube *ca, *cb;
    ca = (MkCube *)this;
    cb = (MkCube *)&ms;
    return *ca != *cb;
  }
  else if (isCircle() && ms.isCircle())
  {
    MkCircle *ca, *cb;
    ca = (MkCircle *)this;
    cb = (MkCircle *)&ms;
    return *ca != *cb;
  }
  else if (isSphere() && ms.isSphere())
  {
    MkSphere *sa, *sb;
    sa = (MkSphere *)this;
    sb = (MkSphere *)&ms;
    return *sa != *sb;
  }
  else if (isCylinder() && ms.isCylinder())
  {
    MkCylinder *ca, *cb;
    ca = (MkCylinder *)this;
    cb = (MkCylinder *)&ms;
    return *ca != *cb;
  }

  MkDebug("MkShape::operator!=() is called.\n");
  return true;
}

MkShapeList::MkShapeList()
{
  FirstMyShapes = NULL;
  CurrentMyShapes = NULL;
  LastMyShapes = NULL;
  NumberOfShape = 0;
}

MkShapeList::MkShapeList(MkShape *ms)
{
  FirstMyShapes = ms;
  FirstMyShapes->NextShape = NULL;
  FirstMyShapes->PrevShape = NULL;
  CurrentMyShapes = ms;
  CurrentMyShapes->NextShape = NULL;
  CurrentMyShapes->PrevShape = NULL;
  LastMyShapes = ms;
  LastMyShapes->NextShape = NULL;
  LastMyShapes->PrevShape = NULL;

  NumberOfShape = 1;
}

MkShapeList::~MkShapeList()
{
  while (CurrentMyShapes)
    Delete(CurrentMyShapes);
  NumberOfShape = 0;
  FirstMyShapes = NULL;
  CurrentMyShapes = NULL;
  LastMyShapes = NULL;
}

bool MkShapeList::Add(MkShape *ms)
{
  if (LastMyShapes == NULL)
  {
    FirstMyShapes = ms;
    FirstMyShapes->NextShape = NULL;
    FirstMyShapes->PrevShape = NULL;
    CurrentMyShapes = ms;
    CurrentMyShapes->NextShape = NULL;
    CurrentMyShapes->PrevShape = NULL;
    LastMyShapes = ms;
    LastMyShapes->NextShape = NULL;
    LastMyShapes->PrevShape = NULL;
    NumberOfShape++;
    return true;
  }
  else
    LastMyShapes->NextShape = ms;

  ms->NextShape = NULL;
  ms->PrevShape = LastMyShapes;
  LastMyShapes = ms;
  NumberOfShape++;
  return true;
}

bool MkShapeList::Insert(MkShape *) // Not implemented yet, but it will be used
{                                   // when new object inserted current position.
  return true;
}

bool MkShapeList::Delete(MkShape *ms)
{
  MkShape *cur = FirstMyShapes;

  for (int i = 0; i < NumberOfShape; i++, cur = cur->Next())
  {
    if (cur == ms)
      break;
  }
  if (cur == ms)
  {
    if (cur->Prev())
    {
      cur->PrevShape->NextShape = cur->Next();
    }
    else
    {
      FirstMyShapes = cur->Next();
      if (FirstMyShapes)
        FirstMyShapes->PrevShape = NULL;
    }
    if (cur->Next())
      cur->NextShape->PrevShape = cur->Prev();
    else
    {
      LastMyShapes = cur->Prev();
      if (LastMyShapes)
        LastMyShapes->NextShape = NULL;
    }
    if (cur == CurrentMyShapes)
    {
      if (cur->Prev())
        CurrentMyShapes = cur->Prev();
      else if (cur->Next())
        CurrentMyShapes = cur->Next();
      else
        CurrentMyShapes = NULL;
    }

    delete ms;
    ms = NULL;
    NumberOfShape--;
    return true;
  }
  else
    return false;
}

bool MkShapeList::Clear()
{
  bool flag = true;
  MkShape *cur = LastMyShapes;

  if (!cur)
    return true;

  while (flag && LastMyShapes != FirstMyShapes)
  {
    flag = Delete(cur);
    cur = LastMyShapes;
  }
  if (flag)
  {
    Delete(FirstMyShapes);
    FirstMyShapes = NULL;
    CurrentMyShapes = NULL;
    LastMyShapes = NULL;
    NumberOfShape = 0;
    return true;
  }
  else
    return false;
}

double MkShapeList::GetArea()
{
  return 0; // ���ľ� �Ѵ�.
}

#ifdef __BCPLUSPLUS__
void MkShapeList::Draw(TObject *Sender)
{
  if (String(Sender->ClassName()) == String("MkPaintBox"))
  {
    MkPaintBox *pb = (MkPaintBox *)Sender;
    for (int i = 0; i < GetNumberOfShape(); i++)
    {
      (*this)[i].Draw(pb);
    }
  }
}
#endif

MkShape &MkShapeList::operator[](int i)
{
  MkShape *ms;
  if (i > NumberOfShape)
    throw;
  ms = FirstMyShapes;
  for (int j = 0; j < i; j++)
    ms = ms->NextShape;
  return *ms;
}

MkShapeList &MkShapeList::operator=(MkShapeList &sl)
{
  //    MkShape *ms;
  this->~MkShapeList();
  for (int i = 0; i < sl.NumberOfShape; i++)
  {
    if (sl[i].ClassName() == "MkLine")
    {
      MkLine *ms = new MkLine;
      *ms = (MkLine &)sl[i];
      Add(ms);
    }
    else if (sl[i].ClassName() == "MkArc")
    {
      MkArc *ms = new MkArc;
      *ms = (MkArc &)sl[i];
      Add(ms);
    }
    else if (sl[i].ClassName() == "MkCircle")
    {
      MkCircle *ms = new MkCircle;
      *ms = (MkCircle &)sl[i];
      Add(ms);
    }
    else
    {
      MkDebug("Error MkShapeList operator=");
      MkShape *ms = new MkShape;
      *ms = sl[i];
      Add(ms);
    }
  }
  return (*this);
}

//---------------------------------------------------------------------------
