//---------------------------------------------------------------------------
#include "MkRect.hpp"
MkRect NullRect(0);

MkRect::MkRect()
{
  Origin = NullPoint;
  Width = 0;
  Height = 0;
  className = "MkRect";
}

MkRect::MkRect(int)
{
  Origin = NullPoint;
  Width = 0;
  Height = 0;
  className = "MkRect";
}

MkRect::MkRect(MkPoint &pnt)
{
  Origin = pnt;
  Width = 0;
  Height = 0;
  className = "MkRect";
}

MkRect::MkRect(float x, float y)
{
  Origin.X = x;
  Origin.Y = y;
  Origin.Z = 0;
  Width = 0;
  Height = 0;
  className = "MkRect";
}

MkPoints &MkRect::GetIntPoints(MkLine &line)
{
  static MkPoints pnts;
  MkPoint p;
  MkLine l;

  pnts.Clear();
  for (int i = 0; i < 4; i++)
  {
    l = (*this)(i);
    l[0].Z = 0;
    l[1].Z = 0;
    l.SetFiniteness(true);
    if (line && l)
    {
      p = line.GetIntPoint(l);
      pnts.Add(p);
    }
  }
  return pnts;
}

bool MkRect::IsIn(MkPoint &pnt)
{
  return min(Origin.X, Origin.X + Width) < pnt.X && pnt.X < max(Origin.X, Origin.X + Width) && min(Origin.Y, Origin.Y + Height) < pnt.Y && pnt.Y < max(Origin.Y, Origin.Y + Height);
}

bool MkRect::IsInSurface(MkPoint &pnt, float thick)
{
  float d = thick + 100;
  MkLine l;
  l.SetLine(Origin.X, Origin.Y, Origin.X + Width, Origin.Y);
  d = l.CalDist(pnt);
  l.SetLine(Origin.X, Origin.Y, Origin.X, Origin.Y + Height);
  d = min(d, l.CalDist(pnt));
  l.SetLine(Origin.X, Origin.Y + Height, Origin.X + Width, Origin.Y + Height);
  d = l.CalDist(pnt);
  l.SetLine(Origin.X + Width, Origin.Y, Origin.X + Width, Origin.Y + Height);
  d = min(d, l.CalDist(pnt));
  if (d < thick)
    return true;
  else
    return false;
}

bool MkRect::IsInSpace(MkPoint &pnt)
{
  return IsIn(pnt);
}

void MkRect::GetCross(MkLine &l, MkPoints &p)
{
  MkLine line[4];
  MkPoint pnt[4];
  int i, cnt;
  pnt[0].SetPoint(Origin.X, Origin.Y, Origin.Z);
  pnt[1].SetPoint(Origin.X + Width, Origin.Y, Origin.Z);
  pnt[2].SetPoint(Origin.X + Width, Origin.Y + Height, Origin.Z);
  pnt[3].SetPoint(Origin.X, Origin.Y + Height, Origin.Z);
  line[0].SetLine(pnt[0], pnt[1]);
  line[1].SetLine(pnt[1], pnt[2]);
  line[2].SetLine(pnt[2], pnt[3]);
  line[3].SetLine(pnt[3], pnt[0]);

  cnt = 0;
  for (i = 0; i < 4; i++)
  {
    if (l.IsIntersect(line[i]))
      cnt++;
  }

  p.Initialize(cnt);

  cnt = 0;
  for (i = 0; i < 4; i++)
  {
    if (l.IsIntersect(line[i]))
    {
      p[cnt] = l & line[i];
      cnt++;
    }
  }
}

MkPoint &MkRect::operator[](int i)
{
  static MkPoint p;
  if (i < 0 || i > 4)
    return NullPoint;

  p.SetPoint(Origin.X + Width * (i == 1 || i == 2 ? 1 : 0), Origin.Y + Height * (i == 2 || i == 3 ? 1 : 0), 0);

  return p;
}

MkLine &MkRect::operator()(int i)
{
  static MkLine l;
  if (i < 0 || i > 4)
    return NullLine;

  l.SetLine((*this)[i], (*this)[(i == 3) ? 0 : i + 1]);

  l.SetFiniteness(true);
  return l;
}

#if defined(_MSC_VER) && defined(_WINDOWS_)
void MkRect::Draw(MkPaint *pb)
{
  MkColor c;
  for (int i = 0; i < 4; i++)
    (*this)(i).Draw(pb);
}
#endif

//---------------------------------------------------------------------------
// MkRects::MkRects()
// {
//   FSize = 0;
//   FRect = (MkRect *)NULL;
// }

// MkRects::MkRects(int size)
// {
//   if (size <= 0)
//   {
//     FSize = 0;
//     FRect = (MkRect *)NULL;
//   }
//   FSize = size;
//   FRect = new MkRect[FSize];
//   if (!FRect)
//   {
//     FSize = 0;
//     return;
//   }
// }

// MkRects::~MkRects()
// {
//   if (FRect)
//   {
//     delete[](MkRect *) FRect;
//     FSize = 0;
//   }
// }

// bool MkRects::Initialize(int size)
// {
//   if (size <= 0)
//   {
//     FSize = 0;
//     FRect = (MkRect *)NULL;
//     return false;
//   }

//   if (!FRect)
//   {
//     delete[](MkRect *) FRect;
//     FRect = (MkRect *)NULL;
//   }

//   FSize = size;
//   FRect = new MkRect[FSize];
//   if (!FRect)
//   {
//     FSize = 0;
//     return false;
//   }
//   return true;
// }

// bool MkRects::Initialize(int size, MkRect *cube)
// {
//   if (size <= 0)
//   {
//     FSize = 0;
//     FRect = (MkRect *)NULL;
//     return false;
//   }

//   if (!FRect)
//   {
//     delete[](MkRect *) FRect;
//     FRect = (MkRect *)NULL;
//   }

//   FSize = size;
//   FRect = new MkRect[FSize];
//   if (!FRect)
//   {
//     FSize = 0;
//     return false;
//   }
//   for (int i = 0; i < FSize; i++)
//     FRect[i] = cube[i];

//   return true;
// }

// MkRect &MkRects::operator()(int i)
// {
//   if (i >= 0 && i < FSize)
//     return FRect[i];
//   else
//     return NullRect;
// }

// MkRect &MkRects::operator[](int i)
// {
//   if (i >= 0 && i < FSize)
//     return FRect[i];
//   else
//     return NullRect;
// }

// MkRects &MkRects::operator=(MkRects &a)
// {
//   Initialize(a.GetSize());
//   for (int i = 0; i < a.GetSize(); i++)
//     FRect[i] = a[i];
//   return *this;
// }

// #if defined(_MSC_VER) && defined(_WINDOWS_)
// void MkRects::Draw(MkPaint *pb)
// {
//   for (int i = 0; i < FSize; i++)
//     FRect[i].Draw(pb);
// }
// #endif

//---------------------------------------------------------------------------
