// All the member functions are assumed to be closed polygon!!!
//---------------------------------------------------------------------------
// This module is general purposed simple graphic class to store, draw,
// manipulate object. It is well suited to VCL component, but not restricted.
// It forms the base for the higher level class, such as tunnel component.
//
// Copyright (c) 1999 Myung Kyu Song, ESCO Consultant Co., Ltd.
#include "MkPolygon.hpp"

MkPolygon NullPolygon(0);

MkPolygon::MkPolygon(int size, MkPoint *rp) : MkPoints(size, rp)
{
  if (rp == NULL)
  {
    PointState.Initialize(size);
    Closeness = false;
    Convexity = false;
    Crossing = false;
    isLengthChanged = true;
    isAreaChanged = true;
    isCrossingChanged = true;
    isClosenessChanged = true;
    isConvexityChanged = true;
  }
  else
  {
    PointState.Initialize(size);
    CheckConvexity();
    CheckCrossing();
    CheckFullness();
    isLengthChanged = true;
    isAreaChanged = true;
  }
  className = "MkPolygon";
}

MkPolygon::MkPolygon(int size) : MkPoints(size)
{
  if (size > 0)
  {
    PointState.Initialize(size);
  }
  Closeness = false;
  Convexity = false;
  Crossing = false;
  isLengthChanged = true;
  isAreaChanged = true;
  isCrossingChanged = true;
  isClosenessChanged = true;
  isConvexityChanged = true;
  className = "MkPolygon";
}

MkPolygon::MkPolygon() : MkPoints()
{
  Closeness = false;
  Convexity = false;
  Crossing = false;
  isLengthChanged = true;
  isAreaChanged = true;
  isCrossingChanged = true;
  isClosenessChanged = true;
  isConvexityChanged = true;
  className = "MkPolygon";
}

MkPolygon::~MkPolygon()
{
}

void MkPolygon::Initialize(int size)
{
  MkPoints::Initialize(size);
  PointState.Initialize(size);
}

void MkPolygon::Initialize(int size, MkPoint *rps)
{
  MkPoints::Initialize(size, rps);
  PointState.Initialize(size);
}

double MkPolygon::CalLength()
{
  if (!isLengthChanged)
    return FLength;
  else
  {
    double TotLen = 0;
    for (int i = 0; i < GetSize(); i++)
      TotLen += (*this)(i).GetLength();
    isLengthChanged = false;
    return FLength = TotLen;
  }
}

double MkPolygon::CalArea()
{
  if (GetSize() < 3)
  {
    FArea = 0;
    return 0;
  }
  if (!isAreaChanged)
    return FArea;
  else
  {
    double TotArea = 0;
    MkPoint rp[3];
    MkTriangle tri;
    for (int i = 0; i < GetSize(); i++)
    {
      PointState(i) = 1;
    }

    SetCurrent(0);
    rp[0] = (*this)[0];
    rp[1] = (*this)[1];
    rp[2] = (*this)[2];

    tri.Reset(rp);
    for (; GetAlivePoint() > 3;)
    {
      if (tri.isValid())
      {
        TotArea += tri.GetArea();
        PointState(AlivedNext(Current())) = 0;
        rp[0] = CurrentPoint();
        rp[1] = AlivedNextPoint(Current());
        rp[2] = AlivedNextPoint(AlivedNext(Current()));
        tri.Reset(rp);
      }
      else if (fabs(tri.CrossProduct()) < 1.0e-6)
      {
        PointState(AlivedNext(Current())) = 0;
        rp[0] = CurrentPoint();
        rp[1] = AlivedNextPoint(Current());
        rp[2] = AlivedNextPoint(AlivedNext(Current()));
        tri.Reset(rp);
      }
      else
      {
        AlivedNext();
        rp[0] = CurrentPoint();
        rp[1] = AlivedNextPoint(Current());
        rp[2] = AlivedNextPoint(AlivedNext(Current()));
        tri.Reset(rp);
      }
    }
    if (GetAlivePoint() == 3)
    {
      rp[0] = CurrentPoint();
      rp[1] = AlivedNextPoint(Current());
      rp[2] = AlivedNextPoint(AlivedNext(Current()));
      tri.Reset(rp);
      TotArea += tri.GetArea();
      PointState(Current()) = 0;
      PointState(AlivedNext(Current())) = 0;
      PointState(AlivedNext(AlivedNext(Current()))) = 0;
    }
    FArea = TotArea;
    if (GetAlivePoint() != 0)
      MkDebug("Error occured when calc polygon area");
  }
  isAreaChanged = false;
  return FArea;
}

double MkPolygon::CalArea2()
{
  if (GetSize() < 3)
  {
    FArea = 0;
    return 0;
  }
  if (!isAreaChanged)
    return FArea;
  else
  {
    double TotArea = 0;
    for (int i = 0; i < GetSize() - 1; i++)
    {
      TotArea = TotArea + ((*this)[i].X - (*this)[i + 1].X) * ((*this)[i].Y + (*this)[i + 1].Y) / 2;
    }
    TotArea = TotArea + ((*this)[GetSize() - 1].X - (*this)[0].X) * ((*this)[GetSize() - 1].Y + (*this)[0].Y) / 2;
    FArea = TotArea;
  }
  isAreaChanged = false;
  return FArea;
}

double MkPolygon::GetLength()
{
  return isLengthChanged ? CalLength() : FLength;
}

double MkPolygon::GetArea()
{
  return isAreaChanged ? CalArea2() : FArea;
}

void MkPolygon::CheckCrossing()
{
  int i, j;
  Crossing = false;
  if (isFilled())
  {
    for (i = 0; i < GetSize() - 1; i++)
    {
      for (j = i; j < GetSize(); j++)
      {
        Crossing = Crossing && ((*this)(i) && (*this)(j));
      }
    }
  }
}

bool MkPolygon::isCrossWith(MkLine &l)
{
  int i;
  bool flag = false;
  for (i = 0; i < FSize - 1; i++)
  {
    flag = flag || (l && (*this)(i));
  }
  return flag;
}

bool MkPolygon::isCrossWithX(double x)
{
  int i;
  bool flag = false;
  for (i = 0; i < FSize - 1; i++)
  {
    if (min(FPoint[i].X, FPoint[i + 1].X) < x &&
        x < max(FPoint[i].X, FPoint[i + 1].X))
      return true;
  }
  return false;
}

bool MkPolygon::isCrossWithY(double y)
{
  int i;
  bool flag = false;
  for (i = 0; i < FSize - 1; i++)
  {
    if (min(FPoint[i].Y, FPoint[i + 1].Y) < y &&
        y < max(FPoint[i].Y, FPoint[i + 1].Y))
      return true;
  }
  return false;
}

void MkPolygon::getCrossWith(MkLine &l, MkPoints &pnts)
{
  int i, cnt = 0;
  for (i = 0; i < FSize - (Closeness ? 0 : 1); i++)
  {
    MkLine line(FPoint[i], FPoint[(i + 1 == FSize) ? 0 : i + 1]);
    line.SetFiniteness(true);
    if (l && line)
      cnt++;
  }

  if (cnt == 0)
  {
    pnts.Clear();
    return;
  }

  pnts.Initialize(cnt);
  cnt = 0;

  for (i = 0; i < FSize - (Closeness ? 0 : 1); i++)
  {
    MkLine line(FPoint[i], FPoint[(i + 1 == FSize) ? 0 : i + 1]);
    line.SetFiniteness(true);
    if (l && line)
    {
      pnts[cnt] = l & line;
      cnt++;
    }
  }
}

void MkPolygon::getCrossWithX(double x, MkPoints &pnts)
{
  int i, cnt = 0;
  for (i = 0; i < FSize - (Closeness ? 0 : 1); i++)
  {
    MkPoint &a = FPoint[i];
    MkPoint &b = FPoint[(i + 1 == FSize) ? 0 : i + 1];
    if (min(a.X, b.X) < x && x < max(a.X, b.X))
      cnt++;
  }

  if (cnt == 0)
  {
    pnts.Clear();
    return;
  }

  pnts.Initialize(cnt);
  cnt = 0;

  for (i = 0; i < FSize; i++)
  {
    MkPoint &a = FPoint[i];
    MkPoint &b = FPoint[(i + 1 == FSize) ? 0 : i + 1];
    if (min(a.X, b.X) < x && x < max(a.X, b.X))
    {
      pnts[cnt].X = x;
      pnts[cnt].Y = a.Y + (b.Y - a.Y) / (b.X - a.X) * (x - a.X);
      cnt++;
    }
  }
}

void MkPolygon::getCrossWithY(double y, MkPoints &pnts)
{
  int i, cnt = 0;
  for (i = 0; i < FSize; i++)
  {
    MkPoint &a = FPoint[i];
    MkPoint &b = FPoint[(i + 1 == FSize) ? 0 : i + 1];
    if (min(a.Y, b.Y) < y && y < max(a.Y, b.Y))
      cnt++;
  }

  if (cnt == 0)
  {
    pnts.Clear();
    return;
  }

  pnts.Initialize(cnt);
  cnt = 0;

  for (i = 0; i < FSize; i++)
  {
    MkPoint &a = FPoint[i];
    MkPoint &b = FPoint[(i + 1 == FSize) ? 0 : i + 1];
    if (min(a.Y, b.Y) < y && y < max(a.Y, b.Y))
    {
      pnts[cnt].Y = y;
      pnts[cnt].X = a.X + (b.X - a.X) / (b.Y - a.Y) * (y - a.Y);
      cnt++;
    }
  }
}

bool MkPolygon::IsIn(MkPoint rp) // valid only for convexity
{
  static MkVector<double> v1(3), v2(3), c1(3), c2(3); //,c3(3),c4(3);
  MkLine rl1, rl2;

  rl1.SetLine(rp, (*this)[0]);
  rl2.SetLine(rp, (*this)[1]);
  v1[0] = rl1.GetL();
  v1[1] = rl1.GetM();
  v1[2] = rl1.GetN();
  v2[0] = rl2.GetL();
  v2[1] = rl2.GetM();
  v2[2] = rl2.GetN();
  v1.Cross(v2, c1);

  for (int i = 1; i < FSize; i++)
  {
    rl1.SetLine(rp, (*this)[i]);
    rl2.SetLine(rp, (*this)[(i + 1 < FSize) ? i + 1 : 0]);
    v1[0] = rl1.GetL();
    v1[1] = rl1.GetM();
    v1[2] = rl1.GetN();
    v2[0] = rl2.GetL();
    v2[1] = rl2.GetM();
    v2[2] = rl2.GetN();
    v1.Cross(v2, c2);
    if (c1 * c2 < 0)
      return false;
    c1 = c2;
  }
  return true;
}

bool MkPolygon::IsIn(MkCircle &rc)
{
  double min_dist;
  MkLine rl;
  if (!IsIn(rc.GetCenter()))
    return false;

  rl = (*this)(0);
  min_dist = rl.CalDist(rc.GetCenter());
  for (int i = 1; i < FSize - 1; i++)
  {
    rl = (*this)(i);
    min_dist = rl.CalDist(rc.GetCenter()) < min_dist ? rl.CalDist(rc.GetCenter()) : min_dist;
  }
  return rc.GetRadius() < min_dist;
}

void MkPolygon::CheckConvexity()
{
  Convexity = false;
}

void MkPolygon::CheckFullness()
{
  Fullness = false;
}

void MkPolygon::Offset(MkPolygon &poly, bool dir, double space) //dir : left=false, right=true
{
  int i;
  bool flag = true;
  MkPoint pnt;
  MkVector<double> vec[2], norm[2];

  poly.Clear();

  if (FSize < 3)
    return;

  vec[0].SetVector(FPoint[1].X - FPoint[0].X, FPoint[1].Y - FPoint[0].Y, FPoint[1].Z - FPoint[0].Z);
  vec[1].SetVector(FPoint[2].X - FPoint[1].X, FPoint[2].Y - FPoint[1].Y, FPoint[2].Z - FPoint[1].Z);
  vec[0].Normalize();
  vec[1].Normalize();

  vec[0].Cross(vec[1], norm[0]);
  norm[0] = norm[0] * 100; // for the precision of normal vector
  for (i = 1; i < FSize - 2; i++)
  {
    vec[0].SetVector(FPoint[i + 1].X - FPoint[i].X, FPoint[i + 1].Y - FPoint[i].Y, FPoint[i + 1].Z - FPoint[i].Z);
    vec[1].SetVector(FPoint[i + 2].X - FPoint[i + 1].X, FPoint[i + 2].Y - FPoint[i + 1].Y, FPoint[i + 2].Z - FPoint[i + 1].Z);
    vec[0].Normalize();
    vec[1].Normalize();
    vec[0].Cross(vec[1], norm[1]);
    norm[1] = norm[1] * 100; // for the precision of normal vector

    flag = flag && (norm[0] == norm[1]);
    if (!flag)
      return;

    norm[0] = norm[1];
  }
  norm[0].Normalize();
  Offset(poly, norm[0], dir, space);
}

void MkPolygon::Offset(MkPolygon &poly, MkVector<double> &up, bool dir /*left=false, right=true */, double space)
{
  int i;
  MkPoint pnt[2], p;
  MkLines lines;
  MkLine l;
  static MkVector<double> v(3), d(3);
  poly.Clear();

  lines.Initialize(FSize - 1);
  for (i = 0; i < FSize - 1; i++)
  {
    pnt[0] = FPoint[i];
    pnt[1] = FPoint[i + 1];
    v[0] = pnt[1].X - pnt[0].X;
    v[1] = pnt[1].Y - pnt[0].Y;
    v[2] = pnt[1].Z - pnt[0].Z;
    v.Normalize();

    ((dir) ? v.Cross(up, d) : up.Cross(v, d));

    pnt[0].SetPoint(FPoint[i].X + space * d[0],
                    FPoint[i].Y + space * d[1],
                    FPoint[i].Z + space * d[2]);

    pnt[1].SetPoint(FPoint[i + 1].X + space * d[0],
                    FPoint[i + 1].Y + space * d[1],
                    FPoint[i + 1].Z + space * d[2]);
    lines[i].SetLine(pnt[0], pnt[1]);
    lines[i].SetFiniteness(false);
  }
  for (i = 0; i < FSize - 2; i++)
  {
    if (lines[i][1] == lines[i + 1][0])
      ;
    else if (lines[i] && lines[i + 1])
    {
      p = lines[i].GetIntPoint(lines[i + 1]);
      //      p = lines[i] & lines[i+1];
      lines[i][1] = p;
      lines[i + 1][0] = p;
    }
    else
    {
      p = (lines[i][1] + lines[i + 1][0]) / 2;
      lines[i][1] = p;
      lines[i + 1][0] = p;
    }
  }
  poly.Initialize(FSize);
  for (i = 0; i < FSize - 1; i++)
  {
    poly[i] = lines[i][0];
  }
  poly[i] = lines[i - 1][1];
}

MkPoint MkPolygon::Measure(double dist)
{
  if (dist > GetLength())
    return NullPoint;
  if (FSize < 2)
    return NullPoint;

  double Len[2] = {0, 0};
  Len[0] = (*this)(0).GetLength();
  if (dist <= Len[0])
  {
    return (*this)(0).GetDivision(dist / Len[0]);
  }

  Len[1] = Len[0];
  for (int i = 1; i < GetSize(); i++)
  {
    Len[1] += (*this)(i).GetLength();
    if (Len[0] < dist && dist <= Len[1] && fabs(Len[1] - Len[0]) > EPS)
      return (*this)(i).GetDivision((dist - Len[0]) / (Len[1] - Len[0]));
    Len[0] = Len[1];
  }
  return NullPoint;
}

double MkPolygon::Measure(MkPoint pnt)
{
  int i, close;
  double len;
  bool flag = false;

  close = Closeness ? 0 : 1;

  for (i = 0; i < FSize - close; i++)
    flag = flag || (*this)(i).IsInLine(pnt);

  if (fabs(FLength) < EPS)
  {
    CalLength();
    if (fabs(FLength) < EPS)
      return -1.0;
  }
  if (!flag)
    return -1.0;

  len = 0;
  for (i = 0; i < FSize - close; i++)
  {
    if ((*this)(i).IsInLine(pnt))
    {
      len += CalDist(FPoint[i], pnt);
      break;
    }
    len += (*this)(i).GetLength();
  }
  return len;
}

bool MkPolygon::InverseDirection()
{
  int i;
  if (FSize <= 0)
    return false;
  for (i = 0; i < FSize / 2; i++)
    Swap(i, FSize - i - 1);
  return true;
}

bool MkPolygon::AddInBetween(MkPoints &pnts)
{
  int i, j, k, n, np = 0;
  MkPolygon &a = (*this);
  MkPoints Points; // new points replace current points
  MkInt lnum;

  lnum.Initialize(pnts.GetSize());
  for (i = 0; i < pnts.GetSize(); i++)
  {
    lnum(i) = -1;
  }

  for (i = 0; i < pnts.GetSize(); i++)
  {
    for (j = 0; j < FSize; j++)
      if (a(j).IsInLine(pnts[i]) && a[j] != pnts[i] && a[(j + 1 == FSize) ? 0 : j + 1] != pnts[i])
        lnum(i) = j;
  }

  np = FSize;
  for (i = 0; i < pnts.GetSize(); i++)
    if (lnum(i) != -1)
    {
      np++;
    }

  Points.Initialize(np); // initial points + intersection points on polygon

  np = 0;
  for (j = 0; j < FSize; j++)
  {
    n = 0;
    for (i = 0; i < pnts.GetSize(); i++)
    {
      if (lnum(i) == j)
        n++;
    }
    if (n == 0)
    {
      Points[np] = FPoint[j == FSize ? 0 : j];
      np++;
      continue;
    }
    MkInt pnum(n);
    MkDouble dnum(n);

    n = 0;
    for (i = 0; i < pnts.GetSize(); i++)
    {
      if (lnum(i) == j)
      {
        pnum(n) = i;
        n++;
      }
    }
    for (i = 0; i < n; i++)
    {
      dnum(i) = CalDist(pnts[pnum(i)], FPoint[j]);
    }

    for (i = 0; i < n - 1; i++)
      for (k = i + 1; k < n; k++)
        if (dnum(i) > dnum(k))
        {
          swap(dnum(i), dnum(k));
          swap(pnum(i), pnum(k));
        }

    Points[np] = FPoint[j == FSize ? 0 : j];
    np++;
    if (np >= Points.GetSize())
      return false;
    for (i = 0; i < n; i++)
    {
      Points[np] = pnts[pnum(i)];
      np++;
      if (np >= Points.GetSize())
      {
        Initialize(Points.GetSize(), Points.GetPoints());
        return false;
      }
    }
  }
  Points[np] = FPoint[FSize - 1];
  Initialize(Points.GetSize(), Points.GetPoints());
  return true;
}

bool MkPolygon::Merge(MkPolygon &poly)
{
  int i, size, nSame = 0;
  bool isSame1, isSame2;

  if (poly.GetSize() <= 0)
    return false;
  isSame1 = (FPoint[0] == poly[poly.GetSize() - 1]);
  isSame2 = (FPoint[FSize - 1] == poly[0]);

  size = FSize + poly.GetSize() - (isSame1 ? 1 : 0) - (isSame2 ? 1 : 0);

  MkPoint *pnt = new MkPoint[size];
  if (!pnt)
    return false;
  for (i = 0; i < FSize; i++)
  {
    pnt[i] = FPoint[i];
  }

  for (i = FSize; i < size; i++)
    pnt[i] = poly[i - FSize + (isSame2 ? 1 : 0)];

  Initialize(size, pnt);
  return true;
}

MkVector<double> &MkPolygon::GetVector()
{
  int i;
  static MkVector<double> Normal(3);
  MkVector<double> avec(3), bvec(3);

  Normal[0] = 0;
  Normal[1] = 0;
  Normal[2] = 0;
  avec = (*this)(0).GetVector();

  for (i = 1; i < FSize; i++)
  {
    bvec = (*this)(i).GetVector();
    if (avec != bvec)
      break;
  }
  if (avec == bvec)
    return Normal;

  avec.Cross(bvec, Normal);
  return Normal;
}

bool MkPolygon::FindInter(MkPolygon &b, MkPoints &pnts)
{
  int i, j, npnt = 0;
  MkDouble dist;
  MkPoint pnt;
  MkPolygon &a = *this;

  for (i = 0; i < a.GetSize(); i++)
  {
    for (j = 0; j < b.GetSize(); j++)
    {
      MkLine al, bl;
      al = a(i);
      al.SetFiniteness(true);
      bl = b(j);
      bl.SetFiniteness(true);
      if (al && bl)
        npnt++;
    }
  }
  if (npnt <= 0 || npnt % 2 != 0)
  {
    return false;
  }

  dist.Initialize(npnt);
  pnts.Initialize(npnt);

  npnt = 0;
  for (i = 0; i < a.GetSize(); i++)
  {
    for (j = 0; j < b.GetSize(); j++)
    {
      if (a(i) && b(j))
      {
        MkLine l = b(j);
        pnts[npnt] = a(i) & l;
        npnt++;
      }
    }
  }

  for (i = 0; i < pnts.GetSize(); i++)
  {
    dist(i) = a.Measure(pnts[i]);
  }

  for (i = 0; i < pnts.GetSize() - 1; i++)
  {
    for (j = i + 1; j < pnts.GetSize(); j++)
    {
      if (dist(i) > dist(j))
      {
        ::Swap(pnts[i], pnts[j]);
        swap(dist(i), dist(j));
      }
    }
  }

  pnt = pnts[0];
  for (i = 1; i < pnts.GetSize(); i++)
  {
    pnts[i - 1] = pnts[i];
  }
  pnts[pnts.GetSize() - 1] = pnt;

  return true;
}

int MkPolygon::FindPoly(MkPolygon &b, MkBoolType type) // find number of polygon
{
  int i, j, ai, bi, npoly = 0;
  bool flag, dir;
  MkPolygon &a = *this;
  MkPolygon at; // a's points + intersection points
  MkPolygon bt; // b's points + intersection points

  MkPoints pnts;
  MkPoint pnt, ap, bp, sp;
  MkVector<double> avec(3), bvec(3), cvec(3);
  MkInt pass;

  if (a.GetSize() < 3 || b.GetSize() < 3)
    return 0;

  avec = a.GetVector();
  bvec = b.GetVector();

  if (avec.GetLength() < EPS || bvec.GetLength() < EPS)
    return 0;

  avec.Normalize();
  bvec.Normalize();

  if (fabs(avec[0] - bvec[0]) > EPS * 1000 ||
      fabs(avec[1] - bvec[1]) > EPS * 1000 ||
      fabs(avec[2] - bvec[2]) > EPS * 1000)
    return 0;

  if (!FindInter(b, pnts))
    return 0;

  pass.Initialize(pnts.GetSize());

  at = a;
  bt = b;
  at.AddInBetween(pnts); // add points between line segment.
  bt.AddInBetween(pnts); // if point is not on any line, point will be ignored
  if (type == btSub)
    bt.InverseDirection();

  flag = true;
  pnt = pnts[0];
  sp = pnt;
  pass(0) = 1;

  while (flag)
  {
    ap = at.NextPoint(sp);
    bp = bt.NextPoint(sp);
    avec.SetVector(ap.X - sp.X, ap.Y - sp.Y, ap.Z - sp.Z);
    bvec.SetVector(bp.X - sp.X, bp.Y - sp.Y, bp.Z - sp.Z);
    avec.Cross(bvec, cvec);

    if (type == btUni)
      dir = cvec * a.GetVector() > 0;
    else if (type == btInt)
      dir = cvec * a.GetVector() < 0;
    else if (type == btSub)
      dir = cvec * a.GetVector() < 0;

    if (dir)
    {
      while (!pnts.hasPoint(ap))
      {
        ap = at.NextPoint(ap);
      }
      sp = ap;
      if (pnts.hasPoint(sp))
        pass(pnts.numPoint(sp)) = 1;
      if (pnt == sp)
      {
        npoly++;
        for (i = 0; i < pnts.GetSize(); i++)
          if (pass(i) == 0)
            break;
        if (i >= pnts.GetSize())
          return npoly;
        pnt = pnts[i];
        sp = pnt;
      }
    }
    else
    {
      while (!pnts.hasPoint(bp))
      {
        bp = bt.NextPoint(bp);
      }
      sp = bp;
      if (pnts.hasPoint(sp))
        pass(pnts.numPoint(sp)) = 1;
      if (pnt == sp)
      {
        npoly++;
        for (i = 0; i < pnts.GetSize(); i++)
          if (pass(i) == 0)
            break;
        if (i >= pnts.GetSize())
          return npoly;
        pnt = pnts[i];
        sp = pnt;
      }
    }
  }
  return 0;
}

bool MkPolygon::BoolSub(MkPolygon &b, MkPolygons &c) // subtraction
{
  int i, j, ai, bi, npnt = 0, npoly = 0, n /*n of pnts */;
  bool flag;
  MkPolygon &a = *this;
  MkPolygon at; // a's points + intersection points
  MkPolygon bt; // b's points + intersection points

  MkPoints pnts, p;
  MkPoint pnt, ap, bp, sp;
  MkVector<double> avec(3), bvec(3), cvec(3);
  MkInt pass;

  c.Clear();
  if (a.GetSize() < 3 || b.GetSize() < 3)
    return false;

  avec = a.GetVector();
  bvec = b.GetVector();

  if (avec.GetLength() < EPS || bvec.GetLength() < EPS)
    return false;

  avec.Normalize();
  bvec.Normalize();

  if (fabs(avec[0] - bvec[0]) > EPS * 1000 ||
      fabs(avec[1] - bvec[1]) > EPS * 1000 ||
      fabs(avec[2] - bvec[2]) > EPS * 1000)
    return false;

  if (!FindInter(b, pnts) && pnts.GetSize() == 0)
  {
    c.Initialize(1);
    c[0] = *this;
    return true;
  }

  if ((npoly = FindPoly(b, btSub)) == 0)
    return false;
  c.Initialize(npoly);

  pass.Initialize(pnts.GetSize());

  at = a;
  bt = b;
  at.AddInBetween(pnts); // add points between line segment.
  bt.AddInBetween(pnts); // if point is not on any line, point will be ignored
  bt.InverseDirection();

  flag = true;
  n = 0;
  pnt = pnts[0];
  sp = pnt;
  pass(0) = 1;
  npoly = 0;

  while (flag)
  {
    c[npoly].Add(sp);
    ap = at.NextPoint(sp);
    bp = bt.NextPoint(sp);
    avec.SetVector(ap.X - sp.X, ap.Y - sp.Y, ap.Z - sp.Z);
    bvec.SetVector(bp.X - sp.X, bp.Y - sp.Y, bp.Z - sp.Z);
    avec.Cross(bvec, cvec);
    if (cvec * a.GetVector() < 0)
    {
      c[npoly].Add(ap);
      while (!pnts.hasPoint(ap))
      {
        ap = at.NextPoint(ap);
        c[npoly].Add(ap);
      }
      sp = ap;
      if (pnts.hasPoint(sp))
        pass(pnts.numPoint(sp)) = 1;
      if (pnt == sp)
      {
        for (i = 0; i < pnts.GetSize(); i++)
          if (pass(i) == 0)
            break;
        if (i >= pnts.GetSize())
        {
          return true;
        }
        pnt = pnts[i];
        sp = pnt;
        npoly++;
      }
    }
    else
    {
      c[npoly].Add(bp);
      while (!pnts.hasPoint(bp))
      {
        bp = bt.NextPoint(bp);
        c[npoly].Add(bp);
      }
      sp = bp;
      if (pnts.hasPoint(sp))
        pass(pnts.numPoint(sp)) = 1;
      if (pnt == sp)
      {
        for (i = 0; i < pnts.GetSize(); i++)
          if (pass(i) == 0)
            break;
        if (i >= pnts.GetSize())
        {
          return true;
        }
        pnt = pnts[i];
        sp = pnt;
        npoly++;
      }
    }
  }
  return true;
}

bool MkPolygon::BoolUni(MkPolygon &b, MkPolygons &c) // union
{
  int i, j, ai, bi, npnt = 0, npoly = 0, n /*n of pnts */;
  bool flag;
  MkPolygon &a = *this;
  MkPolygon at; // a's points + intersection points
  MkPolygon bt; // b's points + intersection points

  MkPoints pnts, p;
  MkPoint pnt, ap, bp, sp;
  MkVector<double> avec(3), bvec(3), cvec(3);
  MkInt pass;

  c.Clear();
  if (a.GetSize() < 3 || b.GetSize() < 3)
    return false;

  avec = a.GetVector();
  bvec = b.GetVector();

  if (avec.GetLength() < EPS || bvec.GetLength() < EPS)
    return false;

  avec.Normalize();
  bvec.Normalize();

  if (fabs(avec[0] - bvec[0]) > EPS * 1000 ||
      fabs(avec[1] - bvec[1]) > EPS * 1000 ||
      fabs(avec[2] - bvec[2]) > EPS * 1000)
    return false;

  if (!FindInter(b, pnts))
    return false;

  if ((npoly = FindPoly(b, btUni)) == 0)
    return false;
  c.Initialize(npoly);

  pass.Initialize(pnts.GetSize());

  at = a;
  bt = b;
  at.AddInBetween(pnts); // add points between line segment.
  bt.AddInBetween(pnts); // if point is not on any line, point will be ignored

  flag = true;
  n = 0;
  pnt = pnts[0];
  sp = pnt;
  pass(0) = 1;
  npoly = 0;

  while (flag)
  {
    c[npoly].Add(sp);
    ap = at.NextPoint(sp);
    bp = bt.NextPoint(sp);
    avec.SetVector(ap.X - sp.X, ap.Y - sp.Y, ap.Z - sp.Z);
    bvec.SetVector(bp.X - sp.X, bp.Y - sp.Y, bp.Z - sp.Z);
    avec.Cross(bvec, cvec);
    if (cvec * a.GetVector() > 0)
    {
      c[npoly].Add(ap);
      while (!pnts.hasPoint(at.NextPoint(ap)))
      {
        ap = at.NextPoint(ap);
        c[npoly].Add(ap);
      }
      sp = at.NextPoint(ap);
      if (pnts.hasPoint(sp))
        pass(pnts.numPoint(sp)) = 1;
      if (pnt == sp)
      {
        for (i = 0; i < pnts.GetSize(); i++)
          if (pass(i) == 0)
            break;
        if (i >= pnts.GetSize())
        {
          return true;
        }
        pnt = pnts[i];
        sp = pnt;
        npoly++;
      }
    }
    else
    {
      c[npoly].Add(bp);
      while (!pnts.hasPoint(bt.NextPoint(bp)))
      {
        bp = bt.NextPoint(bp);
        c[npoly].Add(bp);
      }
      sp = bt.NextPoint(bp);
      if (pnts.hasPoint(sp))
        pass(pnts.numPoint(sp)) = 1;
      if (pnt == sp)
      {
        for (i = 0; i < pnts.GetSize(); i++)
          if (pass(i) == 0)
            break;
        if (i >= pnts.GetSize())
        {
          return true;
        }
        pnt = pnts[i];
        sp = pnt;
        npoly++;
      }
    }
  }
  return true;
}

bool MkPolygon::BoolInt(MkPolygon &b, MkPolygons &c) // intersection
{
  int i, j, ai, bi, npnt = 0, npoly = 0, n /*n of pnts */;
  bool flag;
  MkPolygon &a = *this;
  MkPolygon at; // a's points + intersection points
  MkPolygon bt; // b's points + intersection points

  MkPoints pnts, p;
  MkPoint pnt, ap, bp, sp;
  MkVector<double> avec(3), bvec(3), cvec(3);
  MkInt pass;

  c.Clear();
  if (a.GetSize() < 3 || b.GetSize() < 3)
    return false;

  avec = a.GetVector();
  bvec = b.GetVector();

  if (avec.GetLength() < EPS || bvec.GetLength() < EPS)
    return false;

  avec.Normalize();
  bvec.Normalize();

  if (fabs(avec[0] - bvec[0]) > EPS * 1000 ||
      fabs(avec[1] - bvec[1]) > EPS * 1000 ||
      fabs(avec[2] - bvec[2]) > EPS * 1000)
    return false;

  if (!FindInter(b, pnts))
    return false;

  if ((npoly = FindPoly(b, btInt)) == 0)
    return false;
  c.Initialize(npoly);

  pass.Initialize(pnts.GetSize());

  at = a;
  bt = b;
  at.AddInBetween(pnts); // add points between line segment.
  bt.AddInBetween(pnts); // if point is not on any line, point will be ignored

  flag = true;
  n = 0;
  pnt = pnts[0];
  sp = pnt;
  pass(0) = 1;
  npoly = 0;

  while (flag)
  {
    c[npoly].Add(sp);
    ap = at.NextPoint(sp);
    bp = bt.NextPoint(sp);
    avec.SetVector(ap.X - sp.X, ap.Y - sp.Y, ap.Z - sp.Z);
    bvec.SetVector(bp.X - sp.X, bp.Y - sp.Y, bp.Z - sp.Z);
    avec.Cross(bvec, cvec);
    if (cvec * a.GetVector() < 0)
    {
      c[npoly].Add(ap);
      while (!pnts.hasPoint(at.NextPoint(ap)))
      {
        ap = at.NextPoint(ap);
        c[npoly].Add(ap);
      }
      sp = at.NextPoint(ap);
      if (pnts.hasPoint(sp))
        pass(pnts.numPoint(sp)) = 1;
      if (pnt == sp)
      {
        for (i = 0; i < pnts.GetSize(); i++)
          if (pass(i) == 0)
            break;
        if (i >= pnts.GetSize())
        {
          return true;
        }
        pnt = pnts[i];
        sp = pnt;
        npoly++;
      }
    }
    else
    {
      c[npoly].Add(bp);
      while (!pnts.hasPoint(bt.NextPoint(bp)))
      {
        bp = bt.NextPoint(bp);
        c[npoly].Add(bp);
      }
      sp = bt.NextPoint(bp);
      if (pnts.hasPoint(sp))
        pass(pnts.numPoint(sp)) = 1;
      if (pnt == sp)
      {
        for (i = 0; i < pnts.GetSize(); i++)
          if (pass(i) == 0)
            break;
        if (i >= pnts.GetSize())
        {
          return true;
        }
        pnt = pnts[i];
        sp = pnt;
        npoly++;
      }
    }
  }
  return true;
}

/*
bool MkPolygon::BoolSub(MkPolygon &b, MkPolygons &c)
{
  int i, j, ai, bi, npnt = 0, npoly = 0;
  MkPolygon &a = *this, asub, bsub;
  MkPoints pnts, p;
  MkPoint pnt;
  MkVector<double> avec(3), bvec(3), cvec(3);

  c.Clear();
  if(a.GetSize()<3 || b.GetSize()<3) return false;

  avec = a(0).GetVector() & a(1).GetVector();
  bvec = b(0).GetVector() & b(1).GetVector();
  avec.Normalize();
  bvec.Normalize();

  if (fabs(avec[0]-bvec[0])>EPS*1000 ||
      fabs(avec[1]-bvec[1])>EPS*1000 ||
      fabs(avec[2]-bvec[2])>EPS*1000) return false;

  if(!FindInter(b,pnts) && pnts.GetSize()==0) {
    c.Initialize(1);
    c[0] = *this;
    return true;
  }

  npoly = pnts.GetSize()/2;
  c.Initialize(npoly);

  for(i=0;i<npoly;i++) {
    p.Initialize(2);
    p[0] = pnts[2*i];
    p[1] = pnts[2*i+1];

    asub.Clear();
    bsub.Clear();

    a.Extract(p,asub);
    b.Extract(p,bsub);
    bsub.InverseDirection();

    asub.Merge(bsub);
    asub.SetCloseness(true);

    c[i] = asub;
  }

  a.InverseDirection();
  b.InverseDirection();

  return true;
}

bool MkPolygon::BoolAdd(MkPolygon &b, MkPolygons &c)
{
  int i, j, ai, bi, npnt = 0, npoly = 0;
  MkPolygon &a = *this, asub, bsub;
  MkPoints pnts, p;
  MkPoint pnt;
  MkVector<double> avec(3), bvec(3), cvec(3);

  c.Clear();
  if(a.GetSize()<3 || b.GetSize()<3) return false;

  avec = a(0).GetVector() & a(1).GetVector();
  bvec = b(0).GetVector() & b(1).GetVector();
  avec.Normalize();
  bvec.Normalize();

  if (fabs(avec[0]-bvec[0])>EPS*1000 ||
      fabs(avec[1]-bvec[1])>EPS*1000 ||
      fabs(avec[2]-bvec[2])>EPS*1000) return false;

  if(!FindInter(b,pnts)) return false;

  npoly = (pnts.GetSize()>2)?1:2;
  c.Initialize(npoly);

  if(npoly==2) {
    c[0] = a;
    c[1] = b;
    return true;
  }

  asub.Clear();
  for(i=0;i<pnts.GetSize();i++) {
    p.Initialize(2);
    p[0] = pnts[i];
    p[1] = pnts[(i+1==pnts.GetSize())?0:i+1];

    bsub.Clear();

    if(i%2) a.Extract(p,bsub);
    else b.Extract(p,bsub);

    asub.Merge(bsub);
  }
  c[0] = asub;

  return true;
}

bool MkPolygon::BoolUni(MkPolygon &b, MkPolygons &c)
{
  int i, j, ai, bi, npnt = 0, npoly = 0;
  MkPolygon &a = *this, asub, bsub;
  MkPoints pnts, p;
  MkPoint pnt;
  MkVector<double> avec(3), bvec(3), cvec(3);

  c.Clear();
  if(a.GetSize()<3 || b.GetSize()<3) return false;

  avec = a(0).GetVector() & a(1).GetVector();
  bvec = b(0).GetVector() & b(1).GetVector();
  avec.Normalize();
  bvec.Normalize();

  if (fabs(avec[0]-bvec[0])>EPS*1000 ||
      fabs(avec[1]-bvec[1])>EPS*1000 ||
      fabs(avec[2]-bvec[2])>EPS*1000) return false;

  if(!FindInter(b,pnts)) return false;

  npoly = npnt/2;
  b.InverseDirection();
  c.Initialize(npoly);

  for(i=0;i<npoly;i++) {
    p.Initialize(2);
    p[0] = pnts[2*i];
    p[1] = pnts[2*i+1];

    asub.Clear();
    bsub.Clear();

    b.Extract(p,bsub);
    bsub.InverseDirection();
    a.Extract(p,asub);
    asub.Merge(bsub);

    c[i] = asub;
  }

  b.InverseDirection();

  return true;
}
*/

void MkPolygon::Extract(MkPoints &b, MkPolygon &c)
{
  bool found = false;
  int i, cnt = 0, index[2];
  double meas[2];
  MkPolygon &a = *this;

  c.Clear();
  if (b.GetSize() != 2 || !Closeness)
    return;

  for (i = 0; i < a.GetSize(); i++)
  {
    if (a(i).IsInLine(b[0]))
      cnt++;
    if (a(i).IsInLine(b[1]))
      cnt++;
  }
  if (cnt != 2)
    return; // strange

  for (i = 0; i < a.GetSize(); i++)
  {
    if (a(i).IsInLine(b[0]))
      index[0] = i;
    if (a(i).IsInLine(b[1]))
      index[1] = i;
  }

  meas[0] = Measure(b[0]);
  meas[1] = Measure(b[1]);

  if (index[0] < index[1])
  {
    c.Add(b[0]);
    for (i = index[0]; i < index[1]; i++)
      c.Add(FPoint[i + 1]);
    c.Add(b[1]);
  }
  else if (index[0] == index[1])
  {
    if (meas[0] < meas[1])
    {
      c.Add(b[0]);
      c.Add(b[1]);
      return;
    }
    else
    {
      c.Add(b[0]);
      for (i = index[0]; i < FSize; i++)
        c.Add(FPoint[(i + 1 == FSize) ? 0 : i + 1]);
      for (i = 0; i < index[0]; i++)
        c.Add(FPoint[i + 1]);
      c.Add(b[1]);
      return;
    }
  }
  else
  {
    c.Add(b[0]);
    for (i = index[0]; i < FSize; i++)
      c.Add(FPoint[(i + 1 == FSize) ? 0 : i + 1]);
    for (i = 0; i < index[1]; i++)
      c.Add(FPoint[i + 1]);
    c.Add(b[1]);
  }
}

MkLine &MkPolygon::operator()(int i)
{
  static MkLine ml;
  if (i >= 0 && i < GetSize() - 1)
  {
#ifdef __BCPLUSPLUS__
    ml.SetLine((*this)[i], (*this)[i + 1], Color);
#else
    ml.SetLine((*this)[i], (*this)[i + 1]);
#endif
    return ml;
  }
  else if (Closeness && i == GetSize() - 1)
  {
#ifdef __BCPLUSPLUS__
    ml.SetLine((*this)[i], (*this)[0], Color);
#else
    ml.SetLine((*this)[i], (*this)[0]);
#endif
    return ml;
  }
  else
    return NullLine;
}

MkPoint &MkPolygon::operator[](int i)
{
  isAreaChanged = isLengthChanged = true;
  if (i < FSizeOfArray && FSize < i + 1)
    FSize = i + 1;
  return MkPoints::operator[](i);
}

int MkPolygon::Next()
{
  FCurrent++;
  if (FCurrent >= GetSize())
    FCurrent = 0;
  return FCurrent;
}

int MkPolygon::Prev()
{
  FCurrent--;
  if (FCurrent < 0)
    FCurrent = GetSize() - 1;
  return FCurrent;
}

int MkPolygon::Next(int i)
{
  i++;
  if (i >= GetSize())
    i = 0;
  return i;
}

int MkPolygon::Prev(int i)
{
  i--;
  if (i < 0)
    i = GetSize() - 1;
  return i;
}

int MkPolygon::AlivedNext()
{
  int cnt = 0;
  while (!PointState(Next()) && cnt < GetSize())
    cnt++;
  if (cnt >= GetSize())
  {
    MkDebug("There is no vaild next!");
    return 0;
  }
  else
    return FCurrent;
}

int MkPolygon::AlivedPrev()
{
  int cnt = 0;
  while (!PointState(Prev()) && cnt < GetSize())
    cnt++;
  if (cnt >= GetSize())
  {
    MkDebug("There is no vaild prev!");
    return 0;
  }
  else
    return FCurrent;
}

int MkPolygon::AlivedNext(int i)
{
  int cnt = 0, k = i;
  k = Next(k);
  while (!PointState(k) && cnt < GetSize())
  {
    cnt++;
    k = Next(k);
  }
  if (cnt >= GetSize())
  {
    MkDebug("There is no vaild next of index i!");
    return 0;
  }
  else
    return k;
}

int MkPolygon::AlivedPrev(int i)
{
  int cnt = 0, k = i;
  while (!PointState(k = Prev(k)) && cnt < GetSize())
  {
    cnt++;
  }
  if (cnt >= GetSize())
  {
    MkDebug("There is no vaild prev of index i!");
    return 0;
  }
  else
    return k;
}

MkPoint &MkPolygon::NextPoint()
{
  return FPoint[Next()];
}

MkPoint &MkPolygon::PrevPoint()
{
  return FPoint[Prev()];
}

MkPoint &MkPolygon::NextPoint(int i)
{
  return FPoint[Next(i)];
}

MkPoint &MkPolygon::PrevPoint(int i)
{
  return FPoint[Prev(i)];
}

MkPoint &MkPolygon::NextPoint(MkPoint &pnt)
{
  int i;
  for (i = 0; i < FSize; i++)
    if (FPoint[i] == pnt)
      break;

  if (i == FSize)
    return NullPoint;
  else
    return FPoint[Next(i)];
}

MkPoint &MkPolygon::PrevPoint(MkPoint &pnt)
{
  int i;
  for (i = 0; i < FSize; i++)
    if (FPoint[i] == pnt)
      break;

  if (i == FSize)
    return NullPoint;
  else
    return FPoint[Prev(i)];
}

MkPoint &MkPolygon::AlivedNextPoint()
{
  int cnt = 0;
  while (!PointState(Next()) && cnt < GetSize())
    cnt++;
  if (cnt >= GetSize())
  {
    MkDebug("There is no vaild next!");
    return NullPoint;
  }
  else
    return (*this)[FCurrent];
}

MkPoint &MkPolygon::AlivedPrevPoint()
{
  int cnt = 0;
  while (!PointState(Prev()) && cnt < GetSize())
    cnt++;
  if (cnt >= GetSize())
  {
    MkDebug("There is no vaild prev!");
    return NullPoint;
  }
  else
    return (*this)[FCurrent];
}

MkPoint &MkPolygon::AlivedNextPoint(int i)
{
  int cnt = 0, k = i;
  while (!PointState(k = Next(k)) && cnt < GetSize())
  {
    cnt++;
  }
  if (cnt >= GetSize())
  {
    MkDebug("There is no vaild next of index i!");
    return NullPoint;
  }
  else
    return (*this)[k];
}

MkPoint &MkPolygon::AlivedPrevPoint(int i)
{
  int cnt = 0, k = i;
  while (!PointState(k = Prev(k)) && cnt < GetSize())
  {
    cnt++;
  }
  if (cnt >= GetSize())
  {
    MkDebug("There is no vaild prev of index i!");
    return NullPoint;
  }
  else
    return (*this)[k];
}

int MkPolygon::GetAlivePoint()
{
  int k = 0;
  for (int i = 0; i < GetSize(); i++)
    k += PointState(i);
  return k;
}

bool MkPolygon::Out(char *fname)
{
  FILE *fp;
  fp = fopen(fname, "a");
  if (!fp)
  {
    MkDebug(fname);
    MkDebug(" is not found, so fp is null and return false\n");
    return false;
  }

  fprintf(fp, "  .begin of polygon\n");
  for (int i = 0; i < FSize; i++)
    fprintf(fp, "  %15.8f, %15.8f\n", FPoint[i].X, FPoint[i].Y);
  fprintf(fp, "  .end of polygon\n");
  fclose(fp);
}

MkPolygon &MkPolygon::operator=(MkPolygon &polygon)
{
  MkPoints::operator=((MkPoints &)polygon);

  Closeness = polygon.Closeness;
  Convexity = polygon.Convexity;
  Crossing = polygon.Crossing;
  Fullness = polygon.Fullness;
  isLengthChanged = polygon.isLengthChanged;
  isAreaChanged = polygon.isAreaChanged;
  isCrossingChanged = polygon.isCrossingChanged;
  isClosenessChanged = polygon.isClosenessChanged;
  isConvexityChanged = polygon.isConvexityChanged;

  FLength = polygon.FLength;
  FArea = polygon.FArea;
  PointState = polygon.PointState;
  FCurrent = polygon.FCurrent;

  return *this;
}
bool MkPolygon::operator!=(MkPolygon &polygon)
{
  bool flag = true;

  flag = flag && MkPoints::operator==(polygon);

  flag = flag && Closeness == polygon.Closeness;
  flag = flag && Convexity == polygon.Convexity;
  flag = flag && Crossing == polygon.Crossing;
  flag = flag && Fullness == polygon.Fullness;
  flag = flag && isLengthChanged == polygon.isLengthChanged;
  flag = flag && isAreaChanged == polygon.isAreaChanged;
  flag = flag && isCrossingChanged == polygon.isCrossingChanged;
  flag = flag && isClosenessChanged == polygon.isClosenessChanged;
  flag = flag && isConvexityChanged == polygon.isConvexityChanged;

  flag = flag && FLength == polygon.FLength;
  flag = flag && FArea == polygon.FArea;
  flag = flag && PointState == polygon.PointState;
  flag = flag && FCurrent == polygon.FCurrent;

  return flag;
}

#ifdef __BCPLUSPLUS__
void MkPolygon::Draw(TObject *Sender)
{
  int c;
  TColor C;
  c = Closeness ? 0 : -1;
  if (String(Sender->ClassName()) == String("MkPaintBox"))
  {
    MkPoints::Draw(Sender);
    MkPaintBox *pb = (MkPaintBox *)Sender;
    C = pb->Canvas->Pen->Color;
    pb->Canvas->Pen->Color = Color;

    for (int i = 0; i < GetSize() + c; i++)
    {
      (*this)(i).Draw(Sender);
    }
    pb->Canvas->Pen->Color = C;
  }
}
#endif

#if defined(_MSC_VER) && defined(_WINDOWS_)
void MkPolygon::Draw(MkPaint *pb)
{
  int c = Closeness ? 0 : -1;

  for (int i = 0; i < GetSize() + c; i++)
    (*this)(i).Draw(pb);
}
#endif

void GetSubPolygon(double ymin, double ymax, MkPolygon &inpoly, MkPolygon &outpoly) // special procedure for application of load and subreact
{                                                                                   // extract subpolygon between the depth y1 and y2
  int l, k;
  char str[256];
  MkLine line, tl;
  MkPoint pnt[2], isp[2], esp[2]; //pnt is edge point, isp is intersection point, esp is whatever you want it to be...
  MkPoints pnts;
  pnts.Clear();
  outpoly.Clear();
  double pmax = -1.0e20, pmin = 1.0e20;

  for (int i = 0; i < inpoly.GetSize(); i++)
  {
    if (pmax < inpoly[i].X)
      pmax = inpoly[i].X;
    if (pmin > inpoly[i].X)
      pmin = inpoly[i].X;
  }

  pmax = pmax + fabs(pmax * 0.1);
  pmin = pmin - fabs(pmin * 0.1);

  if (ymin > ymax)
    swap(ymin, ymax);

  pnt[0].X = 0;
  pnt[0].Y = ymax;
  pnt[1].X = 0;
  pnt[1].Y = ymin;

  for (k = 0; k < inpoly.GetSize(); k++)
  {
    if (pnt[1].Y < inpoly[k].Y + EPS && inpoly[k].Y < pnt[0].Y + EPS)
      pnts.Add(inpoly[k]);
  }
  for (k = 0; k < pnts.GetSize() - 1; k++) // sort
    for (l = k + 1; l < pnts.GetSize(); l++)
      if (pnts[k].Y < pnts[l].Y)
        Swap(pnts[k], pnts[l]);

  line.SetLine(MkPoint(pmin, pnt[0].Y, 0), MkPoint(pmax, pnt[0].Y, 0));
  line.SetFiniteness(true);
  for (k = 0; k < inpoly.GetSize() - 1; k++)
  {
    tl = inpoly(k);
    tl.Extend(1.01);
    tl.SetFiniteness(true);
    if (line && tl)
      isp[0] = line.GetIntPoint(tl);
  }
  esp[0] = pnts[0];
  esp[0].X = 0;

  line.SetLine(MkPoint(pmin, pnt[1].Y, 0), MkPoint(pmax, pnt[1].Y, 0));
  line.SetFiniteness(true);
  for (k = 0; k < inpoly.GetSize() - 1; k++)
  {
    tl = inpoly(k);
    tl.Extend(1.01);
    tl.SetFiniteness(true);
    if (line && tl)
      isp[1] = line.GetIntPoint(tl);
  }
  esp[1] = pnts[pnts.GetSize() - 1];
  esp[1].X = 0;

  if (pnt[0] == pnts[0])
    ;
  else if (pnt[0] == isp[0])
    outpoly.Add(pnt[0]);
  else if (fabs(pnt[0].Y - isp[0].Y) < EPS)
  {
    outpoly.Add(pnt[0]);
    outpoly.Add(isp[0]);
  }
  else if (isp[0] == NullPoint && esp[0] != pnt[0])
  {
    outpoly.Add(pnt[0]);
    outpoly.Add(esp[0]);
  }
  else
    outpoly.Add(pnt[0]);
  for (k = 0; k < pnts.GetSize(); k++)
    outpoly.Add(pnts[k]);
  if (pnt[1] == pnts[pnts.GetSize() - 1])
    ;
  else if (pnt[1] == isp[1])
    outpoly.Add(pnt[1]);
  else if (fabs(pnt[1].Y - isp[1].Y) < EPS)
  {
    outpoly.Add(isp[1]);
    outpoly.Add(pnt[1]);
  }
  else if (isp[1] == NullPoint && esp[1] != pnt[1])
  {
    outpoly.Add(esp[1]);
    outpoly.Add(pnt[1]);
  }
  else
  {
    outpoly.Add(pnt[1]);
  }

  sprintf(str, "ymax = %10.5f, ymin = %10.5f \n", ymax, ymin);
  MkDebug(str);
  for (k = 0; k < outpoly.GetSize(); k++)
  {
    sprintf(str, "%d-th point of outpolyline is (%f, %f, %f)\n", k, outpoly[k].X, outpoly[k].Y, outpoly[k].Z);
    MkDebug(str);
  }
}

void GetSubParam(int i, MkPolygon &in, double &aj, double &bj, double &lj1, double &lj) // special procedure for application of load and subreact
{                                                                                       // calc pamameter of line segment. x and y axis are inverted.
  if (i > in.GetSize() - 1 || i < 0)
    return;    // if the polygon is not sorted in y direction, then there will
  MkLine line; // be some errors in values.
  MkPoint sp, ep;
  double a, b, c;

  aj = bj = lj1 = lj = 0;

  sp.X = -in[i].Y + in[0].Y;
  sp.Y = in[i].X;
  ep.X = -in[i + 1].Y + in[0].Y;
  ep.Y = in[i + 1].X;

  if (sp.X > ep.X)
    Swap(sp, ep);

  line.SetLine(sp, ep);

  lj1 = line[1].X;
  lj = line[0].X;

  a = line.GetA();
  b = line.GetB();
  c = line.GetC();

  if (fabs(b) < EPS)
    return;
  aj = -a / b;
  bj = c / b;
}
//---------------------------------------------------------------------------
MkPolygons::MkPolygons(int size, MkPolygon *polys)
{

  if (size < 0)
  {
    MkDebug("::MkPolygons - MkPolygons(int size)");
    ;
    return;
  }

  FSizeOfArray = FSize = size;
  if (FSize == 0)
  {
    FPolygon = NULL;
    return;
  }

  FPolygon = new MkPolygon[FSize];
  for (int i = 0; i < FSize; i++)
    (*this)[i] = polys[i];
}

MkPolygons::MkPolygons(int size)
{
  if (size < 0)
  {
    MkDebug("::MkPolygons - MkPolygons(int size)");
    ;
    return;
  }

  FSize = FSizeOfArray = size;

  if (FSizeOfArray == 0)
  {
    FPolygon = NULL;
    return;
  }

  FPolygon = new MkPolygon[FSizeOfArray];
}

MkPolygons::~MkPolygons()
{
  FSizeOfArray = FSize = 0;
  if (FPolygon)
  {
    delete[] FPolygon;
    FPolygon = NULL;
  }
}

void MkPolygons::Initialize(int size)
{
  if (size < 0)
  {
    MkDebug("::MkPolygons - Initialize(int size)");
    ;
    return;
  }
  if (FSizeOfArray == size)
    return;

  FSize = FSizeOfArray = size;

  if (FSizeOfArray == 0)
  {
    if (FPolygon != NULL)
      delete[](MkPolygon *) FPolygon;
    FPolygon = NULL;
    return;
  }

  if (FPolygon != NULL)
    delete[](MkPolygon *) FPolygon;
  FPolygon = new MkPolygon[FSizeOfArray];
}

void MkPolygons::Initialize(int size, MkPolygon *polys)
{

  if (size < 0 || polys == NULL)
  {
    MkDebug("::MkPolygons - Initialize(int size)");
    ;
    return;
  }
  if (FSizeOfArray == size)
    return;
  FSize = FSizeOfArray = size;
  if (FSizeOfArray == 0)
  {
    if (FPolygon != NULL)
      delete[](MkPolygon *) FPolygon;
    FPolygon = NULL;
    return;
  }

  if (FPolygon != NULL)
    delete[](MkPolygon *) FPolygon;
  FPolygon = new MkPolygon[FSizeOfArray];
  for (int i = 0; i < FSizeOfArray; i++)
    FPolygon[i] = polys[i];
}

int MkPolygons::Grow(int delta)
{
  int i;
  MkPolygon *poly = NULL;

  if (!(poly = new MkPolygon[FSizeOfArray + delta]))
    return FSizeOfArray;

  for (i = 0; i < FSize; i++)
    poly[i] = FPolygon[i];
  for (i = FSize; i < FSizeOfArray + delta; i++)
    poly[i] = NullPolygon;
  if (FPolygon)
  {
    delete[](MkPolygon *) FPolygon;
    FPolygon = NULL;
  }
  FPolygon = poly;
  FSizeOfArray = FSizeOfArray + delta;

  return FSizeOfArray;
}

int MkPolygons::Shrink(int delta)
{
  int i;
  MkPolygon *poly = NULL;

  if (!(poly = new MkPolygon[FSizeOfArray - delta]))
    return FSizeOfArray;

  for (i = 0; i < FSize; i++)
    poly[i] = FPolygon[i];
  for (i = FSize; i < FSizeOfArray - delta; i++)
    poly[i] = NullPolygon;
  if (FPolygon)
  {
    delete[](MkPolygon *) FPolygon;
    FPolygon = NULL;
  }
  FPolygon = poly;
  FSizeOfArray = FSizeOfArray - delta;

  return FSizeOfArray;
}

bool MkPolygons::Add(MkPolygon &poly)
{
  int tmp = FSizeOfArray;
  bool flag = false;
  for (int i = 0; i < FSize; i++)
    if (FPolygon[i] == poly)
      flag = true;

  if (flag)
    return false;
  if (FSize >= FSizeOfArray)
  {
    Grow(FSize - FSizeOfArray + 10);
    if (tmp == FSizeOfArray)
      return false;
  }
  FSize++;
  FPolygon[FSize - 1] = poly;

  return true;
}

bool MkPolygons::Add(int index, MkPolygon &poly)
{
  int tmp = FSizeOfArray;

  if (FSize >= FSizeOfArray)
    Grow(FSize - FSizeOfArray + 1);
  if (tmp == FSizeOfArray)
    return false;

  for (int i = FSize - 1; i >= index; i--)
    FPolygon[i + 1] = FPolygon[i];
  FSize++;
  FPolygon[index] = poly;
  return true;
}

bool MkPolygons::Delete(MkPolygon &poly)
{
  int i;
  for (i = 0; i < FSize; i++)
  {
    if (FPolygon[i] == poly)
      break;
  }
  if (i == FSize)
    return false;
  if (FPolygon[i] == poly)
  {
    for (int j = i; j < FSize - 1; j++)
      FPolygon[j] = FPolygon[j + 1];
  }
  FSize--;
  FPolygon[FSize] = NullPolygon;
  return true;
}

bool MkPolygons::Delete(int index)
{
  for (int j = index; j < FSize - 1; j++)
    FPolygon[j] = FPolygon[j + 1];

  FSize--;
  FPolygon[FSize] = NullPolygon;
  return true;
}

bool MkPolygons::Clear()
{
  FSizeOfArray = FSize = 0;
  if (FPolygon)
  {
    delete[] FPolygon;
    FPolygon = NULL;
  }
  return true;
}

MkPolygon &MkPolygons::operator[](int i)
{
  if (FSizeOfArray == 0)
    return NullPolygon;
  if (i >= FSize && i < FSizeOfArray)
    FSize = i + 1;

  if (i >= 0 && i < FSize)
    return FPolygon[i];
  else
    return NullPolygon;
}

MkPolygons &MkPolygons::operator=(MkPolygons &polys)
{
  int i;

  Clear();
  FSize = polys.FSize;
  FSizeOfArray = polys.FSizeOfArray;
  if (FSize == 0)
  {
    FPolygon = NULL;
    return *this;
  }
  this->FPolygon = new MkPolygon[FSizeOfArray];

  for (i = 0; i < FSize; i++)
    FPolygon[i] = polys.FPolygon[i];
  for (i = FSize; i < FSizeOfArray; i++)
    FPolygon[i] = NullPolygon;

  return *this;
}

bool MkPolygons::operator==(MkPolygons &polys)
{
  int i;

  if (FSize != polys.FSize)
    return false;
  for (i = 0; i < FSize; i++)
    if (this->FPolygon[i] != polys.FPolygon[i])
      return false;

  return true;
}

#ifdef __BCPLUSPLUS__
void MkPolygons::Draw(TObject *Sender)
{
  for (int i = 0; i < FSize; i++)
    FPolygon[i].Draw(Sender);
}
#endif

#if defined(_MSC_VER) && defined(_WINDOWS_)
void MkPolygons::Draw(MkPaint *pb)
{
  for (int i = 0; i < FSize; i++)
    FPolygon[i].Draw(pb);
}
#endif
