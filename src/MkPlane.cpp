//---------------------------------------------------------------------------
#include "MkPlane.hpp"

MkPointsPlane NullPointPlane(0);
MkPlane NullPlane(0);
MkJointPlane NullJoint(0, 0);
MkJointPlanes NullJoints(0);
MkPennyJoint NullPenny(0, 0);
MkPennyJoints NullPennys(0);

MkPlane::MkPlane()
{
  ScaleX = ScaleY = 1;
  Alpha = Beta = Gamma = 0;
  LD.X = -1;
  LD.Y = -1;
  LU.X = -1;
  LU.Y = 1;
  RD.X = 1;
  RD.Y = -1;
  RU.X = 1;
  RU.Y = 1;
  FCenter.Set(0, 0, 0);
  CalcABCD();
  className = "MkPlane";
}

MkPlane::MkPlane(int i)
{
  ScaleX = ScaleY = 1;
  Alpha = Beta = Gamma = 0;
  LD.X = -1;
  LD.Y = -1;
  LU.X = -1;
  LU.Y = 1;
  RD.X = 1;
  RD.Y = -1;
  RU.X = 1;
  RU.Y = 1;
  FCenter.Set(0, 0, 0);
  CalcABCD();
  className = "MkPlane";
}

MkPlane::MkPlane(MkPoint rp1, MkPoint rp2, MkPoint rp3)
{
  SetPoints(rp1, rp2, rp3);
  className = "MkPlane";
}

MkPlane::MkPlane(MkPoint *rps)
{
  MkPoint rp1, rp2, rp3;
  rp1 = rps[0];
  rp2 = rps[1];
  rp3 = rps[2];
  SetPoints(rp1, rp2, rp3);
  className = "MkPlane";
}

void MkPlane::ResetAll()
{
  ScaleX = ScaleY = 1;
  Alpha = Beta = Gamma = 0;
  LD.X = -1;
  LD.Y = -1;
  LU.X = -1;
  LU.Y = 1;
  RD.X = 1;
  RD.Y = -1;
  RU.X = 1;
  RU.Y = 1;
  FCenter.Set(0, 0, 0);
  CalcABCD();
}

void MkPlane::ResetScale()
{
  ScaleX = ScaleY = 1;
  ApplyMatrix();
}

void MkPlane::ResetRotate()
{
  Alpha = Beta = Gamma = 0;
  ApplyMatrix();
  CalcABCD();
}

void MkPlane::ResetTranslate()
{
  FCenter.Set(0, 0, 0);
  ApplyMatrix();
}

//    MkPoint LD,RD,LU,RU;
//    MkPoint FCenter;
//    double ScaleX;
//    double ScaleY;
//    double Alpha,Beta,Gamma; // �� x,y,z �������� ȸ��, ������ x,y,z ��
//    TColor Color;
//    double A,B,C,D; // Ax + By + Cz = D

void MkPlane::SetPoints(MkPoint rp1, MkPoint rp2, MkPoint rp3)
{
  double x, y, z;

  ResetAll();

  x = rp2.X - rp1.X;
  y = rp2.Y - rp1.Y;
  z = rp2.Z - rp1.Z;
  MkVector<double> vec1(x, y, z);

  x = rp3.X - rp1.X;
  y = rp3.Y - rp1.Y;
  z = rp3.Z - rp1.Z;
  MkVector<double> vec2(x, y, z);

  vec1.Unify();
  vec2.Unify();
  MkVector<double> vec3;
  /*    vec1.Cross(vec2,vec3);  // cross product need some extra work
    A = vec3(0);
    B = vec3(1);
    C = vec3(2);
    */

  A = vec1(1) * vec2(2) - vec1(2) * vec2(1);
  B = vec1(2) * vec2(0) - vec1(0) * vec2(2);
  C = vec1(0) * vec2(1) - vec1(1) * vec2(0);

  x = (rp1.X + rp2.X + rp3.X) / 3.0;
  y = (rp1.Y + rp2.Y + rp3.Y) / 3.0;
  z = (rp1.Z + rp2.Z + rp3.Z) / 3.0;
  FCenter.Set(x, y, z);
  D = A * x + B * y + C * z;

  //MkDebug("A %f B %f C %f D %f\n",A,B,C,D);
  //MkDebug("v1x %f v1y  %f v1z %f \n",vec1(0),vec1(1),vec1(2));
  //MkDebug("v2x %f v2y  %f v2z %f \n",vec2(0),vec2(1),vec2(2));
}

void MkPlane::SetPoints(MkPoint *rps)
{
  MkPoint rp1, rp2, rp3;
  rp1 = rps[0];
  rp2 = rps[1];
  rp3 = rps[2];
  SetPoints(rp1, rp2, rp3);
}

void MkPlane::SetScale(double scaleX, double scaleY)
{
  ScaleX = scaleX;
  ScaleY = scaleY;
  ApplyMatrix();
}

void MkPlane::SetRotate(double alpha, double beta, double gamma)
{
  Alpha = alpha;
  Beta = beta;
  Gamma = gamma;
  ApplyMatrix();
  CalcABCD();
}

void MkPlane::SetTranslate(double x, double y, double z)
{
  FCenter.Set(x, y, z);
  ApplyMatrix();
}

void MkPlane::SetTranslate(MkPoint trans)
{
  FCenter.Set(trans);
  ApplyMatrix();
}

//������ ���� ȸ��
void MkPlane::Rotate(double alpha, double beta, double gamma)
{
  MkMatrix4<double> rm;
  rm.LoadIdentity();
  rm.Rotate(alpha, beta, gamma);

  LD = LD * rm;
  RD = RD * rm;
  LU = LU * rm;
  RU = RU * rm;
  FCenter = FCenter * rm;

  MkPoint rp;
  double x, y, z, len;
  rp = RU - LU;

  x = rp.X;
  y = rp.Y;
  len = sqrt(x * x + y * y);

  x /= len;
  y /= len;

  Gamma = acos(x);
  if (fabs(sin(Gamma) - y) > 0.001)
    Gamma = -Gamma;
  Gamma = Gamma * 180 / M_PI;

  rp.RotateInZ(-Gamma);

  x = rp.X;
  z = rp.Z;
  len = sqrt(x * x + z * z);

  x /= len;
  z /= len;
  Beta = acos(x);
  if (fabs(sin(Beta) - z) > 0.001)
    Beta = -Beta;
  Beta = Beta * 180 / M_PI;

  rp = RU - RD;
  rp.RotateInZ(-Gamma);
  rp.RotateInY(-Beta);

  y = rp.Y;
  z = rp.Z;
  len = sqrt(y * y + z * z);

  y /= len;
  z /= len;
  Alpha = acos(y);
  if (fabs(sin(Alpha) - z) > 0.001)
    Alpha = -Alpha;
  Alpha = Alpha * 180 / M_PI;
}

//������ ���� ȸ��
void MkPlane::RotateInX(double ang)
{
  MkMatrix4<double> rm;
  rm.LoadIdentity();
  rm.RotateInX(ang);

  LD = LD * rm;
  RD = RD * rm;
  LU = LU * rm;
  RU = RU * rm;
  FCenter = FCenter * rm;

  MkPoint rp;
  double x, y, z, len;
  rp = RU - LU;

  x = rp.X;
  y = rp.Y;
  len = sqrt(x * x + y * y);

  x /= len;
  y /= len;

  Gamma = acos(x);
  if (fabs(sin(Gamma) - y) > 0.001)
    Gamma = -Gamma;
  Gamma = Gamma * 180 / M_PI;

  rp.RotateInZ(-Gamma);

  x = rp.X;
  z = rp.Z;
  len = sqrt(x * x + z * z);

  x /= len;
  z /= len;
  Beta = acos(x);
  if (fabs(sin(Beta) - z) > 0.001)
    Beta = -Beta;
  Beta = Beta * 180 / M_PI;

  rp = RU - RD;
  rp.RotateInZ(-Gamma);
  rp.RotateInY(-Beta);

  y = rp.Y;
  z = rp.Z;
  len = sqrt(y * y + z * z);

  y /= len;
  z /= len;
  Alpha = acos(y);
  if (fabs(sin(Alpha) - z) > 0.001)
    Alpha = -Alpha;
  Alpha = Alpha * 180 / M_PI;
}

//������ ���� ȸ��
void MkPlane::RotateInY(double ang)
{
  MkMatrix4<double> rm;
  rm.LoadIdentity();
  rm.RotateInY(ang);

  LD = LD * rm;
  RD = RD * rm;
  LU = LU * rm;
  RU = RU * rm;
  FCenter = FCenter * rm;

  MkPoint rp;
  double x, y, z, len;
  rp = RU - LU;

  x = rp.X;
  y = rp.Y;
  len = sqrt(x * x + y * y);

  x /= len;
  y /= len;

  Gamma = acos(x);
  if (fabs(sin(Gamma) - y) > 0.001)
    Gamma = -Gamma;
  Gamma = Gamma * 180 / M_PI;

  rp.RotateInZ(-Gamma);

  x = rp.X;
  z = rp.Z;
  len = sqrt(x * x + z * z);

  x /= len;
  z /= len;
  Beta = acos(x);
  if (fabs(sin(Beta) - z) > 0.001)
    Beta = -Beta;
  Beta = Beta * 180 / M_PI;

  rp = RU - RD;
  rp.RotateInZ(-Gamma);
  rp.RotateInY(-Beta);

  y = rp.Y;
  z = rp.Z;
  len = sqrt(y * y + z * z);

  y /= len;
  z /= len;
  Alpha = acos(y);
  if (fabs(sin(Alpha) - z) > 0.001)
    Alpha = -Alpha;
  Alpha = Alpha * 180 / M_PI;
}

//������ ���� ȸ��
void MkPlane::RotateInZ(double ang)
{
  MkMatrix4<double> rm;
  rm.LoadIdentity();
  rm.RotateInZ(ang);

  LD = LD * rm;
  RD = RD * rm;
  LU = LU * rm;
  RU = RU * rm;
  FCenter = FCenter * rm;

  MkPoint rp;
  double x, y, z, len;
  rp = RU - LU;

  x = rp.X;
  y = rp.Y;
  len = sqrt(x * x + y * y);

  x /= len;
  y /= len;

  Gamma = acos(x);
  if (fabs(sin(Gamma) - y) > 0.001)
    Gamma = -Gamma;
  Gamma = Gamma * 180 / M_PI;

  rp.RotateInZ(-Gamma);

  x = rp.X;
  z = rp.Z;
  len = sqrt(x * x + z * z);

  x /= len;
  z /= len;
  Beta = acos(x);
  if (fabs(sin(Beta) - z) > 0.001)
    Beta = -Beta;
  Beta = Beta * 180 / M_PI;

  rp = RU - RD;
  rp.RotateInZ(-Gamma);
  rp.RotateInY(-Beta);

  y = rp.Y;
  z = rp.Z;
  len = sqrt(y * y + z * z);

  y /= len;
  z /= len;
  Alpha = acos(y);
  if (fabs(sin(Alpha) - z) > 0.001)
    Alpha = -Alpha;
  Alpha = Alpha * 180 / M_PI;
}

//������ ���� ȸ��
void MkPlane::RotateInA(double ang, double l, double m, double n)
{
  MkMatrix4<double> rm;
  rm.LoadIdentity();
  rm.RotateInA(ang, l, m, n);

  LD = LD * rm;
  RD = RD * rm;
  LU = LU * rm;
  RU = RU * rm;
  FCenter = FCenter * rm;

  MkPoint rp;
  double x, y, z, len;
  rp = RU - LU;

  x = rp.X;
  y = rp.Y;
  len = sqrt(x * x + y * y);

  x /= len;
  y /= len;

  Gamma = acos(x);
  if (fabs(sin(Gamma) - y) > 0.001)
    Gamma = -Gamma;
  Gamma = Gamma * 180 / M_PI;

  rp.RotateInZ(-Gamma);

  x = rp.X;
  z = rp.Z;
  len = sqrt(x * x + z * z);

  x /= len;
  z /= len;
  Beta = acos(x);
  if (fabs(sin(Beta) - z) > 0.001)
    Beta = -Beta;
  Beta = Beta * 180 / M_PI;

  rp = RU - RD;
  rp.RotateInZ(-Gamma);
  rp.RotateInY(-Beta);

  y = rp.Y;
  z = rp.Z;
  len = sqrt(y * y + z * z);

  y /= len;
  z /= len;
  Alpha = acos(y);
  if (fabs(sin(Alpha) - z) > 0.001)
    Alpha = -Alpha;
  Alpha = Alpha * 180 / M_PI;
}

//Trans�� �߽����� �� ȸ�� Trans�� �̵����� ����.
void MkPlane::RotateCen(double alpha, double beta, double gamma)
{
  MkPoint rp(A, B, C);
  MkMatrix4<double> rm;
  rm.Identity();
  rm.Rotate(alpha, beta, gamma);

  LD = LD * rm;
  RD = RD * rm;
  LU = LU * rm;
  RU = RU * rm;
  FCenter = FCenter * rm;

  rp = rp * rm;
  rp.Unify();
  A = rp.X;
  B = rp.Y;
  C = rp.Z;
  D = A * FCenter.X + B * FCenter.Y + C * FCenter.Z;

  rp.GetAng(Alpha, Beta, Gamma);
}

//Trans�� �߽����� �� ȸ�� Trans�� �̵����� ����.
void MkPlane::RotateInXCen(double ang)
{
  MkPoint rp(A, B, C);
  MkMatrix4<double> rm;
  rm.Identity();
  rm.RotateInX(ang);

  LD = LD * rm;
  RD = RD * rm;
  LU = LU * rm;
  RU = RU * rm;
  FCenter = FCenter * rm;

  rp = rp * rm;
  rp.Unify();
  A = rp.X;
  B = rp.Y;
  C = rp.Z;
  D = A * FCenter.X + B * FCenter.Y + C * FCenter.Z;

  rp.GetAng(Alpha, Beta, Gamma);
}

//Trans�� �߽����� �� ȸ�� Trans�� �̵����� ����.
void MkPlane::RotateInYCen(double ang)
{
  MkPoint rp(A, B, C);
  MkMatrix4<double> rm;
  rm.Identity();
  rm.RotateInY(ang);

  LD = LD * rm;
  RD = RD * rm;
  LU = LU * rm;
  RU = RU * rm;
  FCenter = FCenter * rm;

  rp = rp * rm;
  rp.Unify();
  A = rp.X;
  B = rp.Y;
  C = rp.Z;
  D = A * FCenter.X + B * FCenter.Y + C * FCenter.Z;

  rp.GetAng(Alpha, Beta, Gamma);
}

//Trans�� �߽����� �� ȸ�� Trans�� �̵����� ����.
void MkPlane::RotateInZCen(double ang)
{
  MkPoint rp(A, B, C);
  MkMatrix4<double> rm;
  rm.Identity();
  rm.RotateInZ(ang);

  LD = LD * rm;
  RD = RD * rm;
  LU = LU * rm;
  RU = RU * rm;
  FCenter = FCenter * rm;

  rp = rp * rm;
  rp.Unify();
  A = rp.X;
  B = rp.Y;
  C = rp.Z;
  D = A * FCenter.X + B * FCenter.Y + C * FCenter.Z;

  rp.GetAng(Alpha, Beta, Gamma);
}

//Trans�� �߽����� �� ȸ�� Trans�� �̵����� ����.
void MkPlane::RotateInACen(double ang, double l, double m, double n)
{
  MkPoint rp(A, B, C);
  MkMatrix4<double> rm;
  rm.Identity();
  rm.RotateInA(ang, l, m, n);

  LD = LD * rm;
  RD = RD * rm;
  LU = LU * rm;
  RU = RU * rm;
  FCenter = FCenter * rm;

  rp = rp * rm;
  rp.Unify();
  A = rp.X;
  B = rp.Y;
  C = rp.Z;
  D = A * FCenter.X + B * FCenter.Y + C * FCenter.Z;

  rp.GetAng(Alpha, Beta, Gamma);
}

bool MkPlane::GetIntParam(MkLine &rl, double &t1, double &t2)
{
  double l, m, n;
  double len;
  double t[2] = {0, 0};
  int cnt = 0;
  MkPoint sp, ep, rp;
  sp = rl[0];
  ep = rl[1];
  bool isPerpend;
  bool flag;

  l = ep.X - sp.X;
  m = ep.Y - sp.Y;
  n = ep.Z - sp.Z;

  len = sqrt(l * l + m * m + n * n);

  l = l / len;
  m = m / len;
  n = n / len;

  CalcABCD();

  isPerpend = false;
  if (fabs(l * A + m * B + n * C) < 0.001)
    isPerpend = true;

  if (isPerpend)
  {
    for (int i = 0; i < 4; i++)
      if (rl && (*this)(i))
      {
        if (rl.GetIntParam((*this)(i), t[cnt]))
          cnt++;
      }

    if (cnt == 0)
      return false;

    t1 = t[0];
    t2 = t[1];
    flag = true;
  }
  else
  {
    MkPoint rp;
    cnt = 1;
    t1 = (D - A * sp.X - B * sp.Y - C * sp.Z) / (A * rl.DeltaX() + B * rl.DeltaY() + C * rl.DeltaZ());
    t2 = 0;
    flag = true;
  }

  if (flag && rl.GetFiniteness())
  {
    flag = t1 > 0 && t1 < 1 && t2 > 0 && t2 < 1;
  }

  if (flag && isFinite && !isPerpend)
  {
    rp.Set(sp.X + A * t1, sp.Y + B * t1, sp.Z + C * t1);
    flag = IsIn(rp);
  }
  return flag;
}

void MkPlane::CalcABCD() //Ax+By+Cz=D
{
  MkPoint rp(0, 0, 1);
  MkMatrix4<double> rm;

  rm.Identity();
  rm.Rotate(Alpha, Beta, Gamma);
  rp = rp * rm;
  rp.Unify();

  A = rp.X;
  B = rp.Y;
  C = rp.Z;

  D = A * FCenter.X + B * FCenter.Y + C * FCenter.Z;
}

bool MkPlane::IsIntersect(MkLine rl)
{
  double t, dot;
  CalcABCD();
  if (rl.GetLength() < 0.0001)
    return false;

  dot = (A * rl.DeltaX() + B * rl.DeltaY() + C * rl.DeltaZ());
  if (fabs(dot) < 0.001)
    return false;

  t = (D - A * rl[0].X - B * rl[0].Y - C * rl[0].Z) / dot;
  if (t < 1 && t > 0)
    return true;
  else
    return false;
}

bool MkPlane::IsIntersect(MkTriangle &rt)
{
  int cnt = 0;
  for (int i = 0; i < 3; i++)
    if (IsIntersect(rt(i)))
      cnt++;
  if (cnt == 1)
    return false;
  return (cnt == 2);
}

bool MkPlane::IsInPlane(MkLine rl)
{
  for (int i = 0; i < 4; i++)
    if (!(*this)(i).IsInSamePlane(rl))
      return false;
  return true;
}

bool MkPlane::IsIn(MkPoint rp) // test is qued...^^
{
  CalcABCD();
  if (fabs(A * rp.X + B * rp.Y + C * rp.Z - D) > 0.001)
    return false;
  MkVector<double> v1, v2, c1, c2, c3, c4;
  MkLine rl1, rl2;

  rl1.SetLine(rp, (*this)[0]);
  rl2.SetLine(rp, (*this)[1]);
  v1.SetVector(rl1.GetL(), rl1.GetM(), rl1.GetN());
  v2.SetVector(rl2.GetL(), rl2.GetM(), rl2.GetN());
  v1.Cross(v2, c1);

  rl1.SetLine(rp, (*this)[1]);
  rl2.SetLine(rp, (*this)[2]);
  v1.SetVector(rl1.GetL(), rl1.GetM(), rl1.GetN());
  v2.SetVector(rl2.GetL(), rl2.GetM(), rl2.GetN());
  v1.Cross(v2, c2);

  rl1.SetLine(rp, (*this)[2]);
  rl2.SetLine(rp, (*this)[3]);
  v1.SetVector(rl1.GetL(), rl1.GetM(), rl1.GetN());
  v2.SetVector(rl2.GetL(), rl2.GetM(), rl2.GetN());
  v1.Cross(v2, c3);

  rl1.SetLine(rp, (*this)[3]);
  rl2.SetLine(rp, (*this)[0]);
  v1.SetVector(rl1.GetL(), rl1.GetM(), rl1.GetN());
  v2.SetVector(rl2.GetL(), rl2.GetM(), rl2.GetN());
  v1.Cross(v2, c4);

  if ((c1 * c2) * (c2 * c3) * (c3 * c4) * (c4 * c1) < 0)
    return false;
  else
    return true;
}

bool MkPlane::IsIn(MkLine rl)
{
  return IsIn(rl[0]) && IsIn(rl[1]);
}

bool MkPlane::IsInSurface(MkPoint &pnt, double thick) // infinitive plane
{
  return fabs(GetDistance(pnt)) < thick;
}

bool MkPlane::IsInnerSpace(MkPoint &pnt) // infinitive plane
{
  return GetDistance(pnt) > 0;
}

MkPoint &MkPlane::CalcIntPnt(MkLine rl)
{
  double t;
  static MkPoint rp;
  CalcABCD();
  if (rl.GetLength() < 0.0001)
    return NullPoint;
  t = (D - A * rl[0].X - B * rl[0].Y - C * rl[0].Z) / (A * rl.DeltaX() + B * rl.DeltaY() + C * rl.DeltaZ());
  rp.X = t * rl.DeltaX() + rl[0].X;
  rp.Y = t * rl.DeltaY() + rl[0].Y;
  rp.Z = t * rl.DeltaZ() + rl[0].Z;
  return rp;
}

MkLine &MkPlane::CalcIntLine(MkTriangle &rt)
{
  MkPoint rp[2];
  static MkLine rl;
  int cnt = 0;
  for (int i = 0; i < 3; i++)
    if (IsIntersect(rt(i)))
    {
      rp[cnt] = CalcIntPnt(rt(i));
      cnt++;
    }
  if (cnt == 2)
  {
    rl.SetLine(rp[0], rp[1]);
    return rl;
  }
  else
    return NullLine;
}

MkPoint &MkPlane::CalcIntPnt(MkTriangle &rt) // Shortest point
{
  static MkPoint rp;
  MkLine rl;
  rl = CalcIntLine(rt);
  rl.SetFiniteness(true);
  if (rl != NullLine)
  {
    rp = rl.GetNearestPnt(FCenter);
    return rp;
  }
  return NullPoint;
}

MkPoint &MkPlane::CalcIntPnt(MkTriangles &rts) // Shortest point
{
  static MkPoint sp;
  MkPoint rp;
  double sdis, rdis;
  bool is_first;
  is_first = true;

  for (int i = 0; i < rts.GetSize(); i++)
  {
    rp = CalcIntPnt(rts[i]);
    if (rp != NullPoint)
    {
      if (is_first)
      {
        sp = rp;
        sdis = CalDist(FCenter, sp);
        is_first = false;
      }
      else
      {
        rdis = CalDist(FCenter, rp);
        if (rdis < sdis)
        {
          sp = rp;
          sdis = rdis;
        }
      }
    }
  }
  return sp;
}

void MkPlane::ApplyMatrix()
{
  MkPoint LD(-1, -1, 0);
  MkPoint LU(-1, 1, 0);
  MkPoint RD(1, -1, 0);
  MkPoint RU(1, 1, 0);

  MkMatrix4<double> rm;
  rm.Identity();
  rm.Translate(FCenter.X, FCenter.Y, FCenter.Z);
  rm.Rotate(Alpha, Beta, Gamma);
  rm.Scale(ScaleX, ScaleY, 1);

  LD = LD * rm;
  LU = LU * rm;
  RD = RD * rm;
  RU = RU * rm;
}

#ifdef __BCPLUSPLUS__
void MkPlane::Draw(TObject *Sender)
{
  for (int i = 0; i < 4; i++)
    (*this)(i).Draw(Sender);
}
#endif

#if defined(_MSC_VER) && defined(_WINDOWS_)
void MkPlane::Draw(MkPaint *pb)
{
  for (int i = 0; i < 4; i++)
    (*this)(i).Draw(pb);
}
#endif

void MkPlane::operator*(MkMatrix4<double> &rm)
{
  // not yet deceided
}

MkPoint &MkPlane::operator[](int i)
{
  if (i == 0)
    return LD;
  else if (i == 1)
    return RD;
  else if (i == 2)
    return RU;
  else if (i == 3)
    return LU;
  else
    return NullPoint;
}

MkLine &MkPlane::operator()(int i)
{
#ifdef __BCPLUSPLUS__
  Edge[i].SetLine((*this)[i % 4], (*this)[(i + 1) % 4], clWhite);
#else
  Edge[i].SetLine((*this)[i % 4], (*this)[(i + 1) % 4]);
#endif
  return Edge[i];
}

bool MkPlane::operator==(MkPlane &rp)
{
  if (LD != rp.LD)
    return false;
  if (RD != rp.RD)
    return false;
  if (LU != rp.LU)
    return false;
  if (RU != rp.RU)
    return false;

  if (FCenter != rp.FCenter)
    return false;
  if (ScaleX != rp.ScaleX)
    return false;
  if (ScaleY != rp.ScaleY)
    return false;

  if (Alpha != rp.Alpha)
    return false;
  if (Beta != rp.Beta)
    return false;
  if (Gamma != rp.Gamma)
    return false;
  return true;
}

bool MkPlane::operator!=(MkPlane &rp)
{
  if (LD == rp.LD)
    return false;
  if (RD == rp.RD)
    return false;
  if (LU == rp.LU)
    return false;
  if (RU == rp.RU)
    return false;

  if (FCenter == rp.FCenter)
    return false;
  if (ScaleX == rp.ScaleX)
    return false;
  if (ScaleY == rp.ScaleY)
    return false;

  if (Alpha == rp.Alpha)
    return false;
  if (Beta == rp.Beta)
    return false;
  if (Gamma == rp.Gamma)
    return false;
  return true;
}

MkPlane &MkPlane::operator=(MkPlane &rp)
{
  MkShape::operator=((MkShape &)rp);
  LD = rp.LD;
  RD = rp.RD;
  LU = rp.LU;
  RU = rp.RU;
  FCenter = rp.FCenter;
  ScaleX = rp.ScaleX;
  ScaleY = rp.ScaleY;
  Alpha = rp.Alpha;
  Beta = rp.Beta;
  Gamma = rp.Gamma;
  return *this;
}
//---------------------------------------------------------------------------

MkPointsPlane::MkPointsPlane() : MkPoints() {}

MkPointsPlane::MkPointsPlane(int size) : MkPoints(size) {}

MkPointsPlane::MkPointsPlane(int size, MkPoint *rps) : MkPoints(size, rps) {}

MkPointsPlane::MkPointsPlane(int size, MkPoint *rps, MkOrient fo) : MkPoints(size, rps)
{
  MkPoint rp(0, 0, 1);
  MkMatrix4<double> rm;
  FOrient = fo;
  Alpha = 0;
  Beta = -fo.Dip;
  Gamma = -fo.DipDir;

  rm.Identity();
  rm.Rotate(Alpha, Beta, Gamma);
  rp = rp * rm;
  rp.Unify();

  A = rp.X;
  B = rp.Y;
  C = rp.Z;
  D = A * FPoint[0].X + B * FPoint[0].Y + C * FPoint[0].Z;
}

void MkPointsPlane::SetOrient(MkOrient fo)
{
  MkPoint rp(0, 0, 1);
  MkMatrix4<double> rm;
  FOrient = fo;

  Alpha = 0;
  Beta = -fo.Dip;
  Gamma = -fo.DipDir;

  rm.Identity();
  rm.Rotate(Alpha, Beta, Gamma);
  rp = rp * rm;
  rp.Unify();

  A = rp.X;
  B = rp.Y;
  C = rp.Z;

  D = A * FPoint[0].X + B * FPoint[0].Y + C * FPoint[0].Z;
}

bool MkPointsPlane::operator==(MkPointsPlane &pp)
{
  bool flag;
  flag = true;
  flag = flag && FOrient == pp.FOrient;
  flag = flag && FHeightMode == pp.FHeightMode;
  flag = flag && Alpha == pp.Alpha;
  flag = flag && Beta == pp.Beta;
  flag = flag && Gamma == pp.Gamma;
  flag = flag && A == pp.A;
  flag = flag && B == pp.B;
  flag = flag && C == pp.C;
  flag = flag && D == pp.D;
  return flag;
}

bool MkPointsPlane::operator!=(MkPointsPlane &pp)
{
  bool flag;
  flag = true;
  flag = flag && FOrient == pp.FOrient;
  flag = flag && FHeightMode == pp.FHeightMode;
  flag = flag && Alpha == pp.Alpha;
  flag = flag && Beta == pp.Beta;
  flag = flag && Gamma == pp.Gamma;
  flag = flag && A == pp.A;
  flag = flag && B == pp.B;
  flag = flag && C == pp.C;
  flag = flag && D == pp.D;
  return !flag;
}

bool MkPointsPlane::IsIn(MkPoint rp)
{
  MkPoint rps[2];
  MkLine rls[2];
  bool flag;
  flag = true;
  for (int i = 0; i < FSize; i++)
  {
    rps[0] = FPoint[i >= FSize ? i - FSize : i];
    rps[1] = FPoint[i + 1 >= FSize ? i + 1 - FSize : i + 1];
    rps[0].Z = rps[1].Z = rp.Z = 0;
    rls[0].SetLine(rps[0], rps[1]);
    rls[1].SetLine(rps[0], rp);
    flag = flag && (rls[0] * rls[1] > -0.01);
  }
  return flag;
}

double MkPointsPlane::CalcHeight(MkPoint rp)
{
  return (fabs(C) > 0.0001) ? (D - A * rp.X - B * rp.Y) / C : MAX_VAL;
}

#if defined(__GL_H__)
void MkPointsPlane::Draw()
{
}
#endif

#ifdef __BCPLUSPLUS__
void MkPointsPlane::Draw(TObject *)
{
}
#endif

#if defined(_MSC_VER) && defined(_WINDOWS_)
void MkPointsPlane::Draw(MkPaint *)
{
}
#endif

//---------------------------------------------------------------------------

MkJointPlane::MkJointPlane() : MkPlane()
{
}

MkJointPlane::MkJointPlane(double dip, double dipdir)
{
  SetRotate(-dip, 0, -dipdir);
  Dip = dip;
  DipDir = dipdir;
}

void MkJointPlane::SetDipNDir(double dip, double dipdir)
{
  SetRotate(-dip, 0, -dipdir);
  Dip = dip;
  DipDir = dipdir;
}

void MkJointPlane::CalcABCD() //Ax+By+Cz=D
{
  MkPoint rp(0, 0, 1);
  MkMatrix4<double> rm;

  rm.Identity();
  rm.Rotate(Alpha, Beta, Gamma);
  rp = rp * rm;
  rp.Unify();

  A = rp.X;
  B = rp.Y;
  C = rp.Z;
  D = A * FCenter.X + B * FCenter.Y + C * FCenter.Z;
}

bool MkJointPlane::operator==(MkJointPlane &j)
{
  bool flag;
  flag = MkPlane::operator==((MkPlane &)j);
  flag = flag && fabs(Dip - j.Dip) < EPS && fabs(DipDir - j.DipDir) < EPS &&
         fabs(Aperture - j.Aperture) < EPS && Sort == j.Sort && Form == j.Form &&
         Condition == j.Condition;
  return flag;
}

bool MkJointPlane::operator!=(MkJointPlane &j)
{
  return !operator==(j);
}

//---------------------------------------------------------------------------

MkPennyJoint::MkPennyJoint() : MkJointPlane()
{
}

MkPennyJoint::MkPennyJoint(double dip, double dipdir) : MkJointPlane(dip, dipdir)
{
}

bool MkPennyJoint::GetOutline(MkPolygon &poly, int div)
{
  double ang, x, y, z;
  poly.Initialize(div);
  for (int i = 0; i < div; i++)
  {
    ang = i / (div - 1) * 3.141592 / 180.0;
    x = FCenter.X + FRadius * cos(ang);
    y = FCenter.Y + FRadius * sin(ang);
    z = FCenter.Z;
    poly[i].SetPoint(x, y, z);
    poly[i].Rotate(Alpha, Beta, Gamma);
  }
  return true;
}

bool MkPennyJoint::IsIntersect(MkLine rl) // need to be tested
{
  double t;
  CalcABCD();
  if (rl.GetLength() < 0.0001)
    return false;
  t = (D - A * rl[0].X - B * rl[0].Y - C * rl[0].Z) / (A * rl.DeltaX() + B * rl.DeltaY() + C * rl.DeltaZ());
  if (t < 1 && t > 0)
  {
    double x, y, z, len;
    x = rl[0].X + t * rl.DeltaX();
    y = rl[0].Y + t * rl.DeltaY();
    z = rl[0].Z + t * rl.DeltaZ();

    len = (FCenter.X - x) * (FCenter.X - x) + (FCenter.Y - y) * (FCenter.Y - y) + (FCenter.Z - z) * (FCenter.Z - z);
    len = sqrt(len);

    if (len > FRadius)
      return false;
    else
      return true;
  }
  else
    return false;
}

bool MkPennyJoint::IsIntersect(MkTriangle &rt)
{
  int cnt = 0;
  for (int i = 0; i < 3; i++)
    if (IsIntersect(rt(i)))
      cnt++;
  if (cnt == 1)
    return false;
  return (cnt == 2);
}

bool MkPennyJoint::IsIntersect(MkPennyJoint &pj)
{
  double A1, B1, C1, D1;
  double A2, B2, C2, D2;
  double A3, B3, C3, D3;
  double x0, y0, z0;

  CalcABCD();
  pj.CalcABCD();

  A1 = A;
  B1 = B;
  C1 = C;
  D1 = D;
  A2 = pj.A;
  B2 = pj.B;
  C2 = pj.C;
  A3 = B1 * C2 - B2 * C1;
  B3 = C1 * A2 - C2 * A1;
  C3 = A1 * B2 - A2 * B1;
  D2 = pj.D;
  D3 = 0;

  double len = sqrt(A3 * A3 + B3 * B3 + C3 * C3);
  A3 = A3 / len;
  B3 = B3 / len;
  C3 = C3 / len;

  MkMatrix<double> mat(3, 3);
  MkVector<double> vec(3);
  MkVector<double> res(3);

  mat(0, 0) = A1;
  mat(0, 1) = B1;
  mat(0, 2) = C1;
  mat(1, 0) = A2;
  mat(1, 1) = B2;
  mat(1, 2) = C2;
  mat(2, 0) = A3;
  mat(2, 1) = B3;
  mat(2, 2) = C3;
  vec(0) = D1;
  vec(1) = D2;
  vec(2) = D3;

  res = operator!(mat) * vec; //check this out!

  x0 = res(0);
  y0 = res(1);
  z0 = res(2);

  double t[4] = {10, 0, 0, 0};
  MkPoint sp(x0, y0, z0);
  MkPoint ep = MkPoint(x0 + t[0] * A3, y0 + t[0] * B3, z0 + t[0] * C3);

  MkLine rl(sp, ep);
  rl.SetFiniteness(false);

  MkSphere sp1(FCenter, FRadius), sp2(pj.FCenter, pj.FRadius);
  if (!((sp1 && rl) && (sp2 && rl)))
    return false;
  sp1.GetIntParam(rl, t[0], t[1]);
  sp2.GetIntParam(rl, t[2], t[3]);
  if (t[0] > t[1])
    swap(t[0], t[1]);
  if (t[2] > t[3])
    swap(t[2], t[3]);

  if (t[0] > t[3])
    return false;
  else if (t[2] > t[1])
    return false;
  else
    return true;
}

// ����� ������ ���� ���ǰ˻��� ���� �ʰ� ���� ����̶�� �����Ͽ���
// ������ ����� �������θ� �����Ѵ�.
bool MkPennyJoint::IsIntersect(MkPlane &rp)
{
  double A1, B1, C1, D1;
  double A2, B2, C2, D2;
  double A3, B3, C3, D3;
  double x0, y0, z0;

  CalcABCD();
  rp.CalcABCD();

  A1 = A;
  B1 = B;
  C1 = C;
  D1 = D;
  A2 = rp.GetA();
  B2 = rp.GetB();
  C2 = rp.GetC();
  D2 = rp.GetD();
  A3 = B1 * C2 - B2 * C1;
  B3 = C1 * A2 - C2 * A1;
  C3 = A1 * B2 - A2 * B1;
  D3 = 0;

  double len = sqrt(A3 * A3 + B3 * B3 + C3 * C3);
  A3 = A3 / len;
  B3 = B3 / len;
  C3 = C3 / len;

  MkMatrix<double> mat(3, 3);
  MkVector<double> vec(3);
  MkVector<double> res(3);

  mat(0, 0) = A1;
  mat(0, 1) = B1;
  mat(0, 2) = C1;
  mat(1, 0) = A2;
  mat(1, 1) = B2;
  mat(1, 2) = C2;
  mat(2, 0) = A3;
  mat(2, 1) = B3;
  mat(2, 2) = C3;
  vec(0) = D1;
  vec(1) = D2;
  vec(2) = D3;

  res = operator!(mat) * vec; //check this out!

  x0 = res(0);
  y0 = res(1);
  z0 = res(2);

  double t[4] = {10, 0, 0, 0};
  MkPoint sp(x0, y0, z0);
  MkPoint ep = MkPoint(x0 + t[0] * A3, y0 + t[0] * B3, z0 + t[0] * C3);

  MkLine rl(sp, ep);
  rl.SetFiniteness(false);

  MkSphere sp1(FCenter, FRadius);
  if (!((sp1 && rl) && (rp.IsInPlane(rl))))
    return false;

  sp1.GetIntParam(rl, t[0], t[1]);
  if (!rp.GetIntParam(rl, t[2], t[3]))
    return false;

  if (t[0] > t[1])
    swap(t[0], t[1]);
  if (t[2] > t[3])
    swap(t[2], t[3]);

  if (t[0] > t[3])
    return false;
  else if (t[2] > t[1])
    return false;
  else
    return true;
}

bool MkPennyJoint::IsIntersect(MkCylinder &rc)
{
  double x0, y0, z0;
  double x1, y1, z1;
  double len;
  double A, B, C, D;

  MkPennyJoint pj = *this;
  MkPoint cen = FCenter;
  rc.RotateSpace(pj);
  rc.RotateSpace(cen);

  x0 = cen.X;
  y0 = cen.Y;
  z0 = cen.Z;

  len = sqrt(x0 * x0 + y0 * y0);
  x1 = x0 * rc.GetRadius() / len;
  y1 = y0 * rc.GetRadius() / len;

  pj.CalcABCD();
  A = pj.GetA();
  B = pj.GetB();
  C = pj.GetC();
  D = pj.GetD();
  if (fabs(C) < 0.001)
    z1 = z0;
  else
    z1 = (D - A * x1 - B * y1) / C;

  len = sqrt((x0 - x1) * (x0 - x1) + (y0 - y1) * (y0 - y1) + (z0 - z1) * (z0 - z1));
  if (len < pj.FRadius)
    return true;
  else
    return false;
}

MkPoint &MkPennyJoint::CalcIntPnt(MkLine &rl)
{
  double t;
  static MkPoint mp;
  CalcABCD();
  if (rl.GetLength() < 0.0001)
    return NullPoint;
  t = (D - A * rl[0].X - B * rl[0].Y - C * rl[0].Z) / (A * rl.DeltaX() + B * rl.DeltaY() + C * rl.DeltaZ());
  if (t < 1 && t > 0)
  {
    double x, y, z, len;
    x = rl[0].X + t * rl.DeltaX();
    y = rl[0].Y + t * rl.DeltaY();
    z = rl[0].Z + t * rl.DeltaZ();

    len = (FCenter.X - x) * (FCenter.X - x) + (FCenter.Y - y) * (FCenter.Y - y) + (FCenter.Z - z) * (FCenter.Z - z);
    len = sqrt(len);

    if (len > FRadius)
      return NullPoint;
    else
    {
      mp.SetPoint(x, y, z);
      return mp;
    }
  }
  else
    return NullPoint;
}

MkLine &MkPennyJoint::CalcIntLine(MkTriangle &rt)
{
  MkPoint rp[2];
  MkLine rl;
  int cnt = 0;
  for (int i = 0; i < 3; i++)
    if (IsIntersect(rt(i)))
    {
      rp[cnt] = CalcIntPnt(rt(i));
      cnt++;
    }
  if (cnt == 2)
  {
    FLine.SetLine(rp[0], rp[1]);
    return FLine;
  }
  else
    return NullLine;
}

MkPoint &MkPennyJoint::CalcIntPnt(MkTriangle &rt) //Shortest point
{
  double dist;
  static MkPoint rp;
  rp = MkPlane::CalcIntPnt(rt);
  dist = CalDist(FCenter, rp);
  if (dist > FRadius)
    return NullPoint;
  else
    return rp;
}

MkPoint &MkPennyJoint::CalcIntPnt(MkTriangles &rts) //Shortest point
{
  double dist;
  static MkPoint rp;
  rp = MkPlane::CalcIntPnt(rts);
  dist = CalDist(FCenter, rp);
  if (dist > FRadius)
    return NullPoint;
  else
    return rp;
}

MkLine &MkPennyJoint::CalcIntLine(MkPennyJoint &pj)
{
  double A1, B1, C1, D1;
  double A2, B2, C2, D2;
  double A3, B3, C3, D3;
  double x0, y0, z0;

  CalcABCD();
  pj.CalcABCD();

  A1 = A;
  B1 = B;
  C1 = C;
  D1 = D;
  A2 = pj.A;
  B2 = pj.B;
  C2 = pj.C;
  D2 = pj.D;
  A3 = B1 * C2 - B2 * C1;
  B3 = C1 * A2 - C2 * A1;
  C3 = A1 * B2 - A2 * B1;
  D3 = 0;

  double len = sqrt(A3 * A3 + B3 * B3 + C3 * C3);
  A3 = A3 / len;
  B3 = B3 / len;
  C3 = C3 / len;

  MkMatrix<double> mat(3, 3);
  MkVector<double> vec(3);
  MkVector<double> res(3);

  mat(0, 0) = A1;
  mat(0, 1) = B1;
  mat(0, 2) = C1;
  mat(1, 0) = A2;
  mat(1, 1) = B2;
  mat(1, 2) = C2;
  mat(2, 0) = A3;
  mat(2, 1) = B3;
  mat(2, 2) = C3;
  vec(0) = D1;
  vec(1) = D2;
  vec(2) = D3;

  res = operator!(mat) * vec; //check this out!

  x0 = res(0);
  y0 = res(1);
  z0 = res(2);

  double t[4] = {10, 0, 0, 0};
  MkPoint sp(x0, y0, z0);
  MkPoint ep = MkPoint(x0 + t[0] * A3, y0 + t[0] * B3, z0 + t[0] * C3);

  MkLine rl(sp, ep);
  rl.SetFiniteness(false);

  MkSphere sp1(FCenter, FRadius), sp2(pj.FCenter, pj.FRadius);
  if (!((sp1 && rl) && (sp2 && rl)))
    return NullLine;
  sp1.GetIntParam(rl, t[0], t[1]);
  sp2.GetIntParam(rl, t[2], t[3]);
  if (t[0] > t[1])
    swap(t[0], t[1]);
  if (t[2] > t[3])
    swap(t[2], t[3]);

  if (t[0] > t[3])
    return NullLine;
  else if (t[2] > t[1])
    return NullLine;

  for (int i = 0; i < 4; i++)
  {
    for (int j = i + 1; j < 4; j++)
      if (t[i] > t[j])
        swap(t[i], t[j]);
  } // sorting...^^

  sp.SetPoint(x0 + t[1] * A3, y0 + t[1] * B3, z0 + t[1] * C3);
  ep.SetPoint(x0 + t[2] * A3, y0 + t[2] * B3, z0 + t[2] * C3);
  FLine.SetLine(sp, ep);
  return FLine;
}

MkLine &MkPennyJoint::CalcIntLine(MkPlane &rp)
{
  double A1, B1, C1, D1;
  double A2, B2, C2, D2;
  double A3, B3, C3, D3;
  double x0, y0, z0;
  double x1, y1, z1;
  double x2, y2, z2;
  double t[4] = {10, 0, 0, 0};

  CalcABCD();
  rp.CalcABCD();

  A1 = A;
  B1 = B;
  C1 = C;
  D1 = D;
  A2 = rp.GetA();
  B2 = rp.GetB();
  C2 = rp.GetC();
  D2 = rp.GetD();
  A3 = B1 * C2 - B2 * C1;
  B3 = C1 * A2 - C2 * A1;
  C3 = A1 * B2 - A2 * B1;
  D3 = 0;

  double len = sqrt(A3 * A3 + B3 * B3 + C3 * C3);
  A3 = A3 / len;
  B3 = B3 / len;
  C3 = C3 / len;

  MkMatrix<double> mat(3, 3);
  MkVector<double> vec(3);
  MkVector<double> res(3);

  mat(0, 0) = A1;
  mat(0, 1) = B1;
  mat(0, 2) = C1;
  mat(1, 0) = A2;
  mat(1, 1) = B2;
  mat(1, 2) = C2;
  mat(2, 0) = A3;
  mat(2, 1) = B3;
  mat(2, 2) = C3;
  vec(0) = D1;
  vec(1) = D2;
  vec(2) = D3;

  res = operator!(mat) * vec; //check this out!

  x0 = res(0);
  y0 = res(1);
  z0 = res(2);

  x1 = x0 - t[0] * A3;
  y1 = y0 - t[0] * B3;
  z1 = z0 - t[0] * C3;

  x2 = x0 + t[0] * A3;
  y2 = y0 + t[0] * B3;
  z2 = z0 + t[0] * C3;

  MkPoint sp(x1, y1, z1);
  MkPoint ep(x2, y2, z2);

  MkLine rl(sp, ep);
  rl.SetFiniteness(false);

  MkSphere sp1(FCenter, FRadius);
  if (!((sp1 && rl) && (rp.IsInPlane(rl))))
    return NullLine;

  sp1.GetIntParam(rl, t[0], t[1]);
  rp.GetIntParam(rl, t[2], t[3]);
  if (t[0] > t[1])
    swap(t[0], t[1]);
  if (t[2] > t[3])
    swap(t[2], t[3]);

  if (t[0] > t[3])
    return NullLine;
  else if (t[2] > t[1])
    return NullLine;

  for (int i = 0; i < 4; i++)
  {
    for (int j = i + 1; j < 4; j++)
      if (t[i] > t[j])
        swap(t[i], t[j]);
  } // sorting...^^

  sp.SetPoint(x1 + t[1] * A3, y1 + t[1] * B3, z1 + t[1] * C3);
  ep.SetPoint(x1 + t[2] * A3, y1 + t[2] * B3, z1 + t[2] * C3);
  FLine.SetLine(sp, ep);
  return FLine;
}

MkPoint &MkPennyJoint::CalcIntPnt(MkCylinder &rc)
{
  double x0, y0, z0;
  double x1, y1, z1;
  double len;
  double A, B, C, D;

  MkPennyJoint pj = *this;
  MkPoint cen = FCenter;
  rc.RotateSpace(pj);
  rc.RotateSpace(cen);

  x0 = cen.X;
  y0 = cen.Y;
  z0 = cen.Z;

  len = sqrt(x0 * x0 + y0 * y0);
  x1 = x0 * rc.GetRadius() / len;
  y1 = y0 * rc.GetRadius() / len;

  pj.CalcABCD();
  A = pj.GetA();
  B = pj.GetB();
  C = pj.GetC();
  D = pj.GetD();
  if (fabs(C) < 0.001)
    z1 = z0;
  else
    z1 = (D - A * x1 - B * y1) / C;

  static MkPoint rp;
  rp.SetPoint(x1, y1, z1);
  rc.UnRotateSpace(rp);

  len = sqrt((x0 - x1) * (x0 - x1) + (y0 - y1) * (y0 - y1) + (z0 - z1) * (z0 - z1));
  if (len < pj.FRadius)
    return rp;
  else
    return NullPoint;
}

bool MkPennyJoint::operator==(MkPennyJoint &j)
{
  bool flag;
  flag = MkJointPlane::operator==((MkJointPlane &)j);
  flag = flag && fabs(FRadius - j.FRadius) < EPS;
  return flag;
}

bool MkPennyJoint::operator!=(MkPennyJoint &j)
{
  return !operator==(j);
}

#ifdef __BCPLUSPLUS__
void MkPennyJoint::Draw(TObject *Sender)
{
  TColor C;
  TPenStyle PS;
  MkPolygon poly;
  if (String(Sender->ClassName()) == String("MkPaintBox"))
  {
    MkPaintBox *pb = (MkPaintBox *)Sender;
    C = pb->Canvas->Pen->Color;
    PS = pb->Canvas->Pen->Style;
    pb->Canvas->Pen->Color = Color;
    //      pb->Canvas->Pen->Style = PenStyle;

    GetOutline(poly, 12);
    poly.Draw(Sender);

    pb->Canvas->Pen->Color = C;
    pb->Canvas->Pen->Style = PS;
  }
}
#endif

#if defined(_MSC_VER) && defined(_WINDOWS_)
void MkPennyJoint::Draw(MkPaint *pb)
{
  MkPolygon poly;

  GetOutline(poly, 12);
  poly.Draw(pb);
}
#endif

//---------------------------------------------------------------------------

// MkPlanes::MkPlanes(int size, MkPlane *jp)
// {
//   if (size < 0)
//   {
//     MkDebug("::MkPlanes - MkPlanes(int size)");
//     return;
//   }

//   FSize = size;
//   if (FSize == 0)
//   {
//     FPlane = NULL;
//     return;
//   }

//   FPlane = new MkPlane[FSize];
//   for (int i = 0; i < FSize; i++)
//     FPlane[i] = jp[i];
// }

// MkPlanes::MkPlanes(int size)
// {
//   if (size < 0)
//   {
//     MkDebug("::MkPlanes - MkPlanes(int size)");
//     return;
//   }

//   FSize = size;
//   if (FSize == 0)
//   {
//     FPlane = NULL;
//     return;
//   }

//   FPlane = new MkPlane[FSize];
// }

// MkPlanes::~MkPlanes()
// {
//   FSize = 0;
//   if (FPlane)
//   {
//     delete[](MkPlane *) FPlane;
//     FPlane = NULL;
//   }
// }

// void MkPlanes::Initialize(int size)
// {
//   if (size < 0)
//   {
//     MkDebug("::MkPlanes - Initialize(int size)");
//     ;
//     return;
//   }
//   if (FSize == size)
//     return;

//   FSize = size;
//   if (FSize == 0)
//   {
//     if (FPlane != NULL)
//       delete[](MkPlane *) FPlane;
//     FPlane = NULL;
//     return;
//   }

//   if (FPlane != NULL)
//     delete[](MkPlane *) FPlane;
//   FPlane = new MkPlane[FSize];
// }

// void MkPlanes::Initialize(int size, MkPlane *jp)
// {
//   if (size < 0)
//   {
//     MkDebug("::MkPlanes - Initialize(int size)");
//     ;
//     return;
//   }

//   FSize = size;
//   if (FSize == 0)
//   {
//     if (FPlane != NULL)
//       delete[](MkPlane *) FPlane;
//     FPlane = NULL;
//     return;
//   }

//   if (FPlane != NULL)
//     delete[](MkPlane *) FPlane;
//   FPlane = new MkPlane[FSize];
//   for (int i = 0; i < FSize; i++)
//     FPlane[i] = jp[i];
// }

// bool MkPlanes::Clear()
// {
//   FSize = 0;
//   if (FPlane)
//   {
//     delete[](MkPlane *) FPlane;
//     FPlane = NULL;
//   }
//   return true;
// }

// MkPlane &MkPlanes::operator[](int i)
// {
//   if (FSize == 0)
//     return NullPlane;
//   else if (i >= 0 && i < FSize)
//     return FPlane[i];
//   else
//     return NullPlane;
// }

// MkPlanes &MkPlanes::operator=(MkPlanes &jps)
// {
//   int i;

//   Clear();
//   FSize = jps.FSize;
//   if (FSize == 0)
//   {
//     this->FPlane = NULL;
//     return *this;
//   }
//   this->FPlane = new MkPlane[FSize];

//   for (i = 0; i < FSize; i++)
//     this->FPlane[i] = jps.FPlane[i];

//   return *this;
// }

// bool MkPlanes::operator==(MkPlanes &Reals)
// {
//   int i;

//   if (FSize != Reals.FSize)
//     return false;
//   for (i = 0; i < FSize; i++)
//     if (this->FPlane[i] != Reals.FPlane[i])
//       return false;

//   return true;
// }

// void MkPlanes::Translate(double x, double y, double z)
// {
//   for (int i = 0; i < GetSize(); i++)
//     FPlane[i].Translate(x, y, z);
// }

// void MkPlanes::Translate(MkPoint rp)
// {
//   for (int i = 0; i < GetSize(); i++)
//     FPlane[i].Translate(rp);
// }

// void MkPlanes::Rotate(double alpha, double beta, double gamma)
// {
//   for (int i = 0; i < GetSize(); i++)
//     FPlane[i].Rotate(alpha, beta, gamma);
// }

// void MkPlanes::RotateInX(double ang)
// {
//   for (int i = 0; i < GetSize(); i++)
//     FPlane[i].RotateInX(ang);
// }

// void MkPlanes::RotateInY(double ang)
// {
//   for (int i = 0; i < GetSize(); i++)
//     FPlane[i].RotateInY(ang);
// }

// void MkPlanes::RotateInZ(double ang)
// {
//   for (int i = 0; i < GetSize(); i++)
//     FPlane[i].RotateInZ(ang);
// }

// void MkPlanes::RotateInA(double ang, double l, double m, double n)
// {
//   for (int i = 0; i < GetSize(); i++)
//     FPlane[i].RotateInA(ang, l, m, n);
// }

// #ifdef __BCPLUSPLUS__
// void MkPlanes::Draw(TObject *Sender)
// {
//   for (int i = 0; i < GetSize(); i++)
//     FPlane[i].Draw(Sender);
// }
// #endif

// #if defined(_MSC_VER) && defined(_WINDOWS_)
// void MkPlanes::Draw(MkPaint *pb)
// {
//   for (int i = 0; i < GetSize(); i++)
//     FPlane[i].Draw(pb);
// }
// #endif

// //----------------------------------------------------------------
// MkPointsPlanes::MkPointsPlanes(int size, MkPointsPlane *jp)
// {
//   if (size < 0)
//   {
//     MkDebug("::MkPointsPlanes - MkPointsPlanes(int size)");
//     return;
//   }

//   FSize = size;
//   if (FSize == 0)
//   {
//     FPoints = NULL;
//     return;
//   }

//   FPoints = new MkPointsPlane[FSize];
//   for (int i = 0; i < FSize; i++)
//     FPoints[i] = jp[i];
// }

// MkPointsPlanes::MkPointsPlanes(int size)
// {
//   if (size < 0)
//   {
//     MkDebug("::MkPointsPlanes - MkPointsPlanes(int size)");
//     return;
//   }

//   FSize = size;
//   if (FSize == 0)
//   {
//     FPoints = NULL;
//     return;
//   }

//   FPoints = new MkPointsPlane[FSize];
// }

// MkPointsPlanes::~MkPointsPlanes()
// {
//   FSize = 0;
//   if (FPoints)
//   {
//     delete[](MkPointsPlane *) FPoints;
//     FPoints = NULL;
//   }
// }

// void MkPointsPlanes::Initialize(int size)
// {
//   if (size < 0)
//   {
//     MkDebug("::MkPointsPlanes - Initialize(int size)");
//     ;
//     return;
//   }
//   if (FSize == size)
//     return;

//   FSize = size;
//   if (FSize == 0)
//   {
//     if (FPoints != NULL)
//       delete[](MkPointsPlane *) FPoints;
//     FPoints = NULL;
//     return;
//   }

//   if (FPoints != NULL)
//     delete[](MkPointsPlane *) FPoints;
//   FPoints = new MkPointsPlane[FSize];
// }

// void MkPointsPlanes::Initialize(int size, MkPointsPlane *jp)
// {
//   if (size < 0)
//   {
//     MkDebug("::MkPointsPlanes - Initialize(int size)");
//     ;
//     return;
//   }

//   FSize = size;
//   if (FSize == 0)
//   {
//     if (FPoints != NULL)
//       delete[](MkPointsPlane *) FPoints;
//     FPoints = NULL;
//     return;
//   }

//   if (FPoints != NULL)
//     delete[](MkPointsPlane *) FPoints;
//   FPoints = new MkPointsPlane[FSize];
//   for (int i = 0; i < FSize; i++)
//     FPoints[i] = jp[i];
// }

// bool MkPointsPlanes::Clear()
// {
//   FSize = 0;
//   if (FPoints)
//   {
//     delete[](MkPointsPlane *) FPoints;
//     FPoints = NULL;
//   }
//   return true;
// }

// MkPointsPlane &MkPointsPlanes::operator[](int i)
// {
//   if (FSize == 0)
//     return NullPointPlane;
//   else if (i >= 0 && i < FSize)
//     return FPoints[i];
//   else
//     return NullPointPlane;
// }

// MkPointsPlanes &MkPointsPlanes::operator=(MkPointsPlanes &pps)
// {
//   int i;

//   Clear();
//   FSize = pps.FSize;
//   if (FSize == 0)
//   {
//     this->FPoints = NULL;
//     return *this;
//   }
//   this->FPoints = new MkPointsPlane[FSize];

//   for (i = 0; i < FSize; i++)
//     this->FPoints[i] = pps.FPoints[i];

//   return *this;
// }

// bool MkPointsPlanes::operator==(MkPointsPlanes &planes)
// {
//   int i;

//   if (FSize != planes.FSize)
//     return false;
//   for (i = 0; i < FSize; i++)
//     if (this->FPoints[i] != planes.FPoints[i])
//       return false;

//   return true;
// }

// #ifdef __BCPLUSPLUS__
// void MkPointsPlanes::Draw(TObject *Sender)
// {
//   for (int i = 0; i < FSize; i++)
//     FPoints[i].Draw(Sender);
// }
// #endif

// #if defined(_MSC_VER) && defined(_WINDOWS_)
// void MkPointsPlanes::Draw(MkPaint *pb)
// {
//   for (int i = 0; i < FSize; i++)
//     FPoints[i].Draw(pb);
// }
// #endif
// //--------------------------------------------------------------
// MkJointPlanes::MkJointPlanes(int size, MkJointPlane *jp)
// {
//   if (size < 0)
//   {
//     MkDebug("::MkJointPlanes - MkJointPlanes(int size)");
//     return;
//   }

//   FSize = size;
//   if (FSize == 0)
//   {
//     FJoint = NULL;
//     return;
//   }

//   FJoint = new MkJointPlane[FSize];
//   for (int i = 0; i < FSize; i++)
//     FJoint[i] = jp[i];
// }

// MkJointPlanes::MkJointPlanes(int size)
// {
//   if (size < 0)
//   {
//     MkDebug("::MkJointPlanes - MkJointPlanes(int size)");
//     return;
//   }

//   FSize = size;
//   if (FSize == 0)
//   {
//     FJoint = NULL;
//     return;
//   }

//   FJoint = new MkJointPlane[FSize];
// }

// MkJointPlanes::~MkJointPlanes()
// {
//   FSize = 0;
//   if (FJoint)
//   {
//     delete[](MkJointPlane *) FJoint;
//     FJoint = NULL;
//   }
// }

// void MkJointPlanes::Initialize(int size)
// {
//   if (size < 0)
//   {
//     MkDebug("::MkJointPlanes - Initialize(int size)");
//     ;
//     return;
//   }
//   if (FSize == size)
//     return;

//   FSize = size;
//   if (FSize == 0)
//   {
//     if (FJoint != NULL)
//       delete[](MkJointPlane *) FJoint;
//     FJoint = NULL;
//     return;
//   }

//   if (FJoint != NULL)
//     delete[](MkJointPlane *) FJoint;
//   FJoint = new MkJointPlane[FSize];
// }

// void MkJointPlanes::Initialize(int size, MkJointPlane *jp)
// {
//   if (size < 0)
//   {
//     MkDebug("::MkJointPlanes - Initialize(int size)");
//     ;
//     return;
//   }

//   FSize = size;
//   if (FSize == 0)
//   {
//     if (FJoint != NULL)
//       delete[](MkJointPlane *) FJoint;
//     FJoint = NULL;
//     return;
//   }

//   if (FJoint != NULL)
//     delete[](MkJointPlane *) FJoint;
//   FJoint = new MkJointPlane[FSize];
//   for (int i = 0; i < FSize; i++)
//     FJoint[i] = jp[i];
// }

// bool MkJointPlanes::Clear()
// {
//   FSize = 0;
//   if (FJoint)
//   {
//     delete[](MkJointPlane *) FJoint;
//     FJoint = NULL;
//   }
//   return true;
// }

// MkJointPlane &MkJointPlanes::operator[](int i)
// {
//   if (FSize == 0)
//     return NullJoint;
//   else if (i >= 0 && i < FSize)
//     return FJoint[i];
//   else
//     return NullJoint;
// }

// MkJointPlanes &MkJointPlanes::operator=(MkJointPlanes &jps)
// {
//   int i;

//   Clear();
//   FSize = jps.FSize;
//   if (FSize == 0)
//   {
//     this->FJoint = NULL;
//     return *this;
//   }
//   this->FJoint = new MkJointPlane[FSize];

//   for (i = 0; i < FSize; i++)
//     this->FJoint[i] = jps.FJoint[i];

//   return *this;
// }

// bool MkJointPlanes::operator==(MkJointPlanes &joints)
// {
//   int i;

//   if (FSize != joints.FSize)
//     return false;
//   for (i = 0; i < FSize; i++)
//     if (this->FJoint[i] != joints.FJoint[i])
//       return false;

//   return true;
// }
// void MkJointPlanes::Translate(double x, double y, double z)
// {
//   for (int i = 0; i < GetSize(); i++)
//     FJoint[i].Translate(x, y, z);
// }

// void MkJointPlanes::Translate(MkPoint rp)
// {
//   for (int i = 0; i < GetSize(); i++)
//     FJoint[i].Translate(rp);
// }

// #ifdef __BCPLUSPLUS__
// void MkJointPlanes::Draw(TObject *Sender)
// {
//   for (int i = 0; i < GetSize(); i++)
//     FJoint[i].Draw(Sender);
// }
// #endif

// #if defined(_MSC_VER) && defined(_WINDOWS_)
// void MkJointPlanes::Draw(MkPaint *pb)
// {
//   for (int i = 0; i < GetSize(); i++)
//     FJoint[i].Draw(pb);
// }
// #endif

// //---------------------------------------------------------------------------
// MkPennyJoints::MkPennyJoints(int size, MkPennyJoint *jp)
// {
//   if (size < 0)
//   {
//     MkDebug("::MkPennyJoints - MkPennyJoints(int size)");
//     return;
//   }

//   FSize = size;
//   if (FSize == 0)
//   {
//     FPenny = NULL;
//     return;
//   }

//   FPenny = new MkPennyJoint[FSize];
//   for (int i = 0; i < FSize; i++)
//     FPenny[i] = jp[i];
// }

// MkPennyJoints::MkPennyJoints(int size)
// {
//   if (size < 0)
//   {
//     MkDebug("::MkPennyJoints - MkPennyJoints(int size)");
//     return;
//   }

//   FSize = size;
//   if (FSize == 0)
//   {
//     FPenny = NULL;
//     return;
//   }

//   FPenny = new MkPennyJoint[FSize];
// }

// MkPennyJoints::~MkPennyJoints()
// {
//   FSize = 0;
//   if (FPenny)
//   {
//     delete[](MkPennyJoint *) FPenny;
//     FPenny = NULL;
//   }
// }

// void MkPennyJoints::Initialize(int size)
// {
//   if (size < 0)
//   {
//     MkDebug("::MkPennyJoints - Initialize(int size)");
//     ;
//     return;
//   }
//   if (FSize == size)
//     return;

//   FSize = size;
//   if (FSize == 0)
//   {
//     if (FPenny != NULL)
//       delete[](MkPennyJoint *) FPenny;
//     FPenny = NULL;
//     return;
//   }

//   if (FPenny != NULL)
//     delete[](MkPennyJoint *) FPenny;
//   FPenny = new MkPennyJoint[FSize];
// }

// void MkPennyJoints::Initialize(int size, MkPennyJoint *jp)
// {
//   if (size < 0)
//   {
//     MkDebug("::MkPennyJoints - Initialize(int size)");
//     ;
//     return;
//   }

//   FSize = size;
//   if (FSize == 0)
//   {
//     if (FPenny != NULL)
//       delete[](MkPennyJoint *) FPenny;
//     FPenny = NULL;
//     return;
//   }

//   if (FPenny != NULL)
//     delete[](MkPennyJoint *) FPenny;
//   FPenny = new MkPennyJoint[FSize];
//   for (int i = 0; i < FSize; i++)
//     FPenny[i] = jp[i];
// }

// bool MkPennyJoints::Clear()
// {
//   FSize = 0;
//   if (FPenny)
//   {
//     delete[](MkPennyJoint *) FPenny;
//     FPenny = NULL;
//   }
//   return true;
// }

// MkPennyJoint &MkPennyJoints::operator[](int i)
// {
//   if (FSize == 0)
//     return NullPenny;
//   else if (i >= 0 && i < FSize)
//     return FPenny[i];
//   else
//     return NullPenny;
// }

// MkPennyJoints &MkPennyJoints::operator=(MkPennyJoints &jps)
// {
//   int i;

//   Clear();
//   FSize = jps.FSize;
//   if (FSize == 0)
//   {
//     FPenny = NULL;
//     return *this;
//   }
//   FPenny = new MkPennyJoint[FSize];

//   for (i = 0; i < FSize; i++)
//     FPenny[i] = jps.FPenny[i];

//   return *this;
// }

// bool MkPennyJoints::operator==(MkPennyJoints &joints)
// {
//   int i;

//   if (FSize != joints.FSize)
//     return false;
//   for (i = 0; i < FSize; i++)
//     if (this->FPenny[i] != joints.FPenny[i])
//       return false;

//   return true;
// }

// #ifdef __BCPLUSPLUS__
// void MkPennyJoints::Draw(TObject *Sender)
// {
//   for (int i = 0; i < GetSize(); i++)
//     FPenny[i].Draw(Sender);
// }
// #endif

// #if defined(_MSC_VER) && defined(_WINDOWS_)
// void MkPennyJoints::Draw(MkPaint *pb)
// {
//   for (int i = 0; i < GetSize(); i++)
//     FPenny[i].Draw(pb);
// }
// #endif
//---------------------------------------------------------------------------
