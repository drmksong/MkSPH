//---------------------------------------------------------------------------
// Author : M. K. Song. Seoul, Korea. 1999.2.20
//---------------------------------------------------------------------------
#include "MkParticle.hpp"
//---------------------------------------------------------------------------

MkParticle NullParticle(0, 0, 0);
MkParticles NullParticles(0);
//double  Zero = 0;

void MkParticle::Initialize()
{
    R = G = B = 0;
    Radius = 0;
    Density = Density1 = 0; //??
    Mass = 0;
    Viscosity = 0;
    Temp = Temp1 = 0; //??
    Pressure = Pressure1 = 0;
    XForceInt = YForceInt = 0;
    XForceVis = YForceVis = 0;
    XForceArtVis = YForceArtVis = 0;
    XForceExt = YForceExt = 0;
    XAccel[0] = YAccel[0] = XAccel[1] = YAccel[1] = XAccel[2] = YAccel[2] = 0;                   // 1 : t-dt,  2 : t, 3 : t+dt
    XVelocity[0] = YVelocity[0] = XVelocity[1] = YVelocity[1] = XVelocity[2] = YVelocity[2] = 0; // 1 : t-dt, 2 : t, 3 : t+dt
    XPosition[3] = YPosition[3] = XPosition[1] = YPosition[1] = XPosition[2] = YPosition[2] = 0; // 1 : t-dt, 2 : t, 3 : t+dt
    GridNum = 0;
    isFixed = false;
}

void MkParticle::SetupTempColor()
{
    R = (Temp1 > 50) ? (Temp1 - 50) / 50.0 : 0;
    B = (Temp1 < 50) ? 1 - Temp1 / 50.0 : 0;
    G = (Temp1 < 50) ? Temp1 / 50.0 : 1 - (Temp1 - 50.0) / 50.0;
}

//#if defined(__GL_H__) && defined(_WINDOWS_)
#if defined(__GL_H__)
void MkParticle::Draw()
{
    int c;
    //  double  r,g,b;
    //  c = Color || 0x0000FF;
    //  r = double (c)/0xFF;
    //  c = Color || 0x00FF00;
    //  c >>= 8;
    //  g = double (c)/0xFF;
    //  c = Color || 0xFF0000;
    //  c >>= 16;
    //  b = double (c)/0xFF;
    //  printf("R:%f, G:%f, B:%f\n",R,G,B);
    glPushMatrix();
    glTranslatef(X, Y, Z);
    glColor3f(R, G, B);
    glutSolidSphere(0.05, 16, 16);
    //auxSoildSphere(2.0f);
    glPopMatrix();
    //  printf("particle (%10.3f,%10.3f,%10.3f) is drawing\n",X,Y,Z);
}
#endif

MkParticle &operator*(MkParticle &rp, MkMatrix4<double> &rm)
{
    static MkParticle rp_t;
    rp_t = rp;

    rp_t.X = rp.X * rm(0, 0) + rp.Y * rm(0, 1) + rp.Z * rm(0, 2) + 1 * rm(0, 3);
    rp_t.Y = rp.X * rm(1, 0) + rp.Y * rm(1, 1) + rp.Z * rm(1, 2) + 1 * rm(1, 3);
    rp_t.Z = rp.X * rm(2, 0) + rp.Y * rm(2, 1) + rp.Z * rm(2, 2) + 1 * rm(2, 3);

    return rp_t;
}

MkParticle &operator*(MkParticle &rp, double f)
{
    static MkParticle rp_t;
    rp_t = rp;

    rp_t.X *= f;
    rp_t.Y *= f;
    rp_t.Z *= f;
    return rp_t;
}

MkParticles::MkParticles(int size, MkParticle *rps)
{

    if (size < 0)
    {
        MkDebug("::MkParticles - MkParticles(int size)");
        ;
        return;
    }

    FSizeOfArray = FSize = size;
    if (FSize == 0)
    {
        FParticle = NULL;
        return;
    }

    FParticle = new MkParticle[FSize];
    for (int i = 0; i < FSize; i++)
        (*this)[i] = rps[i];
}

MkParticles::MkParticles(int size)
{
    if (size < 0)
    {
        MkDebug("::MkParticles - MkParticles(int size)");
        return;
    }

    FSize = FSizeOfArray = size;

    if (FSizeOfArray == 0)
    {
        FParticle = NULL;
        return;
    }

    FParticle = new MkParticle[FSizeOfArray];
}

MkParticles::~MkParticles()
{
    FSizeOfArray = FSize = 0;
    if (FParticle)
    {
        delete[] FParticle;
        FParticle = NULL;
    }
}

void MkParticles::Initialize(int size)
{
    if (size < 0)
    {
        MkDebug("::MkParticles - Initialize(int size)");
        ;
        return;
    }
    if (FSizeOfArray == size)
        return;

    FSize = FSizeOfArray = size;

    if (FSizeOfArray == 0)
    {
        if (FParticle != NULL)
            delete[](MkParticle *) FParticle;
        FParticle = NULL;
        return;
    }

    if (FParticle != NULL)
        delete[](MkParticle *) FParticle;
    FParticle = new MkParticle[FSizeOfArray];
}

void MkParticles::Initialize(int size, MkParticle *rps)
{

    if (size < 0 || rps == NULL)
    {
        MkDebug("::MkParticles - Initialize(int size)");
        ;
        return;
    }
    if (FSizeOfArray == size)
    {
        FSize = size;
        for (int i = 0; i < FSizeOfArray; i++)
            FParticle[i] = rps[i];
        return;
    }

    FSize = FSizeOfArray = size;
    if (FSizeOfArray == 0)
    {
        if (FParticle != NULL)
            delete[](MkParticle *) FParticle;
        FParticle = NULL;
        return;
    }

    if (FParticle != NULL)
        delete[](MkParticle *) FParticle;
    FParticle = new MkParticle[FSizeOfArray];
    for (int i = 0; i < FSizeOfArray; i++)
        FParticle[i] = rps[i];
}

void MkParticles::FindCenter() // find maximum and minimum for each axis and average it, the center
{
    if (FSize == 0)
    {
        FCenter = NullParticle;
        return;
    }

    double x1, y1, z1, x2, y2, z2;
    x1 = (*this)[0].X;
    y1 = (*this)[0].Y;
    z1 = (*this)[0].Z;
    x2 = (*this)[0].X;
    y2 = (*this)[0].Y;
    z2 = (*this)[0].Z;

    for (int i = 1; i < GetSize(); i++)
    {
        x1 = x1 < (*this)[i].X ? x1 : (*this)[i].X;
        y1 = y1 < (*this)[i].Y ? y1 : (*this)[i].Y;
        z1 = z1 < (*this)[i].Z ? z1 : (*this)[i].Z;
        x2 = x2 > (*this)[i].X ? x2 : (*this)[i].X;
        y2 = y2 > (*this)[i].Y ? y2 : (*this)[i].Y;
        z2 = z2 > (*this)[i].Z ? z2 : (*this)[i].Z;
    }
    FCenter.X = (x1 + x2) / 2.0;
    FCenter.Y = (y1 + y2) / 2.0;
    FCenter.Z = (z1 + z2) / 2.0;
}

int MkParticles::Grow(int delta)
{
    int i;
    MkParticle *rp = NULL;

    if (!(rp = new MkParticle[FSizeOfArray + delta]))
        return FSizeOfArray;

    for (i = 0; i < FSize; i++)
        rp[i] = FParticle[i];
    for (i = FSize; i < FSizeOfArray + delta; i++)
        rp[i] = NullParticle;
    if (FParticle)
    {
        delete[](MkParticle *) FParticle;
        FParticle = NULL;
    }
    FParticle = rp;
    FSizeOfArray = FSizeOfArray + delta;
    return FSizeOfArray;
}

int MkParticles::Shrink(int delta)
{
    int i;
    MkParticle *rp = NULL;

    if (!(rp = new MkParticle[FSizeOfArray - delta]))
        return FSizeOfArray;

    for (i = 0; i < FSize; i++)
        rp[i] = FParticle[i];
    for (i = FSize; i < FSizeOfArray - delta; i++)
        rp[i] = NullParticle;
    if (FParticle)
    {
        delete[](MkParticle *) FParticle;
        FParticle = NULL;
    }
    FParticle = rp;
    FSizeOfArray = FSizeOfArray - delta;
    return FSizeOfArray;
}

bool MkParticles::Add(MkParticle point)
{
    if (FSize == FSizeOfArray)
        Grow(FSize - FSizeOfArray + 10);
    if (FSize == FSizeOfArray)
        return false;

    FSize++;
    FParticle[FSize - 1] = point;
    return true;
}

bool MkParticles::Add(int index, MkParticle point)
{
    if (FSize == FSizeOfArray)
        Grow(FSize - FSizeOfArray + 10);
    if (FSize == FSizeOfArray)
        return false;

    for (int i = FSize - 1; i >= index; i--)
        FParticle[i + 1] = FParticle[i];
    FSize++;
    FParticle[index] = point;
    return true;
}

bool MkParticles::Delete(MkParticle point)
{
    int i;
    for (i = 0; i < FSize; i++)
    {
        if (FParticle[i] == point)
            break;
    }
    if (i == FSize)
        return false;
    if (FParticle[i] == point)
    {
        for (int j = i; j < FSize - 1; j++)
            FParticle[j] = FParticle[j + 1];
    }
    FSize--;
    FParticle[FSize] = NullParticle;
    return true;
}

bool MkParticles::Delete(int index)
{
    for (int j = index; j < FSize - 1; j++)
        FParticle[j] = FParticle[j + 1];

    FSize--;
    FParticle[FSize] = NullParticle;
    return true;
}

bool MkParticles::Swap(int i, int j)
{
    MkParticle p;
    if (i >= FSize || j >= FSize)
        return false;
    p = FParticle[i];
    FParticle[i] = FParticle[j];
    FParticle[j] = p;
    return true;
}

bool MkParticles::Clear()
{
    FSizeOfArray = FSize = 0;
    if (FParticle)
    {
        delete[] FParticle;
        FParticle = NULL;
    }
    return true;
}

MkParticle &MkParticles::operator[](int i)
{
    if (FSizeOfArray == 0)
        return NullParticle;
    if (i >= FSize && i < FSizeOfArray)
        FSize = i + 1;

    if (i >= 0 && i < FSize)
        return FParticle[i];
    else
        return NullParticle;
}

MkParticles &MkParticles::operator=(MkParticles &points)
{
    int i;

    Clear();
    FSize = points.FSize;
    FSizeOfArray = points.FSizeOfArray;
    if (FSize == 0)
    {
        FParticle = NULL;
        return *this;
    }
    this->FParticle = new MkParticle[FSizeOfArray];

    for (i = 0; i < FSize; i++)
        FParticle[i] = points.FParticle[i];
    for (i = FSize; i < FSizeOfArray; i++)
        FParticle[i] = NullParticle;

    return *this;
}

bool MkParticles::operator==(MkParticles &points)
{
    int i;

    if (FSize != points.FSize)
        return false;
    for (i = 0; i < FSize; i++)
        if (this->FParticle[i] != points.FParticle[i])
            return false;

    return true;
}

MkParticles &MkParticles::operator*=(MkMatrix4<double> &rm)
{
    for (int i = 0; i < FSize; i++)
        this->FParticle[i] = this->FParticle[i] * rm;
    return *this;
}

MkParticles &operator*(MkParticles &rps, MkMatrix4<double> &rm)
{
    static MkParticles rps_t;
    rps_t = rps;
    for (int i = 0; i < rps.FSize; i++)
        rps_t.FParticle[i] = rps.FParticle[i] * rm;
    return rps_t;
}

MkParticles &MkParticles::Translate(MkParticle rp)
{
    for (int i = 0; i < FSize; i++)
        this->FParticle[i].Translate(rp);
    return *this;
}

MkParticles &MkParticles::Translate(double x, double y, double z)
{
    for (int i = 0; i < FSize; i++)
        this->FParticle[i].Translate(x, y, z);
    return *this;
}

MkParticles &MkParticles::Rotate(double alpha, double beta, double gamma)
{
    for (int i = 0; i < FSize; i++)
        this->FParticle[i].Rotate(alpha, beta, gamma);
    return *this;
}

MkParticles &MkParticles::RotateInX(double ang)
{
    for (int i = 0; i < FSize; i++)
        this->FParticle[i].RotateInX(ang);
    return *this;
}

MkParticles &MkParticles::RotateInY(double ang)
{
    for (int i = 0; i < FSize; i++)
        this->FParticle[i].RotateInY(ang);
    return *this;
}

MkParticles &MkParticles::RotateInZ(double ang)
{
    for (int i = 0; i < FSize; i++)
        this->FParticle[i].RotateInZ(ang);
    return *this;
}

MkParticles &MkParticles::RotateInA(double ang, double l, double m, double n)
{
    for (int i = 0; i < FSize; i++)
        this->FParticle[i].RotateInA(ang, l, m, n);
    return *this;
}

MkParticles &MkParticles::Scale(double sx, double sy, double sz)
{
    for (int i = 0; i < FSize; i++)
        this->FParticle[i].Scale(sx, sy, sz);
    return *this;
}

#ifdef __BCPLUSPLUS__
void MkParticles::Draw(TObject *Sender)
{
    for (int i = 0; i < FSize; i++)
        (*this)[i].Draw(Sender);
}
#endif

#if defined(_MSC_VER) && defined(_WINDOWS_)
void MkParticles::Draw(MkPaint *pb)
{
    for (int i = 0; i < FSize; i++)
        (*this)[i].Draw(pb);
}
#endif

#if defined(__GL_H__)
void MkParticles::Draw()
{
    for (int i = 0; i < FSize; i++)
        (*this)[i].Draw();
}
#endif
//---------------------------------------------------------------------------
