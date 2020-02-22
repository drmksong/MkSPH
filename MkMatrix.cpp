//---------------------------------------------------------------------------
#include "MkMatrix.hpp"

//---------------------------------------------------------------------------
// index
//  [0][0] [0][1] [0][2] [0][3]
//  [1][0] [1][1] [1][2] [1][3]
//  [2][0] [2][1] [2][2] [2][3]
//  [3][0] [3][1] [3][2] [3][3]
//---------------------------------------------------------------------------
// template <class T> MkMatrix NullMatrix;

double Zero = 0;

template <class T>
MkMatrix4<T>::MkMatrix4()
{
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            FMatrix[i][j] = 0;
}

template <class T>
void MkMatrix4<T>::Identity()
{
    int i;
    for (i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            FMatrix[i][j] = 0;

    for (i = 0; i < 4; i++)
        FMatrix[i][i] = 1;
}

template <class T>
void MkMatrix4<T>::Clear()
{
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            FMatrix[i][j] = 0;
}

template <class T>
void MkMatrix4<T>::Translate(T x, T y, T z)
{
    MkMatrix4 rm;
    rm.Identity();
    rm.FMatrix[0][3] = x;
    rm.FMatrix[1][3] = y;
    rm.FMatrix[2][3] = z;
    *this = *this * rm;
}

template <class T>
void MkMatrix4<T>::Rotate(T alpha, T beta, T gamma)
// x,y,z axis each , don't use this. it is obsolete
// use RotateInX, RotateInY, RotateInZ and RotateInA
// alpha is ang for rotation in x
// beta  is ang for rotation in y
// gamma is ang for rotation in z
// the order is x->y->z
// This routine is changed at 2001.10.3
// So, you must check it...
{
    MkMatrix4 rm1, rm2;

    alpha = alpha * M_PI / 180;
    beta = beta * M_PI / 180;
    gamma = gamma * M_PI / 180;

    if (fabs(alpha - int(alpha)) < 0.001)
        alpha = int(alpha);
    if (fabs(beta - int(beta)) < 0.001)
        beta = int(beta);
    if (fabs(gamma - int(gamma)) < 0.001)
        gamma = int(gamma);
    rm2.Identity();
    rm1.Identity();

    rm1.FMatrix[0][0] = cos(gamma);
    rm1.FMatrix[0][1] = -sin(gamma);
    rm1.FMatrix[1][0] = sin(gamma);
    rm1.FMatrix[1][1] = cos(gamma);

    rm2 *= rm1;

    rm1.Identity();

    rm1.FMatrix[0][0] = cos(beta);
    rm1.FMatrix[0][2] = sin(beta);
    rm1.FMatrix[2][0] = -sin(beta);
    rm1.FMatrix[2][2] = cos(beta);

    rm2 *= rm1;

    rm1.Identity();

    rm1.FMatrix[1][1] = cos(alpha);
    rm1.FMatrix[1][2] = -sin(alpha);
    rm1.FMatrix[2][1] = sin(alpha);
    rm1.FMatrix[2][2] = cos(alpha);

    rm2 *= rm1;

    *this *= rm2;
}

template <class T>
void MkMatrix4<T>::RotateInX(T ang)
{
    MkMatrix4 rm1;

    ang = ang * M_PI / 180;
    if (fabs(ang - int(ang)) < 0.001)
        ang = int(ang);

    rm1.Identity();

    rm1.FMatrix[1][1] = cos(ang);
    rm1.FMatrix[1][2] = -sin(ang);
    rm1.FMatrix[2][1] = sin(ang);
    rm1.FMatrix[2][2] = cos(ang);

    *this *= rm1;
}

template <class T>
void MkMatrix4<T>::RotateInY(T ang)
{
    MkMatrix4 rm1;

    ang = ang * M_PI / 180;
    if (fabs(ang - int(ang)) < 0.001)
        ang = int(ang);

    rm1.Identity();

    rm1.FMatrix[0][0] = cos(ang);
    rm1.FMatrix[0][2] = sin(ang);
    rm1.FMatrix[2][0] = -sin(ang);
    rm1.FMatrix[2][2] = cos(ang);

    *this *= rm1;
}

template <class T>
void MkMatrix4<T>::RotateInZ(T ang)
{
    MkMatrix4 rm1, rm2;

    ang = ang * M_PI / 180;
    if (fabs(ang - int(ang)) < 0.001)
        ang = int(ang);

    rm1.Identity();

    rm1.FMatrix[0][0] = cos(ang);
    rm1.FMatrix[0][1] = -sin(ang);
    rm1.FMatrix[1][0] = sin(ang);
    rm1.FMatrix[1][1] = cos(ang);

    *this *= rm1;
}

template <class T>
void MkMatrix4<T>::RotateInA(T theta, T l, T m, T n)
{
    T len;
    T beta, gamma;
    len = sqrt(l * l + m * m + n * n);
    l /= len;
    m /= len;
    n /= len;
    beta = acos(n);

    if (fabs(l) > 0.0001)
    {
        gamma = acos(sin(beta) / l);
        if (l != sin(beta) * cos(gamma))
            gamma = -gamma;
        beta = beta * 180 / M_PI;
        gamma = gamma * 180 / M_PI;

        RotateInZ(-gamma);
        RotateInY(-beta);
        RotateInZ(theta);
        RotateInY(beta);
        RotateInZ(gamma);
    }
    else if (fabs(m) > 0.0001)
    {
        gamma = asin(sin(beta) / m);
        if (l != sin(beta) * sin(gamma))
            gamma = M_PI - gamma;
        beta = beta * 180 / M_PI;
        gamma = gamma * 180 / M_PI;

        RotateInZ(-gamma);
        RotateInY(-beta);
        RotateInZ(theta);
        RotateInY(beta);
        RotateInZ(gamma);
    }
    else
        RotateInZ(theta);
}

template <class T>
void MkMatrix4<T>::Scale(T sx, T sy, T sz)
{
    MkMatrix4 rm;
    rm.Identity();
    rm.FMatrix[0][0] = sx;
    rm.FMatrix[1][1] = sy;
    rm.FMatrix[2][2] = sz;
    *this *= rm;
}

template <class T>
MkMatrix4<T> &MkMatrix4<T>::operator=(MkMatrix4<T> &rm)
{
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            FMatrix[i][j] = rm.FMatrix[i][j];
    FI = rm.FI;
    FJ = rm.FJ;
    return *this;
}

template <class T>
MkMatrix4<T> &MkMatrix4<T>::operator*=(T f)
{
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            FMatrix[i][j] *= f;

    return *this;
}

template <class T>
MkMatrix4<T> &MkMatrix4<T>::operator*=(MkMatrix4<T> &rm)
{
    static MkMatrix4<T> rm_t;
    rm_t.Clear();

    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            for (int k = 0; k < 4; k++)
            {
                rm_t.FMatrix[i][j] += FMatrix[i][k] * rm.FMatrix[k][j];
            }

    return *this = rm_t;
}

template <class T>
MkMatrix4<T> &MkMatrix4<T>::operator+=(T f)
{
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            FMatrix[i][j] += f;

    return *this;
}

template <class T>
MkMatrix4<T> &MkMatrix4<T>::operator+=(MkMatrix4<T> &rm)
{
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            FMatrix[i][j] += rm.FMatrix[i][j];

    return *this;
}

template <class T>
MkMatrix4<T> &MkMatrix4<T>::operator-=(T f)
{
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            FMatrix[i][j] -= f;

    return *this;
}

template <class T>
MkMatrix4<T> &MkMatrix4<T>::operator-=(MkMatrix4<T> &rm)
{
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            FMatrix[i][j] -= rm.FMatrix[i][j];

    return *this;
}

#ifdef __BCPLUSPLUS__
template <class T>
void MkMatrix4<T>::Out(TMemo *memo)
{
    char str[256], s[256];
    memo->Lines->Add("Output of Matrix ");
    sprintf(str, "A-th row [%10.3f %10.3f %10.3f %10.3f ]", 0.0, 1.0, 2.0, 3.0);
    memo->Lines->Add(str);
    for (int i = 0; i < 4; i++)
    {
        sprintf(str, "%d-th row [", i);
        for (int j = 0; j < 4; j++)
        {
            sprintf(s, "%10.3f ", FMatrix[i][j]);
            strcat(str, s);
        }
        sprintf(s, "]");
        strcat(str, s);
        memo->Lines->Add(str);
    }
}
#endif

template <class T>
void MkMatrix4<T>::Out(char *fname)
{
    FILE *fp;
    char str[256], s[256];

    fp = fopen(fname, "a");
    fprintf(fp, "Output of Matrix ");
    sprintf(str, "A-th row [%10.3f %10.3f %10.3f %10.3f ]", 0.0, 1.0, 2.0, 3.0);
    fprintf(fp, str);
    for (int i = 0; i < 4; i++)
    {
        sprintf(str, "%d-th row [", i);
        for (int j = 0; j < 4; j++)
        {
            sprintf(s, "%10.3f ", FMatrix[i][j]);
            strcat(str, s);
        }
        sprintf(s, "]");
        strcat(str, s);
        fprintf(fp, str);
    }
    fclose(fp);
}

//---------------------------------------------------------------------------
template <class T>
MkVector<T>::MkVector()
{
    FSize = 0;
    FVector.Clear();
    FVectType = vtNone;
}

template <class T>
MkVector<T>::MkVector(int sz)
{
    char str[256];
    FSize = sz;
    try
    {
        FVector.Initialize(sz);
    }
    catch (Size s)
    {
        sprintf(str, "Vector::Constructor Size error of MkArray<T> size(%d)\n", s.X);
        MkDebug(str);
        throw Size(0);
    }
    catch (Alloc a)
    {
        sprintf(str, "Vector::Allocation Error MkArray<T>");
        MkDebug(str);
        throw Alloc();
    }

    FVectType = vtNone;
}

template <class T>
MkVector<T>::MkVector(int sz, VectType vt)
{
    char str[256];
    FSize = sz;
    try
    {
        FVector.Initialize(sz);
    }
    catch (Size s)
    {
        sprintf(str, "Vector::Constructor Size error of MkArray<T> size(%d)\n", s.X);
        MkDebug(str);
        throw Size(0);
    }
    catch (Alloc a)
    {
        sprintf(str, "Vector::Allocation Error MkArray<T>");
        MkDebug(str);
        throw Alloc();
    }

    FVectType = vt;
}

template <class T>
MkVector<T>::MkVector(MkArray<T> &b)
{
    FSize = b.getSzX();
    FVector.CopyFrom(b);
    FVectType = vtNone;
}

template <class T>
MkVector<T>::MkVector(T l, T m, T n)
{
    FVector.Initialize(3);
    FSize = 3;
    FVector(0) = l;
    FVector(1) = m;
    FVector(2) = n;
    FVectType = vtNone;
}

template <class T>
MkVector<T>::~MkVector()
{
}

template <class T>
void MkVector<T>::SetVector(int sz)
{
    FSize = sz;
    FVector.Initialize(sz);
    //   FMatrix.Clear();
}

template <class T>
void MkVector<T>::SetVector(T l, T m, T n)
{
    FVector.Initialize(3);
    FSize = 3;
    FVector(0) = l;
    FVector(1) = m;
    FVector(2) = n;
    FVectType = vtNone;
}

template <class T>
void MkVector<T>::Clear()
{
    FSize = 0;
    FVector.Clear();
    //    FMatrix.Clear();
    FVectType = vtNone;
}

template <class T>
void MkVector<T>::Unify()
{
    T sum = 0;
    for (int i = 0; i < FSize; i++)
        sum += FVector(i) * FVector(i);
    sum = sqrt(sum);
    if (sum > 0.0001)
        for (int i = 0; i < FSize; i++)
            FVector(i) = FVector(i) / sum;
}

template <class T>
T MkVector<T>::Dot(MkVector &vect)
{
    T dot = 0;
    //    if (FSize <= 0 || FSize != vect.FSize || FVectType != vtRow || vect.FVectType != vtCol) return 0;
    if (FSize <= 0 || FSize != vect.FSize)
        return 0;
    for (int i = 0; i < FSize; i++)
        dot += FVector(i) * vect(i);
    return dot;
}

template <class T>
void MkVector<T>::Cross(MkVector &vect, MkVector &target) // not yet finished, because I do not find general cross product
{
    char str[256];
    try
    {
        if (FSize != target.GetSize())
            target.Initialize(FSize);
    }
    catch (Size s)
    {
        sprintf(str, "Vector::Cross Size error of MkArray<T> size(%d)\n", s.X);
        MkDebug(str);
        throw Size(0);
    }
    catch (Alloc a)
    {
        sprintf(str, "Vector::Allocation Error MkArray<T>");
        MkDebug(str);
        throw Alloc();
    }

    //    l3 = m1*n2 - m2*n1;
    //    m3 = n1*l2 - n2*l1;
    //    n3 = l1*m2 - l2*m1;  �̰ź��� �����ϱ� �ٶ�...-.-;

    for (int i = 0; i < FSize; i++)
    {
        target.FVector(i) = FVector(i + 1 >= FSize ? i + 1 - FSize : i + 1) * vect(i + 2 >= FSize ? i + 2 - FSize : i + 2) - vect(i + 1 >= FSize ? i + 1 - FSize : i + 1) * FVector(i + 2 >= FSize ? i + 2 - FSize : i + 2); //�򰥸�...^^
    }
}

template <class T>
T &MkVector<T>::operator()(int i)
{
    return FVector(i);
}

template <class T>
T &MkVector<T>::operator[](int i)
{
    return FVector(i);
}

template <class T>
MkVector<T> &MkVector<T>::operator=(MkVector<T> &vect)
{
    FVector.CopyFrom(vect.FVector);
    FSize = vect.FSize;
    FVectType = vect.FVectType;
    return *this;
}

template <class T>
bool MkVector<T>::operator==(MkVector<T> &v)
{
    return FSize == v.FSize && FVectType == v.FVectType && FVector == v.FVector;
}

template <class T>
bool MkVector<T>::operator!=(MkVector &v)
{
    return !operator==(v);
}

template <class T>
MkVector<T> &MkVector<T>::operator+=(MkVector &vect)
{
    static MkVector<T> NullVector(0);
    if (FSize != vect.FSize)
        return NullVector;
    for (int i = 0; i < FSize; i++)
        FVector(i) += vect.FVector(i);
    return *this;
}

template <class T>
MkVector<T> &MkVector<T>::operator-=(MkVector &vect)
{
    static MkVector<T> NullVector(0);
    if (FSize != vect.FSize)
        return NullVector;
    for (int i = 0; i < FSize; i++)
        FVector(i) -= vect.FVector(i);
    return *this;
}

template <class T>
MkVector<T> &MkVector<T>::operator*=(T a)
{
    for (int i = 0; i < FSize; i++)
        FVector(i) *= a;
    return *this;
}

template <class T>
MkVector<T> &MkVector<T>::operator/=(T a)
{
    if (fabs(a) > EPS)
        for (int i = 0; i < FSize; i++)
            FVector(i) /= a;

    return *this;
}

template <class T>
MkVector<T> &MkVector<T>::operator+=(T a)
{
    for (int i = 0; i < FSize; i++)
        FVector(i) += a;
    return *this;
}

template <class T>
MkVector<T> &MkVector<T>::operator-=(T a)
{
    for (int i = 0; i < FSize; i++)
        FVector(i) -= a;
    return *this;
}

#ifdef __BCPLUSPLUS__
template <class T>
void MkVector<T>::Out(TMemo *memo)
{
    char str[256], s[256];
    memo->Lines->Add("Output of Vector ");
    sprintf(str, "vector [");
    for (int i = 0; i < FSize; i++)
    {
        sprintf(s, "%10.6f ", FVector(i));
        strcat(str, s);
        if (strlen(str) > 200)
        {
            strcat(str, "...");
            break;
        }
    }
    sprintf(s, "]\n");
    strcat(str, s);
    memo->Lines->Add(str);
}
#endif
template <class T>
void MkVector<T>::Out(char *fname)
{
    FILE *fp;
    char str[256], s[256];

    fp = fopen(fname, "a");
    fprintf(fp, "Output of Vector ");
    sprintf(str, "vector [");
    for (int i = 0; i < FSize; i++)
    {
        sprintf(s, "%10.6f ", FVector(i));
        strcat(str, s);
        if (strlen(str) > 200)
        {
            strcat(str, "...");
            break;
        }
    }
    sprintf(s, "]\n");
    strcat(str, s);
    fprintf(fp, str);
    fclose(fp);
}

//---------------------------------------------------------------------------
template <class T>
MkMatrix<T>::MkMatrix()
{
    FMatrix.Clear();
    FIndex.Clear();
    FD = 0;
    FI = FJ = 0;
    FMatType = mtNormal;
}

template <class T>
MkMatrix<T>::MkMatrix(int sz_x, int sz_y)
{
    char str[256];
    try
    {
        FMatrix.Initialize(sz_x, sz_y);
    }
    catch (Size s)
    {
        sprintf(str, "Matrix::Constructor Size error of MkArray<T> size(%d,%d)\n", s.X, s.Y);
        MkDebug(str);
        throw Size(0);
    }
    catch (Alloc a)
    {
        sprintf(str, "Vector::Allocation Error MkArray<T>");
        MkDebug(str);
        throw Alloc();
    }

    try
    {
        if (sz_x == sz_y)
            FIndex.Initialize(sz_x);
    }
    catch (Size s)
    {
        sprintf(str, "Matrix::Constructor Size error of MkArray<T> size(%d)\n", s.X);
        MkDebug(str);
        throw Size(0);
    }
    catch (Alloc a)
    {
        sprintf(str, "Vector::Allocation Error MkArray<T>");
        MkDebug(str);
        throw Alloc();
    }

    FD = 0;
    FI = sz_x;
    FJ = sz_y;
    FMatType = mtNormal;
}

template <class T>
MkMatrix<T>::MkMatrix(MkMatrix<T> &tm)
{
    *this = tm;
}

template <class T>
MkMatrix<T>::MkMatrix(MkArray<T> &fm)
{
    if (fm.getSzZ() != 1)
    {
        FMatrix.Clear();
        FIndex.Clear();
        FD = 0;
        FI = FJ = 0;
        FMatType = mtNormal;
    }
    else
    {
        FI = fm.getSzX();
        FJ = fm.getSzY();
        FMatrix.CopyFrom(fm);
        FIndex.Initialize(FI);
        FD = 0;
        FMatType = mtNormal;
    }
}

template <class T>
MkMatrix<T>::~MkMatrix()
{
}

template <class T>
void MkMatrix<T>::SetMatrix(int sz_x, int sz_y)
{
    FMatrix.Initialize(sz_x, sz_y);
    if (sz_x == sz_y)
        FIndex.Initialize(sz_x);
    FI = sz_x;
    FJ = sz_y;
    FMatType = mtNormal;
}

template <class T>
void MkMatrix<T>::Identify()
{
    if (FI != FJ)
        return;
    for (int i = 0; i < FI; i++)
        for (int j = 0; j < FJ; j++)
        {
            if (i == j)
                FMatrix(i, j) = 1;
            else
                FMatrix(i, j) = 0;
        }
    FMatType = mtNormal;
}

template <class T>
void MkMatrix<T>::Clear()
{
    FMatrix.Clear();
    FIndex.Clear();
    FI = FJ = 0;
    FMatType = mtNormal;
}

template <class T>
bool MkMatrix<T>::Transpose()
{
    MkArray<T> A(FJ, FI);
    for (int i = 0; i < FI; i++)
        for (int j = 0; j < FJ; j++)
            A(j, i) = FMatrix(i, j);

    FMatrix.CopyFrom(A);

    int tmp;
    tmp = FI;
    FI = FJ;
    FJ = tmp;

    FMatType = FMatType == mtNormal ? mtTransposed : mtNormal;
    return true;
}

template <class T>
bool MkMatrix<T>::Invert()
{
    int i, j, n;
    if (FI != FJ)
        return false;
    n = FI;

    MkArray<T> y(n, n);
    MkVector<T> TempVect(n);
    LUDecompose();

    for (j = 0; j < n; j++)
    {
        for (i = 0; i < n; i++)
            TempVect(i) = 0;
        TempVect(j) = 1.0;
        LUBackSubstitute(TempVect);
        for (i = 0; i < n; i++)
            y(i, j) = TempVect(i);
    }

    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            FMatrix(i, j) = y(i, j);

    FMatType = mtInverted;
    return true;
}

//  T **a, int n, int *indx, T FD

template <class T>
bool MkMatrix<T>::LUDecompose()
{
    int i, imax, j, k;
    T big, dum, sum, temp;
    MkArray<T> vv;

    if (FMatType == mtLUDecomposed)
        return true;
    if (FI != FJ)
        return false;
    int n = FI;

    vv.Initialize(n);
    FD = 1;

    for (i = 0; i < n; i++)
    {
        big = 0.0;
        for (j = 0; j < n; j++)
            if ((temp = fabs(FMatrix(i, j))) > big)
                big = temp;
        if (fabs(big) < EPS)
            MkDebug("Singular matrix in routine MkMatrix<T>::LUDecompose");
        vv(i) = 1.0 / big;
    }

    for (j = 0; j < n; j++)
    {
        for (i = 0; i < j; i++)
        {
            sum = FMatrix(i, j);
            for (k = 0; k < i; k++)
                sum -= FMatrix(i, k) * FMatrix(k, j);
            FMatrix(i, j) = sum;
        }
        big = 0.0;
        for (i = j; i < n; i++)
        {
            sum = FMatrix(i, j);
            for (k = 0; k < j; k++)
                sum -= FMatrix(i, k) * FMatrix(k, j);
            FMatrix(i, j) = sum;
            if ((dum = vv(i) * fabs(sum)) >= big)
            {
                big = dum;
                imax = i;
            }
        }
        if (j != imax)
        {
            for (k = 0; k < n; k++)
            {
                dum = FMatrix(imax, k);
                FMatrix(imax, k) = FMatrix(j, k);
                FMatrix(j, k) = dum;
            }
            FD = -(FD);
            vv(imax) = vv(j);
        }
        FIndex(j) = imax;
        if (FMatrix(j, j) == 0.0)
            FMatrix(j, j) = TINY;
        if (j != n)
        {
            dum = 1.0 / (FMatrix(j, j));
            for (i = j + 1; i < n; i++)
                FMatrix(i, j) *= dum;
        }
    }
    FMatType = mtLUDecomposed;
    return true;
}

//  b is input and output. Solving equations AX = B.
//  b as input B
//  b as output X

template <class T>
bool MkMatrix<T>::LUBackSubstitute(MkVector<T> &b)
{
    int i, ii = -1, ip, j;
    T sum;
    int n;

    if (FI != FJ)
        return false;
    n = FI;
    if (FMatType != mtLUDecomposed)
        return false;

    for (i = 0; i < n; i++)
    {
        ip = FIndex(i);
        sum = b(ip);
        b(ip) = b(i);
        if (ii >= 0)
            for (j = ii; j <= i - 1; j++)
                sum -= FMatrix(i, j) * b(j);
        else if (sum)
            ii = i;
        b(i) = sum;
    }

    for (i = n - 1; i >= 0; i--)
    {
        sum = b(i);
        for (j = i + 1; j < n; j++)
            sum -= FMatrix(i, j) * b(j);
        b(i) = sum / FMatrix(i, i);
    }
    return true;
}

template <class T>
bool MkMatrix<T>::GaussSeidel(MkVector<T> &X0, MkVector<T> &B) // X0 is initial and return X
{
    int m = 0;
    int NN = FI;
    int MaxIter = 100, iter = 0;
    //    T max_x,x,x0;
    //    bool flag;
    MkArray<T> U(NN);
    MkArray<T> X(NN);
    MkArray<T> A;

    if (FMatType == mtLUDecomposed)
        return false;
    if (FI != FJ)
        return false;

    A.CopyFrom(FMatrix);

    for (int i = 0; i < NN; i++)
    {
        U(i) = X0(i);
        B(i) /= A(i, i);
        for (int j = 0; j < NN; j++)
            A(i, j) = (i == j) ? A(i, j) : A(i, j) / A(i, i);
        A(i, i) = 0;
    }

    while (m < NN && iter < MaxIter)
    {
        for (int i = 0; i < NN; i++)
        {
            X(i) = 0;
            for (int j = 0; j < NN; j++)
            {
                X(i) = X(i) - A(i, j) * U(j);
            }
            X(i) += B(i);
            U(i) = X(i);
        }
        m = 0;
        //in case x and x0 is too small to compare, then normalize the value to its bigger one.
        /*
	x = X(m);
	x0 = X0(m);
	max_x = (fabs(x)>fabs(x0))?fabs(x):fabs(x0);

	if(EPS*EPS<max_x&&max_x<EPS) flag = flag && (fabs(x-x0)/max_x < EPS*10);
	else if(max_x>EPS) flag = flag && (fabs(x-x0) < EPS);
	else flag = flag && true;
*/
        while ((m < NN) && fabs(X(m) - X0(m)) / max((T)fabs(X(m)), (T)fabs(X0(m))) < 1e-10)
            m++;
        if (m < NN)
            for (int i = 0; i < NN; i++)
            {
                X0(i) = X(i);
            }
        iter++;
    }
    return true;
}

template <class T>
bool MkMatrix<T>::GaussSeidelDouble(MkVector<T> &X0, MkVector<T> &B) // X0 is initial and return X
{
    int m = 0;
    int NN = FI;
    int MaxIter = 100000, iter = 0;
    //    T max_x,x,x0;
    //    bool flag;
    MkArray<T> U(NN);
    MkArray<T> X(NN);
    MkArray<T> A;

    if (FMatType == mtLUDecomposed)
        return false;
    if (FI != FJ)
        return false;

    A.CopyFrom(FMatrix);

    for (int i = 0; i < NN; i++)
    {
        U(i) = X0(i);
        B(i) /= A(i, i);
        for (int j = 0; j < NN; j++)
            A(i, j) = (i == j) ? A(i, j) : A(i, j) / A(i, i);
        A(i, i) = 0;
    }

    while (m < NN && iter < MaxIter)
    {
        for (int i = 0; i < NN; i++)
        {
            X(i) = 0;
            for (int j = 0; j < NN; j++)
            {
                X(i) = X(i) - A(i, j) * U(j);
            }
            X(i) += B(i);
            U(i) = X(i);
        }
        m = 0;
        //in case x and x0 is too small to compare, then normalize the value to its bigger one.
        /*
	x = X(m);
	x0 = X0(m);
	max_x = (fabs(x)>fabs(x0))?fabs(x):fabs(x0);

	if(EPS*EPS<max_x&&max_x<EPS) flag = flag && (fabs(x-x0)/max_x < EPS*10);
	else if(max_x>EPS) flag = flag && (fabs(x-x0) < EPS);
	else flag = flag && true;
*/
        while ((m < NN) && fabs(X(m) - X0(m)) / max((T)fabs(X(m)), (T)fabs(X0(m))) < 1e-10)
            m++;
        if (m < NN)
            for (int i = 0; i < NN; i++)
            {
                X0(i) = X(i);
            }
        iter++;
    }
    return true;
}

template <class T>
bool MkMatrix<T>::isSingular()
{
    int i, j;
    T big, temp;

    if (FI != FJ)
        return true;

    FD = 1;

    for (i = 0; i < FI; i++)
    {
        big = 0.0;
        for (j = 0; j < FI; j++)
            if ((temp = fabs(FMatrix(i, j))) > big)
                big = temp;
        if (fabs(big) <= EPS)
            return true;
    }
    return false;
}

template <class T>
bool MkMatrix<T>::Solve(MkVector<T> &B)
{
    if (isSingular())
        return false;
    MkVector<T> X;
    MkMatrix A(FMatrix);
    X = B;
    A.LUDecompose();
    A.LUBackSubstitute(X);
    MkMatrix A1(FMatrix);
    A1.GaussSeidel(X, B);
    B = X;
    return true;
}

template <class T>
bool MkMatrix<T>::Solve(MkVector<T> &B, SolveType solve_type)
{
    if (isSingular())
        return false;

    if (solve_type == stLUD)
    {
        MkVector<T> X;
        MkMatrix A(FMatrix);
        X = B;
        A.LUDecompose();
        A.LUBackSubstitute(X);
        B = X;
        return true;
    }
    else if (solve_type == stGauss)
    {
        MkVector<T> X;
        MkMatrix A(FMatrix);
        X = B;
        A.GaussSeidel(X, B);
        B = X;
        return true;
    }
    else if (solve_type == stHybrid)
        return Solve(B);
    else
        return false;
}

template <class T>
MkMatrix<T> &MkMatrix<T>::GetTranspose()
{
    static MkMatrix<T> NullMatrix;
    if (Transpose())
        return *this;
    else
        return NullMatrix;
}

template <class T>
MkMatrix<T> &MkMatrix<T>::GetInvert()
{
    static MkMatrix<T> NullMatrix;
    if (Invert())
        return *this;
    else
        return NullMatrix;
}

template <class T>
MkMatrix<T> &MkMatrix<T>::operator!() // matrix inversion
{
    return GetInvert();
}

template <class T>
MkMatrix<T> &MkMatrix<T>::operator~() // matrix tranpose
{
    return GetTranspose();
}

template <class T>
MkMatrix<T> &MkMatrix<T>::operator*=(MkMatrix<T> &m)
{
    static MkMatrix<T> NullMatrix;
    if (FJ != m.FI)
        return NullMatrix;
    static MkMatrix m_t(FI, m.FJ);
    T sum;
    for (int i = 0; i < FI; i++)
    {
        for (int j = 0; j < m.FJ; j++)
        {
            sum = 0;
            for (int k = 0; k < FJ; k++)
                sum += FMatrix(i, k) * m(k, j);
            m_t(i, j) = sum;
        }
    }
    return *this = m_t;
}

template <class T>
MkVector<T> &MkMatrix<T>::operator*=(MkVector<T> &v)
{
    static MkVector<T> NullVector;
    static MkVector<T> FVector;
    if (FJ != v.GetFI())
        return NullVector;
    FVector.Initialize(FI);
    T sum;
    for (int i = 0; i < FI; i++)
    {
        sum = 0;
        for (int j = 0; j < FJ; j++)
            sum += FMatrix(i, j) * v(j);
        FVector(i) = sum;
    }
    return FVector;
}

template <class T>
MkMatrix<T> &MkMatrix<T>::operator=(MkMatrix<T> &m)
{
    FMatrix.CopyFrom(m.FMatrix);
    FIndex.CopyFrom(m.FIndex);
    FD = m.FD;
    FI = m.FI;
    FJ = m.FJ;
    FMatType = m.FMatType;
    return *this;
}

template <class T>
T &MkMatrix<T>::operator()(int i, int j)
{
    return FMatrix(i, j);
}

#ifdef __BCPLUSPLUS__
template <class T>
void MkMatrix<T>::Out(TMemo *memo)
{
    char str[256], s[256];
    memo->Lines->Add("Output of Matrix");
    sprintf(str, " n [%10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f]",
            0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);
    memo->Lines->Add(str);
    for (int i = 0; i < GetFI(); i++)
    {
        sprintf(str, "%2d [", i);
        for (int j = 0; j < GetFJ(); j++)
        {
            sprintf(s, "%10.3f ", FMatrix(i, j));
            strcat(str, s);
            if (strlen(str) > 200)
            {
                strcat(str, "...");
                break;
            }
        }
        sprintf(s, "]");
        strcat(str, s);
        memo->Lines->Add(str);
    }
}
#endif

template <class T>
void MkMatrix<T>::Out(char *fname)
{
    FILE *fp;
    char str[256], s[256];

    fp = fopen(fname, "a");
    fprintf(fp, "Output of Matrix\n");
    sprintf(str, " n [%10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f]\n",
            0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);
    fprintf(fp, str);
    for (int i = 0; i < GetFI(); i++)
    {
        sprintf(str, "%2d [", i);
        for (int j = 0; j < GetFJ(); j++)
        {
            sprintf(s, "%10.3f ", FMatrix(i, j));
            strcat(str, s);
            if (strlen(str) > 200)
            {
                strcat(str, "...");
                break;
            }
        }
        sprintf(s, "]\n");
        strcat(str, s);
        fprintf(fp, str);
    }
    fclose(fp);
}
//---------------------------------------------------------------------------

template class MkMatrix4<double>;
template class MkMatrix4<float>;
template class MkVector<float>;
template class MkVector<double>;
template class MkMatrix<double>;
template class MkMatrix<float>;
