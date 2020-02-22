//---------------------------------------------------------------------------
#ifndef MkFaultH
#define MkFaultH

#include "MkPoint.hpp"

class MkFault
{
private:
   double Dip, DipDir, Thickness, Length, Station;
   MkPoint Center;
   MkPoint Location; // �ܺο��� �׻� �������־�� �Ѵ�.
public:
   MkFault();
   MkFault(double station, double dip, double dip_dir, double thickness, double length);
   void SetDip(double dip) { Dip = dip; }
   void SetDipDir(double dip_dir) { DipDir = dip_dir; }
   void SetThickness(double thickness) { Thickness = thickness; }
   void SetLength(double length) { Length = length; }
   void SetStation(double station) { Station = station; }
   void SetStation(int i, double j) { Station = i * 1000 + j; }
   void SetCenter(MkPoint cen) { Center = cen; }
   double GetDip() { return Dip; }
   double GetDipDir() { return DipDir; }
   double GetThickness() { return Thickness; }
   double GetLength() { return Length; }
   double GetStation() { return Station; }
   MkPoint GetCenter() { return Center; }
   MkFault &operator=(MkFault &fault)
   {
      Dip = fault.Dip;
      DipDir = fault.DipDir;
      Thickness = fault.Thickness;
      Length = fault.Length;
      Station = fault.Station;
      Center = fault.Center;
      return *this;
   };
};

class MkFaults
{
private:
   int FSize;
   MkFault *FFault;

public:
   MkFaults();
   MkFaults(int);
   ~MkFaults();
   bool Initialize(int size);
   bool Initialize(int size, MkFault *fault);
   void Clear();

   MkFault &operator()(int);
   MkFault &operator[](int);
   MkFaults &operator=(MkFaults &a);

   int GetSize() { return FSize; };
};
extern MkFault NullFault;
//---------------------------------------------------------------------------
#endif
