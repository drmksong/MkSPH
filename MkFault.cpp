//---------------------------------------------------------------------------
#include "MkFault.h"

MkFault NullFault(0, 0, 0, 0, 0);

MkFault::MkFault()
{
   Station = 0;
   Dip = 0;
   DipDir = 0;
   Thickness = 0;
   Length = 0;
}

MkFault::MkFault(double station, double dip, double dip_dir, double thickness, double length)
{
   Station = station;
   Dip = dip;
   DipDir = dip_dir;
   Thickness = thickness;
   Length = length;
}

MkFaults::MkFaults()
{
   FSize = 0;
   FFault = (MkFault *)NULL;
}

MkFaults::MkFaults(int size)
{
   if (size <= 0)
   {
      FSize = 0;
      FFault = (MkFault *)NULL;
   }
   FSize = size;
   FFault = new MkFault[FSize];
   if (!FFault)
   {
      FSize = 0;
      return;
   }
}

MkFaults::~MkFaults()
{
   if (FFault)
   {
      delete[](MkFault *) FFault;
      FSize = 0;
   }
}

bool MkFaults::Initialize(int size)
{
   if (size <= 0)
   {
      FSize = 0;
      FFault = (MkFault *)NULL;
      return false;
   }

   if (!FFault)
   {
      delete[](MkFault *) FFault;
      FFault = (MkFault *)NULL;
   }

   FSize = size;
   FFault = new MkFault[FSize];
   if (!FFault)
   {
      FSize = 0;
      return false;
   }
   return true;
}

bool MkFaults::Initialize(int size, MkFault *fault)
{
   if (size <= 0)
   {
      FSize = 0;
      FFault = (MkFault *)NULL;
      return false;
   }

   if (!FFault)
   {
      delete[](MkFault *) FFault;
      FFault = (MkFault *)NULL;
   }

   FSize = size;
   FFault = new MkFault[FSize];
   if (!FFault)
   {
      FSize = 0;
      return false;
   }
   for (int i = 0; i < FSize; i++)
      FFault[i] = fault[i];

   return true;
}

MkFault &MkFaults::operator()(int i)
{
   if (i >= 0 && i < FSize)
      return FFault[i];
   else
      return NullFault;
}

MkFault &MkFaults::operator[](int i)
{
   if (i >= 0 && i < FSize)
      return FFault[i];
   else
      return NullFault;
}

MkFaults &MkFaults::operator=(MkFaults &a)
{
   Initialize(a.GetSize());
   for (int i = 0; i < a.GetSize(); i++)
      FFault[i] = a[i];
}
//---------------------------------------------------------------------------
#pragma package(smart_init)
