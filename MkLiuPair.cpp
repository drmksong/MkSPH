//---------------------------------------------------------------------------
#include "MkLiuPair.hpp"

MkLiuPair NullLiuPair(0);

void MkLiuPair::Clear()
{
   I = 0;
   J = 0;
   Dist = 0;
   dX = 0;
   dY = 0;
   dZ = 0;
   W = 0;
   ;
   dWdX = 0;
   dWdY = 0;
   dWdZ = 0;
   Zero = 0;
}

double &MkLiuPair::DX(int i)
{
   if (i < 0 || i > 2)
   {
      MkDebug("MkLiuPair::dX index outof range\n");
      exit(-1);
      return Zero;
   }
   return (i == 0) ? dX : ((i == 1) ? dY : dZ);
};

double &MkLiuPair::DWDX(int i)
{
   if (i < 0 || i > 2)
   {
      MkDebug("MkLiuPair::dWdX index outof range\n");
      exit(-1);
      return Zero;
   }
   return (i == 0) ? dWdX : ((i == 1) ? dWdY : dWdZ);
}

//-----------------------------------------------------------------------------
MkLiuPairs::MkLiuPairs()
{
   FSize = 0;
   FLiuPair.reset();
}

MkLiuPairs::MkLiuPairs(int size)
{
   if (size <= 0)
   {
      MkDebug("MkLiuPairs::Constructor negative Size error ");
      throw Size(std::string("MkLiuPairs::Constructor Throws Size Exception "), size);
   }
   FSize = size;
   try
   {
      FLiuPair.reset(new MkLiuPair[FSize]);
   }
   catch (std::bad_alloc &a)
   {
      MkDebug("MkLiuPairs::Constructor mem allocation error");
      throw Alloc(a.what());
   }
}

MkLiuPairs::~MkLiuPairs()
{
   FLiuPair.reset();
   FSize = 0;
}

void MkLiuPairs::Clear()
{
   FLiuPair.reset();
   FSize = 0;
}

bool MkLiuPairs::Initialize(int size)
{
   Clear();
   if (size <= 0)
   {
      MkDebug("MkLiuPairs::Initialize negative Size error ");
      throw Size(std::string("MkLiuPairs::Initialize Throws Size Exception "), size);
   }

   FSize = size;

   try
   {
      FLiuPair.reset(new MkLiuPair[FSize]);
   }
   catch (std::bad_alloc &a)
   {
      MkDebug("MkLiuPairs::Constructor mem allocation error");
      throw Alloc(a.what());
   }

   return true;
}

bool MkLiuPairs::Initialize(int size, MkLiuPair *fault)
{
   Clear();
   if (size <= 0)
   {
      MkDebug("MkLiuPairs::Initialize negative Size error ");
      throw Size(std::string("MkLiuPairs::Initialize Throws Size Exception "), size);
   }
   if (fault == NULL)
   {
      MkDebug("MkLiuPairs::Initialize null parameter error ");
      throw Size(std::string("MkLiuPairs::Initialize Throws null parameter Size Exception "), size);
   }

   FSize = size;

   try
   {
      FLiuPair.reset(new MkLiuPair[FSize]);
   }
   catch (std::bad_alloc &a)
   {
      MkDebug("MkLiuPairs::Constructor mem allocation error");
      throw Alloc(a.what());
   }

   for (int i = 0; i < FSize; i++)
      FLiuPair[i] = fault[i];

   return true;
}

MkLiuPair &MkLiuPairs::operator()(int i)
{
   if (i >= 0 && i < FSize)
      return FLiuPair[i];
   else
   {
      MkDebug("MkLiuPair::operator() range exception");
      throw Range(std::string("MkLiuPair::operator() throw range exception"), i);
   }
}

MkLiuPair &MkLiuPairs::operator[](int i)
{
   if (i >= 0 && i < FSize)
      return FLiuPair[i];
   else
   {
      MkDebug("MkLiuPair::operator() range exception");
      throw Range(std::string("MkLiuPair::operator() throw range exception"), i);
   }
}

MkLiuPairs &MkLiuPairs::operator=(MkLiuPairs &a)
{
   try
   {
      Initialize(a.GetSize());
   }
   catch (Alloc &a)
   {
      MkDebug("MkLiuPairs::operator= mem allocation error while preparing copying");
      throw Alloc(a.what());
   }

   for (int i = 0; i < a.GetSize(); i++)
      FLiuPair[i] = a[i];
}
//---------------------------------------------------------------------------
#pragma package(smart_init)
