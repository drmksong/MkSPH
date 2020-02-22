//---------------------------------------------------------------------------
#ifndef MkContainerHPP
#define MkContainerHPP

#include <boost/shared_array.hpp>
#include "MkMisc.hpp"

// #include "MkCube.hpp"
// #include "MkLine.hpp"
// #include "MkRect.hpp"
// #include "MkSphere.hpp"
// #include "MkTriangle.hpp"
// #include "MkCylinder.hpp"

template <class T>
class MkContainer
{
protected:
  boost::shared_array<T> F;
  int FSize;

public:
  MkContainer(int size)
  {

    //  try {
    if (size <= 0)
    {
      MkDebug("::MkContainer<T> - MkContainer(int size)");
      throw Size(std::string("MkContainer<T>::size below zero strange error"), size);
    }

    FSize = size;

    try
    {
      F.reset(new T[FSize]);
    }
    catch (std::bad_alloc &a)
    {
      MkDebug("MkContainer Contructor memory allocation std::bad_alloc thrown\n");
      throw Alloc(std::string("MkContainer<T>::bad_alloc error"));
    }
  }

  MkContainer()
  {
    FSize = 0;
    F.reset();
  }
  ~MkContainer()
  {
    F.reset();
  }
  void Initialize(int size)
  {
    Clear();
    FSize = size;
    if (FSize == 0)
    {
      F.reset();
      return;
    }

    try
    {
      F.reset(new T[FSize]);
    }
    catch (std::bad_alloc &a)
    {
      MkDebug("MkContainer Contructor memory allocation std::bad_alloc thrown\n");
      throw Alloc(std::string("MkContainer<T>::bad_alloc error"));
    }
  }
  void Clear()
  {
    FSize = 0;
    F.reset();
  }
  int GetSize() { return FSize; }
  T &operator[](int i)
  {
    if (i >= 0 && i < FSize)
    {
      return F[i];
    }
    else
    {
      MkDebug("MkContainer<T>::operator [%d] out of bound\n", i);
      throw Range(std::string("MkContainer operator [] Range thrown"), i);
    }
  }
  MkContainer<T> &operator=(MkContainer<T> &spheres)
  {
    int i;
    FSize = spheres.FSize;

    try
    {
      F.reset(new T[FSize]);
    }
    catch (std::bad_alloc &a)
    {
      MkDebug("MkContainer Contructor memory allocation std::bad_alloc thrown\n");
      throw Alloc(std::string("MkContainer<T>::operator=(MkContainer<T>)"));
    }

    for (i = 0; i < FSize; i++)
      F[i] = spheres.F[i];

    return *this;
  }

#ifdef __BCPLUSPLUS__
  void Draw(TObject *)
  {
    for (int i = 0; i < FSize; i++)
      F[i].Draw(Sender);
  }
#endif

#if defined(_MSC_VER) && defined(_WINDOWS_)
  void Draw(MkPaint *)
  {
    for (int i = 0; i < FSize; i++)
      F[i].Draw(pb);
  }
#endif

  class Alloc
  {
  public:
    std::string What;
    Alloc(std::string what) : What(what) {}
  };
  class Size
  {
  public:
    std::string What;
    int N;
    Size(std::string what, int n) : What(what), N(n) {}
  };
  class Range
  {
  public:
    std::string What;
    int N;
    Range(std::string what, int n) : What(what), N(n) {}
  };
};

#endif