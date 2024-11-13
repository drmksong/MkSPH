//--------------------------------------------------------------------------------------------
#ifndef LiuPair_H
#define LiuPair_H
//--------------------------------------------------------------------------------------------
#include <stdlib.h>
#include <stddef.h>
#include <boost/shared_array.hpp>
#include "MkMisc.hpp"

class MkLiuPair
{
public:
  int I, J;
  double Dist;
  double dX;
  double dY;
  double dZ;
  double W;
  double dWdX;
  double dWdY;
  double dWdZ;
  double Zero;

public:
  MkLiuPair() { Clear(); }
  MkLiuPair(int i) { Clear(); }
  ~MkLiuPair(){};
  void Clear();
  double &DX(int i);
  double &DWDX(int i);
};

class MkLiuPairs
{
private:
  int FSize;
  boost::shared_array<MkLiuPair> FLiuPair;

public:
  MkLiuPairs();
  MkLiuPairs(int);
  ~MkLiuPairs();
  bool Initialize(int size);
  bool Initialize(int size, MkLiuPair *pair);
  void Clear();

  MkLiuPair &operator()(int);
  MkLiuPair &operator[](int);
  MkLiuPairs &operator=(MkLiuPairs &a);

  int GetSize() { return FSize; };

  class Alloc
  {
  public:
    std::string What;
    Alloc(std::string what) : What(what) {}
    std::string what() { return What; }
  };
  class Size
  {
  public:
    std::string What;
    int N;
    Size(std::string what, int n) : What(what), N(n) {}
    std::string what() { return What; }
  };
  class Range
  {
  public:
    std::string What;
    int N;
    Range(std::string what, int n) : What(what), N(n) {}
    std::string what() { return What; }
  };
};
extern MkLiuPair NullLiuPair;

#endif
