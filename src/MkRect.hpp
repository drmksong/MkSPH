//---------------------------------------------------------------------------
#ifndef MkRectHPP
#define MkRectHPP
#include <string>
#include "MkContainer.hpp"
#include "MkShape.hpp"
#include "MkLine.hpp"
//---------------------------------------------------------------------------
class MkRect : public MkShape
{
protected:
  MkPoint Center;
  MkPoint Origin;
  float Width, Height;
  std::string className;

public:
  MkRect();
  MkRect(int);
  MkRect(MkPoint &p);
  MkRect(float x, float y);

public:
  void SetOrigin(MkPoint &pnt) { Origin = pnt; }
  void SetHeight(float h) { Height = h; }
  void SetWidth(float w) { Width = w; }

public:
  MkPoint &GetOrigin() { return Origin; }
  float GetHeight() { return Height; }
  float GetWidth() { return Width; }
  void GetCross(MkLine &, MkPoints &);
  float GetTop() { return max(Origin.Y + Height, Origin.Y); }
  float GetBot() { return min(Origin.Y + Height, Origin.Y); }
  float GetLeft() { return min(Origin.X + Width, Origin.X); }
  float GetRight() { return max(Origin.X + Width, Origin.X); }
  MkPoints &GetIntPoints(MkLine &line);
  MkPoint GetCenter()
  {
    Center = Origin;
    Center.X += Width / 2;
    Center.Y += Height / 2;
    return Center;
  }

  bool isRect() { return true; }
  bool IsIn(MkPoint &pnt);
  bool IsInSurface(MkPoint &pnt, float thick);
  bool IsInSpace(MkPoint &pnt);
  MkRect &operator=(MkRect &rect)
  {
    Origin = rect.Origin;
    Width = rect.Width;
    Height = rect.Height;
    return *this;
  }
  bool operator==(MkRect &rect) { return Origin == rect.Origin && Width == rect.Width && Height == rect.Height; }
  bool operator!=(MkRect &rect) { return Origin != rect.Origin || Width != rect.Width || Height != rect.Height; }
  MkPoint &operator[](int i);
  MkLine &operator()(int i);

#ifdef __BCPLUSPLUS__
  AnsiString ClassName()
  {
    return AnsiString("MkRect");
  };
#else
  std::string ClassName()
  {
    return className;
  }
#endif

#ifdef __BCPLUSPLUS__
  void Draw(TObject *);
#endif

#if defined(_MSC_VER) && defined(_WINDOWS_)
  void Draw(MkPaint *);
#endif
};

//---------------------------------------------------------------------------
// class MkRects
// {
// private:
//   int FSize;
//   MkRect *FRect;

// public:
//   MkRects();
//   MkRects(int);
//   ~MkRects();
//   bool Initialize(int size);
//   bool Initialize(int size, MkRect *fault);
//   void Clear();

//   MkRect &operator()(int);
//   MkRect &operator[](int);
//   MkRects &operator=(MkRects &a);

//   int GetSize() { return FSize; };

// #ifdef __BCPLUSPLUS__
//   void Draw(TObject *);
// #endif

// #if defined(_MSC_VER) && defined(_WINDOWS_)
//   void Draw(MkPaint *);
// #endif
// };
// extern MkRect NullRect;
//---------------------------------------------------------------------------

typedef MkContainer<MkRect> MkRects;
template class MkContainer<MkRect>;

extern MkRect NullRect;
extern MkRects NullRects;

#endif
