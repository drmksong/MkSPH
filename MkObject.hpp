//---------------------------------------------------------------------------
#ifndef MkObjectH
#define MkObjectH
//---------------------------------------------------------------------------
class MkObject { // use it if needed
public:
  MkObject(){};
  MkObject(int){};
  ~MkObject(){};
  virtual void Clear(){}
  virtual bool operator==(MkObject &){return false;}
  virtual bool operator!=(MkObject &){return true;}
};
#endif
