
#include "MkContainer.hpp"
// template <class T>
// MkContainer<T>::MkContainer(int size)
// {
//   //  try {
//   if (size <= 0)
//   {
//     MkDebug("::MkContainer<T> - MkContainer(int size)");
//     throw Size(std::string("MkContainer<T>::size below zero strange error"), size);
//   }

//   FSize = size;

//   try
//   {
//     F.reset(new T[FSize]);
//   }
//   catch (std::bad_alloc &a)
//   {
//     MkDebug("MkContainer Contructor memory allocation std::bad_alloc thrown\n");
//     throw Alloc(std::string("MkContainer<T>::bad_alloc error"));
//   }
// }

// template <class T>
// void MkContainer<T>::Initialize(int size)
// {
//   Clear();
//   FSize = size;
//   if (FSize == 0)
//   {
//     F.reset();
//     return;
//   }

//   try
//   {
//     F.reset(new T[FSize]);
//   }
//   catch (std::bad_alloc &a)
//   {
//     MkDebug("MkContainer Contructor memory allocation std::bad_alloc thrown\n");
//     throw Alloc(std::string("MkContainer<T>::bad_alloc error"));
//   }
// }

// template <class T>
// void MkContainer<T>::Clear()
// {
//   FSize = 0;
//   F.reset();
// }

// template <class T>
// T &MkContainer<T>::operator[](int i)
// {
//   if (i >= 0 && i < FSize)
//   {
//     return F[i];
//   }
//   else
//   {
//     MkDebug("MkContainer<T>::operator [%d] out of bound\n", i);
//     throw Range(std::string("MkContainer operator [] Range thrown"), i);
//   }
// }

// template <class T>
// MkContainer<T> &MkContainer<T>::operator=(MkContainer<T> &spheres)
// {
//   int i;
//   FSize = spheres.FSize;

//   try
//   {
//     F.reset(new T[FSize]);
//   }
//   catch (std::bad_alloc &a)
//   {
//     MkDebug("MkContainer Contructor memory allocation std::bad_alloc thrown\n");
//     throw Alloc(std::string("MkContainer<T>::operator=(MkContainer<T>)"));
//   }

//   for (i = 0; i < FSize; i++)
//     F[i] = spheres.F[i];

//   return *this;
// }

// #ifdef __BCPLUSPLUS__
// template <class T>
// void MkContainer<T>::Draw(TObject *Sender)
// {
//   for (int i = 0; i < FSize; i++)
//     F[i].Draw(Sender);
// }
// #endif

// #if defined(_MSC_VER) && defined(_WINDOWS_)
// template <class T>
// void MkContainer<T>::Draw(MkPaint *pb)
// {
//   for (int i = 0; i < FSize; i++)
//     F[i].Draw(pb);
// }
// #endif
