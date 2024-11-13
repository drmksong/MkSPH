//---------------------------------------------------------------------------
#ifndef MkColorH
#define MkColorH

enum MkColor {clMin=-0x7fffffff-1, clMax=0x7fffffff};
enum CType {ctNone, ctMass, ctEta, ctRho, ctPress, ctXVel, ctYVel, ctZVel, ctForce, ctTemp};
typedef MkColor *PColor;

static const MkColor clBlack = (MkColor)0x0;
static const MkColor clMaroon = (MkColor)0x80;
static const MkColor clGreen = (MkColor)0x8000;
static const MkColor clOlive = (MkColor)0x8080;
static const MkColor clNavy = (MkColor)0x800000;
static const MkColor clPurple = (MkColor)0x800080;
static const MkColor clTeal = (MkColor)0x808000;
static const MkColor clGray = (MkColor)0x808080;
static const MkColor clSilver = (MkColor)0xc0c0c0;
static const MkColor clRed = (MkColor)0xff;
static const MkColor clLime = (MkColor)0xff00;
static const MkColor clYellow = (MkColor)0xffff;
static const MkColor clBlue = (MkColor)0xff0000;
static const MkColor clFuchsia = (MkColor)0xff00ff;
static const MkColor clAqua = (MkColor)0xffff00;
static const MkColor clLtGray = (MkColor)0xc0c0c0;
static const MkColor clDkGray = (MkColor)0x808080;
static const MkColor clWhite = (MkColor)0xffffff;
static const MkColor clNone = (MkColor)0x1fffffff;
static const MkColor clDefault = (MkColor)0x20000000;
#endif
