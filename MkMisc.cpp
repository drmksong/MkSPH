//---------------------------------------------------------------------------
#include "MkMisc.hpp"

#ifdef __BCPLUSPLUS__
#include <dialogs.hpp>
#include "AdvGrid.hpp"
#endif

bool is_eq(double i, double j)
{
  return (fabs(i - j) < FTOL) ? true : false;
}

#ifndef min
int min(int x, int y)
{
  return (x >= y) ? y : x;
}

double min(double x, double y)
{
  return (x >= y) ? y : x;
}
#endif

#ifndef max
int max(int x, int y)
{
  return (x >= y) ? x : y;
}
double max(double x, double y)
{
  return (x >= y) ? x : y;
}
#endif

#if !defined(_MSC_VER) && !defined(_WINDOWS_)
void dull(const char *__format, ...)
{
  printf("dull is called\n");
  return;
}
#endif

void dumprintf(const char *__format, ...)
{
  return;
}

void swap(int &a, int &b)
{
  int tmp;
  tmp = a;
  a = b;
  b = tmp;
  return;
}

void swap(double &a, double &b)
{
  double tmp;
  tmp = a;
  a = b;
  b = tmp;
  return;
}

#ifdef __BCPLUSPLUS__
void swap(TColor a, TColor b)
{
  TColor tmp;
  tmp = a;
  a = b;
  b = tmp;
  return;
}

void Swap(TObject *Sender, int i, int j)
{
  int k;
  AnsiString str;

  TAdvStringGrid *grid = dynamic_cast<TAdvStringGrid *>(Sender);
  if (!grid)
    return;

  for (k = 1; k < grid->ColCount; k++)
  {
    str = grid->Cells[k][i];
    grid->Cells[k][i] = grid->Cells[k][j];
    grid->Cells[k][j] = str;
  }
}
#endif

int delta(int a, int b)
{
  return (a == b) ? 1 : 0;
}

bool ExtractFileExt(char *ext, char *str)
{
  char *e, s[256];
  e = strchr(str, '.');
  strcpy(ext, e);
  sprintf(s, "Extension of %s is %s\n", str, ext);
  MkDebug(s);
  if (e)
    return true;
  else
    return false;
}

bool TrimLeft(char *&dest, char *src)
{
  if (!src)
    return false;
  for (dest = src; *dest == ' ' || *dest == '\t'; dest++)
    ;
  if (dest)
    return true;
  else
    return false;
}

bool ExtractStr(char *des, int n, char *src)
{
  char *s, *e;
  int len, c = 0;

  *des = '\0';
  s = e = src;
  while (n > c && e)
  {
    s = e;
    TrimLeft(e, s);
    len = strlen(e);
    s = e;
    e = strchr(s, ' ');
    c++;
  }

  if (!e)
    return false;

  TrimLeft(e, e);
  strcpy(des, e);
  e = strchr(des, ' ');

  if (e == NULL)
    return false;
  else if (*e == '\0')
    return false;
  else if (*e == ' ')
    *e = '\0';
  return true;
}

bool ToLower(char *str)
{
  int len;
  char c;
  len = strlen(str);
  if (len <= 0)
    return false;

  for (int i = 0; i < len && str[i]; i++)
  {
    c = str[i];
    str[i] = tolower(c);
  }

  return true;
}

bool ToUpper(char *str)
{
  int len;
  char c;
  len = strlen(str);
  if (len <= 0)
    return false;

  for (int i = 0; i < len && str[i]; i++)
  {
    c = str[i];
    str[i] = toupper(c);
  }

  return true;
}

#ifdef __BCPLUSPLUS__
bool ToOnlyAlpha(AnsiString &dest, AnsiString &src)
{
  int i;
  dest = "";
  if (src.Length() <= 0)
    return false;
  for (i = 1; i <= src.Length(); i++)
  {
    if (isalpha(src[i]) || isdigit(src[i]))
      dest = dest + src[i];
  }
  return true;
}

bool RemoveAnd(AnsiString &dest, AnsiString &src)
{
  int i;
  dest = "";
  if (src.Length() <= 0)
    return false;
  for (i = 1; i <= src.Length(); i++)
  {
    if (src[i] != '\&')
      dest = dest + src[i];
  }
  return true;
}

#else
bool ToOnlyAlpha(char *&dest, char *src)
{
  return false;
}
bool RemoveAnd(char *&dest, char *src)
{
  return false;
}

#endif

bool CompSub(char *str, char *txt)
{
  char *a, *b;
  int len;

  if (!strlen(str))
    return false;
  if (!strlen(txt))
    return false;

  len = strlen(str);
  len = min(len, strlen(txt));

  a = new char[len + 1];
  if (!a)
    return false;
  b = new char[len + 1];
  if (!b)
  {
    delete a;
    return false;
  }

  strncpy(a, str, len);
  strncpy(b, txt, len);

  ToLower(a);
  ToLower(b);

  return !strcmp(a, b);

  delete a;
  delete b;
}

int NumOfParam(char *str)
{
  int n = 0;
  char *s, *e;
  int len;

  s = e = str;
  while (e)
  {
    s = e;
    TrimLeft(e, s);
    len = strlen(e);
    s = e;
    e = strchr(s, ' ');
    n++;
  }
  return n;
}

int NumOfParen(char *str)
{
  int i, nl = 0, nr = 0;
  for (i = 0; i < strlen(str); i++)
  {
    if (str[i] == '(')
      nl++;
    if (str[i] == ')')
      nr++;
  }
  if (nl == nr)
    return nl;
  else
    return -1;
}

bool ExtractFromParen(char *str, int n, double &x, double &y)
{
  char *sp = NULL, *ep = NULL;
  int i, j = 0, nl = 0, nr = 0;

  x = 0;
  y = 0;
  for (i = 0; i < strlen(str); i++)
  {
    if (str[i] == '(')
    {
      if (n == nl)
        sp = str + i + 1;
      nl++;
    }
    if (str[i] == ')')
    {
      if (n == nr)
      {
        while (str[i - j] != ' ' && str[i - j] != ',')
          j++;
        if (str[i - j] == ' ' || str[i - j] == ',')
          ep = str + i - j;
      }
      nr++;
    }
  }
  if (!sp || !ep)
    return false;
  sscanf(sp, "%f ", &x);
  sscanf(ep, "%f ", &y);
  return true;
}

#ifdef __BCPLUSPLUS__
bool IsNumber(AnsiString str)
{
  if (!str.Length())
    return false;
  if (str == "-")
    return false;
  for (int i = 1; i <= str.Length(); i++)
  {
    if (!isdigit(str[i]) && str[i] != '-' && str[i] != '.')
      return false;
  }
  return true;
}

AnsiString ShortSteelName(AnsiString str)
{
  int i, pos, len;
  AnsiString sub, sub2, sub3, res;
  pos = str.Pos("��");
  len = str.Length();

  sub = str.SubString(1, pos + 1);
  sub2 = str.SubString(pos + 2, len - pos - 1);

  pos = sub2.Pos("��");
  len = sub2.Length();
  sub3 = sub2.SubString(1, pos - 1);
  res = sub + sub3;
  return res;
}

#else
char *ShortSteelName(char *str)
{
  int i, pos, len;
  static char res[256];
  // need to be implemented
  return res;
}
#endif

bool IsNumber(char *str)
{
  if (!strlen(str))
    return false;
  if (!strcmp(str, "-"))
    return false;
  for (int i = 0; i < strlen(str); i++)
  {
    if (!isdigit(str[i]) && str[i] != '-' && str[i] != '.')
      return false;
  }
  return true;
}

double ShapeFun1(double x, double l)
{
  double p;
  if (fabs(l) < EPS)
    return 0;
  p = x / l;
  return 1 - (3 - 2 * p) * p * p;
}

double ShapeFun2(double x, double l)
{
  double p;
  if (fabs(l) < EPS)
    return 0;
  p = x / l;
  return -x * (1 - p) * (1 - p);
}

double ShapeFun3(double x, double l)
{
  double p;
  if (fabs(l) < EPS)
    return 0;
  p = x / l;
  return (3 - 2 * p) * p * p;
}

double ShapeFun4(double x, double l)
{
  double p;
  if (fabs(l) < EPS)
    return 0;
  p = x / l;
  return x * p * (1 - p);
}

double ShapeFunInteg1(double aj, double bj, double l, double lj1, double lj)
{
  double lj1_2, lj1_3, lj1_4, lj1_5, lj_2, lj_3, lj_4, lj_5, l_2, l_3;
  if (fabs(l) < EPS)
    return 0;
  if (lj1 < lj)
    swap(lj1, lj);

  lj1_2 = lj1 * lj1;
  lj1_3 = lj1_2 * lj1;
  lj1_4 = lj1_2 * lj1_2;
  lj1_5 = lj1_4 * lj1;
  lj_2 = lj * lj;
  lj_3 = lj_2 * lj;
  lj_4 = lj_2 * lj_2;
  lj_5 = lj_4 * lj;
  l_2 = l * l;
  l_3 = l_2 * l;

  return aj * ((lj1_2 - lj_2) / 2 - 3 * (lj1_4 - lj_4) / 4 / l_2 + 2 * (lj1_5 - lj_5) / 5 / l_3) + bj * ((lj1 - lj) - (lj1_3 - lj_3) / l_2 + (lj1_4 - lj_4) / 2 / l_3);
}

double ShapeFunInteg2(double aj, double bj, double l, double lj1, double lj)
{
  double lj1_2, lj1_3, lj1_4, lj1_5, lj_2, lj_3, lj_4, lj_5, l_2, l_3;
  if (fabs(l) < EPS)
    return 0;
  if (lj1 < lj)
    swap(lj1, lj);

  lj1_2 = lj1 * lj1;
  lj1_3 = lj1_2 * lj1;
  lj1_4 = lj1_2 * lj1_2;
  lj1_5 = lj1_4 * lj1;
  lj_2 = lj * lj;
  lj_3 = lj_2 * lj;
  lj_4 = lj_2 * lj_2;
  lj_5 = lj_4 * lj;
  l_2 = l * l;
  l_3 = l_2 * l;

  return aj * (-(lj1_3 - lj_3) / 3 + (lj1_4 - lj_4) / 2 / l - (lj1_5 - lj_5) / 5 / l_2) + bj * (-(lj1_2 - lj_2) / 2 + 2 * (lj1_3 - lj_3) / 3 / l - (lj1_4 - lj_4) / 4 / l_2);
}

double ShapeFunInteg3(double aj, double bj, double l, double lj1, double lj)
{
  double lj1_2, lj1_3, lj1_4, lj1_5, lj_2, lj_3, lj_4, lj_5, l_2, l_3;
  if (fabs(l) < EPS)
    return 0;
  if (lj1 < lj)
    swap(lj1, lj);

  lj1_2 = lj1 * lj1;
  lj1_3 = lj1_2 * lj1;
  lj1_4 = lj1_2 * lj1_2;
  lj1_5 = lj1_4 * lj1;
  lj_2 = lj * lj;
  lj_3 = lj_2 * lj;
  lj_4 = lj_2 * lj_2;
  lj_5 = lj_4 * lj;
  l_2 = l * l;
  l_3 = l_2 * l;

  return aj * (3 * (lj1_4 - lj_4) / 4 / l_2 - 2 * (lj1_5 - lj_5) / 5 / l_3) + bj * ((lj1_3 - lj_3) / l_2 - (lj1_4 - lj_4) / 2 / l_3);
}

double ShapeFunInteg4(double aj, double bj, double l, double lj1, double lj)
{
  double lj1_2, lj1_3, lj1_4, lj1_5, lj_2, lj_3, lj_4, lj_5, l_2, l_3;
  if (fabs(l) < EPS)
    return 0;
  if (lj1 < lj)
    swap(lj1, lj);

  lj1_2 = lj1 * lj1;
  lj1_3 = lj1_2 * lj1;
  lj1_4 = lj1_2 * lj1_2;
  lj1_5 = lj1_4 * lj1;
  lj_2 = lj * lj;
  lj_3 = lj_2 * lj;
  lj_4 = lj_2 * lj_2;
  lj_5 = lj_4 * lj;
  l_2 = l * l;
  l_3 = l_2 * l;

  return aj * ((lj1_4 - lj_4) / 4 / l - (lj1_5 - lj_5) / 5 / l_2) + bj * ((lj1_3 - lj_3) / 3 / l - (lj1_4 - lj_4) / 4 / l_2);
}

double sosu1(double a)
{
  double b, n, m;

  m = a * 10.0;
  b = modf(m, &n);

  if (b >= 0.4999)
    n = n + 1.0;

  return n / 10.0;
}

double sosu2(double a)
{
  double b, n, m;

  m = a * 100.0;
  b = modf(m, &n);

  if (b >= 0.4999)
    n = n + 1.0;

  return n / 100.0;
}

double sosu3(double a)
{
  double b, n, m;

  m = a * 1000.0;
  b = modf(m, &n);

  if (b >= 0.4999)
    n = n + 1.0;

  return n / 1000.0;
}

double blength(double tlen, double alen)
{
  double b;

  b = tlen - alen;
  if (alen > b)
    ;
  else
  {
    for (int i = 1;; i++)
    {
      b = (tlen - 2.0 * alen) / i;
      if (alen > b)
        break;
    }
  }

  return b;
}

//---------------------------------------------------------------------------
