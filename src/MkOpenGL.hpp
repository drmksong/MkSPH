//---------------------------------------------------------------------------
#ifndef MkOpenGLH
#define MkOpenGLH
#include <GL/gl.h>
#include "MkPoint.hpp"
//---------------------------------------------------------------------------
void accFrustum(GLdouble left, GLdouble right, GLdouble bottom,
                GLdouble top, GLdouble near, GLdouble far, GLdouble pixdx,
                GLdouble pixdy, GLdouble eyedx, GLdouble eyedy, GLdouble focus);
void accPerspective(GLdouble fovy, GLdouble aspect, GLdouble near, GLdouble far,
                    GLdouble pixdx, GLdouble pixdy, GLdouble eyedx, GLdouble eyedy,
                    GLdouble focus);

#endif
