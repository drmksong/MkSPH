//---------------------------------------------------------------------------
#include "MkOpenGL.hpp"
//---------------------------------------------------------------------------
#pragma package(smart_init)

#define PI_ M_PI
void accFrustum(GLdouble left, GLdouble right, GLdouble bottom,
                GLdouble top, GLdouble near_, GLdouble far_, GLdouble pixdx,
                GLdouble pixdy, GLdouble eyedx, GLdouble eyedy, GLdouble focus)
{
  GLdouble xwsize, ywsize;
  GLdouble dx, dy;
  GLint viewport[4];

  glGetIntegerv(GL_VIEWPORT, viewport);
  xwsize = right - left;
  ywsize = top - bottom;
  dx = -(pixdx * xwsize / (GLdouble)viewport[2] + eyedx * near_ / focus);
  dy = -(pixdy * ywsize / (GLdouble)viewport[3] + eyedy * near_ / focus);

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glFrustum(left + dx, right + dx, bottom + dy, top + dy, near_, far_);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glTranslatef(-eyedx, -eyedy, 0.0);
}

void accPerspective(GLdouble fovy, GLdouble aspect, GLdouble near_, GLdouble far_,
                    GLdouble pixdx, GLdouble pixdy, GLdouble eyedx, GLdouble eyedy,
                    GLdouble focus)
{
  GLdouble fov2, left, right, bottom, top;
  fov2 = ((fovy * PI_) / 180.0 / 2.0);

  top = near_ / (cos(fov2) / sin(fov2));
  bottom = -top;
  right = top * aspect;
  left = -right;

  accFrustum(left, right, bottom, top, near_, far_, pixdx, pixdy, eyedx, eyedy, focus);
}
