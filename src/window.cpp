/**
 * 
 * OpenGL window for blazar.
 * author: Nick Layden
 * 
 * 
*/

// Graphics libraries


#include "window.hpp"

void initGL() {
    glClearDepth(1.f);
    glClearColor(1., 1., 1., 1.);
    glEnable(GL_DEPTH_TEST);
    glDepthMask(GL_TRUE);
    glEnable(GL_BLEND);
    glHint(GL_LINE_SMOOTH, GL_NICEST);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45,1,0.1,100);
    glPointSize(5);
    glLineWidth(2);
}
