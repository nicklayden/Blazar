#pragma once
#include <math.h>
/*
	A Host of functions for testing the root finder functionality

*/

template <class T>
T gaussian(T x, T y) { 
	return exp(-(x*x + y*y)) - 0.5;
}

template <class T>
T paraboloid(T x, T y) {
    return x*x + y*y;
}

template <class T>
T plane_para_ellipse(T x, T y) {
    return paraboloid(x,y) - plane(x,y);
}

template <class T>
T wave_pattern(T x, T y) {
    return sin(x)*cos(y) - 0.6;
}