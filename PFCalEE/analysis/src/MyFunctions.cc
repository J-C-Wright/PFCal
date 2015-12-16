#define _USE_MATH_DEFINES

#include <utility>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <string>
#include "MyFunctions.h"

using namespace std;

std::pair<float,float> rPhiCoordinates(float x, float y) {

    float phi, r;
    std::pair<float,float> rPhiPoint;

    r = sqrt(pow(x,2)+pow(y,2));
    if (x < 0 && y >= 0) {
        phi = M_PI - atan(fabs(y/x));
    }else if (x < 0 && y < 0) {
        phi = atan(fabs(y/x)) - M_PI;
    }else{
        phi = atan(y/x);
    }

    rPhiPoint.first  = r;
    rPhiPoint.second = r*phi;

    return rPhiPoint;
}

