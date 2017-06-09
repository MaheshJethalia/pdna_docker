#include "config.h"

RotationalAngle::RotationalAngle() {
    alpha = beta = gamma = 0;
}

RotationalAngle::RotationalAngle( int a, int b, int g){
    alpha = a; beta = b; gamma = g;
}
