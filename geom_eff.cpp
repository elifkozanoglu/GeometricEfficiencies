#include <iostream>
#include <cmath>
using namespace std;


//VARIABLES
//z distance between Bismuth source and OM
z0 = 435;

//dimensions of the optical module (OM)
double xLength = 23;
double yLength = 24;

//Coordinates of bottom left corner of OM
xShift = 3;
yShift = 4;


double solidAngle(a, b, z){
    double alpha = a/(2*z);
    double beta = b/(2*z);
    return 4*atan(alpha*beta/(sqrt(1+alpha*alpha+beta*beta)));
}
    
double main(){
    double sAngle;
    if(xShift>0 && yShift>0){
        sAngle =  solidAngle(2*(xLength+xShift),2*(yLength+yShift),z0)
        -solidAngle(2*(xShift),2*(yLength+yShift),z0)
        -solidAngle(2*(xLength+xShift),2*(yShift),z0)
        +solidAngle(2*(xShift),2*(yShift),z0);
    }

    if(xShift<0 && yShift<0){
        sAngle =  solidAngle(2*xLength+xShift,2*yLength+yShift,z0)
        +solidAngle(-2*(xShift),2*yLength+yShift,z0)
        +solidAngle(2*xLength+xShift,-2*(yShift),z0)
        -solidAngle(-2*(xShift),-2*(yShift),z0);
    }

    if(xShift<0 && yShift>0){
        sAngle =  solidAngle(2*xLength+xShift,2*(yLength+yShift),z0)
        +solidAngle(-2*(xShift),2*(yLength+yShift),z0)
        -solidAngle(2*xLength+xShift,2*(yShift),z0)
        -solidAngle(-2*(xShift),-2*(yShift),z0); ;
    }

    if(xShift>0 && yShift<0){
        sAngle =  solidAngle(2*(xLength+xShift),2*yLength+yShift,z0)
        -solidAngle(2*(xShift),2*yLength+yShift,z0)
        +solidAngle(2*(xLength+xShift),-2*(yShift),z0)
        -solidAngle(-2*(xShift),-2*(yShift),z0); ;
    }

    return sAngle/4;
    
}