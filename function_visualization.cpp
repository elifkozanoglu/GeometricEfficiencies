#include <iostream>
#include <cmath>
#include<array>
#include "TH2D.h"
#include "TCanvas.h"

const double z0 = 435.0;


double centerGeometricEfficiency(double a, double b, double z){
    double alpha = abs(a/(2.0*z));
    double beta = abs(b/(2.0*z));
    double solidAngle = 4.0*atan(alpha*beta/(sqrt(1.0+alpha*alpha+beta*beta)));
    return solidAngle/(4.0*M_PI);
}

double before_integration(double x, double y, double z){
    return z/pow(x*x+y*y+z*z, 3/2);
    
}


void function_visualization(){
    int weight = 10000;
    int nbinsx = 2000;
    int nbinsy = 2000;

    double xlow = -2500;
    double xup = 2500;
    double ylow = -2000;
    double yup = 2000;
    
    TH2D *hist = new TH2D("hist", "Efficiency", nbinsx, xlow, xup, nbinsy, ylow, yup);

    for(int i = 0; i < nbinsx; i++){
        for(int j = 0; j < nbinsy; j++){
            double x = hist->GetXaxis()->GetBinCenter(i);
            double y = hist->GetYaxis()->GetBinCenter(j);
            double f = before_integration(2.0*x,2.0*y,z0);
            hist->SetBinContent(i, j, int(100*log10(f)));
        }
    }

    TCanvas *c2 = new TCanvas("c2", "Geometric Efficiency Visualization", 1200, 1000);
    hist->Draw("COLZ");
    c2->Update();
    
}