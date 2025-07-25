/*This code visualizes total or source-wise geometric efficiencies of the optical module source pairs. 

Note on numbering: The OMs are numbered 0-259 starting from bottom left corner and column-wise. (convention of TKEvent by Tomas)
Sources are numbered 0-41 starting from top left corner and row-wise. (Convention of CalibrationModule by Filip)

Also it should be noted that half of the optical modules that are in the top and bottom row are covered by other walls and this code DOES NOT account for that yet.
*/
#include <iostream>
#include <cmath>
#include<array>
#include "TH2D.h"
#include "TCanvas.h"

using namespace std;

//OM dimensions in mm (from TKEvent) (not all of these are used)
const double mw_sizex = 194.0;
const double mw_sizey = 256.0;
const double mw_sizez = 256.0;

const double gv_sizex = 308.0;
const double gv_sizey = 310.0;
const double gv_sizez = 150.0;

const double xw_sizex = 200.0;
const double xw_sizey = 150.0;
const double xw_sizez = 208.5;

const double yLength = 256.0;
const double zLength = 256.0;

//z distance between Bismuth sources and OMs
const double x0 = 435.0;

//reads values of source positions from .txt.in file
std::pair<double, double> get_values_on_line(const std::string& filename, int target_line) {
    std::ifstream infile(filename);
    if (!infile) {
        throw std::runtime_error("Could not open file: " + filename);
    }

    std::string line;
    int current_line = 0; //zero-based indexing 

    while (std::getline(infile, line)) {
        if (current_line == target_line) {
            std::istringstream ss(line);
            std::string token;
            double val1, val2;

            if (std::getline(ss, token, ';')) {
                val1 = std::stod(token);
            } else {
                throw std::runtime_error("Missing first value on line " + std::to_string(target_line));
            }

            if (std::getline(ss, token)) {
                val2 = std::stod(token);
            } else {
                throw std::runtime_error("Missing second value on line " + std::to_string(target_line));
            }

            return {val1, val2};
        }
        ++current_line;
    }

    throw std::out_of_range("Line number " + std::to_string(target_line) + " exceeds file length");
}

//turns OM_number to OM positions
std::array<double,3> OMnum_to_position(int OM_num){
    array<double,4> SWCR;
    array<double,3> xyz;
    //mainwall IT
	if(OM_num < 260) 
	{
		SWCR[0] = 0;
		SWCR[1] = -1;
		SWCR[2] = OM_num / 13;
		SWCR[3] = OM_num % 13;
	}
	//mainwall FR
	else if(OM_num < 520)
	{
		SWCR[0] = 1;
		SWCR[1] = -1;
		SWCR[2] = (OM_num - 260) / 13;
		SWCR[3] = (OM_num - 260) % 13;
	}
	//Xcalo IT	
	else if(OM_num < 584)
	{
		SWCR[0] = 0;
		SWCR[1] = (OM_num - 520) / 32;
		SWCR[2] = ((OM_num - 520) / 16) % 2;
		SWCR[3] = (OM_num -520) % 16;
	}
	//Xcalo FR
	else if(OM_num < 648)
	{
		SWCR[0] = 1;
		SWCR[1] = (OM_num - 520 - 64) / 32;
		SWCR[2] = ((OM_num - 520 - 64) / 16) % 2;
		SWCR[3] = (OM_num -520 - 64) % 16;
	}
	//GVeto IT
	else if(OM_num < 680)
	{
		SWCR[0] = 0;
		SWCR[1] = (OM_num - 520 - 128) / 16;
		SWCR[2] = (OM_num - 520 - 128) % 16;
		SWCR[3] = -1;
	}
	//GVeto FR
	else if(OM_num < 712)
	{
		SWCR[0] = 1;
		SWCR[1] = (OM_num - 520 - 128 - 32) / 16;
		SWCR[2] = (OM_num - 520 - 128 - 32) % 16;
		SWCR[3] = -1;
	}
	
	int OM_type;
	
	if(OM_num < 520)
	{
		OM_type = 1302;
	}
	else if(OM_num < 648)
	{
		OM_type = 1232;
	}
	else
	{
		OM_type = 1252;
	}

	switch(OM_type)
	{
		case 1302: //MW
			if(SWCR[0] == 1)
				xyz[0] = 532.0;
			else
				xyz[0] = -532.0;
				xyz[1] = ((double)SWCR[2]- 9.5) * 259.0;
				xyz[2] = ((double)SWCR[3] - 6) * 259.0;
				
			break;
			
		case 1232: //XW
			if(SWCR[1] == 1)
				xyz[1] = 2580.5;
			else
				xyz[1] = -2580.5;
				
			if(SWCR[0] == 1)
			{
				if(SWCR[2] == 1)
					xyz[0] = 333.0;
				else
					xyz[0] = 130.0;
			}
			else
			{
				if(SWCR[2] == 1)
					xyz[0] = -333.0;

				else
					xyz[0] = -130.0;
			}
			
			xyz[2] = ((double)SWCR[3] - 7.5) * 212.0;
			
			break;
			
		case 1252: //GV
			if(SWCR[0] == 1)
				xyz[0] = 213.5;
			else
				xyz[0] = -213.5;
			if(SWCR[1] == 1)
				xyz[2] = 1625.0;
			else
				xyz[2] = -1625.0;
			if(SWCR[2] > 7)
				xyz[1] = 161.0 + (((double)SWCR[2]-8) * 311.5);
			else
				xyz[1] = -161.0 + (((double)SWCR[2]-7) * 311.5);
			break;	

	}
    return xyz;
}

//calculates the solid angle seen at the origin for a rectangle of dimenions a x b and z away from origin
double centerSolidAngle(double a, double b, double z){
    double alpha = a/(2.0*z);
    double beta = b/(2.0*z);
    return 4.0*atan(alpha*beta/(sqrt(1.0+alpha*alpha+beta*beta)));
}

double solidAngle(double yShift, double zShift, double x){
    double sAngle =  centerSolidAngle(2.0*(yLength+yShift),2.0*(zLength+zShift),x)
        -centerSolidAngle(2.0*(yShift),2.0*(zLength+zShift),x)
        -centerSolidAngle(2.0*(yLength+yShift),2.0*(zShift),x)
        +centerSolidAngle(2.0*(yShift),2.0*(zShift),x);
    
    return sAngle/4.0;
}

double geometricEfficiency_OMS(int OM_number, int source_number){
    //OM number to OM position
    std::array<double,3> OM_pos;
    OM_pos = OMnum_to_position(OM_number);
    
    //bottom left corner positions
    double blc_y = OM_pos[1]-mw_sizey/2.0;
    double blc_z = OM_pos[2]-mw_sizez/2.0;
    
    //Source number to source position (source_number: 0-41)
    std::pair<double, double> positions = get_values_on_line("source_positions.txt.in", source_number);
    
    //Distances between source and bottom left corner of OM
    double yShift = blc_y-positions.first;
    double zShift = blc_z-positions.second;

    double geometric_eff = solidAngle(yShift, zShift, x0)/(4.0*M_PI);

    return geometric_eff;
}

std::array<double, 260> geom_eff(){
    std::array<double, 260> total_efficiencies;
    //activities of sources as measured by Miro on July 1 2018
    //
    double activities[42] = {129.2, 127.9, 138.5, 125.8, 134.3, 128.4, 
                            129.5, 131.0, 125.0, 129.8, 124.7, 132.0, 
                            133.8, 126.6, 131.3, 130.1, 133.0, 122.8, 
                            122.2, 139.1, 130.9, 140.9, 119.8, 135.0, 
                            132.5, 127.4, 131.9, 130.0, 132.6, 128.4, 
                            123.7, 131.1, 125.5, 130.7, 124.9, 132.2, 
                            131.4, 128.0, 137.3, 128.3, 136.3, 127.1};
    
    for(int i = 0; i < 260; i++){
        
        total_efficiencies[i] = 0;
        
        for(int j = 0; j < 42; j++){
            total_efficiencies[i] += activities[j]*geometricEfficiency_OMS(i, j);
            std::cout << "OM " << i << "src " << j << " srcact " << activities[j] << " geomeff " << geometricEfficiency_OMS(i, j) << std::endl;
        }
    }
    return total_efficiencies;
}

void total_eff_visualization(){
    //PICK WHETHER YOU WANT TOTAL(true) OR SOURCE-WISE(false) VISUALIZATION
    bool total_vis = true;
    
    int weight = 1;
    int nbinsx = 20;
    int nbinsy = 13;

    double xlow = -2560;
    double xup = 2560;
    double ylow = -1664;
    double yup = 1664;
    
    TH2D *hist = new TH2D("hist", "Total Efficiency", nbinsx, xlow, xup, nbinsy, ylow, yup);
    
    //TOTAL VISUALIZIATION
    if(total_vis){
        std::array<double, 260> total_eff = geom_eff();
    
        for (int i = 0; i < 260; ++i) {
            int y = i / 13;      // 0 to 19 (column)
            int z = i % 13;      // 0 to 13 (row)
            double a = weight*total_eff[i];
            hist->SetBinContent(y + 1, z + 1, a);
        }        
    }

    //VISUALIZATION FOR EACH SOURCE
    if(!total_vis){
        int source_num = 15;
        double threshold = 0.0001;
        for (int i = 0; i < 260; ++i) {
            int y = i / 13;      // 0 to 19 (column)
            int z = i % 13;      // 0 to 13 (row)
            double a = geometricEfficiency_OMS(i, source_num);
            if(a>threshold){
                hist->SetBinContent(y + 1, z + 1, a);
            } 
        }
    }

    TCanvas *c2 = new TCanvas("c2", "Geometric Efficiency Visualization", 1200, 1000);
    hist->Draw("COLZ");
    c2->Update();
    
}