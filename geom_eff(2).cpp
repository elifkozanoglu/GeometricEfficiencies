#include <iostream>
#include <cmath>
#include<array>
//#include "/pbs/home/e/ekozanoglu/ekozanoglu/GeometricalEfficiencies/Used_TK/TKOMhit.h"

//R__LOAD_LIBRARY(/pbs/home/e/ekozanoglu/ekozanoglu/GeometricalEfficiencies/TKEvent/TKEvent/lib/libTKEvent.so);
using namespace std;

//OM dimensions in mm (from TKEvent)
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
    int current_line = 0;

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
    //std::cout << std::setprecision (17)<< std::endl;   
    //std::cout<< "x position:" << OM_pos[2] << "y position:" << OM_pos[1] << std::endl;
    
    //bottom left corner positions
    double blc_y = OM_pos[2]-mw_sizey/2.0;
    double blc_z = OM_pos[1]-mw_sizez/2.0;
    
    //Source number to source position (source_number: 0-41)
    std::pair<double, double> positions = get_values_on_line("source_positions.txt.in", source_number);
    
    //Distances between source and bottom left corner of OM
    double xShift = blc_x-positions.first;
    double yShift = blc_y-positions.second;

    //std::cout<< "x shift:" << xShift << "y shift:" << yShift << std::endl;
    double geometric_eff = solidAngle(yShift, zShift, z0)/(4.0*M_PI);

    return geometric_eff;
}

void geom_eff(){
    
    std::cout << geometricEfficiency_OMS(170, 7) << std::endl;
  
    
}








