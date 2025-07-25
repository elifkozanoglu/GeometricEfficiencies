#include <TCanvas.h>
#include <TGraph.h>
#include <TLatex.h>
#include <vector>
#include <utility>
#include <iostream>
#include <cmath>
#include<array>

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
            std::cout << "x pos: "<< val1 <<" y pos: "<< val2 << std::endl;
            
            return {val1, val2};
        }
        ++current_line;
    }

    throw std::out_of_range("Line number " + std::to_string(target_line) + " exceeds file length");
}



void plotting() {
    TCanvas *c1 = new TCanvas("c1", "OM Dot Plot", 1200, 1000);
    c1->SetLeftMargin(0.08);
    c1->SetBottomMargin(0.14);
    c1->SetRightMargin(0.08);
    c1->SetTopMargin(0.14);

    const int n = 260;
    const int s = 42;
    double x_pos[n], y_pos[n];
    double x_source[s], y_source[s];

    // Collect coordinates
    for (int i = 0; i < n; i++) {
        int OM_num = i;
        std::array<double, 3> pos = OMnum_to_position(OM_num);
        x_pos[i] = pos[1];  // X coordinate
        y_pos[i] = pos[2];  // Y coordinate
        //std::cout<<"x pos:"<<x_pos[i]<<"y pos:"<<y_pos[i]<<std::endl;
    }

    for(int i = 0; i < s; i++){
        std::pair<double, double> positions_source = get_values_on_line("source_positions.txt.in", i);
        x_source[i] = positions_source.first;
        y_source[i] = positions_source.second;
        std::cout << "SECOND x pos: "<< x_source[i] <<" y pos: "<< y_source[i] << std::endl;
    }


    // Plot dots
    TGraph *graph = new TGraph(n, x_pos, y_pos);
    graph->SetMarkerStyle(1);
    graph->SetMarkerSize(1.2);
    graph->SetMarkerColor(kBlue);
    graph->SetTitle("OM and Source Positions;X;Y");
    graph->Draw("AP");

    TGraph *crossGraph = new TGraph(s, x_source, y_source);
    crossGraph->SetMarkerStyle(20);      
    crossGraph->SetMarkerSize(1);      
    crossGraph->SetMarkerColor(kGreen+2);
    crossGraph->Draw("P same");

    std::vector<TBox*> boxes;
    for (int i = 0; i < n; ++i) {
        double x_min = x_pos[i] - 128;
        double x_max = x_pos[i] + 128;
        double y_min = y_pos[i] - 128;
        double y_max = y_pos[i] + 128;

        TBox *box = new TBox(x_min, y_min, x_max, y_max);
        box->SetLineColor(kRed);     // Outline color
        box->SetFillStyle(0);        // Transparent fill
        box->Draw("same");
        boxes.push_back(box);        // prevent deletion
    }

    
    // Add labels
    TLatex latex;
    latex.SetTextSize(0.015);
    latex.SetTextAlign(12); // left align

    for (int i = 0; i < n; i++) {
        int OM_num = i;
        std::string label = std::to_string(OM_num);
        latex.DrawLatex(x_pos[i] + 0.3, y_pos[i], label.c_str());
    }
    c1->Update();
}
