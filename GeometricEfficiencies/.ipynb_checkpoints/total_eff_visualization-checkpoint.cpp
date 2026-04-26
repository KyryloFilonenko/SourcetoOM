// -*- coding: utf-8 -*-
/*This code visualizes total or source-wise geometric efficiencies of the optical module source pairs. 

Note on numbering: The OMs are numbered 0-259 starting from bottom left corner and column-wise. (convention of TKEvent by Tomas)
Sources are numbered 0-41 starting from top left corner and row-wise. (Convention of CalibrationModule by Filip)

Also it should be noted that half of the optical modules that are in the top and bottom row are covered by other walls and this code DOES NOT account for that yet.

ENHANCED: Added ROOT file output for data storage and further analysis
*/

#include <iostream>
#include <cmath>
#include <array>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <cctype>

#include "TH2F.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TLatex.h"
#include "TBox.h"
#include "TLine.h"
#include "TGraph.h"
#include "TPaveText.h"
#include "TColor.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TSystem.h"

using namespace std;

enum class DisplayMode {
    PERCENT,   // відсотки (як раніше)
    LOG10,     // log10 від відсотків
    PPM        // parts per million (= eps_percent * 1e4)
};

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

// Trim whitespace from both ends
std::string trim(const std::string& str) {
    size_t first = str.find_first_not_of(" \t\r\n");
    if (std::string::npos == first) {
        return str;
    }
    size_t last = str.find_last_not_of(" \t\r\n");
    return str.substr(first, (last - first + 1));
}

std::string displayLabel(double eps_percent, DisplayMode mode) {
    char buf[64];
    switch (mode) {
        case DisplayMode::PERCENT:
            snprintf(buf, sizeof(buf), "%.2f", eps_percent);
            break;
        case DisplayMode::LOG10:
            if (eps_percent > 0.0)
                snprintf(buf, sizeof(buf), "%.1f", std::log10(eps_percent));
            else
                snprintf(buf, sizeof(buf), "-#infty");
            break;
        case DisplayMode::PPM:
            snprintf(buf, sizeof(buf), "%.0f", eps_percent * 10000.0);
            break;
    }
    return std::string(buf);
}

double displayValue(double eps_percent, DisplayMode mode) {
    switch (mode) {
        case DisplayMode::PERCENT: return eps_percent;
        case DisplayMode::LOG10:   return (eps_percent > 0.0) ? std::log10(eps_percent) : -6.0;
        case DisplayMode::PPM:     return eps_percent * 10000.0;
    }
    return eps_percent;
}

const char* modeSuffix(DisplayMode mode) {
    switch (mode) {
        case DisplayMode::PERCENT: return "pct";
        case DisplayMode::LOG10:   return "log";
        case DisplayMode::PPM:     return "ppm";
    }
    return "pct";
}

struct SourcePosition {
    double y;
    double z;
};

std::vector<SourcePosition> read_source_positions(const std::string& filename) {
    std::vector<SourcePosition> positions;
    std::ifstream infile(filename);
    
    if (!infile) {
        throw std::runtime_error("Could not open file: " + filename);
    }
    
    std::string line;
    int lineNum = 0;
    
    while (std::getline(infile, line)) {
        // Skip empty lines and comments
        if (line.empty() || line[0] == '#' || line[0] == ';') {
            continue;
        }
        
        std::istringstream iss(line);
        int id;
        double y, z;
        
        if (!(iss >> id >> y >> z)) {
            throw std::runtime_error("Cannot parse line " + std::to_string(lineNum) + 
                                   ": expected 'id y z'. Line: [" + line + "]");
        }
        
        positions.push_back({y, z});
        lineNum++;
    }
    
    if (positions.size() != 42) {
        throw std::runtime_error("Expected 42 sources, got " + std::to_string(positions.size()));
    }
    
    infile.close();
    return positions;
}

//turns OM_number to OM positions
std::array<double,3> OMnum_to_position(int OM_num){
    array<double,4> SWCR;
    array<double,3> xyz;
    
	if(OM_num < 260) 
	{
		SWCR[0] = 0;
		SWCR[1] = -1;
		SWCR[2] = OM_num / 13;
		SWCR[3] = OM_num % 13;
	}
	else if(OM_num < 520)
	{
		SWCR[0] = 1;
		SWCR[1] = -1;
		SWCR[2] = (OM_num - 260) / 13;
		SWCR[3] = (OM_num - 260) % 13;
	}
	else if(OM_num < 584)
	{
		SWCR[0] = 0;
		SWCR[1] = (OM_num - 520) / 32;
		SWCR[2] = ((OM_num - 520) / 16) % 2;
		SWCR[3] = (OM_num -520) % 16;
	}
	else if(OM_num < 648)
	{
		SWCR[0] = 1;
		SWCR[1] = (OM_num - 520 - 64) / 32;
		SWCR[2] = ((OM_num - 520 - 64) / 16) % 2;
		SWCR[3] = (OM_num -520 - 64) % 16;
	}
	else if(OM_num < 680)
	{
		SWCR[0] = 0;
		SWCR[1] = (OM_num - 520 - 128) / 16;
		SWCR[2] = (OM_num - 520 - 128) % 16;
		SWCR[3] = -1;
	}
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

//calculates the solid angle seen at the origin for a rectangle of dimensions a x b and z away from origin
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

double geometricEfficiency_OMS(int OM_number, int source_number, 
                               const std::vector<SourcePosition>& sourcePositions){
    std::array<double,3> OM_pos = OMnum_to_position(OM_number);
    
    double blc_y = OM_pos[1] - mw_sizey/2.0;
    double blc_z = OM_pos[2] - mw_sizez/2.0;
    
    double yShift = blc_y - sourcePositions[source_number].y;
    double zShift = blc_z - sourcePositions[source_number].z;

    double geometric_eff = solidAngle(yShift, zShift, x0) / (4.0 * M_PI);

    return geometric_eff;
}

void draw_single_source_efficiency(
    int source_number,
    const std::vector<SourcePosition>& sourcePositions,
    DisplayMode dispMode = DisplayMode::PERCENT)
{
    gROOT->SetBatch(true);
    gStyle->SetPalette(kBird);
    
    gSystem->mkdir("plots", true);
    gSystem->mkdir("plots/efficiencies", true);

    std::vector<double> efficiencies(260);
    double eps_min = 1e30, eps_max = -1e30;
    
    std::cout << "Calculating efficiencies for source " << source_number << "...\n";
    
    for (int om = 0; om < 260; om++) {
        efficiencies[om] = geometricEfficiency_OMS(om, source_number, sourcePositions) * 100.0;  // у %
        if (efficiencies[om] > 0.0) {
            if (efficiencies[om] < eps_min) eps_min = efficiencies[om];
            if (efficiencies[om] > eps_max) eps_max = efficiencies[om];
        }
    }

    if (eps_min > eps_max) {
        eps_min = 0.0; eps_max = 1.0;
    }
    if (std::abs(eps_max - eps_min) < 1e-12) eps_max = eps_min + 1.0;

    const int nCols = 20;
    const int nRows = 13;
    const int colMin = 0;
    const int colMax = nCols - 1;
    const int rowMin = 0;
    const int rowMax = nRows - 1;

    const double mL      = 0.08;
    const double mR      = 0.25;
    const double mT      = 0.08;
    const double mB      = 0.08;
    const int    baseH   = 600;
    const double areaH   = baseH * (1.0 - mT - mB);
    const double areaW   = areaH * (double)nCols / (double)nRows;
    const int    canvasW = (int)(areaW / (1.0 - mL - mR)) + 1;
    const int    canvasH = baseH;

    TCanvas* c = new TCanvas(
        Form("cEff_src_%d", source_number), "", canvasW, canvasH);
    c->SetLeftMargin(mL);
    c->SetRightMargin(mR);
    c->SetTopMargin(mT);
    c->SetBottomMargin(mB);

    TH2F* hFrame = new TH2F(
        Form("hFrame_src_%d", source_number), "",
        nCols, colMin - 0.5, colMax + 0.5,
        nRows, rowMin - 0.5, rowMax + 0.5);
    hFrame->SetDirectory(0);
    hFrame->SetStats(0);
    hFrame->GetXaxis()->SetNdivisions(nCols, false);
    hFrame->GetYaxis()->SetNdivisions(nRows, false);
    hFrame->GetXaxis()->SetLabelSize(0.048);
    hFrame->GetYaxis()->SetLabelSize(0.048);
    hFrame->GetXaxis()->SetTickLength(0.01);
    hFrame->GetYaxis()->SetTickLength(0.01);
    
    for (int ic = 0; ic < nCols; ic++)
        hFrame->GetXaxis()->SetBinLabel(ic + 1, Form("%d", colMin + ic));
    for (int ir = 0; ir < nRows; ir++)
        hFrame->GetYaxis()->SetBinLabel(ir + 1, Form("%d", rowMin + ir));
    
    hFrame->GetXaxis()->SetTitle("Column");
    hFrame->GetYaxis()->SetTitle("Row");
    hFrame->Draw("AXIS");

    const double cellW = 0.46;
    const double cellH = 0.46;
    const int nColors = TColor::GetNumberOfColors();

    // Малюємо комірки з кольорами
    for (int om = 0; om < 260; om++) {
        int y = om / 13;      // 0 to 19 (column)
        int z = om % 13;      // 0 to 12 (row)
        
        if (efficiencies[om] <= 0.0) continue;

        double v    = displayValue(efficiencies[om], dispMode);
        double frac = (v - displayValue(eps_min, dispMode)) / 
                     (displayValue(eps_max, dispMode) - displayValue(eps_min, dispMode));
        frac = std::max(0.0, std::min(1.0, frac));

        int colorIdx = TColor::GetColorPalette((int)(frac * (nColors - 1)));

        double bx = colMin + y;
        double by = z;

        TBox* b = new TBox(bx - cellW, by - cellH, bx + cellW, by + cellH);
        b->SetFillColor(colorIdx);
        b->SetLineColor(kGray + 2);
        b->SetLineWidth(1);
        b->Draw("l same");
    }

    // Сітка
    for (int ic = 0; ic <= nCols; ic++) {
        double x = colMin - 0.5 + ic;
        TLine* l = new TLine(x, rowMin - 0.5, x, rowMax + 0.5);
        l->SetLineColor(kBlack); l->SetLineWidth(1); l->SetLineStyle(2);
        l->Draw("same");
    }
    for (int ir = 0; ir <= nRows; ir++) {
        double y = rowMin - 0.5 + ir;
        TLine* l = new TLine(colMin - 0.5, y, colMax + 0.5, y);
        l->SetLineColor(kBlack); l->SetLineWidth(1); l->SetLineStyle(2);
        l->Draw("same");
    }

    // Числові підписи в комірках
    for (int om = 0; om < 260; om++) {
        int y = om / 13;
        int z = om % 13;
        
        if (efficiencies[om] <= 0.0) continue;

        double bx = colMin + y;
        double by = z;

        double v    = displayValue(efficiencies[om], dispMode);
        double frac = (v - displayValue(eps_min, dispMode)) / 
                     (displayValue(eps_max, dispMode) - displayValue(eps_min, dispMode));
        frac = std::max(0.0, std::min(1.0, frac));

        std::string lbl = displayLabel(efficiencies[om], dispMode);
        TLatex* lt = new TLatex(bx, by, lbl.c_str());
        lt->SetTextAlign(22);
        lt->SetTextSize(0.025);
        lt->SetTextColor(frac > 0.55 ? kBlack : kWhite);
        lt->Draw("same");
    }

    // ── Червона крапка — істинне положення джерела ──────────────────────────
    // Перетворення координат джерела (мм) в систему осей графіку:
    //   column = y / 259.0 + 9.5   (колонки 0..19)
    //   row    = z / 259.0 + 6.0   (рядки  0..12)
    {
        double src_y   = sourcePositions[source_number].y;
        double src_z   = sourcePositions[source_number].z;
        double col_src = src_y / 259.0 + 9.5;
        double row_src = src_z / 259.0 + 6.0;

        // Маркер — заповнене коло
        TGraph* gSrc = new TGraph(1);
        gSrc->SetPoint(0, col_src, row_src);
        gSrc->SetMarkerStyle(20);   // заповнене коло
        gSrc->SetMarkerColor(kRed);
        gSrc->SetMarkerSize(0.7);
        gSrc->Draw("P same");
    }
    // ────────────────────────────────────────────────────────────────────────

    // Заголовок
    TLatex ttl;
    ttl.SetNDC();
    ttl.SetTextAlign(22);
    ttl.SetTextSize(0.038);
    ttl.DrawLatex(0.5 * (1.0 - mR + mL), 0.965,
        Form("Source %d - Geometric Efficiency (Main Wall)", source_number));

    // Інформаційна панель
    TPaveText* pvInfo = new TPaveText(1.0 - mR + 0.01, 0.70, 0.995, 0.92, "NDC");
    pvInfo->SetFillColor(kWhite);
    pvInfo->SetBorderSize(1);
    pvInfo->SetTextAlign(12);
    pvInfo->SetTextSize(0.025);
    
    double total_eff = 0.0;
    for (int om = 0; om < 260; om++) total_eff += efficiencies[om];
    
    pvInfo->AddText(Form("Total #varepsilon = %.2f%%", total_eff));
    pvInfo->AddText(Form("Min = %.2e %s", eps_min, 
        dispMode == DisplayMode::PERCENT ? "%" : 
        dispMode == DisplayMode::LOG10 ? "" : "ppm"));
    pvInfo->AddText(Form("Max = %.2e %s", eps_max, 
        dispMode == DisplayMode::PERCENT ? "%" : 
        dispMode == DisplayMode::LOG10 ? "" : "ppm"));
    
    const char* modeStr = (dispMode == DisplayMode::PERCENT) ? "mode: %"
                        : (dispMode == DisplayMode::LOG10)    ? "mode: log_{10}"
                                                              : "mode: ppm";
    pvInfo->AddText(modeStr);
    pvInfo->Draw("same");

    c->SaveAs(Form("plots/efficiencies/source_%02d_%s.png",
                   source_number, modeSuffix(dispMode)));
    std::cout << "Saved: source_" << std::setfill('0') << std::setw(2) << source_number 
              << "_" << modeSuffix(dispMode) << ".png\n";
    
    delete hFrame;
    delete c;
}


void total_eff_visualization(){
    bool total_vis = true;

     std::vector<SourcePosition> sourcePositions = read_source_positions("source_positions.txt");
    
    // Якщо total_vis = false, будуть виведені карти для кожного джерела
    // Якщо total_vis = true, буде одна карта з сумою всіх джерел
    
    //activities of sources as measured by Miro on July 1 2018
    double activities[42] = {129.2, 127.9, 138.5, 125.8, 134.3, 128.4, 
                            129.5, 131.0, 125.0, 129.8, 124.7, 132.0, 
                            133.8, 126.6, 131.3, 130.1, 133.0, 122.8, 
                            122.2, 139.1, 130.9, 140.9, 119.8, 135.0, 
                            132.5, 127.4, 131.9, 130.0, 132.6, 128.4, 
                            123.7, 131.1, 125.5, 130.7, 124.9, 132.2, 
                            131.4, 128.0, 137.3, 128.3, 136.3, 127.1};

    gROOT->SetBatch(true);
    gStyle->SetPalette(kBird);
    
    gSystem->mkdir("plots", true);
    gSystem->mkdir("plots/efficiencies", true);

    const int nCols = 20;
    const int nRows = 13;
    const int colMin = 0;
    const int colMax = nCols - 1;
    const int rowMin = 0;
    const int rowMax = nRows - 1;

    TFile* rootFile = new TFile("plots/efficiencies/geometric_efficiencies.root", "RECREATE");

    // РЕЖИМ 1
    if(total_vis){
        std::vector<double> total_eff(260, 0.0);
        
        std::cout << "Computing total efficiency (all sources)...\n";
        
        for (int i = 0; i < 260; ++i) {
            for (int j = 0; j < 42; ++j) {
                total_eff[i] += activities[j] * geometricEfficiency_OMS(i, j, sourcePositions) * 100.0;
            }
        }
        
        double eps_min = 1e30, eps_max = -1e30;
        for (int i = 0; i < 260; i++) {
            if (total_eff[i] > 0.0) {
                if (total_eff[i] < eps_min) eps_min = total_eff[i];
                if (total_eff[i] > eps_max) eps_max = total_eff[i];
            }
        }
        
        if (eps_min > eps_max) { eps_min = 0.0; eps_max = 1.0; }
        if (std::abs(eps_max - eps_min) < 1e-12) eps_max = eps_min + 1.0;

        rootFile->cd();
        TH2F* hTotalEff = new TH2F("total_efficiency", "Total Geometric Efficiency",
                                   nCols, colMin - 0.5, colMax + 0.5,
                                   nRows, rowMin - 0.5, rowMax + 0.5);
        hTotalEff->SetDirectory(rootFile);
        
        for (int om = 0; om < 260; om++) {
            int col = om / 13;
            int row = om % 13;
            hTotalEff->SetBinContent(col + 1, row + 1, total_eff[om]);
        }
        hTotalEff->Write();

        const double mL = 0.08, mR = 0.25, mT = 0.08, mB = 0.08;
        const int baseH = 600;
        const double areaH = baseH * (1.0 - mT - mB);
        const double areaW = areaH * (double)nCols / (double)nRows;
        const int canvasW = (int)(areaW / (1.0 - mL - mR)) + 1;
        const int canvasH = baseH;

        TCanvas* c = new TCanvas("cTotal", "", canvasW, canvasH);
        c->SetLeftMargin(mL);
        c->SetRightMargin(mR);
        c->SetTopMargin(mT);
        c->SetBottomMargin(mB);

        TH2F* hFrame = new TH2F("hFrame_total", "",
            nCols, colMin - 0.5, colMax + 0.5,
            nRows, rowMin - 0.5, rowMax + 0.5);
        hFrame->SetDirectory(0);
        hFrame->SetStats(0);
        hFrame->GetXaxis()->SetNdivisions(nCols, false);
        hFrame->GetYaxis()->SetNdivisions(nRows, false);
        hFrame->GetXaxis()->SetLabelSize(0.048);
        hFrame->GetYaxis()->SetLabelSize(0.048);
        
        for (int ic = 0; ic < nCols; ic++)
            hFrame->GetXaxis()->SetBinLabel(ic + 1, Form("%d", ic));
        for (int ir = 0; ir < nRows; ir++)
            hFrame->GetYaxis()->SetBinLabel(ir + 1, Form("%d", ir));
        
        hFrame->GetXaxis()->SetTitle("Column");
        hFrame->GetYaxis()->SetTitle("Row");
        hFrame->Draw("AXIS");

        const double cellW = 0.46, cellH = 0.46;
        const int nColors = TColor::GetNumberOfColors();

        for (int om = 0; om < 260; om++) {
            int y = om / 13;
            int z = om % 13;
            
            if (total_eff[om] <= 0.0) continue;

            double v = total_eff[om];
            double frac = (v - eps_min) / (eps_max - eps_min);
            frac = std::max(0.0, std::min(1.0, frac));

            int colorIdx = TColor::GetColorPalette((int)(frac * (nColors - 1)));

            double bx = colMin + y;
            double by = z;

            TBox* b = new TBox(bx - cellW, by - cellH, bx + cellW, by + cellH);
            b->SetFillColor(colorIdx);
            b->SetLineColor(kGray + 2);
            b->SetLineWidth(1);
            b->Draw("l same");
        }

        // Сітка
        for (int ic = 0; ic <= nCols; ic++) {
            double x = colMin - 0.5 + ic;
            TLine* l = new TLine(x, rowMin - 0.5, x, rowMax + 0.5);
            l->SetLineColor(kBlack); l->SetLineWidth(1); l->SetLineStyle(2);
            l->Draw("same");
        }
        for (int ir = 0; ir <= nRows; ir++) {
            double y = rowMin - 0.5 + ir;
            TLine* l = new TLine(colMin - 0.5, y, colMax + 0.5, y);
            l->SetLineColor(kBlack); l->SetLineWidth(1); l->SetLineStyle(2);
            l->Draw("same");
        }

        // Числові підписи
        for (int om = 0; om < 260; om++) {
            int y = om / 13;
            int z = om % 13;
            
            if (total_eff[om] <= 0.0) continue;

            double bx = colMin + y;
            double by = z;
            double v = total_eff[om];
            double frac = (v - eps_min) / (eps_max - eps_min);
            frac = std::max(0.0, std::min(1.0, frac));

            TLatex* lt = new TLatex(bx, by, Form("%.1f", v));
            lt->SetTextAlign(22);
            lt->SetTextSize(0.020);
            lt->SetTextColor(frac > 0.55 ? kBlack : kWhite);
            lt->Draw("same");
        }

        TLatex ttl;
        ttl.SetNDC();
        ttl.SetTextAlign(22);
        ttl.SetTextSize(0.038);
        ttl.DrawLatex(0.5, 0.965, "Total Geometric Efficiency (All Sources)");

        c->SaveAs("plots/efficiencies/total_efficiency.png");
        std::cout << "Saved: total_efficiency.png\n";
        
        delete hFrame;
        delete c;
    }

    // РЕЖИМ 2
    if(total_vis){
        std::cout << "\nGenerating efficiency maps for each source...\n";
        for (int src = 0; src < 42; src++) {
            std::cout << "Processing source " << src << "...\n";
            
            std::vector<double> efficiencies(260);
            double eps_min = 1e30, eps_max = -1e30;
            
            for (int om = 0; om < 260; om++) {
                efficiencies[om] = geometricEfficiency_OMS(om, src, sourcePositions) * 100.0;
                if (efficiencies[om] > 0.0) {
                    if (efficiencies[om] < eps_min) eps_min = efficiencies[om];
                    if (efficiencies[om] > eps_max) eps_max = efficiencies[om];
                }
            }
            
            if (eps_min > eps_max) { eps_min = 0.0; eps_max = 1.0; }
            if (std::abs(eps_max - eps_min) < 1e-12) eps_max = eps_min + 1.0;
            
            rootFile->cd();
            TH2F* hSourceEff = new TH2F(Form("source_%d", src), 
                                        Form("Source %d Geometric Efficiency", src),
                                        nCols, colMin - 0.5, colMax + 0.5,
                                        nRows, rowMin - 0.5, rowMax + 0.5);
            hSourceEff->SetDirectory(rootFile);
            
            for (int om = 0; om < 260; om++) {
                int col = om / 13;
                int row = om % 13;
                hSourceEff->SetBinContent(col + 1, row + 1, efficiencies[om]);
            }
            hSourceEff->Write();
            
            draw_single_source_efficiency(src, sourcePositions, DisplayMode::PERCENT);
        }
    }

    rootFile->Close();
    delete rootFile;

    std::cout << "\n✓ Visualization complete!\n";
    std::cout << "  PNG output saved to: plots/efficiencies/*.png\n";
    std::cout << "  ROOT output saved to: plots/efficiencies/geometric_efficiencies.root\n";
}


int main() {
    try {
        total_eff_visualization();
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
