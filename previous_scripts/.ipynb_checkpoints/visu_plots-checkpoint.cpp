#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TROOT.h>

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>

#include <TGraph.h>
#include <TBox.h>
#include <TLatex.h>
#include <TSystem.h>
#include <TStyle.h>


std::vector<std::pair<double,double>>
read_positions_from_file(const std::string& filename)
{
    std::vector<std::pair<double,double>> positions;

    std::ifstream infile(filename);

    if(!infile.is_open())
    {
        std::cerr<<"ERROR: cannot open file "<<filename<<std::endl;
        return positions;
    }

    std::string line;
    int auto_index = 0;

    while(std::getline(infile,line))
    {
        if(line.empty()) continue;

        double y,z;

        if(line.find(';') != std::string::npos)
        {
            std::replace(line.begin(), line.end(), ';', ' ');
            std::istringstream ss(line);

            if(!(ss >> y >> z))
                continue;

            positions.emplace_back(y,z);
            auto_index++;
        }

        else
        {
            std::istringstream ss(line);

            int index;

            if(!(ss >> index >> y >> z))
                continue;

            positions.emplace_back(y,z);
        }
    }

    infile.close();

    std::cout<<"Loaded "<<positions.size()
             <<" positions from "<<filename<<std::endl;

    return positions;
}


void draw_source_OM_distributions(
        TH2F* hSourceOM_pos,
        TH2F* hSourceOM_neg,
        const std::vector<std::pair<double,double>>& src_pos,
        const std::vector<std::pair<double,double>>& om_pos)
{

    gSystem->mkdir("plots/source_maps", true);

    const int nSources = src_pos.size();
    const int nOMs = om_pos.size();

    const double mw_sizey = 256.0;
    const double mw_sizez = 256.0;

    int px = 1000;
    int py = 800;

    double rl=0.12, rr=0.12, rt=0.08, rb=0.12;

    for(int side=0; side<2; side++)
    {

        TH2F* matrix = (side==0) ? hSourceOM_pos : hSourceOM_neg;

        for(int s=0; s<nSources; s++)
        {

            TCanvas* c_map = new TCanvas(
                Form("c_map_%d_%d",side,s),
                "",
                px,py);

            c_map->SetLeftMargin(rl);
            c_map->SetRightMargin(rr);
            c_map->SetTopMargin(rt);
            c_map->SetBottomMargin(rb);

            c_map->SetFillColor(0);
            c_map->SetFrameBorderMode(0);
            c_map->SetBorderMode(0);
            c_map->SetFrameFillColor(0);

            TH2D* dummy = new TH2D(
                Form("dummy_%d_%d",side,s),
                Form("Source %d side %s; y [mm]; z [mm]",
                s,(side==0?"X>0":"X<0")),
                200,-2500,2500,
                200,-2000,2000);

            dummy->SetDirectory(0);
            dummy->SetStats(0);
            dummy->Draw();

            dummy->GetXaxis()->SetLabelSize(0);
            dummy->GetYaxis()->SetLabelSize(0);
            dummy->GetXaxis()->SetTitle("");
            dummy->GetYaxis()->SetTitle("");
            dummy->GetXaxis()->SetTickLength(0);
            dummy->GetYaxis()->SetTickLength(0);

            double maxHits = 0;

            for(int om=0; om<nOMs; om++)
            {
                double hits = matrix->GetBinContent(s+1,om+1);
                if(hits > maxHits) maxHits = hits;
            }

            if(maxHits == 0) maxHits = 1;

            TH2D* paletteHist = new TH2D(
                Form("palette_%d_%d",side,s),
                "",
                1,0,1,
                100,0,maxHits);

            paletteHist->SetDirectory(0);
            paletteHist->SetStats(0);
            paletteHist->Draw("COLZ same");


            for(int om=0; om<nOMs; om++)
            {

                double y = om_pos[om].first;
                double z = om_pos[om].second;

                double hits = matrix->GetBinContent(s+1,om+1);

                int color = kWhite;

                if(hits>0)
                {
                    double frac = hits/maxHits;
                    color = gStyle->GetColorPalette(int(frac*50));
                }

                TBox* box = new TBox(
                    y-mw_sizey/2,
                    z-mw_sizez/2,
                    y+mw_sizey/2,
                    z+mw_sizez/2);

                box->SetFillColor(color);
                box->SetLineColor(kRed);
                box->Draw("same");

                TLatex* label = new TLatex(
                    y-90,
                    z-90,
                    Form("%d",om));

                label->SetTextSize(0.02);
                label->SetTextColor(kRed+1);
                label->Draw("same");
            }

            for(int i=0;i<nSources;i++)
            {

                double y = src_pos[i].first;
                double z = src_pos[i].second;

                int color = kGreen+2;
                int size = 1.0;

                if(i == s)
                {
                    color = kBlue+1;
                    size = 2.0;
                }

                TGraph* src_g = new TGraph(1,&y,&z);

                src_g->SetMarkerStyle(20);
                src_g->SetMarkerSize(size);
                src_g->SetMarkerColor(color);

                src_g->Draw("P same");

                TLatex* label = new TLatex(
                    y+40,
                    z+40,
                    Form("%d",i));

                label->SetTextSize(0.025);
                label->SetTextColor(color);

                label->Draw("same");
            }

            TString sideName = (side==0) ? "pos" : "neg";

            c_map->SaveAs(
                Form("plots/source_maps/source_%s_%02d.png",
                sideName.Data(),s));

            delete c_map;
        }
    }
}


void visu_plots()
{

gROOT->SetBatch(true);

gSystem->mkdir("plots", true);
gSystem->mkdir("plots/source_maps", true);

TFile* fin = new TFile("vertex_energy_results.root","READ");

if(!fin || fin->IsZombie())
{
    std::cerr<<"ERROR: cannot open vertex_energy_results.root"<<std::endl;
    return;
}

TH1F* hEnergy      = (TH1F*)fin->Get("Energy");
TH2F* hStart_Xpos  = (TH2F*)fin->Get("hStart_Xpos");
TH2F* hStart_Xneg  = (TH2F*)fin->Get("hStart_Xneg");
TH2F* hEnd_Xpos    = (TH2F*)fin->Get("hEnd_Xpos");
TH2F* hEnd_Xneg    = (TH2F*)fin->Get("hEnd_Xneg");
TH2F* hSourceOM_pos = (TH2F*)fin->Get("SourceOMpos");
TH2F* hSourceOM_neg = (TH2F*)fin->Get("SourceOMneg");

if(!hEnergy || !hStart_Xpos || !hStart_Xneg ||
   !hEnd_Xpos || !hEnd_Xneg || !hSourceOM_pos || !hSourceOM_neg)
{
    std::cerr<<"ERROR: one or more histograms missing from ROOT file"<<std::endl;
    return;
}

// Detach from file so they survive after fin is closed
hEnergy->SetDirectory(0);
hStart_Xpos->SetDirectory(0);
hStart_Xneg->SetDirectory(0);
hEnd_Xpos->SetDirectory(0);
hEnd_Xneg->SetDirectory(0);
hSourceOM_pos->SetDirectory(0);
hSourceOM_neg->SetDirectory(0);

fin->Close();


TCanvas c1;
hStart_Xpos->Draw("COLZ");
c1.SaveAs("plots/start_Xpos.png");

TCanvas c2;
hStart_Xneg->Draw("COLZ");
c2.SaveAs("plots/start_Xneg.png");

TCanvas c3;
hEnd_Xpos->Draw("COLZ");
c3.SaveAs("plots/end_Xpos.png");

TCanvas c4;
hEnd_Xneg->Draw("COLZ");
c4.SaveAs("plots/end_Xneg.png");

TCanvas c5;
hEnergy->Draw("SAME");
c5.SaveAs("plots/energy_spectrum.png");

TCanvas c6;
hSourceOM_pos->Draw("COLZ");
c6.SaveAs("plots/source_OM_matrix_fr.png");   // fr is positive

TCanvas c7;
hSourceOM_neg->Draw("COLZ");
c7.SaveAs("plots/source_OM_matrix_ital.png"); // ital is negative


std::vector<std::pair<double,double>> source_positions =
    read_positions_from_file("default_values/source_positions.txt");

std::vector<std::pair<double,double>> om_positions =
    read_positions_from_file("default_values/om_positions.txt");

draw_source_OM_distributions(
    hSourceOM_pos,
    hSourceOM_neg,
    source_positions,
    om_positions);

std::cout<<"Finished"<<std::endl;

}
