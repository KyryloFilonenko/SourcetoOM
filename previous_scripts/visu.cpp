#include <TFile.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TH2D.h>
#include <TGraph.h>
#include <TBox.h>
#include <TLatex.h>
#include <TLegend.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <string>

std::vector<std::pair<double, double>> read_positions_from_file(const std::string& filename)
{
	std::vector<std::pair<double, double>> positions;
	std::ifstream infile(filename);
	std::string line;

	while (std::getline(infile, line)) 
	{
        	std::istringstream ss(line);
        	std::string index_str;
        	double x, y;

        	if (!(ss >> index_str >> x >> y)) 
        	{
                	continue;
                }

        	positions.emplace_back(x, y);
	}

	return positions;
}

std::vector<double> read_activities(const std::string& filename)
{
	std::vector<double> activities;
	std::ifstream infile(filename);
	std::string line;

	while (std::getline(infile, line)) 
	{
        	std::istringstream ss(line);
        	std::string token;
        	double activity;
        	std::getline(ss, token, ';');
        	
        	if (std::getline(ss, token)) 
        	{
                	activity = std::stod(token);
                	activities.push_back(activity);
                }
	}

	return activities;
}

	const int nSources = 42;
	const int nOMs = 260;
	const int nbinsx = 20;
 	const int nbinsy = 13;
	const double mw_sizey = 256.0;
	const double mw_sizez = 256.0;

void visu()
{
	double rl = 0.09, rr = 0.35, rt = 0.1, rb = 0.1;
	double px = 1000;
	double py = (px / (double(nbinsx)/nbinsy)) * (1 - rl - rr) / (1 - rb - rt);

	TFile* file = TFile::Open("merged_week_2025_X-BETE.root");
	if (!file || file->IsZombie()) 
	{
        	std::cerr << "Error opening merged_week_2025_X-BETE.root\n";
        	return;
        }

	std::vector<std::pair<double, double>> src_pos = read_positions_from_file("source_positions.txt");	
	std::vector<std::pair<double, double>> om_pos = read_positions_from_file("om_positions.txt");
	std::vector<double> activities = read_activities("source_activity.txt");

	if (src_pos.size() != nSources || om_pos.size() != nOMs) 
	{
        	std::cerr << "Mismatch in expected source or OM count!\n";
        	return;
        }

	if (activities.size() != nSources) 
        {
        	std::cerr << "Mismatch in activities count!\n";
        	return;
        }

	gROOT->SetBatch(kTRUE);  // No GUI

	// === Plot each individual eps_G[i] map ===
	for (int i = 0; i < nSources; ++i) 
	{
        	TH2D* hist = (TH2D*)file->Get(Form("hist%d", i));
        	if (!hist) 
        	{
                	std::cerr << "Missing histogram: hist" << i << "\n";
                	continue;
        	}

        	TCanvas* c = new TCanvas(Form("c_src_%02d", i), Form("Source %d", i), px, py);
        	c->SetLeftMargin(rl);
       		c->SetTopMargin(rt);
        	c->SetBottomMargin(rb);
        	c->SetRightMargin(rr);
        	c->cd();

        	hist->GetXaxis()->SetTitle("y [mm]");
        	hist->GetYaxis()->SetTitle("z [mm]");
        	hist->Draw("COLZ");

        	double x = src_pos[i].first;
        	double y = src_pos[i].second;
        	TGraph* src_graph = new TGraph(1, &x, &y);
        	src_graph->SetMarkerStyle(20);
        	src_graph->SetMarkerSize(0.6);
        	src_graph->SetMarkerColor(kRed);
        	src_graph->Draw("P same");

        	c->SaveAs(Form("Darya/eps_G_source_%02d.png", i));
        	delete c;
	}

	// === Plot source & OM positions with labels ===
	rr=rr/3;
	TCanvas* c_map = new TCanvas("c_map", "Source & OM positions", px, py);
	c_map->SetLeftMargin(rl);
	c_map->SetTopMargin(rt);
	c_map->SetBottomMargin(rb);
	c_map->SetRightMargin(rr);
	c_map->cd();
	TH2D* dummy = new TH2D("dummy", "Source and OM Positions", nbinsx, -double(nbinsx)/2*mw_sizey, double(nbinsx)/2*mw_sizey, nbinsy, -double(nbinsy)/2*mw_sizey, double(nbinsy)/2*mw_sizey);
	dummy->GetXaxis()->SetTitle("y [mm]");
	dummy->GetYaxis()->SetTitle("z [mm]");
	dummy->Draw();
	dummy->SetStats(0);


	// Sources: green markers with labels
	for (int i = 0; i < nSources; ++i) 
	{
        	double x = src_pos[i].first;
        	double y = src_pos[i].second;

        	TGraph* src_g = new TGraph(1, &x, &y);
        	src_g->SetMarkerStyle(20);
        	src_g->SetMarkerSize(1.0);
        	src_g->SetMarkerColor(kGreen + 2);
        	src_g->Draw("P same");

		TLatex* label = new TLatex(x + 50, y + 50, Form("%d", i));
        	label->SetTextSize(0.025);
        	label->SetTextColor(kGreen + 3);
        	label->Draw("same");
        }

	// OMs: red boxes with labels
	for (int i = 0; i < nOMs; ++i) 
{
    double x = om_pos[i].first;
    double y = om_pos[i].second;

    TBox* box = new TBox(
        x - mw_sizey / 2,
        y - mw_sizez / 2,
        x + mw_sizey / 2,
        y + mw_sizez / 2
    );
    box->SetLineColor(kRed);
    box->SetFillStyle(0);
    box->Draw("same");

    TLatex* label = new TLatex(x - 100, y - 100, Form("%d", i));
    label->SetTextSize(0.025);
    label->SetTextColor(kRed + 1);
    label->Draw("same");
}


	c_map->SaveAs("Darya/source_OM_map.png");
	delete c_map;
	
	rr=rr*3;

	// === Weighted sum of all eps_G[i] histograms by activity ===
	TH2D* weighted_sum = nullptr;

	for (int i = 0; i < nSources; ++i) 
	{
        	TH2D* hist = (TH2D*)file->Get(Form("hist%d", i));
        	if (!hist) 
        		continue;

        	if (!weighted_sum)
        		weighted_sum = (TH2D*)hist->Clone("weighted_sum");
        	else
        		weighted_sum->Add(hist, activities[i]);
	}

	if (weighted_sum) 
	{
        	TCanvas* c_sum = new TCanvas("c_sum", "Weighted Sum eps_G", px, py);
        	c_sum->SetLeftMargin(rl);
        	c_sum->SetTopMargin(rt);
        	c_sum->SetBottomMargin(rb);
        	c_sum->SetRightMargin(rr);
        	c_sum->cd();

        	weighted_sum->GetXaxis()->SetTitle("y [mm]");
        	weighted_sum->GetYaxis()->SetTitle("z [mm]");
        	weighted_sum->SetTitle("eps_G summed OM by OM and weighted by the activities of the sources");
        	weighted_sum->Draw("COLZ");

        	c_sum->SaveAs("Darya/eps_G_weighted_sum.png");
        	delete c_sum;
	}
	weighted_sum->SetStats(false);        // Removes stats box (not legend)


	file->Close();
}