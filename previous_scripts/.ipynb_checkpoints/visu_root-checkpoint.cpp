// -*- coding: utf-8 -*-
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TBox.h>
#include <TLatex.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TROOT.h>

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <string>


// ============================================================
//  Читання позицій з файлу
// ============================================================

std::vector<std::pair<double,double>>
read_positions_from_file(const std::string& filename)
{
    std::vector<std::pair<double,double>> positions;
    std::ifstream infile(filename);

    if(!infile.is_open())
    {
        std::cerr << "ERROR: cannot open file " << filename << std::endl;
        return positions;
    }

    std::string line;
    while(std::getline(infile, line))
    {
        if(line.empty()) continue;

        double y, z;

        if(line.find(';') != std::string::npos)
        {
            std::replace(line.begin(), line.end(), ';', ' ');
            std::istringstream ss(line);
            if(!(ss >> y >> z)) continue;
        }
        else
        {
            std::istringstream ss(line);
            int index;
            if(!(ss >> index >> y >> z)) continue;
        }

        positions.emplace_back(y, z);
    }

    infile.close();
    std::cout << "Loaded " << positions.size()
              << " positions from " << filename << std::endl;
    return positions;
}


// ============================================================
//  Карти точок початку для кожного джерела окремо
//  Область відображення: ±zoom мм навколо центру джерела
// ============================================================

void draw_source_start_vertices(
    const std::vector<TH2F*>& hStart_pos,
    const std::vector<TH2F*>& hStart_neg,
    const std::vector<std::pair<double,double>>& src_pos)
{
    gSystem->mkdir("plots/bin_sources", true);

    const int nSources = (int)src_pos.size();

    for(int s = 0; s < nSources; s++)
    {
        double yc = src_pos[s].first;
        double zc = src_pos[s].second;

        for(int side = 0; side < 2; side++)
        {
            TH2F* h = (side == 0) ? hStart_pos[s] : hStart_neg[s];

            TCanvas* c = new TCanvas(
                Form("c_sv_%d_%d", side, s), "", 900, 800);

            c->SetLeftMargin(0.12);
            c->SetRightMargin(0.15);
            c->SetTopMargin(0.09);
            c->SetBottomMargin(0.11);

            h->SetStats(0);
            h->Draw("COLZ");

            // Червона точка — істинний центр джерела
            TGraph* gSrc = new TGraph(1, &yc, &zc);
            gSrc->SetMarkerStyle(20);
            gSrc->SetMarkerSize(1.8);
            gSrc->SetMarkerColor(kRed);
            gSrc->Draw("P same");

            // Підпис джерела
            double zoom = (h->GetXaxis()->GetXmax() - h->GetXaxis()->GetXmin()) / 2.0;
            TLatex* lbl = new TLatex(
                yc + zoom * 0.08,
                zc + zoom * 0.08,
                Form("src %d  (%.0f, %.0f)", s, yc, zc));
            lbl->SetTextSize(0.032);
            lbl->SetTextColor(kRed);
            lbl->Draw("same");

            TString sideName = (side == 0) ? "pos" : "neg";
            c->SaveAs(Form("plots/bin_sources/start_%s_src%02d.png",
                           sideName.Data(), s));

            delete c;
        }
    }
}



// ============================================================
//  Малювання карт джерело -> OM
//  Легенда праворуч: до 10 кольорових ступенів з підписом кількості подій
// ============================================================

void draw_source_OM_distributions(
    TH2F* hSourceOM_pos,
    TH2F* hSourceOM_neg,
    const std::vector<std::pair<double,double>>& src_pos,
    const std::vector<std::pair<double,double>>& om_pos)
{
    gSystem->mkdir("plots/source_maps", true);

    const int nSources = (int)src_pos.size();
    const int nOMs     = (int)om_pos.size();

    const double mw_sizey = 256.0;
    const double mw_sizez = 256.0;

    // 10 фіксованих кольорів для ступенів (від холодного до теплого)
    // індекс 0 — 0 хітів (білий), 1..10 — ступені від 1 до maxHits
    // const int    nSteps   = 10;
    // const int    stepColors[nSteps] = {
    //     kCyan - 9,    // дуже світло-блакитний  (найменше)
    //     kCyan,
    //     kAzure - 4,
    //     kGreen,
    //     kGreen + 2,
    //     kYellow,
    //     kOrange,
    //     kOrange + 7,
    //     kRed,
    //     kRed + 2      // темно-червоний (найбільше)
    // };

    const int    nSteps   = 20;
    const int    stepColors[nSteps] = {
        kCyan - 9,    // дуже світло-блакитний  (найменше)
        kCyan - 5,
        kCyan - 3,
        kCyan,
        kAzure - 4,
        kAzure -2,
        kGreen,
        kGreen + 1,
        kGreen + 2,
        kGreen + 3,
        kYellow,
        kYellow + 2,
        kYellow +5,
        kOrange,
        kOrange + 3,
        kOrange + 5,
        kOrange + 7,
        kRed,
        kRed + 2,      // темно-червоний (найбільше)
        kRed + 5,
    };

    int    px = 1200, py = 800;
    // правий відступ збільшений для легенди
    double rl = 0.07, rr = 0.20, rt = 0.08, rb = 0.05;

    for(int side = 0; side < 2; side++)
    {
        TH2F* matrix = (side == 0) ? hSourceOM_pos : hSourceOM_neg;

        for(int s = 0; s < nSources; s++)
        {
            TCanvas* c_map = new TCanvas(
                Form("c_map_%d_%d", side, s), "", px, py);

            c_map->SetLeftMargin(rl);
            c_map->SetRightMargin(rr);
            c_map->SetTopMargin(rt);
            c_map->SetBottomMargin(rb);
            c_map->SetFillColor(0);
            c_map->SetFrameBorderMode(0);
            c_map->SetBorderMode(0);
            c_map->SetFrameFillColor(0);

            // --- Порожня рамка — осі приховані ---
            // X розширено до 3600 щоб легенда праворуч поміщалась у frame
            TH2D* dummy = new TH2D(
                Form("dummy_%d_%d", side, s),
                Form("Source %d side %s",
                     s, (side == 0 ? "X>0" : "X<0")),
                200, -2500, 3600,
                200, -2000, 2000);

            dummy->SetDirectory(0);
            dummy->SetStats(0);
            dummy->GetXaxis()->SetLabelSize(0);
            dummy->GetYaxis()->SetLabelSize(0);
            dummy->GetXaxis()->SetTitle("");
            dummy->GetYaxis()->SetTitle("");
            dummy->GetXaxis()->SetTickLength(0);
            dummy->GetYaxis()->SetTickLength(0);
            dummy->Draw("AXIS");

            // --- Максимум хітів по всіх OM для цього джерела ---
            double maxHits = 0;
            for(int om = 0; om < nOMs; om++)
            {
                double hits = matrix->GetBinContent(s + 1, om + 1);
                if(hits > maxHits) maxHits = hits;
            }
            if(maxHits == 0) maxHits = 1;

            // --- Функція: кількість хітів -> індекс ступеня 0..nSteps-1 ---
            // ступінь 0 відповідає 1 хіту (найменше), nSteps-1 = maxHits
            auto hitsToStep = [&](double hits) -> int {
                if(hits <= 0) return -1;  // білий (порожній)
                int step = (int)( (hits - 1) / maxHits * nSteps );
                if(step >= nSteps) step = nSteps - 1;
                return step;
            };

            // --- Заливка OM ---
            for(int om = 0; om < nOMs; om++)
            {
                double y    = om_pos[om].first;
                double z    = om_pos[om].second;
                double hits = matrix->GetBinContent(s + 1, om + 1);

                int step = hitsToStep(hits);
                int fillColor = (step < 0) ? kWhite : stepColors[step];

                TBox* fill = new TBox(
                    y - mw_sizey / 2, z - mw_sizez / 2,
                    y + mw_sizey / 2, z + mw_sizez / 2);
                fill->SetFillColor(fillColor);
                fill->SetLineStyle(0);
                fill->SetLineWidth(0);
                fill->Draw("f same");
            }

            // --- Тонкі червоні границі OM ---
            for(int om = 0; om < nOMs; om++)
            {
                double y = om_pos[om].first;
                double z = om_pos[om].second;

                TBox* border = new TBox(
                    y - mw_sizey / 2, z - mw_sizez / 2,
                    y + mw_sizey / 2, z + mw_sizez / 2);
                border->SetFillStyle(0);
                border->SetLineColor(kRed);
                border->SetLineWidth(1);
                border->Draw("same");
            }

            // --- Номери OM ---
            for(int om = 0; om < nOMs; om++)
            {
                double y = om_pos[om].first;
                double z = om_pos[om].second;

                TLatex* lom = new TLatex(y - 90, z - 90, Form("%d", om));
                lom->SetTextSize(0.018);
                lom->SetTextColor(kRed + 1);
                lom->Draw("same");
            }

            // --- Точки джерел ---
            for(int i = 0; i < nSources; i++)
            {
                double yy = src_pos[i].first;
                double zz = src_pos[i].second;

                int    color = kGreen + 2;
                double size  = 1.0;
                if(i == s){ color = kBlue + 1; size = 2.0; }

                TGraph* src_g = new TGraph(1, &yy, &zz);
                src_g->SetMarkerStyle(20);
                src_g->SetMarkerSize(size);
                src_g->SetMarkerColor(color);
                src_g->Draw("P same");

                TLatex* lsrc = new TLatex(yy + 40, zz + 40, Form("%d", i));
                lsrc->SetTextSize(0.022);
                lsrc->SetTextColor(color);
                lsrc->Draw("same");
            }

            // --------------------------------------------------------
            //  Ручна легенда праворуч від поля OM
            //  TBox не підтримує NDC — конвертуємо в world coordinates
            //  Область осей: X в [-2500, 2500], Z в [-2000, 2000]
            // --------------------------------------------------------
            // Параметри легенди у world coordinates
            // X: одразу праворуч від правої границі OM-поля
            const double legWX0  =  2600.0;   // ліво колірного квадрата
            const double legWX1  =  3200.0;   // право колірного квадрата
            const double legWtxt =  3280.0;   // x підпису
            const double legWtop =  1900.0;   // верх легенди
            const double legWbot = -1900.0;   // низ легенди
            const int    nLegRows = nSteps + 1;  // +1 для "0 подій"
            const double rowH_w  = (legWtop - legWbot) / nLegRows;

            // Заголовок "Events"
            TLatex* legTitle = new TLatex(
                (legWX0 + legWX1) / 2.0, legWtop + 60.0, "Events");
            legTitle->SetTextAlign(21);
            legTitle->SetTextSize(0.028);
            legTitle->SetTextColor(kBlack);
            legTitle->Draw("same");

            // Рядок 0: білий = 0 подій (знизу)
            {
                double wy0 = legWbot;
                double wy1 = legWbot + rowH_w;

                TBox* bleg = new TBox(legWX0, wy0, legWX1, wy1);
                bleg->SetFillColor(kWhite);
                bleg->SetLineColor(kGray + 1);
                bleg->SetLineWidth(1);
                bleg->Draw("same");

                TLatex* tleg = new TLatex(legWtxt, wy0 + rowH_w * 0.25, "0");
                tleg->SetTextSize(0.022);
                tleg->SetTextColor(kBlack);
                tleg->Draw("same");
            }

            // Рядки 1..nSteps: знизу (мало) -> вгору (багато)
            for(int st = 0; st < nSteps; st++)
            {
                // st=0 найменше — рядок одразу над "0", st=nSteps-1 — верх
                double wy0 = legWbot + (st + 1) * rowH_w;
                double wy1 = legWbot + (st + 2) * rowH_w;

                TBox* bleg = new TBox(legWX0, wy0, legWX1, wy1);
                bleg->SetFillColor(stepColors[st]);
                bleg->SetLineColor(kGray + 1);
                bleg->SetLineWidth(1);
                bleg->Draw("same");

                // Діапазон хітів для цього ступеня
                double loHits = 1.0 + (double)st       * maxHits / nSteps;
                double hiHits =       (double)(st + 1) * maxHits / nSteps;
                if(st == nSteps - 1) hiHits = maxHits;

                TString label;
                if((int)loHits == (int)hiHits)
                    label = Form("%d", (int)loHits);
                else
                    label = Form("%d-%d", (int)loHits, (int)hiHits);

                TLatex* tleg = new TLatex(
                    legWtxt, wy0 + rowH_w * 0.25, label.Data());
                tleg->SetTextSize(0.022);
                tleg->SetTextColor(kBlack);
                tleg->Draw("same");
            }

            TString sideName = (side == 0) ? "pos" : "neg";
            c_map->SaveAs(
                Form("plots/source_maps/source_%s_%02d.png",
                     sideName.Data(), s));

            delete c_map;
        }
    }
}


// ============================================================
//  Головна функція
// ============================================================

void visu_root()
{
    gROOT->SetBatch(true);
    gSystem->mkdir("plots", true);

    // --- Відкриваємо оброблений файл ---
    TFile* f = TFile::Open("processed_data.root", "READ");
    if(!f || f->IsZombie())
    {
        std::cerr << "ERROR: cannot open processed_data.root!" << std::endl;
        return;
    }

    TTree* evTree = (TTree*)f->Get("EventTree");
    if(!evTree)
    {
        std::cerr << "ERROR: cannot find EventTree!" << std::endl;
        return;
    }

    // --- Підключаємо гілки ---
    int    event_id, source_id, om_id, side;
    double x_start, y_start, z_start;
    double x_end,   y_end,   z_end;
    double energy;

    evTree->SetBranchAddress("event_id",  &event_id);
    evTree->SetBranchAddress("source_id", &source_id);
    evTree->SetBranchAddress("om_id",     &om_id);
    evTree->SetBranchAddress("side",      &side);
    evTree->SetBranchAddress("x_start",   &x_start);
    evTree->SetBranchAddress("y_start",   &y_start);
    evTree->SetBranchAddress("z_start",   &z_start);
    evTree->SetBranchAddress("x_end",     &x_end);
    evTree->SetBranchAddress("y_end",     &y_end);
    evTree->SetBranchAddress("z_end",     &z_end);
    evTree->SetBranchAddress("energy",    &energy);

    // --- Позиції джерел і OM ---
    std::vector<std::pair<double,double>> source_positions =
        read_positions_from_file("default_values/source_positions.txt");

    std::vector<std::pair<double,double>> om_positions =
        read_positions_from_file("default_values/om_positions.txt");

    const int nSources = 42;
    const int nOM      = 260;

    // --- Загальні гістограми ---
    TH1F* hEnergy = new TH1F(
        "Energy", "Energy distribution;Energy [MeV];Counts",
        2000, 0, 3000);
    hEnergy->SetDirectory(0);

    TH2F* hStart_Xpos = new TH2F("hStart_Xpos", "", 500, -2500, 2500, 500, -2000, 2000);
    TH2F* hStart_Xneg = new TH2F("hStart_Xneg", "", 500, -2500, 2500, 500, -2000, 2000);
    TH2F* hEnd_Xpos   = new TH2F("hEnd_Xpos",   "", 500, -2500, 2500, 500, -2000, 2000);
    TH2F* hEnd_Xneg   = new TH2F("hEnd_Xneg",   "", 500, -2500, 2500, 500, -2000, 2000);

    hStart_Xpos->SetDirectory(0);
    hStart_Xneg->SetDirectory(0);
    hEnd_Xpos->SetDirectory(0);
    hEnd_Xneg->SetDirectory(0);

    TH2F* hSourceOM_pos = new TH2F(
        "SourceOMpos", "Source -> OM matrix;Source ID;OM ID",
        nSources, 0, nSources, nOM, 0, nOM);
    TH2F* hSourceOM_neg = new TH2F(
        "SourceOMneg", "Source -> OM matrix;Source ID;OM ID",
        nSources, 0, nSources, nOM, 0, nOM);

    hSourceOM_pos->SetDirectory(0);
    hSourceOM_neg->SetDirectory(0);

    const double zoom = 25.0;
    const int    nBinsZoom = 50; // 1 бін = 1 мм
    
    std::vector<TH2F*> hSrcStart_pos(nSources);
    std::vector<TH2F*> hSrcStart_neg(nSources);
    
    for(int s = 0; s < nSources; s++)
    {
        double yc = source_positions[s].first;
        double zc = source_positions[s].second;
    
        hSrcStart_pos[s] = new TH2F(
            Form("hSrcStart_pos_%02d", s),
            Form("Start vertices src %d (X>0);Y [mm];Z [mm]", s),
            nBinsZoom, yc - zoom, yc + zoom,
            nBinsZoom, zc - zoom, zc + zoom);
        hSrcStart_pos[s]->SetDirectory(0);
    
        hSrcStart_neg[s] = new TH2F(
            Form("hSrcStart_neg_%02d", s),
            Form("Start vertices src %d (X<0);Y [mm];Z [mm]", s),
            nBinsZoom, yc - zoom, yc + zoom,
            nBinsZoom, zc - zoom, zc + zoom);
        hSrcStart_neg[s]->SetDirectory(0);
    }


    
    // --- Запам'ятовуємо X для заголовків ---
    double x_pos_value = 0, x_neg_value = 0;
    bool   pos_found   = false, neg_found = false;

    // --- Головний цикл ---
    Long64_t nentries = evTree->GetEntries();
    std::cout << "Reading " << nentries << " entries..." << std::endl;

    for(Long64_t ie = 0; ie < nentries; ie++)
    {
        evTree->GetEntry(ie);

        // Start vertex — загальна гістограма
        if(side == 1)
        {
            hStart_Xpos->Fill(y_start, z_start);
            if(!pos_found){ x_pos_value = x_start; pos_found = true; }
        }
        else
        {
            hStart_Xneg->Fill(y_start, z_start);
            if(!neg_found){ x_neg_value = x_start; neg_found = true; }
        }

        // Start vertex — гістограма для конкретного джерела
        if(source_id >= 0 && source_id < nSources)
        {
            if(side == 1) hSrcStart_pos[source_id]->Fill(y_start, z_start);
            else          hSrcStart_neg[source_id]->Fill(y_start, z_start);
        }

        // End vertex
        if(x_end > 0) hEnd_Xpos->Fill(y_end, z_end);
        else          hEnd_Xneg->Fill(y_end, z_end);

        // Енергія
        if(energy > 0) hEnergy->Fill(energy);

        // Матриця джерело -> OM
        if(source_id >= 0 && om_id >= 0)
        {
            if(side == 1) hSourceOM_pos->Fill(source_id, om_id);
            else          hSourceOM_neg->Fill(source_id, om_id);
        }
    }

    f->Close();

    // --- Оновлюємо заголовки загальних гістограм ---
    hStart_Xpos->SetTitle(Form("Start vertices (X = %.1f mm);Y [mm];Z [mm]", x_pos_value));
    hStart_Xneg->SetTitle(Form("Start vertices (X = %.1f mm);Y [mm];Z [mm]", x_neg_value));
    hEnd_Xpos->SetTitle(  Form("End vertices (X #approx %.1f mm);Y [mm];Z [mm]", x_pos_value));
    hEnd_Xneg->SetTitle(  Form("End vertices (X #approx %.1f mm);Y [mm];Z [mm]", x_neg_value));

    // --- Малюємо загальні PNG ---
    // TCanvas c1; hStart_Xpos->Draw("COLZ"); c1.SaveAs("plots/start_Xpos.png");
    // TCanvas c2; hStart_Xneg->Draw("COLZ"); c2.SaveAs("plots/start_Xneg.png");
    // TCanvas c3; hEnd_Xpos->Draw("COLZ");   c3.SaveAs("plots/end_Xpos.png");
    // TCanvas c4; hEnd_Xneg->Draw("COLZ");   c4.SaveAs("plots/end_Xneg.png");
    TCanvas c5; hEnergy->Draw("PE");            c5.SaveAs("plots/energy_spectrum.png");
    // TCanvas c6; hSourceOM_pos->Draw("COLZ"); c6.SaveAs("plots/source_OM_matrix_fr.png");
    // TCanvas c7; hSourceOM_neg->Draw("COLZ"); c7.SaveAs("plots/source_OM_matrix_ital.png");

    // --- Зберігаємо гістограми у ROOT файл ---
    TFile* fout = new TFile("vertex_energy_results.root", "RECREATE");
    hStart_Xpos->Write();
    hStart_Xneg->Write();
    hEnd_Xpos->Write();
    hEnd_Xneg->Write();
    hEnergy->Write();
    hSourceOM_pos->Write();
    hSourceOM_neg->Write();
    for(int s = 0; s < nSources; s++)
    {
        hSrcStart_pos[s]->Write();
        hSrcStart_neg[s]->Write();
    }
    fout->Close();

    // // --- Детальні карти джерело -> OM ---
    // draw_source_OM_distributions(
    //     hSourceOM_pos,
    //     hSourceOM_neg,
    //     source_positions,
    //     om_positions);

    // draw_source_start_vertices(
    // hSrcStart_pos,
    // hSrcStart_neg,
    // source_positions);


    std::cout << "Finished" << std::endl;
}
