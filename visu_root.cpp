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
#include <TLine.h>
#include <TLegend.h>
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
//  Допоміжна функція: знайти глобальний Z-range для набору гістограм
//  Повертає {zMin, zMax} по непустих бінах
// ============================================================

std::pair<double,double>
get_global_zrange(const std::vector<TH2F*>& hists)
{
    double zMin = std::numeric_limits<double>::max();
    double zMax = 0.0;
    for(auto* h : hists)
    {
        if(!h) continue;
        double mn = h->GetMinimum(0);   // мін > 0 (ігноруємо порожні біни)
        double mx = h->GetMaximum();
        if(mn < zMin) zMin = mn;
        if(mx > zMax) zMax = mx;
    }
    if(zMin == std::numeric_limits<double>::max()) zMin = 0;
    return {zMin, zMax};
}


// ============================================================
//  Побудова енергетичного спектру:
//  три криві на одному canvas —
//  French side (X>0), Italian side (X<0), разом
//  У легенді — кількість подій для кожної лінії
//
//  [PARAM] unifiedYScale — якщо true, всі три гістограми отримують
//          однаковий діапазон по осі Y (від 0 до спільного максимуму).
//          За замовчуванням: true.
// ============================================================

void draw_energy_spectra(TH1F* hE_pos, TH1F* hE_neg, TH1F* hE_all,
                         bool unifiedYScale = true)
{
    gSystem->mkdir("plots", true);

    TCanvas* c = new TCanvas("cEnergy", "Energy spectra", 900, 650);
    c->SetLeftMargin(0.12);
    c->SetRightMargin(0.05);
    c->SetTopMargin(0.09);
    c->SetBottomMargin(0.11);

    hE_all->SetLineColor(kBlack);
    hE_all->SetLineWidth(2);
    hE_pos->SetLineColor(kBlue + 1);
    hE_pos->SetLineWidth(2);
    hE_neg->SetLineColor(kRed  + 1);
    hE_neg->SetLineWidth(2);

    hE_all->SetStats(0);
    hE_pos->SetStats(0);
    hE_neg->SetStats(0);

    if(unifiedYScale)
    {
        // Знаходимо спільний максимум по осі Y серед усіх трьох гістограм
        double yMax = std::max({ hE_all->GetMaximum(),
                                 hE_pos->GetMaximum(),
                                 hE_neg->GetMaximum() });
        hE_all->SetMaximum(yMax * 1.05);   // +5% поля зверху
        hE_pos->SetMaximum(yMax * 1.05);
        hE_neg->SetMaximum(yMax * 1.05);
    }
    else
    {
        // Скидаємо можливі попередні обмеження, дозволяємо ROOT вибирати авто
        hE_all->SetMaximum(-1111);
        hE_pos->SetMaximum(-1111);
        hE_neg->SetMaximum(-1111);
    }

    hE_all->SetTitle("Energy spectrum;Energy [MeV];Counts");
    hE_all->Draw("HIST");
    hE_pos->Draw("HIST same");
    hE_neg->Draw("HIST same");

    long long nAll = (long long)hE_all->GetEntries();
    long long nPos = (long long)hE_pos->GetEntries();
    long long nNeg = (long long)hE_neg->GetEntries();

    TLegend* leg = new TLegend(0.55, 0.68, 0.93, 0.88);
    leg->SetBorderSize(1);
    leg->SetFillStyle(1001);
    leg->SetFillColor(kWhite);
    leg->AddEntry(hE_all,
        Form("Both sides  (N = %lld)", nAll), "l");
    leg->AddEntry(hE_pos,
        Form("French side X>0  (N = %lld)", nPos), "l");
    leg->AddEntry(hE_neg,
        Form("Italian side X<0  (N = %lld)", nNeg), "l");
    leg->Draw();

    c->SaveAs("plots/energy_spectra.png");
    delete c;
}


// ============================================================
//  Карти точок початку для кожного джерела окремо
//  При side==0 (X>0, French) координата Y інвертується: y -> -y
//
//  [PARAM] unifiedScale — якщо true (за замовчуванням), одна шкала colorbar
//          для всіх гістограм однієї сторони (pos окремо, neg окремо).
//          Якщо false — кожна гістограма масштабується індивідуально.
// ============================================================

void draw_source_start_vertices(
    const std::vector<TH2F*>& hStart_pos,
    const std::vector<TH2F*>& hStart_neg,
    const std::vector<std::pair<double,double>>& src_pos,
    bool unifiedScale = true)
{
    gSystem->mkdir("plots/bin_sources", true);
    const int nSources = (int)src_pos.size();

    // Обчислюємо глобальні межі лише якщо потрібна спільна шкала
    double zMin_pos = 0, zMax_pos = 0;
    double zMin_neg = 0, zMax_neg = 0;

    if(unifiedScale)
    {
        std::vector<TH2F*> all_pos(hStart_pos.begin(), hStart_pos.end());
        std::vector<TH2F*> all_neg(hStart_neg.begin(), hStart_neg.end());

        auto [zp_min, zp_max] = get_global_zrange(all_pos);
        auto [zn_min, zn_max] = get_global_zrange(all_neg);
        zMin_pos = zp_min; zMax_pos = zp_max;
        zMin_neg = zn_min; zMax_neg = zn_max;
    }

    for(int s = 0; s < nSources; s++)
    {
        double yc_orig = src_pos[s].first;
        double zc      = src_pos[s].second;

        for(int side = 0; side < 2; side++)
        {
            TH2F* h = (side == 0) ? hStart_pos[s] : hStart_neg[s];
            double yc = (side == 0) ? -yc_orig : yc_orig;

            double nEvents = h->GetEntries();
            TString newTitle = Form("%s  (N = %.0f)", h->GetTitle(), nEvents);
            h->SetTitle(newTitle.Data());

            if(unifiedScale)
            {
                double zMin = (side == 0) ? zMin_pos : zMin_neg;
                double zMax = (side == 0) ? zMax_pos : zMax_neg;
                h->SetMinimum(zMin);
                h->SetMaximum(zMax);
            }
            else
            {
                // Автоматичний масштаб ROOT для кожної гістограми окремо
                h->SetMinimum(-1111);
                h->SetMaximum(-1111);
            }

            TCanvas* c = new TCanvas(
                Form("c_sv_%d_%d", side, s), "", 900, 800);
            c->SetLeftMargin(0.12);
            c->SetRightMargin(0.15);
            c->SetTopMargin(0.09);
            c->SetBottomMargin(0.11);
            h->SetStats(0);
            h->Draw("COLZ");

            TGraph* gSrc = new TGraph(1, &yc, &zc);
            gSrc->SetMarkerStyle(20);
            gSrc->SetMarkerSize(1.8);
            gSrc->SetMarkerColor(kRed);
            gSrc->Draw("P same");

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

            // Відновлюємо оригінальний заголовок
            h->SetTitle(TString(h->GetTitle())
                        .ReplaceAll(Form("  (N = %.0f)", nEvents), "")
                        .Data());
            // Скидаємо встановлені межі, щоб не впливати на інші виклики
            h->SetMinimum(-1111);
            h->SetMaximum(-1111);

            delete c;
        }
    }
}


// ============================================================
//  Сумарний розподіл точок вильоту (pos + neg) для кожного джерела
//  Зберігається у plots/bin_sources_sum/
//
//  [PARAM] unifiedScale — якщо true (за замовчуванням), одна шкала colorbar
//          для всіх sum-гістограм.
//          Якщо false — кожна sum-гістограма масштабується індивідуально.
// ============================================================

void draw_source_start_vertices_sum(
    const std::vector<TH2F*>& hStart_pos,
    const std::vector<TH2F*>& hStart_neg,
    const std::vector<std::pair<double,double>>& src_pos,
    bool unifiedScale = false)
{
    gSystem->mkdir("plots/bin_sources_sum", true);
    const int nSources = (int)src_pos.size();

    std::vector<TH2F*> hSum(nSources, nullptr);

    for(int s = 0; s < nSources; s++)
    {
        TH2F* hP = hStart_pos[s];
        TH2F* hN = hStart_neg[s];

        int    nBinsX = hP->GetNbinsX();
        int    nBinsY = hP->GetNbinsY();

        double yMin = hN->GetXaxis()->GetXmin();
        double yMax = hN->GetXaxis()->GetXmax();
        double zMin = hP->GetYaxis()->GetXmin();
        double zMax = hP->GetYaxis()->GetXmax();

        hSum[s] = new TH2F(
            Form("hSrcStart_sum_%02d", s),
            Form("Start vertices src %d (both sides);Y [mm];Z [mm]", s),
            nBinsX, yMin, yMax,
            nBinsY, zMin, zMax);
        hSum[s]->SetDirectory(0);

        // Додаємо neg напряму
        hSum[s]->Add(hN);

        // Додаємо pos з інверсією X
        for(int ix = 1; ix <= nBinsX; ix++)
        {
            int ix_sum = nBinsX + 1 - ix;
            for(int iy = 1; iy <= nBinsY; iy++)
            {
                double w = hP->GetBinContent(ix, iy);
                hSum[s]->SetBinContent(
                    ix_sum, iy,
                    hSum[s]->GetBinContent(ix_sum, iy) + w);
            }
        }
    }

    // Обчислюємо глобальну шкалу лише якщо потрібна
    double zMin_sum = 0, zMax_sum = 0;
    if(unifiedScale)
    {
        std::vector<TH2F*> all_sum(hSum.begin(), hSum.end());
        auto [zm, zM] = get_global_zrange(all_sum);
        zMin_sum = zm; zMax_sum = zM;
    }

    for(int s = 0; s < nSources; s++)
    {
        double yc = src_pos[s].first;
        double zc = src_pos[s].second;

        double nEvents = hSum[s]->GetEntries();
        hSum[s]->SetTitle(
            Form("Start vertices src %d (both sides)  (N = %.0f);Y [mm];Z [mm]",
                 s, nEvents));

        if(unifiedScale)
        {
            hSum[s]->SetMinimum(zMin_sum);
            hSum[s]->SetMaximum(zMax_sum);
        }
        else
        {
            hSum[s]->SetMinimum(-1111);
            hSum[s]->SetMaximum(-1111);
        }

        TCanvas* c = new TCanvas(Form("c_sum_%d", s), "", 900, 800);
        c->SetLeftMargin(0.12);
        c->SetRightMargin(0.15);
        c->SetTopMargin(0.09);
        c->SetBottomMargin(0.11);
        hSum[s]->SetStats(0);
        hSum[s]->Draw("COLZ");

        TGraph* gSrc = new TGraph(1, &yc, &zc);
        gSrc->SetMarkerStyle(20);
        gSrc->SetMarkerSize(1.8);
        gSrc->SetMarkerColor(kRed);
        gSrc->Draw("P same");

        double zoom = (hSum[s]->GetXaxis()->GetXmax() -
                       hSum[s]->GetXaxis()->GetXmin()) / 2.0;
        TLatex* lbl = new TLatex(
            yc + zoom * 0.08, zc + zoom * 0.08,
            Form("src %d  (%.0f, %.0f)", s, yc, zc));
        lbl->SetTextSize(0.032);
        lbl->SetTextColor(kRed);
        lbl->Draw("same");

        c->SaveAs(Form("plots/bin_sources_sum/start_sum_src%02d.png", s));
        delete c;
    }

    for(int s = 0; s < nSources; s++) delete hSum[s];
}


// ============================================================
//  Карти джерело -> OM
//
//  [PARAM] unifiedScale — якщо true (за замовчуванням), одна шкала colorbar
//          для всіх гістограм однієї сторони (pos окремо, neg окремо).
//          Якщо false — кожна гістограма масштабується індивідуально.
// ============================================================

void draw_source_OM_distributions(
    TH2F* hSourceOM_pos,
    TH2F* hSourceOM_neg,
    const std::vector<std::pair<double,double>>& src_pos,
    const std::vector<std::pair<double,double>>& om_pos,
    bool unifiedScale = true)
{
    gSystem->mkdir("plots/source_maps", true);

    const int nSources = (int)src_pos.size();
    const int nOMs     = (int)om_pos.size();

    std::vector<double> unique_y_neg, unique_z;
    for(auto& p : om_pos)
    {
        unique_y_neg.push_back(p.first);
        unique_z    .push_back(p.second);
    }

    std::vector<double> unique_y_pos;
    for(auto& p : om_pos)
        unique_y_pos.push_back(-p.first);

    auto dedup = [](std::vector<double>& v)
    {
        std::sort(v.begin(), v.end());
        std::vector<double> out;
        out.push_back(v[0]);
        for(size_t i = 1; i < v.size(); i++)
            if(std::abs(v[i] - out.back()) > 1.0)
                out.push_back(v[i]);
        v = out;
    };
    dedup(unique_y_pos);
    dedup(unique_y_neg);
    dedup(unique_z);

    const int nColsTotal = (int)unique_y_neg.size();
    const int nRowsTotal = (int)unique_z.size();

    auto coord_to_idx = [](const std::vector<double>& vals, double val) -> int
    {
        int best = 0;
        double bestDist = std::abs(vals[0] - val);
        for(int i = 1; i < (int)vals.size(); i++)
        {
            double d = std::abs(vals[i] - val);
            if(d < bestDist){ bestDist = d; best = i; }
        }
        return best;
    };

    auto calc_step = [](const std::vector<double>& uv) -> double
    {
        if(uv.size() < 2) return 1.0;
        std::vector<double> diffs;
        for(size_t i = 1; i < uv.size(); i++)
            diffs.push_back(uv[i] - uv[i-1]);
        std::sort(diffs.begin(), diffs.end());
        return diffs[diffs.size()/2];
    };

    double stepY_pos = calc_step(unique_y_pos);
    double stepY_neg = calc_step(unique_y_neg);
    double stepZ     = calc_step(unique_z);

    const int colOffset_pos = 20;
    const int colOffset_neg =  0;

    struct SideParam {
        TH2F*                   matrix;
        int                     colOffset;
        int                     nCols;
        int                     nRows;
        const char*             sideName;
        const char*             sideTitle;
        const std::vector<double>* unique_y;
        double                  stepY;
    };

    SideParam sides[2] = {
        { hSourceOM_pos, colOffset_pos, nColsTotal, nRowsTotal,
          "pos", "French side",  &unique_y_pos, stepY_pos },
        { hSourceOM_neg, colOffset_neg, nColsTotal, nRowsTotal,
          "neg", "Italian side", &unique_y_neg, stepY_neg }
    };

    gStyle->SetPalette(kBird);
    gStyle->SetNumberContours(99);

    // Крок 0: будуємо всі hGrid для обох сторін
    std::vector<std::vector<TH2F*>> allGrids(2, std::vector<TH2F*>(nSources, nullptr));

    for(int side = 0; side < 2; side++)
    {
        SideParam& sp = sides[side];
        const int colMin = sp.colOffset;
        const int colMax = sp.colOffset + sp.nCols - 1;
        const int rowMin = 0;
        const int rowMax = sp.nRows - 1;

        for(int s = 0; s < nSources; s++)
        {
            double totalHits = 0;
            for(int om = 0; om < nOMs; om++)
                totalHits += sp.matrix->GetBinContent(s + 1, om + 1);

            TH2F* hGrid = new TH2F(
                Form("hGrid_%d_%d", side, s),
                Form("Source %d  %s  (total events: %.0f);Column;Row",
                     s, sp.sideTitle, totalHits),
                sp.nCols, colMin - 0.5, colMax + 0.5,
                sp.nRows, rowMin - 0.5, rowMax + 0.5);
            hGrid->SetDirectory(0);
            hGrid->SetStats(0);

            for(int om = 0; om < nOMs; om++)
            {
                double hits = sp.matrix->GetBinContent(s + 1, om + 1);
                if(hits <= 0) continue;

                double om_y_for_search = (side == 0) ? -om_pos[om].first
                                                      :  om_pos[om].first;
                int ci = coord_to_idx(*sp.unique_y, om_y_for_search) + sp.colOffset;
                int ri = coord_to_idx(unique_z, om_pos[om].second);

                hGrid->Fill(ci, ri, hits);
            }

            allGrids[side][s] = hGrid;
        }
    }

    // Крок 1: глобальний zMax для кожної сторони (якщо потрібна спільна шкала)
    double globalZmax_pos = 0, globalZmax_neg = 0;
    if(unifiedScale)
    {
        for(int s = 0; s < nSources; s++)
        {
            if(allGrids[0][s]) globalZmax_pos = std::max(globalZmax_pos, allGrids[0][s]->GetMaximum());
            if(allGrids[1][s]) globalZmax_neg = std::max(globalZmax_neg, allGrids[1][s]->GetMaximum());
        }
    }

    // Крок 2: малюємо
    for(int side = 0; side < 2; side++)
    {
        SideParam& sp = sides[side];

        const int colMin = sp.colOffset;
        const int colMax = sp.colOffset + sp.nCols - 1;
        const int rowMin = 0;
        const int rowMax = sp.nRows - 1;

        const double mL = 0.10, mR = 0.14, mT = 0.09, mB = 0.10;
        const int    baseH  = 600;
        const double areaH  = baseH * (1.0 - mT - mB);
        const double areaW  = areaH * (double)sp.nCols / (double)sp.nRows;
        const int    canvasW = (int)(areaW / (1.0 - mL - mR)) + 1;
        const int    canvasH = baseH;

        for(int s = 0; s < nSources; s++)
        {
            TH2F* hGrid = allGrids[side][s];

            if(unifiedScale)
            {
                double globalZmax = (side == 0) ? globalZmax_pos : globalZmax_neg;
                hGrid->SetMinimum(0);
                hGrid->SetMaximum(globalZmax);
            }
            else
            {
                // Автоматичний масштаб ROOT
                hGrid->SetMinimum(-1111);
                hGrid->SetMaximum(-1111);
            }

            TCanvas* c = new TCanvas(
                Form("cGrid_%d_%d", side, s), "", canvasW, canvasH);
            c->SetLeftMargin(mL);
            c->SetRightMargin(mR);
            c->SetTopMargin(mT);
            c->SetBottomMargin(mB);
            c->SetCanvasSize(canvasW, canvasH);

            hGrid->GetXaxis()->SetNdivisions(sp.nCols, false);
            hGrid->GetYaxis()->SetNdivisions(sp.nRows, false);
            hGrid->GetXaxis()->SetLabelSize(0.048);
            hGrid->GetYaxis()->SetLabelSize(0.048);
            hGrid->GetXaxis()->SetTickLength(0.01);
            hGrid->GetYaxis()->SetTickLength(0.01);

            for(int ic = 0; ic < sp.nCols; ic++)
                hGrid->GetXaxis()->SetBinLabel(ic + 1, Form("%d", colMin + ic));
            for(int ir = 0; ir < sp.nRows; ir++)
                hGrid->GetYaxis()->SetBinLabel(ir + 1, Form("%d", rowMin + ir));

            hGrid->Draw("COLZ");

            for(int ic = 0; ic <= sp.nCols; ic++)
            {
                double x = colMin - 0.5 + ic;
                TLine* l = new TLine(x, rowMin - 0.5, x, rowMax + 0.5);
                l->SetLineColor(kBlack);
                l->SetLineWidth(2);
                l->SetLineStyle(2);
                l->Draw("same");
            }
            for(int ir = 0; ir <= sp.nRows; ir++)
            {
                double y = rowMin - 0.5 + ir;
                TLine* l = new TLine(colMin - 0.5, y, colMax + 0.5, y);
                l->SetLineColor(kBlack);
                l->SetLineWidth(2);
                l->SetLineStyle(2);
                l->Draw("same");
            }

            for(int i = 0; i < nSources; i++)
            {
                double sy_orig = src_pos[i].first;
                double sz      = src_pos[i].second;
                double sy = (side == 0) ? -sy_orig : sy_orig;

                double col_src = sp.colOffset +
                    (sy - (*sp.unique_y)[0]) / sp.stepY;
                double row_src =
                    (sz - unique_z[0]) / stepZ;

                int    markerColor = (i == s) ? kRed  : kGreen + 1;
                double markerSize  = (i == s) ? 1.3    : 0.8;
                double textSize    = (i == s) ? 0.045  : 0.0;

                TGraph* gSrc = new TGraph(1, &col_src, &row_src);
                gSrc->SetMarkerStyle(20);
                gSrc->SetMarkerSize(markerSize);
                gSrc->SetMarkerColor(markerColor);
                gSrc->Draw("P same");

                TLatex* lsrc = new TLatex(
                    col_src + 0.15, row_src + 0.15,
                    Form("%d", i));
                lsrc->SetTextSize(textSize);
                lsrc->SetTextColor(markerColor);
                lsrc->Draw("same");
            }

            c->SaveAs(Form("plots/source_maps/source_%s_%02d.png",
                           sp.sideName, s));

            delete hGrid;
            delete c;
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

    TH1F* hEnergy_pos = new TH1F(
        "Energy_pos", "French side;Energy [MeV];Counts",
        2000, 0, 3000);
    hEnergy_pos->SetDirectory(0);

    TH1F* hEnergy_neg = new TH1F(
        "Energy_neg", "Italian side;Energy [MeV];Counts",
        2000, 0, 3000);
    hEnergy_neg->SetDirectory(0);

    TH2F* hStart_Xpos = new TH2F("hStart_Xpos", "", 500, -2500, 2500, 500, -2000, 2000);
    TH2F* hStart_Xneg = new TH2F("hStart_Xneg", "", 500, -2500, 2500, 500, -2000, 2000);
    TH2F* hEnd_Xpos   = new TH2F("hEnd_Xpos",   "", 520, -2600, 2600, 500, -2000, 2000);
    TH2F* hEnd_Xneg   = new TH2F("hEnd_Xneg",   "", 520, -2600, 2600, 500, -2000, 2000);

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

    const double zoom = 50.0;
    const int    nBinsZoom = 100;
    
    std::vector<TH2F*> hSrcStart_pos(nSources);
    std::vector<TH2F*> hSrcStart_neg(nSources);
    
    for(int s = 0; s < nSources; s++)
    {
        double yc_neg = source_positions[s].first;
        double yc_pos = -source_positions[s].first;
        double zc     = source_positions[s].second;
    
        hSrcStart_pos[s] = new TH2F(
            Form("hSrcStart_pos_%02d", s),
            Form("Start vertices src %d (X>0);-Y [mm];Z [mm]", s),
            nBinsZoom, yc_pos - zoom, yc_pos + zoom,
            nBinsZoom, zc     - zoom, zc     + zoom);
        hSrcStart_pos[s]->SetDirectory(0);

        hSrcStart_neg[s] = new TH2F(
            Form("hSrcStart_neg_%02d", s),
            Form("Start vertices src %d (X<0);Y [mm];Z [mm]", s),
            nBinsZoom, yc_neg - zoom, yc_neg + zoom,
            nBinsZoom, zc     - zoom, zc     + zoom);
        hSrcStart_neg[s]->SetDirectory(0);
    }

    
    // --- Запам'ятовуємо X для заголовків ---
    double x_pos_value = 0, x_neg_value = 0;
    bool   pos_found   = false, neg_found = false;

    // --- Головний цикл ---
    Long64_t nentries = evTree->GetEntries();
    // Long64_t nentries = 1e6;  // для тесту
    std::cout << "Reading " << nentries << " entries..." << std::endl;

    for(Long64_t ie = 0; ie < nentries; ie++)
    {
        evTree->GetEntry(ie);

        if(side == 1)   // French side, X > 0
        {
            hStart_Xpos->Fill(-y_start, z_start);
            if(!pos_found){ x_pos_value = x_start; pos_found = true; }
        }
        else            // Italian side, X < 0
        {
            hStart_Xneg->Fill(y_start, z_start);
            if(!neg_found){ x_neg_value = x_start; neg_found = true; }
        }

        if(source_id >= 0 && source_id < nSources)
        {
            if(side == 1) hSrcStart_pos[source_id]->Fill(-y_start, z_start);
            else          hSrcStart_neg[source_id]->Fill( y_start, z_start);
        }

        if(x_end > 0) hEnd_Xpos->Fill(-y_end, z_end);
        else          hEnd_Xneg->Fill( y_end, z_end);

        if(energy > 0)
        {
            hEnergy->Fill(energy);
            if(side == 1) hEnergy_pos->Fill(energy);
            else          hEnergy_neg->Fill(energy);
        }

        if(source_id >= 0 && om_id >= 0)
        {
            if(side == 1) hSourceOM_pos->Fill(source_id, om_id);
            else          hSourceOM_neg->Fill(source_id, om_id);
        }
    }

    f->Close();

    // --- Оновлюємо заголовки загальних гістограм ---
    hStart_Xpos->SetTitle(Form("Start vertices (X = %.1f mm);-Y [mm];Z [mm]", x_pos_value));
    hStart_Xneg->SetTitle(Form("Start vertices (X = %.1f mm);Y [mm];Z [mm]",  x_neg_value));
    hEnd_Xpos->SetTitle( Form("End vertices (X #approx %.1f mm);-Y [mm];Z [mm]", x_pos_value));
    hEnd_Xneg->SetTitle( Form("End vertices (X #approx %.1f mm);Y [mm];Z [mm]",  x_neg_value));

    // ============================================================
    //  Єдина шкала colorbar для пар: start_Xpos/start_Xneg
    //  та end_Xpos/end_Xneg
    // ============================================================

    // Пара start
    {
        double zmaxStart = std::max(hStart_Xpos->GetMaximum(),
                                    hStart_Xneg->GetMaximum());
        hStart_Xpos->SetMinimum(0); hStart_Xpos->SetMaximum(zmaxStart);
        hStart_Xneg->SetMinimum(0); hStart_Xneg->SetMaximum(zmaxStart);
    }
    // Пара end
    {
        double zmaxEnd = std::max(hEnd_Xpos->GetMaximum(),
                                  hEnd_Xneg->GetMaximum());
        hEnd_Xpos->SetMinimum(0); hEnd_Xpos->SetMaximum(zmaxEnd);
        hEnd_Xneg->SetMinimum(0); hEnd_Xneg->SetMaximum(zmaxEnd);
    }

    // --- Малюємо загальні PNG ---
    TCanvas c1("c1_start_Xpos", "", 1300, 700); c1.SetRightMargin(0.35); hStart_Xpos->Draw("COLZ"); c1.SaveAs("plots/start_Xpos.png");
    TCanvas c2("c2_start_Xpos", "", 1300, 700); c2.SetRightMargin(0.35); hStart_Xneg->Draw("COLZ"); c2.SaveAs("plots/start_Xneg.png");

    // [NEW] ======================================================
    //  Побудова end_Xpos та end_Xneg з сіткою оптичних модулів
    //  (стиль аналогічний draw_source_OM_distributions)
    // ============================================================

    // --- Допоміжні функції для побудови унікальних координатних сіток ---

    // Обчислюємо унікальні відсортовані значення Y і Z з om_positions
    // Для X>0: cx = -y_om; для X<0: cx = +y_om
    auto calc_step_vec = [](std::vector<double>& v) -> double
    {
        if(v.size() < 2) return 1.0;
        std::vector<double> diffs;
        for(size_t i = 1; i < v.size(); i++)
            diffs.push_back(v[i] - v[i-1]);
        std::sort(diffs.begin(), diffs.end());
        return diffs[diffs.size() / 2];
    };

    // Будуємо унікальні Z (спільні для обох сторін)
    std::vector<double> om_z_vals;
    for(auto& p : om_positions)
        om_z_vals.push_back(p.second);
    std::sort(om_z_vals.begin(), om_z_vals.end());
    om_z_vals.erase(std::unique(om_z_vals.begin(), om_z_vals.end()), om_z_vals.end());
    // Залишаємо унікальні з допуском 1 мм
    {
        std::vector<double> tmp;
        tmp.push_back(om_z_vals[0]);
        for(size_t i = 1; i < om_z_vals.size(); i++)
            if(std::abs(om_z_vals[i] - tmp.back()) > 1.0)
                tmp.push_back(om_z_vals[i]);
        om_z_vals = tmp;
    }
    double stepZ_om = calc_step_vec(om_z_vals);

    // Для кожної сторони: унікальні X (= ±y_om)
    auto build_unique_x = [&](bool isXpos) -> std::vector<double>
    {
        std::vector<double> xv;
        for(auto& p : om_positions)
            xv.push_back(isXpos ? -p.first : p.first);
        std::sort(xv.begin(), xv.end());
        std::vector<double> tmp;
        tmp.push_back(xv[0]);
        for(size_t i = 1; i < xv.size(); i++)
            if(std::abs(xv[i] - tmp.back()) > 1.0)
                tmp.push_back(xv[i]);
        return tmp;
    };

    // --- Лямбда: малює червону переривчасту сітку навколо OM-позицій ---
    // halfX, halfZ — половини розміру комірки (крок/2)
    auto draw_om_grid = [&](bool isXpos)
    {
        std::vector<double> unique_x = build_unique_x(isXpos);
        double stepX = calc_step_vec(unique_x);
        double hX    = stepX * 0.5;
        double hZ    = stepZ_om * 0.5;

        for(auto& p : om_positions)
        {
            double cx = isXpos ? -p.first : p.first;
            double cz = p.second;

            // Чотири сторони прямокутника навколо OM
            double x1 = cx - hX, x2 = cx + hX;
            double z1 = cz - hZ, z2 = cz + hZ;

            TLine* lB = new TLine(x1, z1, x2, z1);  // bottom
            TLine* lT = new TLine(x1, z2, x2, z2);  // top
            TLine* lL = new TLine(x1, z1, x1, z2);  // left
            TLine* lR = new TLine(x2, z1, x2, z2);  // right

            for(TLine* l : {lB, lT, lL, lR})
            {
                l->SetLineColor(kRed);
                l->SetLineStyle(2);   // dashed
                l->SetLineWidth(1);
                l->Draw("same");
            }
        }
    };

    // --- Малюємо end_Xpos з сіткою OM ---
    {
        TCanvas c3("c3_end_Xpos", "", 1300, 700);
        c3.SetLeftMargin(0.10);
        c3.SetRightMargin(0.35);   // [CHANGE] ширший правий відступ — colorbar не налазить
        c3.SetTopMargin(0.09);
        c3.SetBottomMargin(0.10);
        // hEnd_Xpos->SetStats(0);
        hEnd_Xpos->Draw("COLZ");
        draw_om_grid(true);
        c3.SaveAs("plots/end_Xpos.png");
    }

    // --- Малюємо end_Xneg з сіткою OM ---
    {
        TCanvas c4("c4_end_Xneg", "", 1300, 700);
        c4.SetLeftMargin(0.10);
        c4.SetRightMargin(0.35);   // [CHANGE] ширший правий відступ — colorbar не налазить
        c4.SetTopMargin(0.09);
        c4.SetBottomMargin(0.10);
        // hEnd_Xneg->SetStats(0);
        hEnd_Xneg->Draw("COLZ");
        draw_om_grid(true);
        c4.SaveAs("plots/end_Xneg.png");
    }
    // [END NEW] ==================================================

    TCanvas c5; hEnergy->Draw("PE");            c5.SaveAs("plots/energy_spectrum.png");
    TCanvas c6; hSourceOM_pos->Draw("COLZ"); c6.SaveAs("plots/source_OM_matrix_fr.png");
    TCanvas c7; hSourceOM_neg->Draw("COLZ"); c7.SaveAs("plots/source_OM_matrix_ital.png");

    // --- Зберігаємо гістограми у ROOT файл ---
    TFile* fout = new TFile("vertex_energy_results.root", "RECREATE");
    hStart_Xpos->Write();
    hStart_Xneg->Write();
    hEnd_Xpos->Write();
    hEnd_Xneg->Write();
    hEnergy->Write();
    hEnergy_pos->Write();
    hEnergy_neg->Write();
    hSourceOM_pos->Write();
    hSourceOM_neg->Write();
    for(int s = 0; s < nSources; s++)
    {
        hSrcStart_pos[s]->Write();
        hSrcStart_neg[s]->Write();
    }
    fout->Close();

    // --- Енергетичний спектр з трьома кривими ---
    draw_energy_spectra(hEnergy_pos, hEnergy_neg, hEnergy);

    // --- Детальні карти джерело -> OM (єдина шкала всередині функції) ---
    draw_source_OM_distributions(
        hSourceOM_pos,
        hSourceOM_neg,
        source_positions,
        om_positions);

    // --- Карти точок початку для кожного джерела (єдина шкала) ---
    draw_source_start_vertices(
        hSrcStart_pos,
        hSrcStart_neg,
        source_positions);

    // --- Сумарний розподіл (pos + neg) для кожного джерела ---
    draw_source_start_vertices_sum(
        hSrcStart_pos,
        hSrcStart_neg,
        source_positions);

    std::cout << "Finished" << std::endl;
}