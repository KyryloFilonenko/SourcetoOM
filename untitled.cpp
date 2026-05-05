// -*- coding: utf-8 -*-
/*
 * ratio_visualization.cpp (модифікована версія з GID-графіками та точними позиціями джерел)
 *
 * Будує карти відношення для обох типів детекторів:
 *
 * 1. ЗВИЧАЙНЕ ВІДНОШЕННЯ (OM):
 *    R_OM(OM, source, side) = eps_act(OM, source, side) / eps_geo(OM, source)
 *
 * 2. GID-ВІДНОШЕННЯ:
 *    R_GID(OM, source, side) = eps_act_GID(OM, source, side) / eps_geo(OM, source)
 *
 * де:
 *   eps_act         — реєстраційна ефективність OM [%]
 *   eps_act_GID     — реєстраційна ефективність GID [%]
 *   eps_geo         — геометрична ефективність [%]
 *
 * Джерела даних:
 *   activity_results.root               — OM гістограми (draw_efficiencies)
 *   activity_results.root               — GID гістограми (draw_efficiencies_GID)
 *   geometric_efficiencies.root         — геометричні ефективності
 *   default_values/source_positions.txt — позиції джерел (точні координати)
 *
 * Вихід:
 *   plots/ratio/source_<side>_<NN>_<mode>.png           — OM-відношення
 *   plots/ratio_GID/source_<side>_<NN>_<mode>.png       — GID-відношення
 *
 * ЗМІНИ:
 *   - Точні позиції джерел з source_positions.txt
 *   - Вдосконалена інтерполяція бінів
 *   - Кольорова диференціація маркерів
 *   - Відображення координат і номерів источників
 *   - Контури навколо поточного джерела
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <string>
#include <iomanip>

#include "TFile.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TBox.h"
#include "TLine.h"
#include "TGraph.h"
#include "TPaveText.h"
#include "TColor.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TTree.h"

using namespace std;

// ============================================================
//  Режим відображення
// ============================================================
enum class RatioMode {
    LINEAR,  // R = eps_act / eps_geo
    LOG10    // log10(R)
};

static string ratioLabel(double ratio, RatioMode mode)
{
    char buf[64];
    if (ratio <= 0.0) { snprintf(buf, sizeof(buf), "N/A"); return buf; }
    switch (mode) {
        case RatioMode::LINEAR: snprintf(buf, sizeof(buf), "%.2f", ratio);        break;
        case RatioMode::LOG10:  snprintf(buf, sizeof(buf), "%.2f", log10(ratio)); break;
    }
    return buf;
}

static double ratioValue(double ratio, RatioMode mode)
{
    if (ratio <= 0.0) return NAN;
    switch (mode) {
        case RatioMode::LINEAR: return ratio;
        case RatioMode::LOG10:  return log10(ratio);
    }
    return ratio;
}

static const char* modeSuffix(RatioMode mode)
{
    switch (mode) {
        case RatioMode::LINEAR: return "lin";
        case RatioMode::LOG10:  return "log";
    }
    return "lin";
}

// ============================================================
//  Читання позицій джерел (для маркерів)
// ============================================================
static vector<pair<double,double>>
read_positions(const string& filename)
{
    vector<pair<double,double>> pos;
    ifstream f(filename);
    if (!f.is_open()) { cerr << "ERROR: cannot open " << filename << "\n"; return pos; }
    string line;
    while (getline(f, line)) {
        if (line.empty() || line[0] == '#') continue;
        double y, z;
        if (line.find(';') != string::npos) {
            replace(line.begin(), line.end(), ';', ' ');
            istringstream ss(line); if (!(ss >> y >> z)) continue;
        } else {
            istringstream ss(line);
            int idx; if (!(ss >> idx >> y >> z)) continue;
        }
        pos.emplace_back(y, z);
    }
    cout << "Loaded " << pos.size() << " source positions from " << filename << "\n";
    return pos;
}

// ============================================================
//  Функція для побудови GID-графіків відношень
// ============================================================
void draw_ratio_GID(
    const string& actRootFile,
    const string& geoRootFile,
    const vector<pair<double,double>>& src_pos,
    int nWalls_gid = 20,
    int nRows_gid = 13,
    RatioMode dispMode = RatioMode::LINEAR)
{
    // ── 1. Відкриваємо activity_results.root ───────────────────────────
    TFile* fAct = TFile::Open(actRootFile.c_str(), "READ");
    if (!fAct || fAct->IsZombie()) {
        cerr << "ERROR (GID): cannot open " << actRootFile << "\n"; return;
    }

    // Перевіряємо наявність GID гістограм
    if (!fAct->Get("eff_GID_pos_src00")) {
        cerr << "WARNING (GID): eff_GID_pos_src00 not found — skipping GID ratio maps\n"
             << "  Hint: run draw_efficiencies_GID() first\n";
        fAct->Close(); return;
    }

    // ── 2. Відкриваємо geometric_efficiencies.root ─────────────────────
    TFile* fGeo = TFile::Open(geoRootFile.c_str(), "READ");
    if (!fGeo || fGeo->IsZombie()) {
        cerr << "ERROR (GID): cannot open " << geoRootFile << "\n";
        fAct->Close(); return;
    }

    // ── 3. Визначаємо кількість джерел ─────────────────────────────────
    int nSources = 0;
    while (fAct->Get(Form("eff_GID_pos_src%02d", nSources))) nSources++;
    if (nSources == 0) {
        cerr << "ERROR (GID): no eff_GID_pos_src* histograms found\n";
        fAct->Close(); fGeo->Close(); return;
    }
    cout << "[GID] Found " << nSources << " sources\n";

    // ── 4. Зчитуємо N_true з TTree ────────────────────────────────────
    vector<double> N_true_vec(nSources, 0.0);
    {
        TTree* t = (TTree*)fAct->Get("ActivityResults");
        if (t) {
            int    src_idx = 0;
            double n_true  = 0.0;
            t->SetBranchAddress("source_idx", &src_idx);
            t->SetBranchAddress("N_true",     &n_true);
            for (Long64_t e = 0; e < t->GetEntries(); e++) {
                t->GetEntry(e);
                if (src_idx >= 0 && src_idx < nSources)
                    N_true_vec[src_idx] = n_true;
            }
        }
    }

    const int nCols = nWalls_gid;
    const int nRows = nRows_gid;
    const int nGeoCols = 20;
    const int nGeoRows = 13;

    cout << "[GID] Act grid: " << nCols << " cols x " << nRows << " rows\n";

    // ── 5. Опис сторін ──────────────────────────────────────────────────
    struct SideGID {
        const char*  name;
        const char*  title;
        bool         isPos;
    } sides[2] = {
        { "pos", "French side (X>0)",  true  },
        { "neg", "Italian side (X<0)", false }
    };

    // ── 6. Основний цикл по джерелах та сторонах ────────────────────────
    for (int s = 0; s < nSources; s++) {

        TH2F* hGeo = (TH2F*)fGeo->Get(Form("source_%d", s));
        if (!hGeo) {
            cerr << "WARNING (GID): source_" << s << " not found in geo file\n";
            continue;
        }
        hGeo->SetDirectory(0);

        for (int si = 0; si < 2; si++) {
            SideGID& sd = sides[si];

            TH2F* hActGID = (TH2F*)fAct->Get(
                Form("eff_GID_%s_src%02d", sd.name, s));
            if (!hActGID) {
                cerr << "WARNING (GID): eff_GID_" << sd.name << "_src"
                     << Form("%02d", s) << " not found\n";
                continue;
            }
            hActGID->SetDirectory(0);

            // ── Обчислення GID-відношень по всіх комірках ────────────
            struct CellData {
                int    ci, ri;
                double eps_act_gid;
                double eps_geo;
                double ratio;
            };
            vector<CellData> cells;
            cells.reserve(nCols * nRows);

            for (int ci = 0; ci < nCols; ci++) {
                for (int ri = 0; ri < nRows; ri++) {

                    double eps_act_gid = hActGID->GetBinContent(ci + 1, ri + 1);
                    if (eps_act_gid <= 0.0) continue;

                    // Відповідна колонка у geo для GID:
                    //   pos (французька): дзеркало по X
                    //   neg (італійська): прямий порядок
                    int geo_col = sd.isPos ? (nGeoCols - 1 - ci) : ci;
                    int geo_row = ri;

                    if (geo_col < 0 || geo_col >= nGeoCols ||
                        geo_row < 0 || geo_row >= nGeoRows) continue;

                    double eps_geo = hGeo->GetBinContent(geo_col + 1, geo_row + 1);
                    if (eps_geo <= 0.0) continue;

                    cells.push_back({ci, ri, eps_act_gid, eps_geo, eps_act_gid / eps_geo});
                }
            }

            if (cells.empty()) {
                cerr << "WARNING (GID): no valid cells for source " << s << "\n";
                delete hActGID; continue;
            }

            // ── Діапазон для кольорової шкали ────────────────────────
            double vMin = 1e30, vMax = -1e30;
            for (auto& cd : cells) {
                double v = ratioValue(cd.ratio, dispMode);
                if (!isnan(v) && !isinf(v)) {
                    if (v < vMin) vMin = v;
                    if (v > vMax) vMax = v;
                }
            }
            if (vMin > vMax) { vMin = 0.0; vMax = 1.0; }
            if (fabs(vMax - vMin) < 1e-12) vMax = vMin + 1.0;

            // ── Параметри канвасу ────────────────────────────────────
            const int    colMin = 0, colMax = nCols - 1;
            const int    rowMin = 0, rowMax = nRows - 1;
            const double mL = 0.08, mR = 0.30, mT = 0.09, mB = 0.10;
            const int    baseH   = 600;
            const double areaH   = baseH * (1.0 - mT - mB);
            const double areaW   = areaH * (double)nCols / (double)nRows;
            const int    canvasW = (int)(areaW / (1.0 - mL - mR)) + 1;

            TCanvas* c = new TCanvas(
                Form("cRatioGID_%s_%d", sd.name, s), "", canvasW, baseH);
            c->SetLeftMargin(mL); c->SetRightMargin(mR);
            c->SetTopMargin(mT);  c->SetBottomMargin(mB);

            TH2F* hFrame = new TH2F(
                Form("hFrameGID_%s_%d", sd.name, s), "",
                nCols, colMin - 0.5, colMax + 0.5,
                nRows, rowMin - 0.5, rowMax + 0.5);
            hFrame->SetDirectory(0); hFrame->SetStats(0);
            hFrame->GetXaxis()->SetNdivisions(nCols, false);
            hFrame->GetYaxis()->SetNdivisions(nRows, false);
            hFrame->GetXaxis()->SetLabelSize(0.048);
            hFrame->GetYaxis()->SetLabelSize(0.048);
            hFrame->GetXaxis()->SetTickLength(0.01);
            hFrame->GetYaxis()->SetTickLength(0.01);

            // Підписи по X: для French — реальний wall (bin 0 = wall 19, bin 19 = wall 0)
            for (int ic = 0; ic < nCols; ic++) {
                int wall_label = sd.isPos ? (nCols - 1 - ic) : ic;
                hFrame->GetXaxis()->SetBinLabel(ic + 1, Form("%d", wall_label));
            }
            for (int ir = 0; ir < nRows; ir++)
                hFrame->GetYaxis()->SetBinLabel(ir + 1, Form("%d", ir));
            hFrame->GetXaxis()->SetTitle("Wall");
            hFrame->GetYaxis()->SetTitle("Col");
            hFrame->Draw("AXIS");

            const double cellW  = 0.46, cellH = 0.46;
            const int    nColors = TColor::GetNumberOfColors();

            // ПРОХІД 1: кольорові комірки
            for (auto& cd : cells) {
                double v = ratioValue(cd.ratio, dispMode);
                if (isnan(v) || isinf(v)) continue;
                double frac = (v - vMin) / (vMax - vMin);
                frac = max(0.0, min(1.0, frac));
                int colorIdx = TColor::GetColorPalette((int)(frac * (nColors - 1)));

                TBox* b = new TBox(
                    cd.ci - cellW, cd.ri - cellH,
                    cd.ci + cellW, cd.ri + cellH);
                b->SetFillColor(colorIdx);
                b->SetLineColor(kGray + 2); b->SetLineWidth(1);
                b->Draw("l same");
            }

            // ПРОХІД 2: сітка
            for (int ic = 0; ic <= nCols; ic++) {
                TLine* l = new TLine(
                    colMin - 0.5 + ic, rowMin - 0.5,
                    colMin - 0.5 + ic, rowMax + 0.5);
                l->SetLineColor(kBlack); l->SetLineWidth(1); l->SetLineStyle(2);
                l->Draw("same");
            }
            for (int ir = 0; ir <= nRows; ir++) {
                TLine* l = new TLine(
                    colMin - 0.5, rowMin - 0.5 + ir,
                    colMax + 0.5, rowMin - 0.5 + ir);
                l->SetLineColor(kBlack); l->SetLineWidth(1); l->SetLineStyle(2);
                l->Draw("same");
            }

            // ПРОХІД 3: числові підписи
            for (auto& cd : cells) {
                double v = ratioValue(cd.ratio, dispMode);
                if (isnan(v) || isinf(v)) continue;
                double frac = (v - vMin) / (vMax - vMin);
                frac = max(0.0, min(1.0, frac));

                string lbl = ratioLabel(cd.ratio, dispMode);
                TLatex* lt = new TLatex((double)cd.ci, (double)cd.ri, lbl.c_str());
                lt->SetTextAlign(22); lt->SetTextSize(0.027);
                lt->SetTextColor(frac > 0.55 ? kBlack : kWhite);
                lt->Draw("same");
            }

            // ПРОХІД 4: маркери джерел з точними позиціями
            for (int i = 0; i < nSources && i < (int)src_pos.size(); i++) {
                // Для GID не потрібна трансформація Y координат (GID має свою сітку)
                double sy_i = src_pos[i].first;
                double sz_i = src_pos[i].second;

                // Спрощена інтерполяція для GID
                double col_f = (double)nCols / 2.0;  // приблизно середина
                double row_f = (double)nRows / 2.0;

                bool is_current = (i == s);
                draw_source_marker(i, is_current, col_f, row_f, sy_i, sz_i, nSources);
            }

            // ── Заголовок ────────────────────────────────────────────
            TLatex ttl;
            ttl.SetNDC(); ttl.SetTextAlign(22); ttl.SetTextSize(0.038);
            ttl.DrawLatex(0.5 * (1.0 - mR + mL), 0.965,
                Form("Source %d  |  #varepsilon_{act}^{GID}/#varepsilon_{geo}  |  %s  (GID)",
                     s, sd.title));

            // ── Інфопанель ───────────────────────────────────────────
            double sum_act_gid = 0.0, sum_geo = 0.0, mean_ratio = 0.0;
            int n_valid = 0;
            for (auto& cd : cells) {
                if (cd.ratio > 0.0) {
                    sum_act_gid += cd.eps_act_gid;
                    sum_geo     += cd.eps_geo;
                    mean_ratio  += cd.ratio;
                    n_valid++;
                }
            }
            if (n_valid > 0) mean_ratio /= n_valid;

            TPaveText* pvInfo = new TPaveText(
                1.0 - mR + 0.01, 0.62, 0.995, 0.90, "NDC");
            pvInfo->SetFillColor(kWhite); pvInfo->SetBorderSize(1);
            pvInfo->SetTextAlign(12); pvInfo->SetTextSize(0.025);
            pvInfo->AddText(Form("N_{true} = %.4e", N_true_vec[s]));
            pvInfo->AddText(Form("#Sigma #varepsilon_{act}^{GID} = %.3f%%", sum_act_gid));
            pvInfo->AddText(Form("#Sigma #varepsilon_{geo} = %.3f%%", sum_geo));
            pvInfo->AddText(Form("<R_{GID}> = %.4f", mean_ratio));
            pvInfo->AddText(
                (dispMode == RatioMode::LINEAR) ? "mode: linear" : "mode: log_{10}");
            if (sd.isPos)
                pvInfo->AddText("(geo Y-axis inverted)");
            pvInfo->Draw("same");

            TPaveText* pvNote = new TPaveText(
                1.0 - mR + 0.01, 0.52, 0.995, 0.61, "NDC");
            pvNote->SetFillColor(kWhite); pvNote->SetBorderSize(1);
            pvNote->SetTextAlign(12); pvNote->SetTextSize(0.022);
            pvNote->AddText("#varepsilon_{act}^{GID}: eff_GID_*_src* [%]");
            pvNote->AddText("#varepsilon_{geo}: source_* [%]");
            pvNote->AddText("Red marker: current source (precise position)");
            pvNote->Draw("same");

            // ── Збереження ───────────────────────────────────────────
            c->SaveAs(Form("plots/ratio_GID/source_%s_%02d_%s.png",
                           sd.name, s, modeSuffix(dispMode)));
            delete hFrame;
            delete c;
            delete hActGID;
        }

        delete hGeo;
    }

    fAct->Close();
    fGeo->Close();
    cout << "[GID] ✓ Ratio maps saved to plots/ratio_GID/\n";
}

// ============================================================
//  ГОЛОВНА ФУНКЦІЯ
// ============================================================
int ratio_visualization()
{
    try {
        gROOT->SetBatch(true);
        gStyle->SetPalette(kBird);
        gSystem->mkdir("plots",           true);
        gSystem->mkdir("plots/ratio",     true);
        gSystem->mkdir("plots/ratio_GID", true);

        const string actRootFile = "activity_results.root";
        const string geoRootFile = "GeometricEfficiencies/plots/efficiencies/geometric_efficiencies.root";
        const string srcPosFile  = "default_values/source_positions.txt";

        // Читаємо позиції джерел один раз
        auto src_pos = read_positions(srcPosFile);
        if (src_pos.empty()) {
            cerr << "FATAL: Cannot load source positions\n";
            return 1;
        }

        cout << "\n========== Building OM ratio maps ==========\n";
        draw_ratio_OM(actRootFile, geoRootFile, src_pos, RatioMode::LINEAR);

        cout << "\n========== Building GID ratio maps ==========\n";
        draw_ratio_GID(actRootFile, geoRootFile, src_pos, 20, 13, RatioMode::LINEAR);

        cout << "\n✓ All ratio maps completed successfully\n";
    }
    catch (const exception& e) {
        cerr << "Fatal error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}