/*
 * ratio_visualization.cpp
 *
 * Будує карти відношення:
 *
 *   R(OM, source, side) = eps_act(OM, source, side) / eps_geo(OM, source)
 *
 * де:
 *   eps_act  — реєстраційна ефективність [%], зчитується напряму з
 *              activity_results.root  →  TH2F "eff_pos_srcNN" / "eff_neg_srcNN"
 *              (записані функцією draw_efficiencies у activity.cpp)
 *
 *   eps_geo  — геометрична ефективність [%], зчитується з
 *              geometric_efficiencies.root  →  TH2F "source_N"
 *              (записані функцією total_eff_visualization)
 *
 * Відповідність осей:
 *   Обидві гістограми індексовані однаково:
 *     X-вісь = column (0-based індекс у сітці OM)
 *     Y-вісь = row    (0-based індекс у сітці OM)
 *   Для французької сторони (pos) колонки у geo дзеркаляться:
 *     geo_col = (nGeoCols - 1) - act_col
 *   Для італійської (neg) відображення пряме:
 *     geo_col = act_col
 *
 * Вхідні файли:
 *   activity_results.root                    — з activity.cpp
 *   geometric_efficiencies.root              — з total_eff_visualization.cpp
 *   default_values/source_positions.txt      — для маркерів джерел
 *
 * Вихід:
 *   plots/ratio/source_<side>_<NN>_<mode>.png
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
//  ГОЛОВНА ФУНКЦІЯ
// ============================================================
void draw_ratio(
    const string& actRootFile  = "activity_results.root",
    const string& geoRootFile  = "GeometricEfficiencies/plots/efficiencies/geometric_efficiencies.root",
    const string& srcPosFile   = "default_values/source_positions.txt",
    RatioMode     dispMode     = RatioMode::LINEAR)
{
    gROOT->SetBatch(true);
    gStyle->SetPalette(kBird);
    gSystem->mkdir("plots",       true);
    gSystem->mkdir("plots/ratio", true);

    // ── 1. Відкриваємо activity_results.root ───────────────────────────
    TFile* fAct = TFile::Open(actRootFile.c_str(), "READ");
    if (!fAct || fAct->IsZombie()) {
        cerr << "ERROR: cannot open " << actRootFile << "\n"; return;
    }

    // ── 2. Відкриваємо geometric_efficiencies.root ─────────────────────
    TFile* fGeo = TFile::Open(geoRootFile.c_str(), "READ");
    if (!fGeo || fGeo->IsZombie()) {
        cerr << "ERROR: cannot open " << geoRootFile << "\n";
        fAct->Close(); return;
    }

    // ── 3. Визначаємо кількість джерел ─────────────────────────────────
    int nSources = 0;
    while (fAct->Get(Form("eff_pos_src%02d", nSources))) nSources++;
    if (nSources == 0) {
        cerr << "ERROR: no eff_pos_src* histograms found in " << actRootFile << "\n";
        fAct->Close(); fGeo->Close(); return;
    }
    cout << "Found " << nSources << " sources in " << actRootFile << "\n";

    // ── 4. Зчитуємо N_true з TTree (для інфопанелі) ────────────────────
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
            cout << "N_true loaded from ActivityResults tree\n";
        } else {
            cerr << "WARNING: ActivityResults tree not found — N_true will show 0\n";
        }
    }

    // ── 5. Позиції джерел (для маркерів) ───────────────────────────────
    auto src_pos = read_positions(srcPosFile);

    // ── 6. Параметри сітки з першої eff-гістограми ─────────────────────
    TH2F* hFirst = (TH2F*)fAct->Get("eff_pos_src00");
    if (!hFirst) {
        cerr << "ERROR: eff_pos_src00 not found\n";
        fAct->Close(); fGeo->Close(); return;
    }
    const int nCols = hFirst->GetNbinsX();
    const int nRows = hFirst->GetNbinsY();

    // Геометрична сітка MW (фіксована)
    const int nGeoCols = 20;
    const int nGeoRows = 13;

    cout << "Act grid: " << nCols << " cols x " << nRows << " rows\n";
    cout << "Geo grid: " << nGeoCols << " cols x " << nGeoRows << " rows\n";

    // Відновлюємо координати бінів для маркерів джерел
    vector<double> uY_neg(nCols), uY_pos(nCols), uZ(nRows);
    for (int ic = 0; ic < nCols; ic++) {
        double center  = hFirst->GetXaxis()->GetBinCenter(ic + 1);
        uY_neg[ic] =  center;
        uY_pos[ic] = -center;
    }
    sort(uY_pos.begin(), uY_pos.end());
    for (int ir = 0; ir < nRows; ir++)
        uZ[ir] = hFirst->GetYaxis()->GetBinCenter(ir + 1);

    // ── 7. Опис сторін ──────────────────────────────────────────────────
    struct Side {
        const char*         name;
        const char*         title;
        bool                isPos;   // true = French (X>0), false = Italian (X<0)
        const vector<double>* uY;
    } sides[2] = {
        { "pos", "French side (X>0)",  true,  &uY_pos },
        { "neg", "Italian side (X<0)", false, &uY_neg }
    };

    auto nearest_idx = [](const vector<double>& v, double val) -> int {
        int best = 0; double bd = fabs(v[0] - val);
        for (int i = 1; i < (int)v.size(); i++) {
            double dd = fabs(v[i] - val);
            if (dd < bd) { bd = dd; best = i; }
        }
        return best;
    };

    // ── 8. Основний цикл по джерелах та сторонах ────────────────────────
    for (int s = 0; s < nSources; s++) {

        // Геометрична ефективність (однакова для обох сторін)
        TH2F* hGeo = (TH2F*)fGeo->Get(Form("source_%d", s));
        if (!hGeo) {
            cerr << "WARNING: source_" << s << " not found in geo file, skipping\n";
            continue;
        }
        hGeo->SetDirectory(0);

        for (int si = 0; si < 2; si++) {
            Side& sd = sides[si];

            // Реєстраційна ефективність для цієї сторони і джерела
            TH2F* hAct = (TH2F*)fAct->Get(
                Form("eff_%s_src%02d", sd.name, s));
            if (!hAct) {
                cerr << "WARNING: eff_" << sd.name << "_src"
                     << Form("%02d", s) << " not found, skipping\n";
                continue;
            }
            hAct->SetDirectory(0);

            // ── Обчислення відношень по всіх комірках ────────────────
            struct CellData {
                int    ci, ri;
                double eps_act;
                double eps_geo;
                double ratio;
            };
            vector<CellData> cells;
            cells.reserve(nCols * nRows);

            for (int ci = 0; ci < nCols; ci++) {
                for (int ri = 0; ri < nRows; ri++) {

                    // eps_act вже збережена у % у гістограмі
                    double eps_act = hAct->GetBinContent(ci + 1, ri + 1);
                    if (eps_act <= 0.0) continue;

                    // Відповідна колонка у geo:
                    //   pos (французька): дзеркало по Y
                    //   neg (італійська): прямий порядок
                    int geo_col = sd.isPos ? (nGeoCols - 1 - ci) : ci;
                    int geo_row = ri;

                    if (geo_col < 0 || geo_col >= nGeoCols ||
                        geo_row < 0 || geo_row >= nGeoRows) continue;

                    // eps_geo теж у %
                    double eps_geo = hGeo->GetBinContent(geo_col + 1, geo_row + 1);
                    if (eps_geo <= 0.0) continue;

                    cells.push_back({ci, ri, eps_act, eps_geo, eps_act / eps_geo});
                }
            }

            if (cells.empty()) {
                cerr << "WARNING: no valid cells for source " << s
                     << " side " << sd.name << "\n";
                delete hAct; continue;
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
                Form("cRatio_%s_%d", sd.name, s), "", canvasW, baseH);
            c->SetLeftMargin(mL); c->SetRightMargin(mR);
            c->SetTopMargin(mT);  c->SetBottomMargin(mB);

            TH2F* hFrame = new TH2F(
                Form("hFrame_%s_%d", sd.name, s), "",
                nCols, colMin - 0.5, colMax + 0.5,
                nRows, rowMin - 0.5, rowMax + 0.5);
            hFrame->SetDirectory(0); hFrame->SetStats(0);
            hFrame->GetXaxis()->SetNdivisions(nCols, false);
            hFrame->GetYaxis()->SetNdivisions(nRows, false);
            hFrame->GetXaxis()->SetLabelSize(0.048);
            hFrame->GetYaxis()->SetLabelSize(0.048);
            hFrame->GetXaxis()->SetTickLength(0.01);
            hFrame->GetYaxis()->SetTickLength(0.01);
            for (int ic = 0; ic < nCols; ic++)
                hFrame->GetXaxis()->SetBinLabel(ic + 1, Form("%d", ic));
            for (int ir = 0; ir < nRows; ir++)
                hFrame->GetYaxis()->SetBinLabel(ir + 1, Form("%d", ir));
            hFrame->GetXaxis()->SetTitle("Column");
            hFrame->GetYaxis()->SetTitle("Row");
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

            // ПРОХІД 4: маркери джерел
            for (int i = 0; i < nSources && i < (int)src_pos.size(); i++) {
                // Для pos: шукаємо по -Y (дзеркало), для neg: по +Y
                double sy_i = sd.isPos ? -src_pos[i].first : src_pos[i].first;
                double sz_i = src_pos[i].second;

                const auto& uY = *sd.uY;
                int ci0 = 0;
                for (int k = 0; k + 1 < (int)uY.size(); k++) {
                    if (sy_i >= uY[k] && sy_i <= uY[k+1]) { ci0 = k; break; }
                    if (sy_i <  uY.front()) { ci0 = 0; break; }
                    if (sy_i >  uY.back())  { ci0 = (int)uY.size() - 2; break; }
                }
                double ay = (uY.size() > 1 && fabs(uY[ci0+1] - uY[ci0]) > 1e-9)
                    ? (sy_i - uY[ci0]) / (uY[ci0+1] - uY[ci0]) : 0.0;
                double col_f = ci0 + ay;

                int ri0 = 0;
                for (int k = 0; k + 1 < (int)uZ.size(); k++) {
                    if (sz_i >= uZ[k] && sz_i <= uZ[k+1]) { ri0 = k; break; }
                    if (sz_i <  uZ.front()) { ri0 = 0; break; }
                    if (sz_i >  uZ.back())  { ri0 = (int)uZ.size() - 2; break; }
                }
                double az = (uZ.size() > 1 && fabs(uZ[ri0+1] - uZ[ri0]) > 1e-9)
                    ? (sz_i - uZ[ri0]) / (uZ[ri0+1] - uZ[ri0]) : 0.0;
                double row_f = ri0 + az;

                double cx = col_f, cy = row_f;
                bool isCurr = (i == s);
                int mColor  = isCurr ? kRed : kGray + 1;

                TGraph* gm = new TGraph(1, &cx, &cy);
                gm->SetMarkerStyle(20);
                gm->SetMarkerSize(isCurr ? 1.2 : 0.0);
                gm->SetMarkerColor(mColor);
                gm->Draw("P same");

                if (isCurr) {
                    TLatex* lbl = new TLatex(cx + 0.15, cy + 0.15, Form("%d", i));
                    lbl->SetTextSize(0.04); lbl->SetTextColor(mColor);
                    lbl->Draw("same");
                }
            }

            // ── Заголовок ────────────────────────────────────────────
            TLatex ttl;
            ttl.SetNDC(); ttl.SetTextAlign(22); ttl.SetTextSize(0.038);
            ttl.DrawLatex(0.5 * (1.0 - mR + mL), 0.965,
                Form("Source %d  |  #varepsilon_{act}/#varepsilon_{geo}  |  %s",
                     s, sd.title));

            // ── Інфопанель ───────────────────────────────────────────
            double sum_act = 0.0, sum_geo = 0.0, mean_ratio = 0.0;
            int n_valid = 0;
            for (auto& cd : cells) {
                if (cd.ratio > 0.0) {
                    sum_act    += cd.eps_act;
                    sum_geo    += cd.eps_geo;
                    mean_ratio += cd.ratio;
                    n_valid++;
                }
            }
            if (n_valid > 0) mean_ratio /= n_valid;

            TPaveText* pvInfo = new TPaveText(
                1.0 - mR + 0.01, 0.62, 0.995, 0.90, "NDC");
            pvInfo->SetFillColor(kWhite); pvInfo->SetBorderSize(1);
            pvInfo->SetTextAlign(12); pvInfo->SetTextSize(0.025);
            pvInfo->AddText(Form("N_{true} = %.4e", N_true_vec[s]));
            pvInfo->AddText(Form("#Sigma #varepsilon_{act} = %.3f%%", sum_act));
            pvInfo->AddText(Form("#Sigma #varepsilon_{geo} = %.3f%%", sum_geo));
            pvInfo->AddText(Form("<R> = %.4f", mean_ratio));
            pvInfo->AddText(
                (dispMode == RatioMode::LINEAR) ? "mode: linear" : "mode: log_{10}");
            if (sd.isPos)
                pvInfo->AddText("(geo Y-axis inverted)");
            pvInfo->Draw("same");

            TPaveText* pvNote = new TPaveText(
                1.0 - mR + 0.01, 0.52, 0.995, 0.61, "NDC");
            pvNote->SetFillColor(kWhite); pvNote->SetBorderSize(1);
            pvNote->SetTextAlign(12); pvNote->SetTextSize(0.022);
            pvNote->AddText("#varepsilon_{act}: eff_*_src* [%]");
            pvNote->AddText("#varepsilon_{geo}: source_* [%]");
            pvNote->Draw("same");

            // ── Збереження ───────────────────────────────────────────
            c->SaveAs(Form("plots/ratio/source_%s_%02d_%s.png",
                           sd.name, s, modeSuffix(dispMode)));
            delete hFrame;
            delete c;
            delete hAct;
        }  // end side loop

        delete hGeo;

        if (s % 5 == 0)
            cout << "Done source " << s << "/" << (nSources - 1) << "\n";
    }  // end source loop

    fAct->Close();
    fGeo->Close();
    cout << "\n✓ Ratio maps saved to plots/ratio/\n";
}

// ============================================================
//  ТОЧКА ВХОДУ
// ============================================================
int ratio_visualization()
{
    try {
        const string actRootFile = "activity_results.root";
        const string geoRootFile = "GeometricEfficiencies/plots/efficiencies/geometric_efficiencies.root";
        const string srcPosFile  = "default_values/source_positions.txt";

        draw_ratio(actRootFile, geoRootFile, srcPosFile, RatioMode::LINEAR);
        //draw_ratio(actRootFile, geoRootFile, srcPosFile, RatioMode::LOG10);
    }
    catch (const exception& e) {
        cerr << "Fatal error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}