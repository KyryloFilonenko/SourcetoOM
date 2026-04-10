#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <ctime>
#include <cmath>
#include <iomanip>

struct Source {
    int    number;
    double activity;
};

bool isLeapYear(int year)
{
    return (year % 4 == 0 && year % 100 != 0) || (year % 400 == 0);
}

bool isValidDate(int year, int month, int day)
{
    if (year < 2018) return false;
    if (month < 1 || month > 12) return false;
    int daysInMonth[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    if (month == 2 && isLeapYear(year)) daysInMonth[1] = 29;
    if (day < 1 || day > daysInMonth[month - 1]) return false;
    return true;
}

long daysBetween(int year, int month, int day)
{
    std::tm inputDate = {};
    inputDate.tm_year = year  - 1900;
    inputDate.tm_mon  = month - 1;
    inputDate.tm_mday = day;

    std::tm fixedDate = {};
    fixedDate.tm_year = 2018 - 1900;
    fixedDate.tm_mon  = 6;
    fixedDate.tm_mday = 1;

    std::time_t timeInput = std::mktime(&inputDate);
    std::time_t timeFixed = std::mktime(&fixedDate);

    return std::abs(static_cast<long>(timeInput - timeFixed)) / 86400L;
}

std::vector<Source> loadSources(const std::string& filename)
{
    std::vector<Source> sources;
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: cannot open the file \"" << filename << "\"\n";
        return sources;
    }

    std::string line;
    int lineNum = 1;
    while (std::getline(file, line)) {
        if (line.empty()) { ++lineNum; continue; }
        std::istringstream ss(line);
        std::string numStr, actStr;
        if (!std::getline(ss, numStr, ';') || !std::getline(ss, actStr)) {
            std::cerr << "Warning: row " << lineNum << " skipped (unacceptable format)\n";
            ++lineNum;
            continue;
        }
        Source src;
        src.number   = std::stoi(numStr);
        src.activity = std::stod(actStr);
        sources.push_back(src);
        ++lineNum;
    }
    return sources;
}

// Стара логіка: інтегруємо від 0 до t, але тепер A0 = активність на початок вимірювань
double particlesEmitted(double A0_Bq, long t_days, double halfLife_years)
{
    if (t_days == 0) return 0.0;
    const double ln2  = std::log(2.0);
    const double hl_s = halfLife_years * 365.0 * 86400.0;
    const double t_s  = static_cast<double>(t_days) * 86400.0;
    return A0_Bq * (hl_s / ln2) * (1.0 - std::pow(2.0, -t_s / hl_s));
}





// ============================================================
//  Карти ефективності реєстрації джерело -> OM
//  Файл: plots/efficiencies/source_XX_side.png
// ============================================================

#include <TFile.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TBox.h>
#include <TLine.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TSystem.h>

// -----------------------------------------------------------
//  Структура з результатами розрахунку активності
//  (заповнюється в activity() і передається сюди)
// -----------------------------------------------------------
struct SourceResult {
    int    number;       // номер джерела (0-based index у матриці)
    double N_true;       // кількість частинок з урахуванням мертвого часу
};


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

// -----------------------------------------------------------
//  Головна функція
// -----------------------------------------------------------
// draw_efficiencies.cpp
//
// Зміни відносно попередньої версії:
//   1. Маркери джерел: червоний (поточне) / сірий (інші)
//   2. Кольорова гама: синьо-зелено-жовта (kBird або kRainBow)
//   3. Colorbar збоку з мінімальним/максимальним підписами
//   4. Позиції джерел — точні (без прив'язки до центрів OM)
//   5. Три режими відображення: PERCENT, LOG10, PPM

#include <TFile.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TBox.h>
#include <TLine.h>
#include <TLatex.h>
#include <TGraph.h>
#include <TPaveText.h>
#include <TColor.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TSystem.h>
#include <iostream>
#include <algorithm>
#include <vector>
#include <string>
#include <cmath>

enum class DisplayMode {
    PERCENT,   // відсотки (як раніше)
    LOG10,     // log10 від відсотків
    PPM        // parts per million (= eps_percent * 1e4)
};

static std::string displayLabel(double eps_percent, DisplayMode mode) {
    char buf[64];
    switch (mode) {
        case DisplayMode::PERCENT:
            snprintf(buf, sizeof(buf), "%.1f", eps_percent);
            break;
        case DisplayMode::LOG10:
            if (eps_percent > 0.0)
                snprintf(buf, sizeof(buf), "%.1f", std::log10(eps_percent / 100));
            else
                snprintf(buf, sizeof(buf), "-#infty");
            break;
        case DisplayMode::PPM:
            // 1% = 10 000 ppm
            snprintf(buf, sizeof(buf), "%.1f", eps_percent * 10000.0);
            break;
    }
    return std::string(buf);
}

static double displayValue(double eps_percent, DisplayMode mode) {
    switch (mode) {
        case DisplayMode::PERCENT: return eps_percent;
        case DisplayMode::LOG10:   return (eps_percent / 100 > 0.0) ? std::log10(eps_percent / 100) : -6.0;
        case DisplayMode::PPM:     return eps_percent * 100.0;
    }
    return eps_percent;
}

static const char* modeSuffix(DisplayMode mode) {
    switch (mode) {
        case DisplayMode::PERCENT: return "pct";
        case DisplayMode::LOG10:   return "log";
        case DisplayMode::PPM:     return "ppm";
    }
    return "pct";
}

// // ─── Побудова colorbar ──────────────────────────────────────────────────────
// // Малює вертикальну кольорову шкалу в NDC-координатах [x1,x2] x [y1,y2]
// // і підписи min/max/title.
// static void DrawColorbar(
//     double x1, double x2, double y1, double y2,
//     double vMin, double vMax,
//     DisplayMode mode,
//     int nBands = 100)
// {
//     const double dY = (y2 - y1) / nBands;

//     for (int ig = 0; ig < nBands; ig++) {
//         double frac = (double)ig / (nBands - 1);
//         // Синьо-зелено-жовта гама: синій (0) → зелений (0.5) → жовтий (1)
//         int colorIdx = TColor::GetColorPalette(
//             (int)(frac * (TColor::GetNumberOfColors() - 1)));
//         TBox* bg = new TBox(x1, y1 + ig * dY, x2, y1 + (ig + 1) * dY);
//         bg->SetFillColor(colorIdx);
//         bg->SetLineColor(colorIdx);
//         bg->Draw("same");
//     }

//     // Рамка
//     TBox* frame = new TBox(x1, y1, x2, y2);
//     frame->SetFillStyle(0);
//     frame->SetLineColor(kBlack);
//     frame->SetLineWidth(1);
//     frame->Draw("same");

//     // Підписи: max (top), середина, min (bottom)
//     const double lx = x2 + 0.004;
//     TLatex ltx;
//     ltx.SetNDC();
//     ltx.SetTextSize(0.021);
//     ltx.SetTextAlign(12);
//     ltx.SetTextColor(kBlack);

//     // Формат числових підписів шкали
//     auto fmtVal = [&](double v) -> std::string {
//         char buf[64];
//         switch (mode) {
//             case DisplayMode::PERCENT: snprintf(buf, sizeof(buf), "%.3f%%", v);      break;
//             case DisplayMode::LOG10:   snprintf(buf, sizeof(buf), "10^{%.2f}", v);   break;
//             case DisplayMode::PPM:     snprintf(buf, sizeof(buf), "%.1f", v);        break;
//         }
//         return buf;
//     };

//     ltx.DrawLatex(lx, y2,              fmtVal(vMax).c_str());
//     ltx.DrawLatex(lx, (y1 + y2) / 2,  fmtVal((vMin + vMax) / 2).c_str());
//     ltx.DrawLatex(lx, y1,              fmtVal(vMin).c_str());

//     // Заголовок шкали
//     const char* title = (mode == DisplayMode::PERCENT) ? "#varepsilon_{i,j}, %"
//                       : (mode == DisplayMode::LOG10)    ? "log_{10}(#varepsilon_{i,j})"
//                                                         : "#varepsilon_{i,j}, ppm";
//     ltx.SetTextAlign(22);
//     ltx.SetTextSize(0.023);
//     ltx.DrawLatex((x1 + x2) / 2, y2 + 0.025, title);
// }

void draw_efficiencies(
    const std::vector<SourceResult>&              results,
    const std::vector<std::pair<double,double>>&  src_pos,
    const std::vector<std::pair<double,double>>&  om_pos,
    const std::string&                            rootFile   = "vertex_energy_results.root",
    bool                                          use5x5     = false,
    DisplayMode                                   dispMode   = DisplayMode::PPM) //PERCENT, LOG10, PPM
{
    gROOT->SetBatch(true);

    gStyle->SetPalette(kBird);

    gSystem->mkdir("plots",              true);
    gSystem->mkdir("plots/efficiencies", true);

    TFile* f = TFile::Open(rootFile.c_str(), "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "draw_efficiencies: cannot open " << rootFile << "\n";
        return;
    }
    TH2F* hSourceOM_pos = (TH2F*)f->Get("SourceOMpos");
    TH2F* hSourceOM_neg = (TH2F*)f->Get("SourceOMneg");
    if (!hSourceOM_pos || !hSourceOM_neg) {
        std::cerr << "draw_efficiencies: SourceOMpos/neg not found\n";
        f->Close(); return;
    }
    hSourceOM_pos->SetDirectory(0);
    hSourceOM_neg->SetDirectory(0);
    f->Close();

    const int nSources = (int)src_pos.size();
    const int nOMs     = (int)om_pos.size();

    auto dedup = [](std::vector<double> v) -> std::vector<double> {
        std::sort(v.begin(), v.end());
        std::vector<double> out;
        out.push_back(v[0]);
        for (size_t i = 1; i < v.size(); i++)
            if (std::abs(v[i] - out.back()) > 1.0)
                out.push_back(v[i]);
        return out;
    };

    auto calc_step = [](const std::vector<double>& uv) -> double {
        if (uv.size() < 2) return 1.0;
        std::vector<double> d;
        for (size_t i = 1; i < uv.size(); i++) d.push_back(uv[i] - uv[i-1]);
        std::sort(d.begin(), d.end());
        return d[d.size()/2];
    };

    auto nearest_idx = [](const std::vector<double>& v, double val) -> int {
        int best = 0; double bd = std::abs(v[0] - val);
        for (int i = 1; i < (int)v.size(); i++) {
            double dd = std::abs(v[i] - val);
            if (dd < bd) { bd = dd; best = i; }
        }
        return best;
    };

    std::vector<double> raw_z;
    for (auto& p : om_pos) raw_z.push_back(p.second);
    std::vector<double> uZ = dedup(raw_z);
    double stepZ = calc_step(uZ);
    const int nRows = (int)uZ.size();

    std::vector<double> raw_y_neg;
    for (auto& p : om_pos) raw_y_neg.push_back(p.first);
    std::vector<double> uY_neg = dedup(raw_y_neg);
    double stepY_neg = calc_step(uY_neg);

    std::vector<double> raw_y_pos;
    for (auto& p : om_pos) raw_y_pos.push_back(-p.first);
    std::vector<double> uY_pos = dedup(raw_y_pos);
    double stepY_pos = calc_step(uY_pos);

    const int nCols = (int)uY_neg.size();

    const int colOffset_pos = 20;
    const int colOffset_neg =  0;

    struct Side {
        TH2F*                      matrix;
        const std::vector<double>* uY;
        double                     stepY;
        int                        colOffset;
        const char*                name;
        const char*                title;
        bool                       isPos;
    } sides[2] = {
        { hSourceOM_pos, &uY_pos, stepY_pos, colOffset_pos, "pos", "French side (X>0)",  true  },
        { hSourceOM_neg, &uY_neg, stepY_neg, colOffset_neg, "neg", "Italian side (X<0)", false }
    };

    for (int si = 0; si < 2; si++) {
        Side& sd = sides[si];

        const int colMin = sd.colOffset;
        const int colMax = sd.colOffset + nCols - 1;
        const int rowMin = 0;
        const int rowMax = nRows - 1;

        for (int s = 0; s < nSources; s++) {

            double N_true = 0.0;
            for (auto& r : results)
                if (r.number == s) { N_true = r.N_true; break; }
            if (N_true <= 0.0) continue;

            std::vector<int>    om_ci(nOMs, -1), om_ri(nOMs, -1);
            std::vector<double> om_eps(nOMs, 0.0);   // завжди в %

            for (int om = 0; om < nOMs; om++) {
                double hits = sd.matrix->GetBinContent(s + 1, om + 1);
                double eps  = (hits > 0) ? (hits / N_true) * 100.0 : 0.0;

                double y_search = sd.isPos ? -om_pos[om].first : om_pos[om].first;
                int ci = nearest_idx(*sd.uY, y_search);
                int ri = nearest_idx(uZ, om_pos[om].second);

                om_ci[om]  = ci;
                om_ri[om]  = ri;
                om_eps[om] = eps;
            }

            double sy_raw  = sd.isPos ? -src_pos[s].first : src_pos[s].first;
            double sz_raw  = src_pos[s].second;

            double src_col_f = 0.0, src_row_f = 0.0;
            {
                const auto& uY = *sd.uY;
                int ci0 = 0;
                for (int i = 0; i + 1 < (int)uY.size(); i++) {
                    if (sy_raw >= uY[i] && sy_raw <= uY[i+1]) { ci0 = i; break; }
                    if (sy_raw < uY[0])  { ci0 = 0; break; }
                    if (sy_raw > uY.back()) { ci0 = (int)uY.size() - 2; break; }
                }
                double alpha = (uY[ci0+1] > uY[ci0])
                    ? (sy_raw - uY[ci0]) / (uY[ci0+1] - uY[ci0]) : 0.0;
                src_col_f = ci0 + alpha;
            }
            {
                int ri0 = 0;
                for (int i = 0; i + 1 < (int)uZ.size(); i++) {
                    if (sz_raw >= uZ[i] && sz_raw <= uZ[i+1]) { ri0 = i; break; }
                    if (sz_raw < uZ[0])  { ri0 = 0; break; }
                    if (sz_raw > uZ.back()) { ri0 = (int)uZ.size() - 2; break; }
                }
                double alpha = (uZ[ri0+1] > uZ[ri0])
                    ? (sz_raw - uZ[ri0]) / (uZ[ri0+1] - uZ[ri0]) : 0.0;
                src_row_f = ri0 + alpha;
            }
            int src_ci = nearest_idx(*sd.uY, sy_raw);
            int src_ri = nearest_idx(uZ, sz_raw);

            auto in5x5 = [&](int ci, int ri) -> bool {
                if (!use5x5) return true;
                return (std::abs(ci - src_ci) <= 2 && std::abs(ri - src_ri) <= 2);
            };

            double eps_total = 0.0;
            for (int om = 0; om < nOMs; om++)
                if (om_ci[om] >= 0 && om_eps[om] > 0.0 && in5x5(om_ci[om], om_ri[om]))
                    eps_total += om_eps[om];

            double vMin =  1e30, vMax = -1e30;
            for (int om = 0; om < nOMs; om++) {
                if (om_eps[om] <= 0.0) continue;
                double v = displayValue(om_eps[om], dispMode);
                if (v < vMin) vMin = v;
                if (v > vMax) vMax = v;
            }
            if (vMin > vMax) { vMin = 0.0; vMax = 1.0; } // fallback
            if (std::abs(vMax - vMin) < 1e-12) vMax = vMin + 1.0;

            const double mL      = 0.08;
            const double mR      = 0.30;
            const double mT      = 0.09;
            const double mB      = 0.10;
            const int    baseH   = 600;
            const double areaH   = baseH * (1.0 - mT - mB);
            const double areaW   = areaH * (double)nCols / (double)nRows;
            const int    canvasW = (int)(areaW / (1.0 - mL - mR)) + 1;
            const int    canvasH = baseH;

            TCanvas* c = new TCanvas(
                Form("cEff_%d_%d", si, s), "", canvasW, canvasH);
            c->SetLeftMargin(mL);
            c->SetRightMargin(mR);
            c->SetTopMargin(mT);
            c->SetBottomMargin(mB);

            TH2F* hFrame = new TH2F(
                Form("hFrame_%d_%d", si, s), "",
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

            if (use5x5) {
                int c0 = std::max(0, src_ci - 2);
                int c1 = std::min(nCols - 1, src_ci + 2);
                int r0 = std::max(0, src_ri - 2);
                int r1 = std::min(nRows - 1, src_ri + 2);
                TBox* box5 = new TBox(
                    colMin + c0 - 0.5, r0 - 0.5,
                    colMin + c1 + 0.5, r1 + 0.5);
                box5->SetFillColorAlpha(kYellow, 0.25);
                box5->SetLineColor(kOrange + 1);
                box5->SetLineWidth(2);
                box5->Draw("same");
            }

            const double cellW = 0.46;
            const double cellH = 0.46;
            const int nColors = TColor::GetNumberOfColors();

            for (int om = 0; om < nOMs; om++) {
                int ci = om_ci[om];
                int ri = om_ri[om];
                if (ci < 0 || ri < 0 || om_eps[om] <= 0.0) continue;

                double v    = displayValue(om_eps[om], dispMode);
                double frac = (v - vMin) / (vMax - vMin);
                frac = std::max(0.0, std::min(1.0, frac));

                int colorIdx = TColor::GetColorPalette((int)(frac * (nColors - 1)));

                double bx = colMin + ci;
                double by = ri;

                TBox* b = new TBox(bx - cellW, by - cellH, bx + cellW, by + cellH);
                b->SetFillColor(colorIdx);
                b->SetLineColor(kGray + 2);
                b->SetLineWidth(1);
                b->Draw("l same");
            }

            // ── Сітка ───────────────────────────────────────────────────────
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

            // ПРОХІД 3: зелені рамки поверх усього
            for (int om = 0; om < nOMs; om++) {
                int ci = om_ci[om];
                int ri = om_ri[om];
                if (ci < 0 || ri < 0 || om_eps[om] <= 0.0) continue;
                if (!in5x5(ci, ri)) continue;

                double bx = colMin + ci;
                double by = ri;

                TBox* border = new TBox(bx - cellW, by - cellH, bx + cellW, by + cellH);
                border->SetFillStyle(0);
                border->SetLineColor(kGreen + 2);
                border->SetLineWidth(3);
                border->Draw("same");
            }

            // ПРОХІД 2: написи і маркери джерел
            for (int i = 0; i < nSources; i++) {
                double sy_i = sd.isPos ? -src_pos[i].first : src_pos[i].first;
                double sz_i = src_pos[i].second;

                double col_f = 0.0, row_f = 0.0;
                {
                    const auto& uY = *sd.uY;
                    int ci0 = 0;
                    for (int k = 0; k + 1 < (int)uY.size(); k++) {
                        if (sy_i >= uY[k] && sy_i <= uY[k+1]) { ci0 = k; break; }
                        if (sy_i < uY[0])  { ci0 = 0; break; }
                        if (sy_i > uY.back()) { ci0 = (int)uY.size() - 2; break; }
                    }
                    double alpha = (uY[ci0+1] > uY[ci0])
                        ? (sy_i - uY[ci0]) / (uY[ci0+1] - uY[ci0]) : 0.0;
                    col_f = ci0 + alpha;
                }
                {
                    int ri0 = 0;
                    for (int k = 0; k + 1 < (int)uZ.size(); k++) {
                        if (sz_i >= uZ[k] && sz_i <= uZ[k+1]) { ri0 = k; break; }
                        if (sz_i < uZ[0])  { ri0 = 0; break; }
                        if (sz_i > uZ.back()) { ri0 = (int)uZ.size() - 2; break; }
                    }
                    double alpha = (uZ[ri0+1] > uZ[ri0])
                        ? (sz_i - uZ[ri0]) / (uZ[ri0+1] - uZ[ri0]) : 0.0;
                    row_f = ri0 + alpha;
                }

                double cx = colMin + col_f;
                double cy = row_f;

                bool isCurrent = (i == s);
                int mColor = isCurrent ? kRed : kGray + 1;

                TGraph* gm = new TGraph(1, &cx, &cy);
                gm->SetMarkerStyle(isCurrent ? 20 : 20);
                gm->SetMarkerSize (isCurrent ? 1.2 : 0.0);
                gm->SetMarkerColor(mColor);
                gm->Draw("P same");

                TLatex* lbl = new TLatex(cx + 0.15, cy + 0.15, Form("%d", i));
                lbl->SetTextSize(isCurrent ? 0.04 : 0.0);
                lbl->SetTextColor(mColor);
                lbl->Draw("same");
            }

            for (int om = 0; om < nOMs; om++) {
                int ci = om_ci[om];
                int ri = om_ri[om];
                if (ci < 0 || ri < 0 || om_eps[om] <= 0.0) continue;

                double v    = displayValue(om_eps[om], dispMode);
                double frac = (v - vMin) / (vMax - vMin);
                frac = std::max(0.0, std::min(1.0, frac));

                double bx = colMin + ci;
                double by = ri;

                std::string lbl = displayLabel(om_eps[om], dispMode);
                TLatex* lt = new TLatex(bx, by, lbl.c_str());
                lt->SetTextAlign(22);
                lt->SetTextSize(0.03);
                lt->SetTextColor(frac > 0.55 ? kBlack : kWhite);
                lt->Draw("same");
            }

            TLatex ttl;
            ttl.SetNDC();
            ttl.SetTextAlign(22);
            ttl.SetTextSize(0.038);
            ttl.DrawLatex(0.5 * (1.0 - mR + mL), 0.965,
                Form("Source %d  |  %s", s, sd.title));

            const double barX1 = 1.0 - mR + 0.015;
            const double barX2 = barX1 + 0.04;    // ширина смуги
            const double barY1 = 0.12;
            const double barY2 = 0.68;

            // DrawColorbar(barX1, barX2, barY1, barY2, vMin, vMax, dispMode);

            const double infoX1 = barX1;
            const double infoX2 = 0.995;
            TPaveText* pvInfo = new TPaveText(infoX1, 0.72, infoX2, 0.91, "NDC");
            pvInfo->SetFillColor(kWhite);
            pvInfo->SetBorderSize(1);
            pvInfo->SetTextAlign(12);
            pvInfo->SetTextSize(0.026);
            pvInfo->AddText(Form("N_{true} = %.4e", N_true));
            pvInfo->AddText(Form("#varepsilon_{%d} = %.3f%%", s, eps_total));
            if (use5x5)
                pvInfo->AddText("(5#times5 window)");
            else
                pvInfo->AddText("(all active OMs)");
            const char* modeStr = (dispMode == DisplayMode::PERCENT) ? "mode: %"
                                : (dispMode == DisplayMode::LOG10)    ? "mode: log_{10}"
                                                                      : "mode: ppm";
            pvInfo->AddText(modeStr);
            pvInfo->Draw("same");

            // // Легенда маркерів (коротка)
            // const double legX1 = infoX1;
            // const double legX2 = infoX2;
            // TPaveText* pvLeg = new TPaveText(legX1, 0.91, legX2, 0.99, "NDC");
            // pvLeg->SetFillColor(kWhite);
            // pvLeg->SetBorderSize(1);
            // pvLeg->SetTextAlign(12);
            // pvLeg->SetTextSize(0.022);
            // pvLeg->AddText(Form("#color[2]{#star} Source %d (this)", s));
            // pvLeg->AddText("#color[920]{#bullet} Other sources");
            // pvLeg->Draw("same");

            c->SaveAs(Form("plots/efficiencies/source_%s_%02d_%s.png",
                           sd.name, s, modeSuffix(dispMode)));
            delete hFrame;
            delete c;
        }
    }

    std::cout << "Efficiency maps saved to plots/efficiencies/\n";
}

void activity()
{
    const std::string dataFile    = "default_values/source_activity.txt";
    const double      halfLife_yr = 31.55;
    const double      dead_time   = 0.0;

    int startYear = 2025, startMonth = 6, startDay = 10;
    int endYear   = 2025, endMonth   = 6, endDay   = 24;

    if (!isValidDate(startYear, startMonth, startDay)) {
        std::cerr << "Invalid start date.\n"; return;
    }
    if (!isValidDate(endYear, endMonth, endDay)) {
        std::cerr << "Invalid end date.\n"; return;
    }

    long t1 = daysBetween(startYear, startMonth, startDay);
    long t2 = daysBetween(endYear,   endMonth,   endDay) - t1;

    if (t2 <= 0) {
        std::cerr << "End date must be after start date.\n"; return;
    }

    double t2_s = static_cast<double>(t2) * 86400.0;

    std::vector<Source> sources = loadSources(dataFile);
    if (sources.empty()) {
        std::cerr << "No data.\n"; return;
    }

    std::cout << "\n";
    std::cout << "Reference date   : 01/07/2018\n";
    std::cout << "Start of measure : "
              << std::setfill('0') << std::setw(2) << startDay << "/"
              << std::setw(2) << startMonth << "/" << startYear
              << "  (" << t1 << " days from ref)\n";
    std::cout << "End of measure   : "
              << std::setw(2) << endDay << "/"
              << std::setw(2) << endMonth << "/" << endYear
              << "  (duration: " << t2 << " days)\n";
    std::cout << "Half-life        : " << halfLife_yr << " y.\n";
    std::cout << "Dead time        : " << std::scientific << dead_time << " s\n\n";

    const int w = 14;
    std::cout << std::left << std::setfill(' ')
              << std::setw(10) << "Source"
              << std::setw(w)  << "A0 (Bq)"
              << std::setw(w)  << "A_start (Bq)"
              << std::setw(w)  << "A_end (Bq)"
              << std::setw(w)  << "N particles"
              << "N_true\n";
    std::cout << std::string(80, '-') << "\n";

    double totalParticles = 0.0;
    double totalTrue      = 0.0;

    std::vector<SourceResult> results;

    for (int idx = 0; idx < (int)sources.size(); idx++) {
        const auto& src = sources[idx];

        double A0_Bq  = src.activity;
        double Astart = A0_Bq * std::pow(2.0,
                        -static_cast<double>(t1) / (halfLife_yr * 365.0));
        double Aend   = Astart * std::pow(2.0,
                        -static_cast<double>(t2) / (halfLife_yr * 365.0));
        double N      = particlesEmitted(Astart, t2, halfLife_yr);
        double dead_fraction = (t2_s > 0.0) ? (N * dead_time) / t2_s : 0.0;
        double N_true = N * (1.0 - dead_fraction);

        totalParticles += N;
        totalTrue      += N_true;

        std::cout << std::left  << std::fixed << std::setprecision(2)
                  << std::setw(10) << src.number
                  << std::setw(w)  << A0_Bq
                  << std::setw(w)  << Astart
                  << std::setw(w)  << Aend
                  << std::scientific << std::setprecision(4)
                  << std::setw(w)  << N
                  << N_true << "\n";

        results.push_back({ idx, N_true });
    }

    std::cout << std::string(80, '-') << "\n";
    std::cout << std::left << std::setfill(' ')
              << std::setw(10) << "Total"
              << std::setw(w)  << ""
              << std::setw(w)  << ""
              << std::setw(w)  << ""
              << std::scientific << std::setprecision(4)
              << std::setw(w)  << totalParticles
              << totalTrue     << "\n\n";

    std::vector<std::pair<double,double>> source_positions =
        read_positions_from_file("default_values/source_positions.txt");
    std::vector<std::pair<double,double>> om_positions =
        read_positions_from_file("default_values/om_positions.txt");

    std::cout << "sources loaded    : " << sources.size()         << "\n";
    std::cout << "results filled    : " << results.size()         << "\n";
    std::cout << "source_positions  : " << source_positions.size()<< "\n";
    std::cout << "om_positions      : " << om_positions.size()    << "\n";

    {
        TFile* test = TFile::Open("vertex_energy_results.root", "READ");
        if (!test || test->IsZombie()) {
            std::cerr << "ERROR: vertex_energy_results.root not found!\n"
                      << "Run visu_root() first.\n";
            if (test) test->Close();
            return;
        }
        TH2F* hp = (TH2F*)test->Get("SourceOMpos");
        TH2F* hn = (TH2F*)test->Get("SourceOMneg");
        if (!hp || !hn) {
            std::cerr << "ERROR: SourceOMpos/neg not found in ROOT file!\n";
            test->Close(); return;
        }
        std::cout << "SourceOMpos bins : "
                  << hp->GetNbinsX() << " x " << hp->GetNbinsY() << "\n";
        std::cout << "SourceOMneg bins : "
                  << hn->GetNbinsX() << " x " << hn->GetNbinsY() << "\n";
        std::cout << "SourceOMpos max  : " << hp->GetMaximum() << "\n";
        std::cout << "SourceOMneg max  : " << hn->GetMaximum() << "\n";
        test->Close();
    }

    draw_efficiencies(
        results,
        source_positions,
        om_positions,
        "vertex_energy_results.root",
        false);
}