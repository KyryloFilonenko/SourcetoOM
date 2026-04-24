#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <ctime>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <string>

#include <TFile.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TBox.h>
#include <TLine.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TColor.h>
#include <TROOT.h>
#include <TTree.h>
#include <TNamed.h>

// ============================================================
//  Структури
// ============================================================

struct Source {
    int    number;
    double activity;
};

struct SourceResult {
    int    number;   // 0-based index
    double N_true;   // кількість частинок з урахуванням мертвого часу
};

// ============================================================
//  Дати і активність
// ============================================================

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
            ++lineNum; continue;
        }
        Source src;
        src.number   = std::stoi(numStr);
        src.activity = std::stod(actStr);
        sources.push_back(src);
        ++lineNum;
    }
    return sources;
}

double particlesEmitted(double A0_Bq, long t_days, double halfLife_years)
{
    if (t_days == 0) return 0.0;
    const double ln2  = std::log(2.0);
    const double hl_s = halfLife_years * 365.0 * 86400.0;
    const double t_s  = static_cast<double>(t_days) * 86400.0;
    return A0_Bq * (hl_s / ln2) * (1.0 - std::pow(2.0, -t_s / hl_s));
}

// ============================================================
//  Читання позицій
// ============================================================

std::vector<std::pair<double,double>>
read_positions_from_file(const std::string& filename)
{
    std::vector<std::pair<double,double>> positions;
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        std::cerr << "ERROR: cannot open file " << filename << std::endl;
        return positions;
    }
    std::string line;
    while (std::getline(infile, line)) {
        if (line.empty()) continue;
        double y, z;
        if (line.find(';') != std::string::npos) {
            std::replace(line.begin(), line.end(), ';', ' ');
            std::istringstream ss(line);
            if (!(ss >> y >> z)) continue;
        } else {
            std::istringstream ss(line);
            int index;
            if (!(ss >> index >> y >> z)) continue;
        }
        positions.emplace_back(y, z);
    }
    infile.close();
    std::cout << "Loaded " << positions.size() << " positions from " << filename << "\n";
    return positions;
}

// ============================================================
//  Режими відображення (для draw_efficiencies)
// ============================================================

enum class DisplayMode {
    PERCENT,
    LOG10,
    PPM
};

static std::string displayLabel(double eps_percent, DisplayMode mode) {
    char buf[64];
    switch (mode) {
        case DisplayMode::PERCENT:
            snprintf(buf, sizeof(buf), "%.1f", eps_percent); break;
        case DisplayMode::LOG10:
            if (eps_percent > 0.0)
                snprintf(buf, sizeof(buf), "%.1f", std::log10(eps_percent / 100));
            else
                snprintf(buf, sizeof(buf), "-#infty");
            break;
        case DisplayMode::PPM:
            snprintf(buf, sizeof(buf), "%.1f", eps_percent * 10000.0); break;
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

// ============================================================
//  draw_efficiencies
//  — малює карти ефективностей реєстрації
//  — НОВІ параметри: effRootFile, effRootMode
//      effRootFile  — ім'я ROOT-файлу для збереження ефективностей
//                     (порожній рядок = не зберігати)
//      effRootMode  — "RECREATE" або "UPDATE"
// ============================================================

void draw_efficiencies(
    const std::vector<SourceResult>&              results,
    const std::vector<std::pair<double,double>>&  src_pos,
    const std::vector<std::pair<double,double>>&  om_pos,
    const std::string&                            rootFile        = "vertex_energy_results.root",
    bool                                          use5x5          = false,
    DisplayMode                                   dispMode        = DisplayMode::LOG10,
    double                                        electron_intens = 1.0,
    const std::string&                            effRootFile     = "activity_results.root",
    const std::string&                            effRootMode     = "UPDATE")
{
    gROOT->SetBatch(true);
    gStyle->SetPalette(kBird);
    gSystem->mkdir("plots",              true);
    gSystem->mkdir("plots/efficiencies", true);

    // ── Відкриваємо вхідний ROOT-файл з hits ───────────────────────────
    TFile* f = TFile::Open(rootFile.c_str(), "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "draw_efficiencies: cannot open " << rootFile << "\n"; return;
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

    // ── Відкриваємо вихідний ROOT-файл для ефективностей ───────────────
    //    Відкриваємо один раз тут, щоб записати всі гістограми разом.
    TFile* fEff = nullptr;
    if (!effRootFile.empty()) {
        fEff = TFile::Open(effRootFile.c_str(), effRootMode.c_str());
        if (!fEff || fEff->IsZombie()) {
            std::cerr << "WARNING: cannot open " << effRootFile
                      << " for writing efficiencies — skipping ROOT output\n";
            fEff = nullptr;
        } else {
            std::cout << "Efficiency ROOT output: " << effRootFile
                      << " (mode=" << effRootMode << ")\n";
        }
    }

    // ── Допоміжні лямбди ───────────────────────────────────────────────
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

    // ── Побудова сіток координат ────────────────────────────────────────
    std::vector<double> raw_z;
    for (auto& p : om_pos) raw_z.push_back(p.second);
    std::vector<double> uZ = dedup(raw_z);
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

    struct Side {
        TH2F*                      matrix;
        const std::vector<double>* uY;
        double                     stepY;
        int                        colOffset;
        const char*                name;
        const char*                title;
        bool                       isPos;
    } sides[2] = {
        { hSourceOM_pos, &uY_pos, stepY_pos, 20, "pos", "French side (X>0)",  true  },
        { hSourceOM_neg, &uY_neg, stepY_neg,  0, "neg", "Italian side (X<0)", false }
    };

    // ── Головний цикл ───────────────────────────────────────────────────
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

            // Обчислюємо ефективності
            std::vector<int>    om_ci(nOMs, -1), om_ri(nOMs, -1);
            std::vector<double> om_eps(nOMs, 0.0);  // у %

            for (int om = 0; om < nOMs; om++) {
                double hits = sd.matrix->GetBinContent(s + 1, om + 1);
                double eps  = (hits > 0) ? (hits / (N_true * electron_intens)) * 100.0 : 0.0;

                double y_search = sd.isPos ? -om_pos[om].first : om_pos[om].first;
                int ci = nearest_idx(*sd.uY, y_search);
                int ri = nearest_idx(uZ, om_pos[om].second);

                om_ci[om]  = ci;
                om_ri[om]  = ri;
                om_eps[om] = eps;
            }

            // ── ЗБЕРЕЖЕННЯ ЕФЕКТИВНОСТЕЙ У ROOT-ФАЙЛ ───────────────────
            //
            //    Назва гістограми: "eff_<side>_src<NN>"
            //    Осі: X = column index у сітці OM (0..nCols-1)
            //         Y = row index у сітці OM   (0..nRows-1)
            //    Вміст bin: eps_registration [%]
            //
            if (fEff) {
                fEff->cd();
                TH2F* hEff = new TH2F(
                    Form("eff_%s_src%02d", sd.name, s),
                    Form("Registration efficiency; %s; source %d; col; row",
                         sd.title, s),
                    nCols, colMin - 0.5, colMax + 0.5,
                    nRows, rowMin - 0.5, rowMax + 0.5);
                hEff->SetDirectory(fEff);

                for (int om = 0; om < nOMs; om++) {
                    int ci = om_ci[om];
                    int ri = om_ri[om];
                    if (ci < 0 || ri < 0 || om_eps[om] <= 0.0) continue;
                    // bin-індекси ROOT: 1-based
                    hEff->SetBinContent(ci + 1, ri + 1, om_eps[om]);
                }
                hEff->Write("", TObject::kOverwrite);
                // не delete — файл забере ownership через SetDirectory(fEff)
            }

            // ── ВІЗУАЛІЗАЦІЯ (без змін відносно оригіналу) ─────────────

            double sy_raw = sd.isPos ? -src_pos[s].first : src_pos[s].first;
            double sz_raw = src_pos[s].second;

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

            double vMin = 1e30, vMax = -1e30;
            for (int om = 0; om < nOMs; om++) {
                if (om_eps[om] <= 0.0) continue;
                double v = displayValue(om_eps[om], dispMode);
                if (v < vMin) vMin = v;
                if (v > vMax) vMax = v;
            }
            if (vMin > vMax) { vMin = 0.0; vMax = 1.0; }
            if (std::abs(vMax - vMin) < 1e-12) vMax = vMin + 1.0;

            const double mL = 0.08, mR = 0.30, mT = 0.09, mB = 0.10;
            const int baseH = 600;
            const double areaH = baseH * (1.0 - mT - mB);
            const double areaW = areaH * (double)nCols / (double)nRows;
            const int canvasW  = (int)(areaW / (1.0 - mL - mR)) + 1;

            TCanvas* c = new TCanvas(Form("cEff_%d_%d", si, s), "", canvasW, baseH);
            c->SetLeftMargin(mL); c->SetRightMargin(mR);
            c->SetTopMargin(mT); c->SetBottomMargin(mB);

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
                int c0 = std::max(0, src_ci - 2), c1 = std::min(nCols-1, src_ci+2);
                int r0 = std::max(0, src_ri - 2), r1 = std::min(nRows-1, src_ri+2);
                TBox* box5 = new TBox(colMin+c0-0.5, r0-0.5, colMin+c1+0.5, r1+0.5);
                box5->SetFillColorAlpha(kYellow, 0.25);
                box5->SetLineColor(kOrange+1); box5->SetLineWidth(2);
                box5->Draw("same");
            }

            const double cellW = 0.46, cellH = 0.46;
            const int nColors = TColor::GetNumberOfColors();

            // ПРОХІД 1: кольорові комірки
            for (int om = 0; om < nOMs; om++) {
                int ci = om_ci[om], ri = om_ri[om];
                if (ci < 0 || ri < 0 || om_eps[om] <= 0.0) continue;

                double v    = displayValue(om_eps[om], dispMode);
                double frac = (v - vMin) / (vMax - vMin);
                frac = std::max(0.0, std::min(1.0, frac));
                int colorIdx = TColor::GetColorPalette((int)(frac * (nColors - 1)));

                TBox* b = new TBox(
                    colMin+ci-cellW, ri-cellH, colMin+ci+cellW, ri+cellH);
                b->SetFillColor(colorIdx);
                b->SetLineColor(kGray+2); b->SetLineWidth(1);
                b->Draw("l same");
            }

            // ПРОХІД 2: сітка
            for (int ic = 0; ic <= nCols; ic++) {
                TLine* l = new TLine(colMin-0.5+ic, rowMin-0.5,
                                     colMin-0.5+ic, rowMax+0.5);
                l->SetLineColor(kBlack); l->SetLineWidth(1); l->SetLineStyle(2);
                l->Draw("same");
            }
            for (int ir = 0; ir <= nRows; ir++) {
                TLine* l = new TLine(colMin-0.5, rowMin-0.5+ir,
                                     colMax+0.5, rowMin-0.5+ir);
                l->SetLineColor(kBlack); l->SetLineWidth(1); l->SetLineStyle(2);
                l->Draw("same");
            }

            // ПРОХІД 3: зелені рамки (5x5)
            for (int om = 0; om < nOMs; om++) {
                int ci = om_ci[om], ri = om_ri[om];
                if (ci < 0 || ri < 0 || om_eps[om] <= 0.0) continue;
                if (!in5x5(ci, ri)) continue;
                TBox* border = new TBox(
                    colMin+ci-cellW, ri-cellH, colMin+ci+cellW, ri+cellH);
                border->SetFillStyle(0);
                border->SetLineColor(kGreen+2); border->SetLineWidth(3);
                border->Draw("same");
            }

            // ПРОХІД 4: маркери джерел
            for (int i = 0; i < nSources; i++) {
                double sy_i = sd.isPos ? -src_pos[i].first : src_pos[i].first;
                double sz_i = src_pos[i].second;

                // інтерполяція позиції
                const auto& uY = *sd.uY;
                int ci0 = 0;
                for (int k = 0; k+1 < (int)uY.size(); k++) {
                    if (sy_i >= uY[k] && sy_i <= uY[k+1]) { ci0 = k; break; }
                    if (sy_i < uY[0])    { ci0 = 0; break; }
                    if (sy_i > uY.back()) { ci0 = (int)uY.size()-2; break; }
                }
                double alpha = (uY[ci0+1] > uY[ci0])
                    ? (sy_i - uY[ci0]) / (uY[ci0+1] - uY[ci0]) : 0.0;
                double col_f = ci0 + alpha;

                int ri0 = 0;
                for (int k = 0; k+1 < (int)uZ.size(); k++) {
                    if (sz_i >= uZ[k] && sz_i <= uZ[k+1]) { ri0 = k; break; }
                    if (sz_i < uZ[0])    { ri0 = 0; break; }
                    if (sz_i > uZ.back()) { ri0 = (int)uZ.size()-2; break; }
                }
                double alpha_z = (uZ[ri0+1] > uZ[ri0])
                    ? (sz_i - uZ[ri0]) / (uZ[ri0+1] - uZ[ri0]) : 0.0;
                double row_f = ri0 + alpha_z;

                double cx = colMin + col_f, cy = row_f;
                bool isCurrent = (i == s);
                int mColor = isCurrent ? kRed : kGray+1;

                TGraph* gm = new TGraph(1, &cx, &cy);
                gm->SetMarkerStyle(20);
                gm->SetMarkerSize(isCurrent ? 1.2 : 0.0);
                gm->SetMarkerColor(mColor);
                gm->Draw("P same");

                TLatex* lbl = new TLatex(cx+0.15, cy+0.15, Form("%d", i));
                lbl->SetTextSize(isCurrent ? 0.04 : 0.0);
                lbl->SetTextColor(mColor);
                lbl->Draw("same");
            }

            // ПРОХІД 5: числові підписи
            for (int om = 0; om < nOMs; om++) {
                int ci = om_ci[om], ri = om_ri[om];
                if (ci < 0 || ri < 0 || om_eps[om] <= 0.0) continue;

                double v    = displayValue(om_eps[om], dispMode);
                double frac = (v - vMin) / (vMax - vMin);
                frac = std::max(0.0, std::min(1.0, frac));

                std::string lbl = displayLabel(om_eps[om], dispMode);
                TLatex* lt = new TLatex(colMin+ci, ri, lbl.c_str());
                lt->SetTextAlign(22); lt->SetTextSize(0.03);
                lt->SetTextColor(frac > 0.55 ? kBlack : kWhite);
                lt->Draw("same");
            }

            // Заголовок
            TLatex ttl;
            ttl.SetNDC(); ttl.SetTextAlign(22); ttl.SetTextSize(0.038);
            ttl.DrawLatex(0.5*(1.0-mR+mL), 0.965,
                Form("Source %d  |  %s", s, sd.title));

            // Інфопанель
            const double infoX1 = 1.0 - mR + 0.015;
            TPaveText* pvInfo = new TPaveText(infoX1, 0.72, 0.995, 0.91, "NDC");
            pvInfo->SetFillColor(kWhite); pvInfo->SetBorderSize(1);
            pvInfo->SetTextAlign(12); pvInfo->SetTextSize(0.026);
            pvInfo->AddText(Form("N_{true} = %.4e", N_true));
            pvInfo->AddText(Form("#varepsilon_{%d} = %.3f%%", s, eps_total));
            pvInfo->AddText(use5x5 ? "(5#times5 window)" : "(all active OMs)");
            pvInfo->AddText(
                (dispMode == DisplayMode::PERCENT) ? "mode: %" :
                (dispMode == DisplayMode::LOG10)   ? "mode: log_{10}" : "mode: ppm");
            pvInfo->Draw("same");

            c->SaveAs(Form("plots/efficiencies/source_%s_%02d_%s.png",
                           sd.name, s, modeSuffix(dispMode)));
            delete hFrame;
            delete c;
        }  // end source loop
    }  // end side loop

    // Закриваємо eff ROOT-файл
    if (fEff) {
        fEff->Close();
        std::cout << "✓ Registration efficiencies saved to: " << effRootFile << "\n";
        std::cout << "  Histograms: eff_pos_src<NN> and eff_neg_src<NN>\n";
        std::cout << "  Content unit: registration efficiency [%]\n";
    }

    std::cout << "Efficiency maps saved to plots/efficiencies/\n";
}

// ============================================================
//  Збереження результатів активності у ROOT файл
// ============================================================

void saveActivityResultsToRoot(
    const std::string&              outputFile,
    const std::vector<SourceResult>& results,
    const std::vector<Source>&      sources,
    int startYear, int startMonth, int startDay,
    int endYear,   int endMonth,   int endDay,
    double halfLife_yr,
    double dead_time,
    double totalParticles,
    double totalTrue,
    double electron_intens)
{
    // Відкриваємо в UPDATE — щоб не затерти eff_*_src* гістограми,
    // які вже записала draw_efficiencies вище
    TFile* outfile = TFile::Open(outputFile.c_str(), "UPDATE");
    if (!outfile || outfile->IsZombie()) {
        std::cerr << "ERROR: cannot open output file " << outputFile << "\n"; return;
    }

    // ── TTree з результатами для кожного джерела ───────────────────────
    TTree* tResults = new TTree("ActivityResults", "Activity calculation results");

    int    source_idx;
    double A0, A_start, A_end, N_particles, N_true_br;
    long   t1_days, t2_days;

    tResults->Branch("source_idx",  &source_idx,  "source_idx/I");
    tResults->Branch("A0",          &A0,          "A0/D");
    tResults->Branch("A_start",     &A_start,     "A_start/D");
    tResults->Branch("A_end",       &A_end,       "A_end/D");
    tResults->Branch("N_particles", &N_particles, "N_particles/D");
    tResults->Branch("N_true",      &N_true_br,   "N_true/D");
    tResults->Branch("t1_days",     &t1_days,     "t1_days/L");
    tResults->Branch("t2_days",     &t2_days,     "t2_days/L");

    t1_days = daysBetween(startYear, startMonth, startDay);
    t2_days = daysBetween(endYear,   endMonth,   endDay) - t1_days;

    for (size_t i = 0; i < sources.size(); i++) {
        const auto& src = sources[i];
        source_idx = src.number;
        A0         = src.activity;
        A_start    = A0    * std::pow(2.0, -static_cast<double>(t1_days) / (halfLife_yr * 365.0));
        A_end      = A_start * std::pow(2.0, -static_cast<double>(t2_days) / (halfLife_yr * 365.0));

        N_true_br   = 0.0;
        N_particles = 0.0;
        for (const auto& r : results) {
            if (r.number == (int)i) {
                N_true_br = r.N_true;
                double t2_s = static_cast<double>(t2_days) * 86400.0;
                double df   = (t2_s > 0.0) ? (N_true_br * dead_time) / t2_s : 0.0;
                N_particles = N_true_br * electron_intens / (1.0 - df + 1e-15);
                break;
            }
        }
        tResults->Fill();
    }

    // ── Метаінформація ─────────────────────────────────────────────────
    TNamed* nHalfLife  = new TNamed("HalfLife_years", Form("%f", halfLife_yr));
    TNamed* nDeadTime  = new TNamed("DeadTime_s",     Form("%e", dead_time));
    TNamed* nStartDate = new TNamed("StartDate",
        Form("%02d/%02d/%d", startDay, startMonth, startYear));
    TNamed* nEndDate   = new TNamed("EndDate",
        Form("%02d/%02d/%d", endDay,   endMonth,   endYear));
    TNamed* nRefDate   = new TNamed("ReferenceDate", "01/07/2018");
    TNamed* nElInt     = new TNamed("ElectronIntensity", Form("%e", electron_intens));
    TNamed* nNote      = new TNamed("EffNote",
        "eff_pos_srcNN / eff_neg_srcNN = registration efficiency [%] per OM cell");

    nHalfLife ->Write("", TObject::kOverwrite);
    nDeadTime ->Write("", TObject::kOverwrite);
    nStartDate->Write("", TObject::kOverwrite);
    nEndDate  ->Write("", TObject::kOverwrite);
    nRefDate  ->Write("", TObject::kOverwrite);
    nElInt    ->Write("", TObject::kOverwrite);
    nNote     ->Write("", TObject::kOverwrite);

    // ── Підсумкова гістограма ──────────────────────────────────────────
    TH1F* hTotals = new TH1F("Totals", "Total particles and true counts", 2, 0, 2);
    hTotals->SetBinContent(1, totalParticles);
    hTotals->SetBinContent(2, totalTrue);
    hTotals->GetXaxis()->SetBinLabel(1, "Total N particles");
    hTotals->GetXaxis()->SetBinLabel(2, "Total N true");
    hTotals->Write("", TObject::kOverwrite);

    tResults->Write("", TObject::kOverwrite);
    outfile->Close();

    std::cout << "\n✓ Activity results saved to " << outputFile << "\n";
    std::cout << "  Contents:\n";
    std::cout << "    TTree  'ActivityResults'     — per-source decay data\n";
    std::cout << "    TH2F   'eff_pos_srcNN'       — registration efficiency, French side [%]\n";
    std::cout << "    TH2F   'eff_neg_srcNN'       — registration efficiency, Italian side [%]\n";
    std::cout << "    TNamed 'HalfLife_years' etc. — calculation parameters\n";
    std::cout << "    TH1F   'Totals'              — total N_particles / N_true\n";
}

// ============================================================
//  ГОЛОВНА ФУНКЦІЯ
// ============================================================

void activity()
{
    const std::string dataFile        = "default_values/source_activity.txt";
    const std::string outputRootFile  = "activity_results.root";
    const double      halfLife_yr     = 31.55;
    const double      dead_time       = 0.0;
    const double      electron_intens = 0.1145;

    const long long t_run_first_start = 1750493319LL;
    const long long t_run_last_start  = 1751285915LL;
    const long long last_run_duration = 1800LL;
    const double    total_duration_s  = 650575.0;

    const long long t_start_unix = t_run_first_start;
    const long long t_end_unix   = t_run_last_start + last_run_duration;

    auto unixToYMD = [](long long ut, int& y, int& m, int& d) {
        std::time_t t = (std::time_t)ut;
        std::tm* ti = std::gmtime(&t);
        y = ti->tm_year + 1900;
        m = ti->tm_mon  + 1;
        d = ti->tm_mday;
    };

    auto printTime = [](const char* label, long long ut) {
        std::time_t t = (std::time_t)ut;
        std::tm* ti = std::gmtime(&t);
        char buf[64];
        std::strftime(buf, sizeof(buf), "%d/%m/%Y %H:%M:%S UTC", ti);
        std::cout << label << buf;
    };

    int startYear, startMonth, startDay;
    int endYear,   endMonth,   endDay;
    unixToYMD(t_start_unix, startYear, startMonth, startDay);
    unixToYMD(t_end_unix,   endYear,   endMonth,   endDay);

    if (!isValidDate(startYear, startMonth, startDay)) {
        std::cerr << "Invalid start date.\n"; return;
    }
    if (!isValidDate(endYear, endMonth, endDay)) {
        std::cerr << "Invalid end date.\n"; return;
    }

    long t1 = daysBetween(startYear, startMonth, startDay);
    long t2 = daysBetween(endYear, endMonth, endDay) - t1;
    if (t2 <= 0) { std::cerr << "End date must be after start date.\n"; return; }

    const double t2_s = total_duration_s;

    std::vector<Source> sources = loadSources(dataFile);
    if (sources.empty()) { std::cerr << "No data.\n"; return; }

    std::cout << "\n";
    std::cout << "Reference date    : 01/07/2018\n";
    printTime("Start of measure  : ", t_start_unix);
    std::cout << "  (run 2204, unix=" << t_start_unix << ", " << t1 << " days from ref)\n";
    printTime("End of measure    : ", t_end_unix);
    std::cout << "  (run 2570 end,  unix=" << t_end_unix << ")\n";
    std::cout << "Data duration     : " << std::fixed << std::setprecision(0)
              << total_duration_s << " s\n";
    std::cout << "Half-life         : " << halfLife_yr << " y.\n";
    std::cout << "Dead time         : " << std::scientific << dead_time << " s\n";
    std::cout << "Electron intensity: " << electron_intens << "\n\n";

    const int w = 14;
    std::cout << std::left << std::setfill(' ')
              << std::setw(10) << "Source"
              << std::setw(w)  << "A0 (Bq)"
              << std::setw(w)  << "A_start (Bq)"
              << std::setw(w)  << "A_end (Bq)"
              << std::setw(w)  << "N particles"
              << "N_true\n";
    std::cout << std::string(80, '-') << "\n";

    double totalParticles = 0.0, totalTrue = 0.0;
    std::vector<SourceResult> results;

    for (int idx = 0; idx < (int)sources.size(); idx++) {
        const auto& src = sources[idx];
        double A0_Bq  = src.activity;
        double Astart = A0_Bq * std::pow(2.0, -(double)t1 / (halfLife_yr * 365.0));
        double Aend   = Astart * std::pow(2.0, -(double)t2 / (halfLife_yr * 365.0));
        double N      = particlesEmitted(Astart, t2, halfLife_yr);
        double df     = (t2_s > 0.0) ? (N * dead_time) / t2_s : 0.0;
        double N_true = N * (1.0 - df);

        totalParticles += N;
        totalTrue      += N_true;

        std::cout << std::left << std::fixed << std::setprecision(2)
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
    std::cout << std::left << std::setw(10) << "Total"
              << std::setw(w) << "" << std::setw(w) << "" << std::setw(w) << ""
              << std::scientific << std::setprecision(4)
              << std::setw(w) << totalParticles << totalTrue << "\n\n";

    auto source_positions = read_positions_from_file("default_values/source_positions.txt");
    auto om_positions     = read_positions_from_file("default_values/om_positions.txt");

    std::cout << "sources loaded    : " << sources.size()          << "\n";
    std::cout << "results filled    : " << results.size()          << "\n";
    std::cout << "source_positions  : " << source_positions.size() << "\n";
    std::cout << "om_positions      : " << om_positions.size()     << "\n";

    {
        TFile* test = TFile::Open("vertex_energy_results.root", "READ");
        if (!test || test->IsZombie()) {
            std::cerr << "ERROR: vertex_energy_results.root not found!\n"
                      << "Run visu_root() first.\n";
            if (test) test->Close(); return;
        }
        TH2F* hp = (TH2F*)test->Get("SourceOMpos");
        TH2F* hn = (TH2F*)test->Get("SourceOMneg");
        if (!hp || !hn) {
            std::cerr << "ERROR: SourceOMpos/neg not found!\n";
            test->Close(); return;
        }
        std::cout << "SourceOMpos bins  : " << hp->GetNbinsX() << " x " << hp->GetNbinsY() << "\n";
        std::cout << "SourceOMneg bins  : " << hn->GetNbinsX() << " x " << hn->GetNbinsY() << "\n";
        std::cout << "SourceOMpos max   : " << hp->GetMaximum() << "\n";
        std::cout << "SourceOMneg max   : " << hn->GetMaximum() << "\n";
        test->Close();
    }

    // ── draw_efficiencies: малює PNG і ОДНОЧАСНО пише eff_*_srcNN у ROOT ──
    //    outputRootFile спочатку створюється тут (RECREATE),
    //    тому передаємо "RECREATE" — saveActivityResultsToRoot потім додасть
    //    решту об'єктів через UPDATE.
    draw_efficiencies(
        results,
        source_positions,
        om_positions,
        "vertex_energy_results.root",
        /*use5x5=*/false,
        DisplayMode::PERCENT,
        electron_intens,
        /*effRootFile=*/outputRootFile,
        /*effRootMode=*/"RECREATE");

    // ── Записуємо decay-дані та метадані (UPDATE — eff_* вже є) ────────
    saveActivityResultsToRoot(
        outputRootFile,
        results,
        sources,
        startYear, startMonth, startDay,
        endYear,   endMonth,   endDay,
        halfLife_yr,
        dead_time,
        totalParticles,
        totalTrue,
        electron_intens);
}

int main() {
    try {
        activity();
    } catch (const std::exception& e) {
        std::cerr << "Fatal: " << e.what() << "\n";
        return 1;
    }
    return 0;
}