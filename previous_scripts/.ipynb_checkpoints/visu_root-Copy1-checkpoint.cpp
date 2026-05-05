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
//  [CHANGE] Побудова енергетичного спектру:
//           три криві на одному canvas —
//           French side (X>0), Italian side (X<0), разом
//           У легенді — кількість подій для кожної лінії
// ============================================================

void draw_energy_spectra(TH1F* hE_pos, TH1F* hE_neg, TH1F* hE_all)
{
    gSystem->mkdir("plots", true);

    TCanvas* c = new TCanvas("cEnergy", "Energy spectra", 900, 650);
    c->SetLeftMargin(0.12);
    c->SetRightMargin(0.05);
    c->SetTopMargin(0.09);
    c->SetBottomMargin(0.11);
    // c->SetLogy();   // логарифмічна шкала по Y — зручніше для спектрів

    // Колір і стиль ліній
    hE_all->SetLineColor(kBlack);
    hE_all->SetLineWidth(2);
    hE_pos->SetLineColor(kBlue + 1);
    hE_pos->SetLineWidth(2);
    hE_neg->SetLineColor(kRed  + 1);
    hE_neg->SetLineWidth(2);

    hE_all->SetStats(0);
    hE_pos->SetStats(0);
    hE_neg->SetStats(0);

    // Малюємо: спочатку той, що має найбільший максимум
    hE_all->SetTitle("Energy spectrum;Energy [MeV];Counts");
    hE_all->Draw("HIST");
    hE_pos->Draw("HIST same");
    hE_neg->Draw("HIST same");

    // Легенда з кількостями подій
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
//  Область відображення: ±zoom мм навколо центру джерела
//
//  [CHANGE] При side==0 (X>0, French) координата Y інвертується: y -> -y
//           Це стосується і центру джерела (yc), і підпису.
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
        double yc_orig = src_pos[s].first;
        double zc      = src_pos[s].second;
        for(int side = 0; side < 2; side++)
        {
            TH2F* h = (side == 0) ? hStart_pos[s] : hStart_neg[s];
            double yc = (side == 0) ? -yc_orig : yc_orig;

            // [CHANGE] Рахуємо кількість подій у гістограмі
            double nEvents = h->GetEntries();

            // [CHANGE] Оновлюємо заголовок: додаємо кількість подій
            // Зберігаємо оригінальний заголовок і дописуємо до нього
            TString newTitle = Form("%s  (N = %.0f)",
                                    h->GetTitle(), nEvents);
            h->SetTitle(newTitle.Data());

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

            // [CHANGE] Відновлюємо оригінальний заголовок після збереження,
            // щоб повторний виклик функції не дублював "(N = ...)" у назві
            h->SetTitle(TString(h->GetTitle())
                        .ReplaceAll(Form("  (N = %.0f)", nEvents), "")
                        .Data());

            delete c;
        }
    }
}


// ============================================================
//  Карти джерело -> OM
//
//  [CHANGE] При X>0 (French, side==0):
//           - координати OM по Y інвертуються: y -> -y
//           - координати джерел по Y інвертуються: sy -> -sy
//           Нумерація стовпців і рядків не змінюється.
//
//  [CHANGE] У заголовок гістограми додається сумарна кількість подій
//           для даного джерела на даній стороні.
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

    // ============================================================
    //  Крок 1. Унікальні Y і Z позиції OM-ів
    //
    //  [CHANGE] Вектори unique_y_pos і unique_y_neg будуються окремо:
    //           для French (pos) береться -y, для Italian (neg) — y як є.
    //           unique_z однаковий для обох сторін.
    // ============================================================

    // Italian (neg): Y без змін
    std::vector<double> unique_y_neg, unique_z;
    for(auto& p : om_pos)
    {
        unique_y_neg.push_back(p.first);
        unique_z    .push_back(p.second);
    }

    // French (pos): Y інвертується
    // [CHANGE] y -> -y для French side
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

    // ============================================================
    //  Крок 2. Пошук найближчого індексу
    // ============================================================
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

    // ============================================================
    //  Крок 3. Крок сітки (медіана різниць)
    // ============================================================
    auto calc_step = [](const std::vector<double>& uv) -> double
    {
        if(uv.size() < 2) return 1.0;
        std::vector<double> diffs;
        for(size_t i = 1; i < uv.size(); i++)
            diffs.push_back(uv[i] - uv[i-1]);
        std::sort(diffs.begin(), diffs.end());
        return diffs[diffs.size()/2];
    };

    // [CHANGE] Окремий крок для pos (інвертований Y) і neg
    double stepY_pos = calc_step(unique_y_pos);
    double stepY_neg = calc_step(unique_y_neg);
    double stepZ     = calc_step(unique_z);

    // ============================================================
    //  Крок 4. Параметри сторін
    // ============================================================
    const int colOffset_pos = 20;   // French:  cols 20..20+nCols-1
    const int colOffset_neg =  0;   // Italian: cols  0..nCols-1

    struct SideParam {
        TH2F*                   matrix;
        int                     colOffset;
        int                     nCols;
        int                     nRows;
        const char*             sideName;
        const char*             sideTitle;
        const std::vector<double>* unique_y;   // [CHANGE] вказівник на потрібний вектор Y
        double                  stepY;          // [CHANGE] крок по Y для даної сторони
    };

    SideParam sides[2] = {
        { hSourceOM_pos, colOffset_pos, nColsTotal, nRowsTotal,
          "pos", "French side",  &unique_y_pos, stepY_pos },   // [CHANGE] pos -> інвертований Y
        { hSourceOM_neg, colOffset_neg, nColsTotal, nRowsTotal,
          "neg", "Italian side", &unique_y_neg, stepY_neg }    // [CHANGE] neg -> оригінальний Y
    };

    gStyle->SetPalette(kBird);
    gStyle->SetNumberContours(99);

    for(int side = 0; side < 2; side++)
    {
        SideParam& sp = sides[side];

        const int colMin = sp.colOffset;
        const int colMax = sp.colOffset + sp.nCols - 1;
        const int rowMin = 0;
        const int rowMax = sp.nRows - 1;

        // Розміри canvas для квадратних клітинок
        const double mL = 0.10, mR = 0.14, mT = 0.09, mB = 0.10;
        const int    baseH  = 600;
        const double areaH  = baseH * (1.0 - mT - mB);
        const double areaW  = areaH * (double)sp.nCols / (double)sp.nRows;
        const int    canvasW = (int)(areaW / (1.0 - mL - mR)) + 1;
        const int    canvasH = baseH;

        for(int s = 0; s < nSources; s++)
        {
            // [CHANGE] Рахуємо сумарну кількість подій для джерела s на цій стороні
            double totalHits = 0;
            for(int om = 0; om < nOMs; om++)
                totalHits += sp.matrix->GetBinContent(s + 1, om + 1);

            // Гістограма: вісь X = колонки, вісь Y = рядки
            // [CHANGE] У заголовку додано сумарну кількість подій
            TH2F* hGrid = new TH2F(
                Form("hGrid_%d_%d", side, s),
                Form("Source %d  %s  (total events: %.0f);Column;Row",
                     s, sp.sideTitle, totalHits),
                sp.nCols, colMin - 0.5, colMax + 0.5,
                sp.nRows, rowMin - 0.5, rowMax + 0.5);
            hGrid->SetDirectory(0);
            hGrid->SetStats(0);

            // Заповнюємо хітами
            // [CHANGE] для French (side==0) використовується -om_pos[om].first
            //          через sp.unique_y, який вже містить інвертовані значення
            for(int om = 0; om < nOMs; om++)
            {
                double hits = sp.matrix->GetBinContent(s + 1, om + 1);
                if(hits <= 0) continue;

                // [CHANGE] Для pos: шукаємо індекс у unique_y_pos (де Y вже -y)
                //          Для neg: шукаємо у unique_y_neg (Y без змін)
                int ci = coord_to_idx(*sp.unique_y, -om_pos[om].first * (side==0 ? 1.0 : -1.0))
                         + sp.colOffset;
                // Спрощено: для side==0 передаємо -y, для side==1 передаємо +y
                // Те саме записати явніше:
                double om_y_for_search = (side == 0) ? -om_pos[om].first
                                                      :  om_pos[om].first;
                ci = coord_to_idx(*sp.unique_y, om_y_for_search) + sp.colOffset;

                int ri = coord_to_idx(unique_z, om_pos[om].second);

                hGrid->Fill(ci, ri, hits);
            }

            // Canvas
            TCanvas* c = new TCanvas(
                Form("cGrid_%d_%d", side, s), "", canvasW, canvasH);
            c->SetLeftMargin(mL);
            c->SetRightMargin(mR);
            c->SetTopMargin(mT);
            c->SetBottomMargin(mB);
            c->SetCanvasSize(canvasW, canvasH);

            // Підписи осей
            hGrid->GetXaxis()->SetNdivisions(sp.nCols, false);
            hGrid->GetYaxis()->SetNdivisions(sp.nRows, false);
            hGrid->GetXaxis()->SetLabelSize(0.028);
            hGrid->GetYaxis()->SetLabelSize(0.028);
            hGrid->GetXaxis()->SetTickLength(0.01);
            hGrid->GetYaxis()->SetTickLength(0.01);

            for(int ic = 0; ic < sp.nCols; ic++)
                hGrid->GetXaxis()->SetBinLabel(ic + 1, Form("%d", colMin + ic));
            for(int ir = 0; ir < sp.nRows; ir++)
                hGrid->GetYaxis()->SetBinLabel(ir + 1, Form("%d", rowMin + ir));

            hGrid->Draw("COLZ");

            // Червона сітка
            for(int ic = 0; ic <= sp.nCols; ic++)
            {
                double x = colMin - 0.5 + ic;
                TLine* l = new TLine(x, rowMin - 0.5, x, rowMax + 0.5);
                l->SetLineColor(kRed);
                l->SetLineWidth(1);
                l->Draw("same");
            }
            for(int ir = 0; ir <= sp.nRows; ir++)
            {
                double y = rowMin - 0.5 + ir;
                TLine* l = new TLine(colMin - 0.5, y, colMax + 0.5, y);
                l->SetLineColor(kRed);
                l->SetLineWidth(1);
                l->Draw("same");
            }

            // ============================================================
            //  Джерела: точне розміщення
            //
            //  [CHANGE] Для French (side==0): sy -> -sy перед перетворенням
            //           Для Italian (side==1): sy без змін
            // ============================================================
            for(int i = 0; i < nSources; i++)
            {
                double sy_orig = src_pos[i].first;
                double sz      = src_pos[i].second;

                // [CHANGE] інвертування Y для French side
                double sy = (side == 0) ? -sy_orig : sy_orig;

                double col_src = sp.colOffset +
                    (sy - (*sp.unique_y)[0]) / sp.stepY;
                double row_src =
                    (sz - unique_z[0]) / stepZ;

                int    markerColor = (i == s) ? kBlue  : kGreen + 2;
                double markerSize  = (i == s) ? 1.8    : 1.0;

                TGraph* gSrc = new TGraph(1, &col_src, &row_src);
                gSrc->SetMarkerStyle(20);
                gSrc->SetMarkerSize(markerSize);
                gSrc->SetMarkerColor(markerColor);
                gSrc->Draw("P same");

                TLatex* lsrc = new TLatex(
                    col_src + 0.15, row_src + 0.15,
                    Form("%d", i));
                lsrc->SetTextSize(0.020);
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

    // [CHANGE] Додаємо окремі гістограми енергій для кожної сторони
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

    const double zoom = 50.0;
    const int    nBinsZoom = 100; // 1 бін = 1 мм
    
    std::vector<TH2F*> hSrcStart_pos(nSources);
    std::vector<TH2F*> hSrcStart_neg(nSources);
    
    for(int s = 0; s < nSources; s++)
    {
        // [CHANGE] Для French side (pos) центр джерела по Y інвертується: yc -> -yc
        double yc_neg = source_positions[s].first;          // Italian: без змін
        double yc_pos = -source_positions[s].first;         // [CHANGE] French: -y
        double zc     = source_positions[s].second;
    
        // [CHANGE] hSrcStart_pos будується навколо (-yc, zc)
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
    std::cout << "Reading " << nentries << " entries..." << std::endl;

    for(Long64_t ie = 0; ie < nentries; ie++)
    {
        evTree->GetEntry(ie);

        if(side == 1)   // French side, X > 0
        {
            // [CHANGE] Для French side заповнюємо з -y_start (інверсія Y)
            hStart_Xpos->Fill(-y_start, z_start);
            if(!pos_found){ x_pos_value = x_start; pos_found = true; }
        }
        else            // Italian side, X < 0
        {
            hStart_Xneg->Fill(y_start, z_start);
            if(!neg_found){ x_neg_value = x_start; neg_found = true; }
        }

        // Start vertex — гістограма для конкретного джерела
        if(source_id >= 0 && source_id < nSources)
        {
            // [CHANGE] French side: -y_start; Italian: y_start
            if(side == 1) hSrcStart_pos[source_id]->Fill(-y_start, z_start);
            else          hSrcStart_neg[source_id]->Fill( y_start, z_start);
        }

        // End vertex
        // [CHANGE] French (x_end > 0): заповнюємо з -y_end
        if(x_end > 0) hEnd_Xpos->Fill(-y_end, z_end);
        else          hEnd_Xneg->Fill( y_end, z_end);

        // Енергія
        if(energy > 0)
        {
            hEnergy->Fill(energy);
            // [CHANGE] Заповнюємо окремі гістограми енергій за стороною
            if(side == 1) hEnergy_pos->Fill(energy);
            else          hEnergy_neg->Fill(energy);
        }

        // Матриця джерело -> OM (координати тут не потрібні, лише ID)
        if(source_id >= 0 && om_id >= 0)
        {
            if(side == 1) hSourceOM_pos->Fill(source_id, om_id);
            else          hSourceOM_neg->Fill(source_id, om_id);
        }
    }

    f->Close();

    // --- Оновлюємо заголовки загальних гістограм ---
    // [CHANGE] Для French side підпис осі "-Y [mm]" щоб відобразити інверсію
    hStart_Xpos->SetTitle(Form("Start vertices (X = %.1f mm);-Y [mm];Z [mm]", x_pos_value));
    hStart_Xneg->SetTitle(Form("Start vertices (X = %.1f mm);Y [mm];Z [mm]",  x_neg_value));
    hEnd_Xpos->SetTitle( Form("End vertices (X #approx %.1f mm);-Y [mm];Z [mm]", x_pos_value));
    hEnd_Xneg->SetTitle( Form("End vertices (X #approx %.1f mm);Y [mm];Z [mm]",  x_neg_value));

    // --- Малюємо загальні PNG ---
    TCanvas c1; hStart_Xpos->Draw("COLZ"); c1.SaveAs("plots/start_Xpos.png");
    TCanvas c2; hStart_Xneg->Draw("COLZ"); c2.SaveAs("plots/start_Xneg.png");
    TCanvas c3; hEnd_Xpos->Draw("COLZ");   c3.SaveAs("plots/end_Xpos.png");
    TCanvas c4; hEnd_Xneg->Draw("COLZ");   c4.SaveAs("plots/end_Xneg.png");
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
    hEnergy_pos->Write();   // [CHANGE]
    hEnergy_neg->Write();   // [CHANGE]
    hSourceOM_pos->Write();
    hSourceOM_neg->Write();
    for(int s = 0; s < nSources; s++)
    {
        hSrcStart_pos[s]->Write();
        hSrcStart_neg[s]->Write();
    }
    fout->Close();

    // [CHANGE] Будуємо графік енергетичного спектру з трьома кривими
    draw_energy_spectra(hEnergy_pos, hEnergy_neg, hEnergy);

    // --- Детальні карти джерело -> OM ---
    draw_source_OM_distributions(
        hSourceOM_pos,
        hSourceOM_neg,
        source_positions,
        om_positions);

    draw_source_start_vertices(
    hSrcStart_pos,
    hSrcStart_neg,
    source_positions);


    std::cout << "Finished" << std::endl;
}
