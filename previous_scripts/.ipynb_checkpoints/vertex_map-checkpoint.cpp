// #include "MiModule/include/MiEvent.h"
// #include "MiModule/include/MiVertex.h"
// #include "MiModule/include/MiVector3D.h"
// #include "MiModule/include/MiPTD.h"
// #include "MiModule/include/MiCD.h"

// #include <TFile.h>
// #include <TTree.h>
// #include <TH1F.h>
// #include <TCanvas.h>
// #include <TGraph2D.h>
// #include <TLegend.h>

// #include <vector>
// #include <iostream>

// R__LOAD_LIBRARY(/sps/nemo/scratch/kfilonen/Falaise/MiModule/lib/libMiModule.so);

// void vertex_map()
// {

// TFile* f = TFile::Open("merged_week_2025_X-BETE.root","READ");
// TTree* t = (TTree*)f->Get("Event");

// MiEvent* Eve = new MiEvent();
// t->SetBranchAddress("Eventdata",&Eve);


// TH1F* hEnergy = new TH1F(
//     "hEnergy",
//     "Energy distribution;Energy [MeV];Counts",
//     200,
//     0,
//     5
// );

// std::vector<double> startX;
// std::vector<double> startY;
// std::vector<double> startZ;

// std::vector<double> endX;
// std::vector<double> endY;
// std::vector<double> endZ;


// Long64_t nentries = t->GetEntries();

// for(Long64_t ie=0; ie<int(nentries / 10000000); ie++)
// {

//     t->GetEntry(ie);
//     Eve->print();

//     MiPTD* ptd = Eve->getPTD();
//     int nParts = Eve->getPTDNoPart();

//     for(int ip=0; ip<nParts; ip++)
//     {

//         MiCDParticle* particle = ptd->getpart(ip);

//         int charge = particle->getcharge();
//         if(charge != 1000) continue;


//         /* =============================================
//          if calohit
//            ============================================= */

//         std::vector<MiCDCaloHit>* caloHits = particle->getcalohitv();
//         if(!caloHits || caloHits->empty()) continue;


//         /* =============================================
//            n of vert
//            ============================================= */

//         int nVerts = Eve->getPTDNoVert(ip);
//         if(nVerts < 2) continue;


//         /* =============================================
//           first ver
//            ============================================= */

//         double x_start = Eve->getPTDverX(ip,0);
//         double y_start = Eve->getPTDverY(ip,0);
//         double z_start = Eve->getPTDverZ(ip,0);

//         startX.push_back(x_start);
//         startY.push_back(y_start);
//         startZ.push_back(z_start);


//         /* =============================================
//            end vert
//            ============================================= */

//         double x_end = Eve->getPTDverX(ip,nVerts-1);
//         double y_end = Eve->getPTDverY(ip,nVerts-1);
//         double z_end = Eve->getPTDverZ(ip,nVerts-1);

//         std::cout << "(" << x_start << "; " << y_start << "; " << z_start << ") -> (" << x_end << "; " << y_end << "; " << z_end << ")" << std::endl;

//         std::cout << nVerts << std::endl;
//         for(int i=0; i<nVerts;i++)
//         {
//             std::cout << Eve->getPTDverX(ip,i) << std::endl;
//         }

//         endX.push_back(x_end);
//         endY.push_back(y_end);
//         endZ.push_back(z_end);


//         /* =============================================
//            energy
//            ============================================= */

//         MiCDCaloHit* hit = particle->getcalohit(0);

//         if(hit)
//         {
//             double energy = hit->getE();
//             hEnergy->Fill(energy);
//             std::cout << energy << std::endl;
//         }

//     }

// }


// TGraph2D* gStart = new TGraph2D(startX.size());
// TGraph2D* gEnd   = new TGraph2D(endX.size());

// for(size_t i=0;i<startX.size();i++)
// {
//     gStart->SetPoint(i,startX[i],startY[i],startZ[i]);
// }

// for(size_t i=0;i<endX.size();i++)
// {
//     gEnd->SetPoint(i,endX[i],endY[i],endZ[i]);
// }


// gStart->SetMarkerColor(kBlue);
// gStart->SetMarkerStyle(20);

// gEnd->SetMarkerColor(kRed);
// gEnd->SetMarkerStyle(20);


// TCanvas* c1 = new TCanvas("c1","3D Vertex Map",1000,800);

// gStart->SetTitle("Track vertices in detector;X [mm];Y [mm];Z [mm]");

// gStart->Draw("P0");
// gEnd->Draw("P0 SAME");


// TLegend* leg = new TLegend(0.75,0.75,0.9,0.9);
// leg->AddEntry(gStart,"Start vertex","p");
// leg->AddEntry(gEnd,"End vertex","p");
// leg->Draw();



// TCanvas* c2 = new TCanvas("c2","Energy spectrum",800,600);

// hEnergy->Draw();


// TFile* fout = new TFile("vertex_energy_results.root","RECREATE");

// gStart->Write("start_vertices");
// gEnd->Write("end_vertices");
// hEnergy->Write();

// fout->Close();


// std::cout << "Finished processing." << std::endl;

// }



#include "MiModule/include/MiEvent.h"
#include "MiModule/include/MiVertex.h"
#include "MiModule/include/MiVector3D.h"
#include "MiModule/include/MiPTD.h"
#include "MiModule/include/MiCD.h"

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TROOT.h>

#include <vector>
#include <iostream>

R__LOAD_LIBRARY(/sps/nemo/scratch/kfilonen/Falaise/MiModule/lib/libMiModule.so);

void vertex_map()
{

/* =============================================
   ROOT batch mode (графіки не відкриваються)
   ============================================= */

gROOT->SetBatch(true);


/* =============================================
   ВІДКРИТТЯ ROOT ФАЙЛУ
   ============================================= */

TFile* f = TFile::Open("merged_week_2025_X-BETE.root","READ");
TTree* t = (TTree*)f->Get("Event");

MiEvent* Eve = new MiEvent();
t->SetBranchAddress("Eventdata",&Eve);


/* =============================================
   ГІСТОГРАМА ЕНЕРГІЙ
   ============================================= */

TH1F* hEnergy = new TH1F(
"Energy",
"Energy distribution;Energy [MeV];Counts",
200,0,5);


/* =============================================
   ГІСТОГРАМИ (Y,Z)
   ============================================= */

TH2F* hStart_Xpos = new TH2F("hStart_Xpos","",500,-2500,2500,500,-2000,2000);
TH2F* hStart_Xneg = new TH2F("hStart_Xneg","",500,-2500,2500,500,-2000,2000);

TH2F* hEnd_Xpos = new TH2F("hEnd_Xpos","",500,-2500,2500,500,-2000,2000);
TH2F* hEnd_Xneg = new TH2F("hEnd_Xneg","",500,-2500,2500,500,-2000,2000);


/* =============================================
   для запису X
   ============================================= */

double x_pos_value = 0;
double x_neg_value = 0;
bool pos_found=false;
bool neg_found=false;


/* =============================================
   ЦИКЛ ПО ПОДІЯХ
   ============================================= */

Long64_t nentries = t->GetEntries();

Long64_t maxEvents = std::min(nentries,(Long64_t)100000);

for(Long64_t ie=0; ie<maxEvents; ie++)
{

    t->GetEntry(ie);

    MiPTD* ptd = Eve->getPTD();
    int nParts = Eve->getPTDNoPart();

    for(int ip=0; ip<nParts; ip++)
    {

        MiCDParticle* particle = ptd->getpart(ip);

        int charge = particle->getcharge();
        if(charge != 1000) continue;

        std::vector<MiCDCaloHit>* caloHits = particle->getcalohitv();
        if(!caloHits || caloHits->empty()) continue;

        int nVerts = Eve->getPTDNoVert(ip);
        if(nVerts < 2) continue;


        /* ========= START VERTEX ========= */

        double x_start = Eve->getPTDverX(ip,0);
        double y_start = Eve->getPTDverY(ip,0);
        double z_start = Eve->getPTDverZ(ip,0);

        if(x_start > 0)
        {
            hStart_Xpos->Fill(y_start,z_start);

            if(!pos_found){
                x_pos_value=x_start;
                pos_found=true;
            }
        }
        else
        {
            hStart_Xneg->Fill(y_start,z_start);

            if(!neg_found){
                x_neg_value=x_start;
                neg_found=true;
            }
        }


        /* ========= END VERTEX ========= */

        double x_end = Eve->getPTDverX(ip,nVerts-1);
        double y_end = Eve->getPTDverY(ip,nVerts-1);
        double z_end = Eve->getPTDverZ(ip,nVerts-1);

        if(x_end > 0)
            hEnd_Xpos->Fill(y_end,z_end);
        else
            hEnd_Xneg->Fill(y_end,z_end);


        /* ========= ENERGY ========= */

        MiCDCaloHit* hit = particle->getcalohit(0);

        if(hit)
        {
            double energy = hit->getE();
            hEnergy->Fill(energy);
        }

    }

}


/* =============================================
   ЗАГОЛОВКИ З X
   ============================================= */

hStart_Xpos->SetTitle(Form("Start vertices (X = %.1f mm);Y [mm];Z [mm]",x_pos_value));
hStart_Xneg->SetTitle(Form("Start vertices (X = %.1f mm);Y [mm];Z [mm]",x_neg_value));

hEnd_Xpos->SetTitle(Form("End vertices (X ≈ %.1f mm);Y [mm];Z [mm]",x_pos_value));
hEnd_Xneg->SetTitle(Form("End vertices (X ≈ %.1f mm);Y [mm];Z [mm]",x_neg_value));


/* =============================================
   ЗБЕРЕЖЕННЯ PNG
   ============================================= */

TCanvas c1;
hStart_Xpos->Draw("COLZ");
c1.SaveAs("start_Xpos.png");

TCanvas c2;
hStart_Xneg->Draw("COLZ");
c2.SaveAs("start_Xneg.png");

TCanvas c3;
hEnd_Xpos->Draw("COLZ");
c3.SaveAs("end_Xpos.png");

TCanvas c4;
hEnd_Xneg->Draw("COLZ");
c4.SaveAs("end_Xneg.png");

TCanvas c5;
hEnergy->Draw();
c5.SaveAs("energy_spectrum.png");


/* =============================================
   ROOT OUTPUT
   ============================================= */

TFile* fout = new TFile("vertex_energy_results.root","RECREATE");

hStart_Xpos->Write();
hStart_Xneg->Write();
hEnd_Xpos->Write();
hEnd_Xneg->Write();
hEnergy->Write();

fout->Close();

std::cout<<"Finished"<<std::endl;

}