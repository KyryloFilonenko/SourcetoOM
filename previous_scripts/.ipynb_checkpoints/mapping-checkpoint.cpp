#include "MiModule/include/MiEvent.h"
#include "MiModule/include/MiVertex.h"
#include "MiModule/include/MiVector3D.h"
#include "MiModule/include/MiPTD.h"
#include "MiModule/include/MiCD.h"

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TROOT.h>

#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>

#include <TSystem.h>

R__LOAD_LIBRARY(/sps/nemo/scratch/kfilonen/Falaise/MiModule/lib/libMiModule.so);


struct EventInfo
{
    int event_id;
    int source_id;
    int om_id;
    int side;

    double x_start,y_start,z_start;
    double x_end,y_end,z_end;
    double energy;
};


int find_closest_source(double y,double z,
                        const std::vector<std::pair<double,double>>& sources)
{
    double min_dist = 1e12;
    int best = -1;

    for(size_t i=0;i<sources.size();i++)
    {
        double dy = y - sources[i].first;
        double dz = z - sources[i].second;

        double dist = sqrt(dy*dy + dz*dz);

        if(dist < min_dist)
        {
            min_dist = dist;
            best = i;
        }
    }

    return best;
}


const double OM_SIZE_Y = 256.0;
const double OM_SIZE_Z = 256.0;

int find_OM(double y,double z,
            const std::vector<std::pair<double,double>>& oms)
{
    for(size_t i=0;i<oms.size();i++)
    {
        double yc = oms[i].first;
        double zc = oms[i].second;

        if( fabs(y-yc) < OM_SIZE_Y/2 &&
            fabs(z-zc) < OM_SIZE_Z/2 )
        {
            return i;
        }
    }

    return -1;
}


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


void mapping()
{

gROOT->SetBatch(true);

TFile* f = TFile::Open("DATA/merged_week_2025_X-BETE.root","READ");
TTree* t = (TTree*)f->Get("Event");

MiEvent* Eve = new MiEvent();
t->SetBranchAddress("Eventdata",&Eve);


TH1F* hEnergy = new TH1F(
"Energy",
"Energy distribution;Energy [MeV];Counts",
200,0,3000);

TH2F* hStart_Xpos = new TH2F("hStart_Xpos","",500,-2500,2500,500,-2000,2000);
TH2F* hStart_Xneg = new TH2F("hStart_Xneg","",500,-2500,2500,500,-2000,2000);

TH2F* hEnd_Xpos = new TH2F("hEnd_Xpos","",500,-2500,2500,500,-2000,2000);
TH2F* hEnd_Xneg = new TH2F("hEnd_Xneg","",500,-2500,2500,500,-2000,2000);


const int nSources = 42;
const int nOM = 260;

TH2F* hSourceOM_pos = new TH2F(
"SourceOMpos",
"Source -> OM matrix;Source ID;OM ID",
nSources,0,nSources,
nOM,0,nOM);

TH2F* hSourceOM_neg = new TH2F(
"SourceOMneg",
"Source -> OM matrix;Source ID;OM ID",
nSources,0,nSources,
nOM,0,nOM);


std::vector<std::pair<double,double>> source_positions =
    read_positions_from_file("default_values/source_positions.txt");

std::vector<std::pair<double,double>> om_positions =
    read_positions_from_file("default_values/om_positions.txt");


std::vector<EventInfo> events;

double x_pos_value = 0;
double x_neg_value = 0;
bool pos_found=false;
bool neg_found=false;


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

        double x_end = Eve->getPTDverX(ip,nVerts-1);
        double y_end = Eve->getPTDverY(ip,nVerts-1);
        double z_end = Eve->getPTDverZ(ip,nVerts-1);

        if(x_end > 0)
            hEnd_Xpos->Fill(y_end,z_end);
        else
            hEnd_Xneg->Fill(y_end,z_end);

        double energy = 0;

        MiCDCaloHit* hit = particle->getcalohit(0);
        if(hit)
        {
            energy = hit->getE();
            if(!std::isnan(energy))
                hEnergy->Fill(energy);
        }


        EventInfo ev;

        ev.event_id = ie;
        ev.side = (x_start>0)?1:-1;

        if(ev.side == 1)
        {
            int source_id_pos = find_closest_source(
                y_start,
                z_start,
                source_positions);

            int om_id_pos = find_OM(
                y_end,
                z_end,
                om_positions);

            if(source_id_pos>=0 && om_id_pos>=0)
                hSourceOM_pos->Fill(source_id_pos,om_id_pos);

            ev.source_id = source_id_pos;
            ev.om_id = om_id_pos;
        }

        if(ev.side == -1)
        {
            int source_id_neg = find_closest_source(
                y_start,
                z_start,
                source_positions);

            int om_id_neg = find_OM(
                y_end,
                z_end,
                om_positions);

            if(source_id_neg>=0 && om_id_neg>=0)
                hSourceOM_neg->Fill(source_id_neg,om_id_neg);

            ev.source_id = source_id_neg;
            ev.om_id = om_id_neg;
        }

        ev.x_start=x_start;
        ev.y_start=y_start;
        ev.z_start=z_start;

        ev.x_end=x_end;
        ev.y_end=y_end;
        ev.z_end=z_end;

        ev.energy=energy;

        events.push_back(ev);

    }

}

hStart_Xpos->SetTitle(Form("Start vertices (X = %.1f mm);Y [mm];Z [mm]",x_pos_value));
hStart_Xneg->SetTitle(Form("Start vertices (X = %.1f mm);Y [mm];Z [mm]",x_neg_value));

hEnd_Xpos->SetTitle(Form("End vertices (X ≈ %.1f mm);Y [mm];Z [mm]",x_pos_value));
hEnd_Xneg->SetTitle(Form("End vertices (X ≈ %.1f mm);Y [mm];Z [mm]",x_neg_value));


TFile* fout = new TFile("vertex_energy_results.root","RECREATE");

hStart_Xpos->Write();
hStart_Xneg->Write();
hEnd_Xpos->Write();
hEnd_Xneg->Write();
hEnergy->Write();
hSourceOM_pos->Write();
hSourceOM_neg->Write();

fout->Close();

std::cout<<"Events stored: "<<events.size()<<std::endl;
std::cout<<"Finished"<<std::endl;

}
