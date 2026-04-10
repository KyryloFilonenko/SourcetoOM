#include "MiModule/include/MiEvent.h"
#include "MiModule/include/MiVertex.h"
#include "MiModule/include/MiVector3D.h"
#include "MiModule/include/MiPTD.h"
#include "MiModule/include/MiCD.h"
#include "MiModule/include/MiGID.h"

#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>

R__LOAD_LIBRARY(/sps/nemo/scratch/kfilonen/Falaise/MiModule/lib/libMiModule.so);


// ============================================================
//  Допоміжні функції
// ============================================================

int find_closest_source(double y, double z,
                        const std::vector<std::pair<double,double>>& sources)
{
    double min_dist = 1e12;
    int    best     = -1;

    for(size_t i = 0; i < sources.size(); i++)
    {
        double dy   = y - sources[i].first;
        double dz   = z - sources[i].second;
        double dist = sqrt(dy*dy + dz*dz);

        if(dist < min_dist){ min_dist = dist; best = (int)i; }
    }

    return best;
}


const double OM_SIZE_Y = 256.0;
const double OM_SIZE_Z = 256.0;

int find_OM(double y, double z,
            const std::vector<std::pair<double,double>>& oms)
{
    for(size_t i = 0; i < oms.size(); i++)
    {
        if( fabs(y - oms[i].first)  < OM_SIZE_Y / 2 &&
            fabs(z - oms[i].second) < OM_SIZE_Z / 2 )
            return (int)i;
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
//  Головна функція
// ============================================================

void making_root()
{
    gROOT->SetBatch(true);

    // --- Вхідний файл ---
    TFile* f = TFile::Open("DATA/merged_week_2025_X.root", "READ");
    if(!f || f->IsZombie())
    {
        std::cerr << "ERROR: cannot open input file!" << std::endl;
        return;
    }

    TTree*   t   = (TTree*)f->Get("Event");
    MiEvent* Eve = new MiEvent();
    t->SetBranchAddress("Eventdata", &Eve);

    // --- Позиції джерел і OM ---
    std::vector<std::pair<double,double>> source_positions =
        read_positions_from_file("default_values/source_positions.txt");

    std::vector<std::pair<double,double>> om_positions =
        read_positions_from_file("default_values/om_positions.txt");

    // --- Вихідний файл ---
    TFile* fout = new TFile("processed_data.root", "RECREATE");
    if(!fout || fout->IsZombie())
    {
        std::cerr << "ERROR: cannot create output file!" << std::endl;
        return;
    }

    // --- Змінні гілок ---
    int    event_id;
    int    source_id;
    int    om_id;
    int    side;        // +1 = X>0,  -1 = X<0
    double x_start, y_start, z_start;
    double x_end,   y_end,   z_end;
    double energy;

    string gid_type;
	string gid_module;
    // GID: зберігаємо як int (стоки з MiGID конвертуємо через stoi)
    int gid_side;
    int gid_wall;
    int gid_column;
    int gid_row;

    // --- Дерево ---
    TTree* evTree = new TTree("EventTree", "Single-electron events");

    evTree->Branch("event_id",  &event_id,  "event_id/I");
    evTree->Branch("source_id", &source_id, "source_id/I");
    evTree->Branch("om_id",     &om_id,     "om_id/I");
    evTree->Branch("side",      &side,      "side/I");
    evTree->Branch("x_start",   &x_start,   "x_start/D");
    evTree->Branch("y_start",   &y_start,   "y_start/D");
    evTree->Branch("z_start",   &z_start,   "z_start/D");
    evTree->Branch("x_end",     &x_end,     "x_end/D");
    evTree->Branch("y_end",     &y_end,     "y_end/D");
    evTree->Branch("z_end",     &z_end,     "z_end/D");
    evTree->Branch("energy",    &energy,    "energy/D");
    
    evTree->Branch("gid_type",  &gid_type,  "gid_type/S");
    evTree->Branch("gid_module",&gid_module,"gid_module/S");
    evTree->Branch("gid_side",  &gid_side,  "gid_side/I");
    evTree->Branch("gid_wall",  &gid_wall,  "gid_wall/I");
    evTree->Branch("gid_column",&gid_column,"gid_column/I");
    evTree->Branch("gid_row",   &gid_row,   "gid_row/I");

    // Лямбда для безпечного перетворення string -> int
    auto safe_stoi = [](const std::string& s) -> int {
        if(s.empty()) return -1;
        try { return std::stoi(s); }
        catch(...) { return -1; }
    };

    
    // --- Головний цикл ---
    Long64_t nentries  = t->GetEntries();
    Long64_t maxEvents = nentries;
    // Long64_t maxEvents = std::min(nentries, (Long64_t)10000);  // для тесту

    std::cout << "Processing " << maxEvents << " entries..." << std::endl;

    for(Long64_t ie = 0; ie < maxEvents; ie++)
    {
        if(ie % 10000 == 0)
            std::cout << "  Entry " << ie << " / " << maxEvents << std::endl;

        t->GetEntry(ie);

        MiPTD* ptd    = Eve->getPTD();
        int    nParts = Eve->getPTDNoPart();

        for(int ip = 0; ip < nParts; ip++)
        {
            MiCDParticle* particle = ptd->getpart(ip);

            // Відбираємо лише електрони
            if(particle->getcharge() != 1000) continue;

            // Потрібен хоча б один удар у калориметрі
            std::vector<MiCDCaloHit>* caloHits = particle->getcalohitv();
            if(!caloHits || caloHits->empty()) continue;

            // Потрібні мінімум 2 вершини (початок і кінець)
            int nVerts = Eve->getPTDNoVert(ip);
            if(nVerts < 2) continue;

            // --- Координати ---
            x_start = Eve->getPTDverX(ip, 0);
            y_start = Eve->getPTDverY(ip, 0);
            z_start = Eve->getPTDverZ(ip, 0);

            x_end = Eve->getPTDverX(ip, nVerts - 1);
            y_end = Eve->getPTDverY(ip, nVerts - 1);
            z_end = Eve->getPTDverZ(ip, nVerts - 1);

            // --- Енергія і GID ---
            MiCDCaloHit* hit = particle->getcalohit(0);

            energy     = -1.0;
            gid_side   = -1;
            gid_wall   = -1;
            gid_column = -1;
            gid_row    = -1;

            if(hit)
            {
                // Енергія
                double e = hit->getE();
                energy = (!std::isnan(e)) ? e : -1.0;

                // GID — getGID() повертає MiGID об'єкт (не pointer)
                MiGID* gid = hit->getGID();

                gid_type   = gid->gettype();
                gid_module   = gid->getmodule();

                gid_side   = safe_stoi(gid->getside());
                gid_wall   = safe_stoi(gid->getwall());
                gid_column = safe_stoi(gid->getcolumn());
                gid_row    = safe_stoi(gid->getrow());
            }

            // --- Метадані ---
            event_id  = (int)ie;
            side      = (x_start > 0) ? 1 : -1;
            source_id = find_closest_source(y_start, z_start, source_positions);
            om_id     = find_OM(y_end, z_end, om_positions);

            evTree->Fill();
        }
    }

    // --- Зберігаємо ---
    Long64_t nStored = evTree->GetEntries();

    fout->cd();
    evTree->Write();
    fout->Close();
    f->Close();

    std::cout << "Done. Events stored: " << nStored << std::endl;
    std::cout << "Output: processed_data.root" << std::endl;
}
