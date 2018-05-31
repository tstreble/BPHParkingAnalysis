// Author: T. Strebler (IC)
// Date:   31 May 2018
//
// Wrapper for NanoAOD tree
// Can be either included in a interpreted macro or compiled in c++
// (use `root-config --glibs --cflags`)
//
// Create the nanoAODTree object from the pointer to the tree, then access the stored objects from it
// Common TTree functions GetEntry (entry), GetEntries() are implemented



#ifndef NANOAODTREE_H
#define NANOAODTREE_H



#include <TROOT.h>
#include <TChain.h>
#include <TTree.h>

const int kMuonMax = 100;
const int kBToKpipiMax = 10000;
const int kGenPartMax = 10000;

using namespace std;

class NanoAODTree {
public :
   TChain          *_tree; 

   //Old NanoAOD branches used
   int run;
   int luminosityBlock;
   int event;

   int nMuon;
   int Muon_charge[kMuonMax];
   float Muon_pt[kMuonMax];
   float Muon_eta[kMuonMax];
   float Muon_phi[kMuonMax];
   float Muon_mass[kMuonMax];
   float Muon_dxy[kMuonMax];
   float Muon_dz[kMuonMax];
   float Muon_pfRelIso04_all[kMuonMax];
   bool Muon_mediumId[kMuonMax];

   int nBToKpipi;
   float BToKpipi_CL_vtx[kBToKpipiMax];
   float BToKpipi_Kpi_CL_vtx[kBToKpipiMax];
   float BToKpipi_Kpi_mass[kBToKpipiMax];
   float BToKpipi_mass[kBToKpipiMax];
   int BToKpipi_piBu_charge[kBToKpipiMax];
   float BToKpipi_piBu_pt[kBToKpipiMax];
   float BToKpipi_piBu_eta[kBToKpipiMax];
   float BToKpipi_piBu_phi[kBToKpipiMax];
   float BToKpipi_kaon_pt[kBToKpipiMax];
   float BToKpipi_kaon_eta[kBToKpipiMax];
   float BToKpipi_kaon_phi[kBToKpipiMax];
   float BToKpipi_piD0_pt[kBToKpipiMax];
   float BToKpipi_piD0_eta[kBToKpipiMax];
   float BToKpipi_piD0_phi[kBToKpipiMax];
 
   int nGenPart;
   int GenPart_pdgId[kGenPartMax];
   int GenPart_genPartIdxMother[kGenPartMax];
   float GenPart_pt[kGenPartMax];
   float GenPart_eta[kGenPartMax];
   float GenPart_phi[kGenPartMax];

   // methods
   NanoAODTree (TChain* tree);
   ~NanoAODTree();
   void Init(TChain* tree);
   Int_t GetEntry(int entry);
   Long64_t GetEntries();
   TChain* GetTree();

};



NanoAODTree::NanoAODTree (TChain* tree)
{
    Init(tree);
}


NanoAODTree::~NanoAODTree() {}



void NanoAODTree::Init(TChain* tree)
{

  // Set branch addresses and branch pointers
  if (!tree) return;
  _tree = tree;  
  _tree->SetMakeClass(1); // needed especially when compiling
  
  _tree->SetBranchAddress("run",&run);
  _tree->SetBranchAddress("luminosityBlock",&luminosityBlock);
  _tree->SetBranchAddress("event",&event);

  _tree->SetBranchAddress("nMuon",&nMuon);  
  _tree->SetBranchAddress("Muon_pt",&Muon_pt);  
  _tree->SetBranchAddress("Muon_eta",&Muon_eta);
  _tree->SetBranchAddress("Muon_phi",&Muon_phi);
  _tree->SetBranchAddress("Muon_mass",&Muon_mass);
  _tree->SetBranchAddress("Muon_dxy",&Muon_dxy);
  _tree->SetBranchAddress("Muon_dz",&Muon_dz);
  _tree->SetBranchAddress("Muon_pfRelIso04_all",&Muon_pfRelIso04_all);
  _tree->SetBranchAddress("Muon_mediumId",&Muon_mediumId);

  _tree->SetBranchAddress("nBToKpipi",&nBToKpipi);
  _tree->SetBranchAddress("BToKpipi_CL_vtx",&BToKpipi_CL_vtx);
  _tree->SetBranchAddress("BToKpipi_Kpi_CL_vtx",&BToKpipi_Kpi_CL_vtx);
  _tree->SetBranchAddress("BToKpipi_Kpi_mass",&BToKpipi_Kpi_mass);
  _tree->SetBranchAddress("BToKpipi_mass",&BToKpipi_mass);
  _tree->SetBranchAddress("BToKpipi_piBu_charge",&BToKpipi_piBu_charge);
  _tree->SetBranchAddress("BToKpipi_piBu_pt",&BToKpipi_piBu_pt);
  _tree->SetBranchAddress("BToKpipi_piBu_eta",&BToKpipi_piBu_eta);
  _tree->SetBranchAddress("BToKpipi_piBu_phi",&BToKpipi_piBu_phi);
  _tree->SetBranchAddress("BToKpipi_kaon_pt",&BToKpipi_kaon_pt);
  _tree->SetBranchAddress("BToKpipi_kaon_eta",&BToKpipi_kaon_eta);
  _tree->SetBranchAddress("BToKpipi_kaon_phi",&BToKpipi_kaon_phi);
  _tree->SetBranchAddress("BToKpipi_piD0_pt",&BToKpipi_piD0_pt);
  _tree->SetBranchAddress("BToKpipi_piD0_eta",&BToKpipi_piD0_eta);
  _tree->SetBranchAddress("BToKpipi_piD0_phi",&BToKpipi_piD0_phi);

  _tree->SetBranchAddress("nGenPart",&nGenPart);
  _tree->SetBranchAddress("GenPart_pdgId",&GenPart_pdgId);
  _tree->SetBranchAddress("GenPart_genPartIdxMother",&GenPart_genPartIdxMother);
  _tree->SetBranchAddress("GenPart_pt",&GenPart_pt);
  _tree->SetBranchAddress("GenPart_eta",&GenPart_phi);

}


Int_t NanoAODTree::GetEntry(int entry)
{
    return _tree->GetEntry(entry);
} 

Long64_t NanoAODTree::GetEntries()
{
    return _tree->GetEntries();
}

TChain* NanoAODTree::GetTree()
{
    return _tree;
}

#endif
