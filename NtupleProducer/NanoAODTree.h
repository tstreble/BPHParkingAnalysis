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
const int kElectronMax = 100;
const int kBToKpipiMax = 1000000;
const int kBToKmumuMax = 50000;
const int kBToKeeMax = 50000;
const int kGenPartMax = 10000;
const int kTrigObjMax = 1000;
const int kPFCandMax = 10000;

using namespace std;

class NanoAODTree {
public :
   TChain          *_tree; 

   //Old NanoAOD branches used
   int run;
   int luminosityBlock;
   long event;

   uint nMuon;
   int Muon_charge[kMuonMax];
   float Muon_pt[kMuonMax];
   float Muon_eta[kMuonMax];
   float Muon_phi[kMuonMax];
   float Muon_mass[kMuonMax];
   float Muon_dxy[kMuonMax];
   float Muon_dz[kMuonMax];
   float Muon_pfRelIso04_all[kMuonMax];
   bool Muon_softId[kMuonMax];
   bool Muon_mediumId[kMuonMax];

   uint nElectron;
   int Electron_charge[kElectronMax];
   float Electron_pt[kElectronMax];
   float Electron_eta[kElectronMax];
   float Electron_phi[kElectronMax];
   float Electron_mass[kElectronMax];
   float Electron_dxy[kElectronMax];
   float Electron_dz[kElectronMax];

   uint nBToKpipi;
   float BToKpipi_CL_vtx[kBToKpipiMax];
   float BToKpipi_cosAlpha[kBToKpipiMax];
   float BToKpipi_Lxy[kBToKpipiMax];
   float BToKpipi_Kpi_CL_vtx[kBToKpipiMax];
   float BToKpipi_Kpi_mass[kBToKpipiMax];
   float BToKpipi_Kpi_pt[kBToKpipiMax];
   float BToKpipi_mass[kBToKpipiMax];
   int BToKpipi_piBu_charge[kBToKpipiMax];
   float BToKpipi_piBu_pt[kBToKpipiMax];
   float BToKpipi_piBu_eta[kBToKpipiMax];
   float BToKpipi_piBu_phi[kBToKpipiMax];
   float BToKpipi_kaon_pt[kBToKpipiMax];
   float BToKpipi_kaon_eta[kBToKpipiMax];
   float BToKpipi_kaon_phi[kBToKpipiMax];
   int BToKpipi_kaon_charge[kBToKpipiMax];
   float BToKpipi_kaon_dz[kBToKpipiMax];
   float BToKpipi_piD0_pt[kBToKpipiMax];
   float BToKpipi_piD0_eta[kBToKpipiMax];
   float BToKpipi_piD0_phi[kBToKpipiMax];
   float BToKpipi_piD0_dz[kBToKpipiMax];
 
   uint nBToKmumu;
   float BToKmumu_CL_vtx[kBToKmumuMax];
   float BToKmumu_mumu_CL_vtx[kBToKmumuMax];
   float BToKmumu_mumu_mass[kBToKmumuMax];
   float BToKmumu_mass[kBToKmumuMax];
   int BToKmumu_kaon_charge[kBToKmumuMax];
   float BToKmumu_kaon_pt[kBToKmumuMax];
   float BToKmumu_kaon_eta[kBToKmumuMax];
   float BToKmumu_kaon_phi[kBToKmumuMax];
   float BToKmumu_mu1_pt[kBToKmumuMax];
   int BToKmumu_mu1_charge[kBToKmumuMax];
   float BToKmumu_mu1_eta[kBToKmumuMax];
   float BToKmumu_mu1_phi[kBToKmumuMax];
   int BToKmumu_mu1_index[kBToKmumuMax];
   int BToKmumu_mu2_charge[kBToKmumuMax];
   float BToKmumu_mu2_pt[kBToKmumuMax];
   float BToKmumu_mu2_eta[kBToKmumuMax];
   float BToKmumu_mu2_phi[kBToKmumuMax];
   int BToKmumu_mu2_index[kBToKmumuMax];


   uint nBToKee;
   float BToKee_CL_vtx[kBToKeeMax];
   float BToKee_ee_CL_vtx[kBToKeeMax];
   float BToKee_ee_mass[kBToKeeMax];
   float BToKee_mass[kBToKeeMax];
   int BToKee_kaon_charge[kBToKeeMax];
   float BToKee_kaon_pt[kBToKeeMax];
   float BToKee_kaon_eta[kBToKeeMax];
   float BToKee_kaon_phi[kBToKeeMax];
   int BToKee_ele1_charge[kBToKeeMax];
   float BToKee_ele1_pt[kBToKeeMax];
   float BToKee_ele1_eta[kBToKeeMax];
   float BToKee_ele1_phi[kBToKeeMax];
   int BToKee_ele2_charge[kBToKeeMax];
   float BToKee_ele2_pt[kBToKeeMax];
   float BToKee_ele2_eta[kBToKeeMax];
   float BToKee_ele2_phi[kBToKeeMax];

   uint nGenPart;
   int GenPart_pdgId[kGenPartMax];
   int GenPart_genPartIdxMother[kGenPartMax];
   float GenPart_pt[kGenPartMax];
   float GenPart_eta[kGenPartMax];
   float GenPart_phi[kGenPartMax];

   bool HLT_Mu8p5_IP3p5_part0;
   bool HLT_Mu8p5_IP3p5_part1;
   bool HLT_Mu8p5_IP3p5_part2;
   bool HLT_Mu8p5_IP3p5_part3;
   bool HLT_Mu8p5_IP3p5_part4;
   bool HLT_Mu8p5_IP3p5_part5;
   bool HLT_Mu10p5_IP3p5_part0;
   bool HLT_Mu10p5_IP3p5_part1;
   bool HLT_Mu10p5_IP3p5_part2;
   bool HLT_Mu10p5_IP3p5_part3;
   bool HLT_Mu10p5_IP3p5_part4;
   bool HLT_Mu10p5_IP3p5_part5;
   bool HLT_Mu9_IP6_part0;
   bool HLT_Mu9_IP6_part1;
   bool HLT_Mu9_IP6_part2;
   bool HLT_Mu9_IP6_part3;
   bool HLT_Mu9_IP6_part4;
   bool HLT_Mu9_IP6_part5;
   bool HLT_Mu8_IP3_part0;
   bool HLT_Mu8_IP3_part1;
   bool HLT_Mu8_IP3_part2;
   bool HLT_Mu8_IP3_part3;
   bool HLT_Mu8_IP3_part4;
   bool HLT_Mu8_IP3_part5;

   uint nTrigObj;
   int TrigObj_id[kTrigObjMax];
   float TrigObj_pt[kTrigObjMax];
   float TrigObj_eta[kTrigObjMax];
   float TrigObj_phi[kTrigObjMax];
   int TrigObj_filterBits[kTrigObjMax];

   uint nPFCand;
   float PFCand_pt[kPFCandMax];
   float PFCand_eta[kPFCandMax];
   float PFCand_phi[kPFCandMax];
   float PFCand_mass[kPFCandMax];
   int PFCand_pdgId[kPFCandMax];
   float PFCand_DCASig[kPFCandMax];
   float PFCand_dz[kPFCandMax];

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
  _tree->GetEntries(); //ROOT issue with TChain which only works after calling GetEntries... (sigh)
  
  _tree->SetBranchAddress("run",&run);
  _tree->SetBranchAddress("luminosityBlock",&luminosityBlock);
  _tree->SetBranchAddress("event",&event);

  _tree->SetBranchAddress("nMuon",&nMuon);  
  _tree->SetBranchAddress("Muon_charge",&Muon_charge);  
  _tree->SetBranchAddress("Muon_pt",&Muon_pt);  
  _tree->SetBranchAddress("Muon_eta",&Muon_eta);
  _tree->SetBranchAddress("Muon_phi",&Muon_phi);
  _tree->SetBranchAddress("Muon_mass",&Muon_mass);
  _tree->SetBranchAddress("Muon_dxy",&Muon_dxy);
  _tree->SetBranchAddress("Muon_dz",&Muon_dz);
  _tree->SetBranchAddress("Muon_pfRelIso04_all",&Muon_pfRelIso04_all);
  _tree->SetBranchAddress("Muon_softId",&Muon_softId);
  _tree->SetBranchAddress("Muon_mediumId",&Muon_mediumId);

  _tree->SetBranchAddress("nElectron",&nElectron);  
  _tree->SetBranchAddress("Electron_charge",&Electron_charge);  
  _tree->SetBranchAddress("Electron_pt",&Electron_pt);  
  _tree->SetBranchAddress("Electron_eta",&Electron_eta);
  _tree->SetBranchAddress("Electron_phi",&Electron_phi);
  _tree->SetBranchAddress("Electron_mass",&Electron_mass);
  _tree->SetBranchAddress("Electron_dxy",&Electron_dxy);
  _tree->SetBranchAddress("Electron_dz",&Electron_dz);

  int BToKpipi_info = _tree->SetBranchAddress("nBToKpipi",&nBToKpipi);
  if(BToKpipi_info>=0){
    _tree->SetBranchAddress("BToKpipi_CL_vtx",&BToKpipi_CL_vtx);
    _tree->SetBranchAddress("BToKpipi_Lxy",&BToKpipi_Lxy);
    _tree->SetBranchAddress("BToKpipi_cosAlpha",&BToKpipi_cosAlpha);
    _tree->SetBranchAddress("BToKpipi_Kpi_CL_vtx",&BToKpipi_Kpi_CL_vtx);
    _tree->SetBranchAddress("BToKpipi_Kpi_mass",&BToKpipi_Kpi_mass);
    _tree->SetBranchAddress("BToKpipi_Kpi_pt",&BToKpipi_Kpi_pt);
    _tree->SetBranchAddress("BToKpipi_mass",&BToKpipi_mass);
    _tree->SetBranchAddress("BToKpipi_piBu_charge",&BToKpipi_piBu_charge);
    _tree->SetBranchAddress("BToKpipi_piBu_pt",&BToKpipi_piBu_pt);
    _tree->SetBranchAddress("BToKpipi_piBu_eta",&BToKpipi_piBu_eta);
    _tree->SetBranchAddress("BToKpipi_piBu_phi",&BToKpipi_piBu_phi);
    _tree->SetBranchAddress("BToKpipi_kaon_charge",&BToKpipi_kaon_charge);
    _tree->SetBranchAddress("BToKpipi_kaon_pt",&BToKpipi_kaon_pt);
    _tree->SetBranchAddress("BToKpipi_kaon_eta",&BToKpipi_kaon_eta);
    _tree->SetBranchAddress("BToKpipi_kaon_phi",&BToKpipi_kaon_phi);
    _tree->SetBranchAddress("BToKpipi_kaon_dz",&BToKpipi_kaon_dz);
    _tree->SetBranchAddress("BToKpipi_piD0_pt",&BToKpipi_piD0_pt);
    _tree->SetBranchAddress("BToKpipi_piD0_eta",&BToKpipi_piD0_eta);
    _tree->SetBranchAddress("BToKpipi_piD0_phi",&BToKpipi_piD0_phi);
    _tree->SetBranchAddress("BToKpipi_piD0_dz",&BToKpipi_piD0_dz);
  }

  int BToKmumu_info = _tree->SetBranchAddress("nBToKmumu",&nBToKmumu);
  if(BToKmumu_info>=0){
    _tree->SetBranchAddress("BToKmumu_CL_vtx",&BToKmumu_CL_vtx);
    _tree->SetBranchAddress("BToKmumu_mumu_CL_vtx",&BToKmumu_mumu_CL_vtx);
    _tree->SetBranchAddress("BToKmumu_mumu_mass",&BToKmumu_mumu_mass);
    _tree->SetBranchAddress("BToKmumu_mass",&BToKmumu_mass);
    _tree->SetBranchAddress("BToKmumu_kaon_charge",&BToKmumu_kaon_charge);
    _tree->SetBranchAddress("BToKmumu_kaon_pt",&BToKmumu_kaon_pt);
    _tree->SetBranchAddress("BToKmumu_kaon_eta",&BToKmumu_kaon_eta);
    _tree->SetBranchAddress("BToKmumu_kaon_phi",&BToKmumu_kaon_phi);
    _tree->SetBranchAddress("BToKmumu_mu1_charge",&BToKmumu_mu1_charge);
    _tree->SetBranchAddress("BToKmumu_mu1_pt",&BToKmumu_mu1_pt);
    _tree->SetBranchAddress("BToKmumu_mu1_eta",&BToKmumu_mu1_eta);
    _tree->SetBranchAddress("BToKmumu_mu1_phi",&BToKmumu_mu1_phi);
    _tree->SetBranchAddress("BToKmumu_mu1_index",&BToKmumu_mu1_index);
    _tree->SetBranchAddress("BToKmumu_mu2_charge",&BToKmumu_mu2_charge);
    _tree->SetBranchAddress("BToKmumu_mu2_pt",&BToKmumu_mu2_pt);
    _tree->SetBranchAddress("BToKmumu_mu2_eta",&BToKmumu_mu2_eta);
    _tree->SetBranchAddress("BToKmumu_mu2_phi",&BToKmumu_mu2_phi);
    _tree->SetBranchAddress("BToKmumu_mu2_index",&BToKmumu_mu2_index);
  }

  int BToKee_info = _tree->SetBranchAddress("nBToKee",&nBToKee);
  if(BToKee_info>=0){
    _tree->SetBranchAddress("BToKee_CL_vtx",&BToKee_CL_vtx);
    _tree->SetBranchAddress("BToKee_ee_CL_vtx",&BToKee_ee_CL_vtx);
    _tree->SetBranchAddress("BToKee_ee_mass",&BToKee_ee_mass);
    _tree->SetBranchAddress("BToKee_mass",&BToKee_mass);
    _tree->SetBranchAddress("BToKee_kaon_charge",&BToKee_kaon_charge);
    _tree->SetBranchAddress("BToKee_kaon_pt",&BToKee_kaon_pt);
    _tree->SetBranchAddress("BToKee_kaon_eta",&BToKee_kaon_eta);
    _tree->SetBranchAddress("BToKee_kaon_phi",&BToKee_kaon_phi);
    _tree->SetBranchAddress("BToKee_ele1_pt",&BToKee_ele1_pt);
    _tree->SetBranchAddress("BToKee_ele1_eta",&BToKee_ele1_eta);
    _tree->SetBranchAddress("BToKee_ele1_phi",&BToKee_ele1_phi);
    _tree->SetBranchAddress("BToKee_ele1_charge",&BToKee_ele1_charge);
    _tree->SetBranchAddress("BToKee_ele2_pt",&BToKee_ele2_pt);
    _tree->SetBranchAddress("BToKee_ele2_eta",&BToKee_ele2_eta);
    _tree->SetBranchAddress("BToKee_ele2_phi",&BToKee_ele2_phi);
    _tree->SetBranchAddress("BToKee_ele2_charge",&BToKee_ele2_charge);
  }

  int isMC = _tree->SetBranchAddress("nGenPart",&nGenPart);
  if(isMC>=0){
    _tree->SetBranchAddress("GenPart_pdgId",&GenPart_pdgId);
    _tree->SetBranchAddress("GenPart_genPartIdxMother",&GenPart_genPartIdxMother);
    _tree->SetBranchAddress("GenPart_pt",&GenPart_pt);
    _tree->SetBranchAddress("GenPart_eta",&GenPart_eta);
    _tree->SetBranchAddress("GenPart_phi",&GenPart_phi);
  }

  int isBPHTrig = _tree->SetBranchAddress("HLT_Mu8p5_IP3p5_part0",&HLT_Mu8p5_IP3p5_part0);
  if(isBPHTrig>=0){
    _tree->SetBranchAddress("HLT_Mu8p5_IP3p5_part1",&HLT_Mu8p5_IP3p5_part1);
    _tree->SetBranchAddress("HLT_Mu8p5_IP3p5_part2",&HLT_Mu8p5_IP3p5_part2);
    _tree->SetBranchAddress("HLT_Mu8p5_IP3p5_part3",&HLT_Mu8p5_IP3p5_part3);
    _tree->SetBranchAddress("HLT_Mu8p5_IP3p5_part4",&HLT_Mu8p5_IP3p5_part4);
    _tree->SetBranchAddress("HLT_Mu8p5_IP3p5_part5",&HLT_Mu8p5_IP3p5_part5);
    _tree->SetBranchAddress("HLT_Mu10p5_IP3p5_part0",&HLT_Mu10p5_IP3p5_part0);
    _tree->SetBranchAddress("HLT_Mu10p5_IP3p5_part1",&HLT_Mu10p5_IP3p5_part1);
    _tree->SetBranchAddress("HLT_Mu10p5_IP3p5_part2",&HLT_Mu10p5_IP3p5_part2);
    _tree->SetBranchAddress("HLT_Mu10p5_IP3p5_part3",&HLT_Mu10p5_IP3p5_part3);
    _tree->SetBranchAddress("HLT_Mu10p5_IP3p5_part4",&HLT_Mu10p5_IP3p5_part4);
    _tree->SetBranchAddress("HLT_Mu10p5_IP3p5_part5",&HLT_Mu10p5_IP3p5_part5);
    _tree->SetBranchAddress("HLT_Mu9_IP6_part0",&HLT_Mu9_IP6_part0);
    _tree->SetBranchAddress("HLT_Mu9_IP6_part1",&HLT_Mu9_IP6_part1);
    _tree->SetBranchAddress("HLT_Mu9_IP6_part2",&HLT_Mu9_IP6_part2);
    _tree->SetBranchAddress("HLT_Mu9_IP6_part3",&HLT_Mu9_IP6_part3);
    _tree->SetBranchAddress("HLT_Mu9_IP6_part4",&HLT_Mu9_IP6_part4);
    _tree->SetBranchAddress("HLT_Mu9_IP6_part5",&HLT_Mu9_IP6_part5);
    _tree->SetBranchAddress("HLT_Mu8_IP3_part0",&HLT_Mu8_IP3_part0);
    _tree->SetBranchAddress("HLT_Mu8_IP3_part1",&HLT_Mu8_IP3_part1);
    _tree->SetBranchAddress("HLT_Mu8_IP3_part2",&HLT_Mu8_IP3_part2);
    _tree->SetBranchAddress("HLT_Mu8_IP3_part3",&HLT_Mu8_IP3_part3);
    _tree->SetBranchAddress("HLT_Mu8_IP3_part4",&HLT_Mu8_IP3_part4);
    _tree->SetBranchAddress("HLT_Mu8_IP3_part5",&HLT_Mu8_IP3_part5);
  }

  _tree->SetBranchAddress("nTrigObj",&nTrigObj);
  _tree->SetBranchAddress("TrigObj_id",&TrigObj_id);
  _tree->SetBranchAddress("TrigObj_pt",&TrigObj_pt);
  _tree->SetBranchAddress("TrigObj_eta",&TrigObj_eta);
  _tree->SetBranchAddress("TrigObj_phi",&TrigObj_phi);
  _tree->SetBranchAddress("TrigObj_filterBits",&TrigObj_filterBits);
  
  int PFCand_info = _tree->SetBranchAddress("nPFCand",&nPFCand);
  if(PFCand_info){
    _tree->SetBranchAddress("PFCand_pt",&PFCand_pt);
    _tree->SetBranchAddress("PFCand_eta",&PFCand_eta);
    _tree->SetBranchAddress("PFCand_phi",&PFCand_phi);
    _tree->SetBranchAddress("PFCand_mass",&PFCand_mass);
    _tree->SetBranchAddress("PFCand_pdgId",&PFCand_pdgId);
    _tree->SetBranchAddress("PFCand_DCASig",&PFCand_DCASig);
    _tree->SetBranchAddress("PFCand_dz",&PFCand_dz);
  }

}


Int_t NanoAODTree::GetEntry(int entry)
{

  int out = _tree->GetEntry(entry);

  if(nMuon>kMuonMax) return -1;
  if(nElectron>kElectronMax) return -1;
  if(nBToKpipi>kBToKpipiMax) return -1;
  if(nBToKmumu>kBToKmumuMax) return -1;
  if(nBToKee>kBToKeeMax) return -1;
  if(nGenPart>kGenPartMax) return -1;
  if(nTrigObj>kTrigObjMax)  return -1;
  if(nPFCand>kPFCandMax) return -1;

  return out;

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
