//compute double ratio nnReso/JPsi  event counts following all levels of selections
//muonTag + 
//gen Acceptance +
//gen Efficiency (reco final state matched to gen level) +
//selections for the analysis in data  

//should be ok for electron final state 
//example to run
//root -l computeAE_charged.C'("/vols/cms/amartell/BParking/ntuPROD/ntu_BToKee_v18_03_22_and_21.root" , "/vols/cms/amartell/BParking/ntuPROD/ntu_BToKJPsiee_v18_06_4_and_5.root")' 


#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>

#include "TROOT.h"
#include "TSystem.h"
#include "TKey.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TString.h"
#include "TCut.h"
#include "TMath.h"
#include "TApplication.h"
#include "TError.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TPad.h"
#include "TMultiGraph.h"


void computeAE_charged(std::string nonResonantFile, std::string ResonantFile){

  gROOT->Reset();
  gROOT->Macro("~/setStyle.C");
  gROOT->Macro("~/setStyle.C");
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  TChain* t1 = new TChain("Events");
  TChain* t2 = new TChain("Events");

  t1->Add(nonResonantFile.c_str());
  t2->Add(ResonantFile.c_str());


  //muon tag with soft ID https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideMuonIdRun2#Soft_Muon
  int nMuonTag_nonReso = t1->GetEntries();
  int nMuonTag_Reso = t2->GetEntries();

  std::cout << " #moun tag events: nnreso = " << nMuonTag_nonReso << " resonant = " << nMuonTag_Reso << std::endl;

  //if false => move to reco ee invariant mass bins

  //selections
  std::string cut_muonTag = "Muon_sel_index != -1";
  // should require that muonTag does not match with gen probes?
  std::string cut_genAcc = "GenPart_e2FromB_index != -1 && GenPart_pt[GenPart_e2FromB_index] > 1. && abs(GenPart_eta[GenPart_e2FromB_index]) < 2.4 && GenPart_e1FromB_index != -1 && GenPart_pt[GenPart_e1FromB_index] > 1. && abs(GenPart_eta[GenPart_e1FromB_index]) < 2.4 && GenPart_KFromB_index != -1 && GenPart_pt[GenPart_KFromB_index] > 1. && abs(GenPart_eta[GenPart_KFromB_index]) < 2.4";
  std::string cut_genEff = "BToKee_gen_index != -1 && BToKee_gen_index == BToKee_sel_index";
  std::string cut_chargeEff = "BToKee_ele1_charge[BToKee_gen_index]*BToKee_ele2_charge[BToKee_gen_index] < 0.";
  std::string cut_alphaEff = "BToKee_cosAlpha[BToKee_gen_index] > 0.99";
  std::string cut_vtxCLEff = "BToKee_CL_vtx[BToKee_gen_index] > 0.1";  
  std::string cut_DCAEff  = "abs(BToKee_kaon_DCASig[BToKee_gen_index]) > 6";

  std::vector<std::string> eeMassCut;
  std::vector<std::string> eeGenMassCut;

  std::vector<float> eeMassBoundary;
  eeMassBoundary.push_back(0.);
  eeMassBoundary.push_back(1.);
  eeMassBoundary.push_back(2.5);
  eeMassBoundary.push_back(2.9);
  eeMassBoundary.push_back(3.3);
  eeMassBoundary.push_back(3.58);

  for(int ij=0; ij<5; ++ij){
    std::string cut = Form("BToKee_eeKFit_ee_mass[BToKee_gen_index] > %.2f && BToKee_eeKFit_ee_mass[BToKee_gen_index] < %.2f",
                           eeMassBoundary.at(ij), eeMassBoundary.at(ij+1));
    eeMassCut.push_back(cut);
    std::string gencut = Form("BToKee_gen_eeMass > %.2f && BToKee_gen_eeMass < %.2f",
                              eeMassBoundary.at(ij), eeMassBoundary.at(ij+1));
    eeGenMassCut.push_back(gencut);
  }

  //non resonant first - JPsi last
  float nEv_muonTag[6] = {0.};
  float nEv_genAcc[6] = {0.};
  float nEv_genEff[6] = {0.};
  float nEv_recoEff[6] = {0.};  /* same as genEff, but computed in reco ee invarinat mass bins*/
  float nEv_chargeEff[6] = {0.};
  float nEv_alphaEff[6] = {0.};
  float nEv_vtxCLEff[6] = {0.};
  float nEv_DCAEff[6] = {0.};

  for(int ij=0; ij<5; ++ij){
    nEv_muonTag[ij] = t1->Draw("Muon_sel_index", (eeGenMassCut.at(ij)).c_str(), "goff");
    nEv_genAcc[ij] = t1->Draw("Muon_sel_index", (cut_genAcc+" && "+eeGenMassCut.at(ij)).c_str(), "goff");
    nEv_genEff[ij] = t1->Draw("Muon_sel_index", (cut_genAcc+" && "+cut_genEff+" && "+eeGenMassCut.at(ij)).c_str(), "goff");
    nEv_recoEff[ij] = t1->Draw("Muon_sel_index", (cut_genAcc+" && "+cut_genEff+" && "+eeMassCut.at(ij)+" && "+eeGenMassCut.at(ij)).c_str(), "goff");
    nEv_chargeEff[ij] = t1->Draw("Muon_sel_index", (cut_genAcc+" && "+cut_genEff+" && "+cut_chargeEff+" && "+eeMassCut.at(ij)+" && "+eeGenMassCut.at(ij)).c_str(), "goff");
    nEv_alphaEff[ij] = t1->Draw("Muon_sel_index", (cut_genAcc+" && "+cut_genEff+" && "+cut_chargeEff+" && "+cut_alphaEff+" && "+eeMassCut.at(ij)+" && "+eeGenMassCut.at(ij)).c_str(), "goff");
    nEv_vtxCLEff[ij] = t1->Draw("Muon_sel_index", (cut_genAcc+" && "+cut_genEff+" && "+cut_chargeEff+" && "+cut_alphaEff+" && "+cut_vtxCLEff+" && "+eeMassCut.at(ij)+" && "+eeGenMassCut.at(ij)).c_str(), "goff");
    nEv_DCAEff[ij] = t1->Draw("Muon_sel_index", (cut_genAcc+" && "+cut_genEff+" && "+cut_chargeEff+" && "+cut_alphaEff+" && "+cut_vtxCLEff+" && "+cut_DCAEff+" && "+eeMassCut.at(ij)+" && "+eeGenMassCut.at(ij)).c_str(), "goff");
  }
  nEv_muonTag[5] = t2->Draw("Muon_sel_index", (eeGenMassCut.at(3)).c_str(), "goff");
  nEv_genAcc[5] = t2->Draw("Muon_sel_index", (cut_genAcc+" && "+eeGenMassCut.at(3)).c_str(), "goff");
  nEv_genEff[5] = t2->Draw("Muon_sel_index", (cut_genAcc+" && "+cut_genEff+" && "+eeGenMassCut.at(3)).c_str(), "goff");
  nEv_recoEff[5] = t2->Draw("Muon_sel_index", (cut_genAcc+" && "+cut_genEff+" && "+eeMassCut.at(3)+" && "+eeGenMassCut.at(3)).c_str(), "goff");
  nEv_chargeEff[5] = t2->Draw("Muon_sel_index", (cut_genAcc+" && "+cut_genEff+" && "+cut_chargeEff+" && "+eeMassCut.at(3)+" && "+eeGenMassCut.at(3)).c_str(), "goff");
  nEv_alphaEff[5] = t2->Draw("Muon_sel_index", (cut_genAcc+" && "+cut_genEff+" && "+cut_chargeEff+" && "+cut_alphaEff+" && "+eeMassCut.at(3)+" && "+eeGenMassCut.at(3)).c_str(), "goff");
  nEv_vtxCLEff[5] = t2->Draw("Muon_sel_index", (cut_genAcc+" && "+cut_genEff+" && "+cut_chargeEff+" && "+cut_alphaEff+" && "+cut_vtxCLEff+" && "+eeMassCut.at(3)+" && "+eeGenMassCut.at(3)).c_str(), "goff");
  nEv_DCAEff[5] = t2->Draw("Muon_sel_index", (cut_genAcc+" && "+cut_genEff+" && "+cut_chargeEff+" && "+cut_alphaEff+" && "+cut_vtxCLEff+" && "+cut_DCAEff+" && "+eeMassCut.at(3)+" && "+eeGenMassCut.at(3)).c_str(), "goff");
  


  //upgrate to 2D plot for future or something better than printout
  for(int ij=0; ij<5; ++ij){
    std::cout <<        " \n \t  mass bin non-Resonant = " << eeGenMassCut.at(ij) << " \t  Resonant (J/Psi)   \t" << std::endl;
    std::cout << " muonTag = \t " << nEv_muonTag[ij] << " \t " << nEv_muonTag[5]  << " \t " << std::endl;
    std::cout << " genAcc = \t " << nEv_genAcc[ij] << " \t " << nEv_genAcc[5] << " \t double ratio = " << (nEv_genAcc[ij]/nEv_muonTag[ij])/(nEv_genAcc[5]/nEv_muonTag[5]) << "\n";
    std::cout << " x genEff = \t " << nEv_genEff[ij] << " \t " << nEv_genEff[5] << " \t double ratio = " << (nEv_genEff[ij]/nEv_muonTag[ij])/(nEv_genEff[5]/nEv_muonTag[5]) << "\n";
    std::cout << " x recoEff = \t " << nEv_recoEff[ij] << " \t " << nEv_recoEff[5] << " \t  double ratio = " << (nEv_recoEff[ij]/nEv_muonTag[ij])/(nEv_recoEff[5]/nEv_muonTag[5]) << "\n";
    std::cout << " x chargeEff = \t " << nEv_chargeEff[ij] << " \t " << nEv_chargeEff[5] << " \t double ratio = " << (nEv_chargeEff[ij]/nEv_muonTag[ij])/(nEv_chargeEff[5]/nEv_muonTag[5]) << "\n";
    std::cout << " x alphaEff = \t " << nEv_alphaEff[ij] << " \t " << nEv_alphaEff[5] << " \t double ratio = " << (nEv_alphaEff[ij]/nEv_muonTag[ij])/(nEv_alphaEff[5]/nEv_muonTag[5]) << "\n";
    std::cout << " x vtxCLEff = \t " << nEv_vtxCLEff[ij] << " \t " << nEv_vtxCLEff[5] << " \t double ratio = " << (nEv_vtxCLEff[ij]/nEv_muonTag[ij])/(nEv_vtxCLEff[5]/nEv_muonTag[5]) << "\n";
    std::cout << " x DCAEff = \t " << nEv_DCAEff[ij] << " \t " << nEv_DCAEff[5] << " \t double ratio = " << (nEv_DCAEff[ij]/nEv_muonTag[ij])/(nEv_DCAEff[5]/nEv_muonTag[5]) << std::endl;
  }


}
