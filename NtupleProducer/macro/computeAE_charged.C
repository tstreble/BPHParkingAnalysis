//compute double ratio nnReso/JPsi  event counts following all levels of selections
//muonTag + 
//gen Acceptance +
//gen Efficiency (reco final state matched to gen level) +
//selections for the analysis in data  

//should be ok for electron final state 
//example to run
//electron final state
//root -l computeAE_charged.C'("/vols/cms/amartell/BParking/ntuPROD/ntu_BToKee.root" , "/vols/cms/amartell/BParking/ntuPROD/ntu_BToKJPsiee.root", 1)' 
//muon final state
//root -l computeAE_charged.C'("/vols/cms/amartell/BParking/ntuPROD/ntu_BToKmumu.root" , "/vols/cms/amartell/BParking/ntuPROD/ntu_BToKJPsimumu.root", 0)' 


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
#include "TStyle.h"
#include "TChain.h"

void computeAE_charged(TString nonResonantFile, TString ResonantFile, int isEleFinalState){

  gROOT->Reset();
  gROOT->Macro("setStyle.C");

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);


  std::cout << " isEleFinalState = " << isEleFinalState << std::endl;

  TChain* t1 = new TChain("Events");
  TChain* t2 = new TChain("Events");

  t1->Add(nonResonantFile);
  t2->Add(ResonantFile);


  //muon tag with soft ID https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideMuonIdRun2#Soft_Muon
  int nMuonTag_nonReso = t1->GetEntries();
  int nMuonTag_Reso = t2->GetEntries();

  std::cout << " #muon tag events: nnreso = " << nMuonTag_nonReso << " resonant = " << nMuonTag_Reso << std::endl;

  //if false => move to reco ee invariant mass bins

  //selections
  TString cut_genAcc = "GenPart_e2FromB_index != -1 && GenPart_pt[GenPart_e2FromB_index] > 1. && abs(GenPart_eta[GenPart_e2FromB_index]) < 2.4 && GenPart_e1FromB_index != -1 && GenPart_pt[GenPart_e1FromB_index] > 1. && abs(GenPart_eta[GenPart_e1FromB_index]) < 2.4 && GenPart_KFromB_index != -1 && GenPart_pt[GenPart_KFromB_index] > 1. && abs(GenPart_eta[GenPart_KFromB_index]) < 2.4";
  TString cut_genEff = cut_genAcc + "&& BToKee_gen_index != -1 && BToKee_gen_index == BToKee_sel_index";
  TString cut_chargeEff = cut_genEff + "&& BToKee_ele1_charge[BToKee_gen_index]*BToKee_ele2_charge[BToKee_gen_index] < 0.";
  TString cut_alphaEff = cut_chargeEff + "&& BToKee_cosAlpha[BToKee_gen_index] > 0.99";
  TString cut_vtxCLEff = cut_alphaEff + "&& BToKee_CL_vtx[BToKee_gen_index] > 0.1";
  TString cut_LxyEff  = cut_vtxCLEff + "&& abs(BToKee_Lxy[BToKee_gen_index]) > 6";

  if(!isEleFinalState){
    cut_genAcc = "GenPart_mu2FromB_index != -1 && GenPart_pt[GenPart_mu2FromB_index] > 1. && abs(GenPart_eta[GenPart_mu2FromB_index]) < 2.4 && GenPart_mu1FromB_index != -1 && GenPart_pt[GenPart_mu1FromB_index] > 1. && abs(GenPart_eta[GenPart_mu1FromB_index]) < 2.4 && GenPart_KFromB_index != -1 && GenPart_pt[GenPart_KFromB_index] > 1. && abs(GenPart_eta[GenPart_KFromB_index]) < 2.4";
    cut_genEff = cut_genAcc + "&& BToKmumu_gen_index != -1 && BToKmumu_gen_index == BToKmumu_sel_index";
    cut_chargeEff = cut_genEff + "&& BToKmumu_mu1_charge[BToKmumu_gen_index]*BToKmumu_mu2_charge[BToKmumu_gen_index] < 0.";
    cut_alphaEff = cut_chargeEff + "&& BToKmumu_cosAlpha[BToKmumu_gen_index] > 0.99";
    cut_vtxCLEff = cut_alphaEff + "&& BToKmumu_CL_vtx[BToKmumu_gen_index] > 0.1";
    cut_LxyEff  = cut_vtxCLEff + "&& abs(BToKmumu_Lxy[BToKmumu_gen_index]) > 6";
  }


  std::vector<TString> llMassCut;
  std::vector<TString> llGenMassCut;

  std::vector<float> llMassBoundary;
  llMassBoundary.push_back(0.);
  llMassBoundary.push_back(1.);
  llMassBoundary.push_back(2.5);
  llMassBoundary.push_back(2.9);
  llMassBoundary.push_back(3.3);
  llMassBoundary.push_back(3.58);

  unsigned JPsi_bin = 3;

  for(int ij=0; ij<5; ++ij){
    TString cut = Form("BToKee_gen_index>=0 && BToKee_eeKFit_ee_mass[BToKee_gen_index] > %.2f && BToKee_eeKFit_ee_mass[BToKee_gen_index] < %.2f",
                           llMassBoundary.at(ij), llMassBoundary.at(ij+1));
    TString gencut = Form("BToKee_gen_eeMass > %.2f && BToKee_gen_eeMass < %.2f",
                              llMassBoundary.at(ij), llMassBoundary.at(ij+1));
    if(!isEleFinalState){
      cut = Form("BToKmumu_gen_index>=0 && BToKmumu_mumuKFit_mumu_mass[BToKmumu_gen_index] > %.2f && BToKmumu_mumuKFit_mumu_mass[BToKmumu_gen_index] < %.2f",
		 llMassBoundary.at(ij), llMassBoundary.at(ij+1));
      gencut = Form("BToKmumu_gen_mumuMass > %.2f && BToKmumu_gen_mumuMass < %.2f",
		    llMassBoundary.at(ij), llMassBoundary.at(ij+1));

    }

    llMassCut.push_back(cut);
    llGenMassCut.push_back(gencut);
  }

  vector<TString> cut_name = {"muonTag","genAcc","genEff","recoEff","chargeEff","alphaEff","vtxCLEff","LxyEff"};

  //non resonant first - JPsi last
  vector<float> nEv;
  nEv.resize(cut_name.size(), 0.);
  vector<vector<float> > nEv_cuts; //nEv_cuts[i][j] = #events passing cut j in mass bin i
  nEv_cuts.resize(llMassCut.size()+1,nEv);

  for(unsigned i=0; i<llMassCut.size()+1; i++){

    TString llGenMCut;
    TString llMCut;
    if(i<llMassCut.size()){
      llGenMCut = llGenMassCut[i];
      llMCut = llMassCut[i];
    }
    else{
      llGenMCut = llGenMassCut[JPsi_bin];
      llMCut = llMassCut[JPsi_bin];
    }

    vector<TString> cuts = {llGenMCut,
			    llGenMCut + "&&" + cut_genAcc,
			    llGenMCut + "&&" + cut_genEff,
			    llGenMCut + "&&" + llMCut + "&&" + cut_genEff,
			    llGenMCut + "&&" + llMCut + "&&" + cut_chargeEff,
			    llGenMCut + "&&" + llMCut + "&&" + cut_alphaEff,
			    llGenMCut + "&&" + llMCut + "&&" + cut_vtxCLEff,
			    llGenMCut + "&&" + llMCut + "&&" + cut_LxyEff};

    for(unsigned j=0; j<cuts.size(); j++){
      if(i<llMassCut.size()) nEv_cuts[i][j] = t1->GetEntries(cuts[j]); //Non-resonant
      else nEv_cuts[i][j] = t2->GetEntries(cuts[j]); //Resonant
    }
  }



  //upgrade to 2D plot for future or something better than printout
  for(unsigned i=0; i<llMassCut.size(); i++){
    std::cout << " \n   \t nonResonant ll mass-bin [" << llMassBoundary.at(i) << "-"<< llMassBoundary.at(i+1) << "] \t       \t Resonant (J/Psi) \t   \t \t" << std::endl;

    for(unsigned j=0; j<cut_name.size(); j++){
      if(j==0){
	std::cout << " " << cut_name[j] << " = \t " << nEv_cuts[i][j] << " \t              \t \t        \t " << nEv_cuts[llMassCut.size()][j]  << " \t " << std::endl;
      }
      else{
	double eff = nEv_cuts[i][j]/nEv_cuts[i][0];
	double eff_JPsi = nEv_cuts[llMassCut.size()][j]/nEv_cuts[llMassCut.size()][0];
	std::cout << " " << cut_name[j] << " = \t " << nEv_cuts[i][j] << " \t /muonTag = " << eff
		  << "               \t " << nEv_cuts[llMassCut.size()][j] << " \t  /muonTag = " << eff_JPsi
		  << " \t double ratio = \t" << eff/eff_JPsi << "\n";
      }

    }

  }

}
