//compute double ratio nnReso/JPsi  event counts following all levels of selections
//muonTag + 
//gen Acceptance +
//gen Efficiency (reco final state matched to gen level) +
//selections for the analysis in data  

//should be ok for electron final state 
//example to run
//electron final state
//root -l computeAE_charged_plotEff.C'("/vols/cms/amartell/BParking/ntuPROD/ntu_BToKee.root" , "/vols/cms/amartell/BParking/ntuPROD/ntu_BToKJPsiee.root", 1)' 
//muon final state
//root -l computeAE_charged_plotEff.C'("/vols/cms/amartell/BParking/ntuPROD/ntu_BToKmumu.root" , "/vols/cms/amartell/BParking/ntuPROD/ntu_BToKJPsimumu.root", 0)' 


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

void computeAE_charged_plotEff(std::string nonResonantFile, std::string ResonantFile, int isEleFinalState){

  gROOT->Reset();
  gROOT->Macro("setStyle.C");

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);


  std::cout << " isEleFinalState = " << isEleFinalState << std::endl;

  TChain* t1 = new TChain("Events");
  TChain* t2 = new TChain("Events");

  t1->Add(nonResonantFile.c_str());
  t2->Add(ResonantFile.c_str());


  //muon tag with soft ID https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideMuonIdRun2#Soft_Muon
  int nMuonTag_nonReso = t1->GetEntries();
  int nMuonTag_Reso = t2->GetEntries();

  std::cout << " #muon tag events: nnreso = " << nMuonTag_nonReso << " resonant = " << nMuonTag_Reso << std::endl;

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

  if(!isEleFinalState){
    cut_muonTag = "Muon_probe_index != -1";
    cut_genAcc = "GenPart_mu2FromB_index != -1 && GenPart_pt[GenPart_mu2FromB_index] > 1. && abs(GenPart_eta[GenPart_mu2FromB_index]) < 2.4 && GenPart_mu1FromB_index != -1 && GenPart_pt[GenPart_mu1FromB_index] > 1. && abs(GenPart_eta[GenPart_mu1FromB_index]) < 2.4 && GenPart_KFromB_index != -1 && GenPart_pt[GenPart_KFromB_index] > 1. && abs(GenPart_eta[GenPart_KFromB_index]) < 2.4";
    cut_genEff = "BToKmumu_gen_index != -1 && BToKmumu_gen_index == BToKmumu_sel_index";
    cut_chargeEff = "BToKmumu_mu1_charge[BToKmumu_gen_index]*BToKmumu_mu2_charge[BToKmumu_gen_index] < 0.";
    cut_alphaEff = "BToKmumu_cosAlpha[BToKmumu_gen_index] > 0.99";
    cut_vtxCLEff = "BToKmumu_CL_vtx[BToKmumu_gen_index] > 0.1";
    cut_DCAEff  = "abs(BToKmumu_kaon_DCASig[BToKmumu_gen_index]) > 6";
  }

  std::vector<std::string> llMassCut;
  std::vector<std::string> llGenMassCut;

  std::vector<float> llMassBoundary;
  llMassBoundary.push_back(0.);
  llMassBoundary.push_back(1.);
  llMassBoundary.push_back(2.5);
  llMassBoundary.push_back(2.9);
  llMassBoundary.push_back(3.3);
  llMassBoundary.push_back(3.58);

  for(int ij=0; ij<5; ++ij){
    std::string cut = Form("BToKee_eeKFit_ee_mass[BToKee_gen_index] > %.2f && BToKee_eeKFit_ee_mass[BToKee_gen_index] < %.2f",
                           llMassBoundary.at(ij), llMassBoundary.at(ij+1));
    std::string gencut = Form("BToKee_gen_eeMass > %.2f && BToKee_gen_eeMass < %.2f",
                              llMassBoundary.at(ij), llMassBoundary.at(ij+1));
    if(!isEleFinalState){
      cut = Form("BToKmumu_mumuKFit_mumu_mass[BToKmumu_gen_index] > %.2f && BToKmumu_mumuKFit_mumu_mass[BToKmumu_gen_index] < %.2f",
		 llMassBoundary.at(ij), llMassBoundary.at(ij+1));
      gencut = Form("BToKmumu_gen_mumuMass > %.2f && BToKmumu_gen_mumuMass < %.2f",
		    llMassBoundary.at(ij), llMassBoundary.at(ij+1));

    }

    llMassCut.push_back(cut);
    llGenMassCut.push_back(gencut);
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


  std::vector<std::string> cutFlow;
  cutFlow.push_back("muonTag");
  cutFlow.push_back("genAcc");
  cutFlow.push_back("genEff");
  cutFlow.push_back("recoEff");
  cutFlow.push_back("chargeEff");
  cutFlow.push_back("alphaEff");
  cutFlow.push_back("vtxCLEff");
  cutFlow.push_back("DCAEff");

  TH2F* efficiency2D[8][6];
  for(int ij=0; ij<5; ++ij){
    efficiency2D[0][ij] = new TH2F(Form((cutFlow.at(0)+"_2D_llmassBin_%.2f-%.2f").c_str(), llMassBoundary.at(ij), llMassBoundary.at(ij+1)), "", 15, -3., 3., 20, 0., 10.);
    efficiency2D[1][ij] = new TH2F(Form((cutFlow.at(1)+"_2D_llmassBin_%.2f-%.2f").c_str(), llMassBoundary.at(ij), llMassBoundary.at(ij+1)), "", 15, -3., 3., 20, 0., 10.);
    efficiency2D[2][ij] = new TH2F(Form((cutFlow.at(2)+"_2D_llmassBin_%.2f-%.2f").c_str(), llMassBoundary.at(ij), llMassBoundary.at(ij+1)), "", 15, -3., 3., 20, 0., 10.);
    efficiency2D[3][ij] = new TH2F(Form((cutFlow.at(3)+"_2D_llmassBin_%.2f-%.2f").c_str(), llMassBoundary.at(ij), llMassBoundary.at(ij+1)), "", 15, -3., 3., 20, 0., 10.);
    efficiency2D[4][ij] = new TH2F(Form((cutFlow.at(4)+"_2D_llmassBin_%.2f-%.2f").c_str(), llMassBoundary.at(ij), llMassBoundary.at(ij+1)), "", 15, -3., 3., 20, 0., 10.);
    efficiency2D[5][ij] = new TH2F(Form((cutFlow.at(5)+"_2D_llmassBin_%.2f-%.2f").c_str(), llMassBoundary.at(ij), llMassBoundary.at(ij+1)), "", 15, -3., 3., 20, 0., 10.);
    efficiency2D[6][ij] = new TH2F(Form((cutFlow.at(6)+"_2D_llmassBin_%.2f-%.2f").c_str(), llMassBoundary.at(ij), llMassBoundary.at(ij+1)), "", 15, -3., 3., 20, 0., 10.);
    efficiency2D[7][ij] = new TH2F(Form((cutFlow.at(7)+"_2D_llmassBin_%.2f-%.2f").c_str(), llMassBoundary.at(ij), llMassBoundary.at(ij+1)), "", 15, -3., 3., 20, 0., 10.);

    efficiency2D[0][ij]->Sumw2();
    efficiency2D[1][ij]->Sumw2();
    efficiency2D[2][ij]->Sumw2();
    efficiency2D[3][ij]->Sumw2();
    efficiency2D[4][ij]->Sumw2();
    efficiency2D[5][ij]->Sumw2();
    efficiency2D[6][ij]->Sumw2();
    efficiency2D[7][ij]->Sumw2();


    if(isEleFinalState){
      nEv_muonTag[ij] = t1->Draw(Form("GenPart_pt[GenPart_e2FromB_index]:GenPart_eta[GenPart_e2FromB_index] >>%s", efficiency2D[0][ij]->GetName()), (llGenMassCut.at(ij)).c_str());
      nEv_genAcc[ij] = t1->Draw(Form("GenPart_pt[GenPart_e2FromB_index]:GenPart_eta[GenPart_e2FromB_index] >>%s", efficiency2D[1][ij]->GetName()), (cut_genAcc+" && "+llGenMassCut.at(ij)).c_str());
      nEv_genEff[ij] = t1->Draw(Form("GenPart_pt[GenPart_e2FromB_index]:GenPart_eta[GenPart_e2FromB_index] >>%s", efficiency2D[2][ij]->GetName()), 
				(cut_genAcc+" && "+cut_genEff+" && "+llGenMassCut.at(ij)).c_str());
      nEv_recoEff[ij] = t1->Draw(Form("GenPart_pt[GenPart_e2FromB_index]:GenPart_eta[GenPart_e2FromB_index] >>%s", efficiency2D[3][ij]->GetName()), 
				 (cut_genAcc+" && "+cut_genEff+" && "+llMassCut.at(ij)+" && "+llGenMassCut.at(ij)).c_str());
      nEv_chargeEff[ij] = t1->Draw(Form("GenPart_pt[GenPart_e2FromB_index]:GenPart_eta[GenPart_e2FromB_index] >>%s", efficiency2D[4][ij]->GetName()), 
				   (cut_genAcc+" && "+cut_genEff+" && "+cut_chargeEff+" && "+llMassCut.at(ij)+" && "+llGenMassCut.at(ij)).c_str());
      nEv_alphaEff[ij] = t1->Draw(Form("GenPart_pt[GenPart_e2FromB_index]:GenPart_eta[GenPart_e2FromB_index] >>%s", efficiency2D[5][ij]->GetName()), 
				  (cut_genAcc+" && "+cut_genEff+" && "+cut_chargeEff+" && "+cut_alphaEff+" && "+llMassCut.at(ij)+" && "+llGenMassCut.at(ij)).c_str());
      nEv_vtxCLEff[ij] = t1->Draw(Form("GenPart_pt[GenPart_e2FromB_index]:GenPart_eta[GenPart_e2FromB_index] >>%s", efficiency2D[6][ij]->GetName()), 
				  (cut_genAcc+" && "+cut_genEff+" && "+cut_chargeEff+" && "+cut_alphaEff+" && "+cut_vtxCLEff+" && "+llMassCut.at(ij)+" && "+llGenMassCut.at(ij)).c_str());
      nEv_DCAEff[ij] = t1->Draw(Form("GenPart_pt[GenPart_e2FromB_index]:GenPart_eta[GenPart_e2FromB_index] >>%s", efficiency2D[7][ij]->GetName()), 
				(cut_genAcc+" && "+cut_genEff+" && "+cut_chargeEff+" && "+cut_alphaEff+" && "+cut_vtxCLEff+" && "+cut_DCAEff+" && "+llMassCut.at(ij)+" && "+llGenMassCut.at(ij)).c_str());
    }
    else{
      nEv_muonTag[ij] = t1->Draw(Form("GenPart_pt[GenPart_mu2FromB_index]:GenPart_eta[GenPart_mu2FromB_index] >>%s", efficiency2D[0][ij]->GetName()), (llGenMassCut.at(ij)).c_str());
      nEv_genAcc[ij] = t1->Draw(Form("GenPart_pt[GenPart_mu2FromB_index]:GenPart_eta[GenPart_mu2FromB_index] >>%s", efficiency2D[1][ij]->GetName()), (cut_genAcc+" && "+llGenMassCut.at(ij)).c_str());
      nEv_genEff[ij] = t1->Draw(Form("GenPart_pt[GenPart_mu2FromB_index]:GenPart_eta[GenPart_mu2FromB_index] >>%s", efficiency2D[2][ij]->GetName()), 
				(cut_genAcc+" && "+cut_genEff+" && "+llGenMassCut.at(ij)).c_str());
      nEv_recoEff[ij] = t1->Draw(Form("GenPart_pt[GenPart_mu2FromB_index]:GenPart_eta[GenPart_mu2FromB_index] >>%s", efficiency2D[3][ij]->GetName()), 
				 (cut_genAcc+" && "+cut_genEff+" && "+llMassCut.at(ij)+" && "+llGenMassCut.at(ij)).c_str());
      nEv_chargeEff[ij] = t1->Draw(Form("GenPart_pt[GenPart_mu2FromB_index]:GenPart_eta[GenPart_mu2FromB_index] >>%s", efficiency2D[4][ij]->GetName()), 
				   (cut_genAcc+" && "+cut_genEff+" && "+cut_chargeEff+" && "+llMassCut.at(ij)+" && "+llGenMassCut.at(ij)).c_str());
      nEv_alphaEff[ij] = t1->Draw(Form("GenPart_pt[GenPart_mu2FromB_index]:GenPart_eta[GenPart_mu2FromB_index] >>%s", efficiency2D[5][ij]->GetName()), 
				  (cut_genAcc+" && "+cut_genEff+" && "+cut_chargeEff+" && "+cut_alphaEff+" && "+llMassCut.at(ij)+" && "+llGenMassCut.at(ij)).c_str());
      nEv_vtxCLEff[ij] = t1->Draw(Form("GenPart_pt[GenPart_mu2FromB_index]:GenPart_eta[GenPart_mu2FromB_index] >>%s", efficiency2D[6][ij]->GetName()), 
				  (cut_genAcc+" && "+cut_genEff+" && "+cut_chargeEff+" && "+cut_alphaEff+" && "+cut_vtxCLEff+" && "+llMassCut.at(ij)+" && "+llGenMassCut.at(ij)).c_str());
      nEv_DCAEff[ij] = t1->Draw(Form("GenPart_pt[GenPart_mu2FromB_index]:GenPart_eta[GenPart_mu2FromB_index] >>%s", efficiency2D[7][ij]->GetName()), 
				(cut_genAcc+" && "+cut_genEff+" && "+cut_chargeEff+" && "+cut_alphaEff+" && "+cut_vtxCLEff+" && "+cut_DCAEff+" && "+llMassCut.at(ij)+" && "+llGenMassCut.at(ij)).c_str());
    }
  }

    efficiency2D[0][5] = new TH2F(Form((cutFlow.at(0)+"_2D_JPsimassBin_%.2f-%.2f").c_str(), llMassBoundary.at(3), llMassBoundary.at(3+1)), "", 15, -3., 3., 20, 0., 10.);
    efficiency2D[1][5] = new TH2F(Form((cutFlow.at(1)+"_2D_JPsimassBin_%.2f-%.2f").c_str(), llMassBoundary.at(3), llMassBoundary.at(3+1)), "", 15, -3., 3., 20, 0., 10.);
    efficiency2D[2][5] = new TH2F(Form((cutFlow.at(2)+"_2D_JPsimassBin_%.2f-%.2f").c_str(), llMassBoundary.at(3), llMassBoundary.at(3+1)), "", 15, -3., 3., 20, 0., 10.);
    efficiency2D[3][5] = new TH2F(Form((cutFlow.at(3)+"_2D_JPsimassBin_%.2f-%.2f").c_str(), llMassBoundary.at(3), llMassBoundary.at(3+1)), "", 15, -3., 3., 20, 0., 10.);
    efficiency2D[4][5] = new TH2F(Form((cutFlow.at(4)+"_2D_JPsimassBin_%.2f-%.2f").c_str(), llMassBoundary.at(3), llMassBoundary.at(3+1)), "", 15, -3., 3., 20, 0., 10.);
    efficiency2D[5][5] = new TH2F(Form((cutFlow.at(5)+"_2D_JPsimassBin_%.2f-%.2f").c_str(), llMassBoundary.at(3), llMassBoundary.at(3+1)), "", 15, -3., 3., 20, 0., 10.);
    efficiency2D[6][5] = new TH2F(Form((cutFlow.at(6)+"_2D_JPsimassBin_%.2f-%.2f").c_str(), llMassBoundary.at(3), llMassBoundary.at(3+1)), "", 15, -3., 3., 20, 0., 10.);
    efficiency2D[7][5] = new TH2F(Form((cutFlow.at(7)+"_2D_JPsimassBin_%.2f-%.2f").c_str(), llMassBoundary.at(3), llMassBoundary.at(3+1)), "", 15, -3., 3., 20, 0., 10.);

    efficiency2D[0][5]->Sumw2();
    efficiency2D[1][5]->Sumw2();
    efficiency2D[2][5]->Sumw2();
    efficiency2D[3][5]->Sumw2();
    efficiency2D[4][5]->Sumw2();
    efficiency2D[5][5]->Sumw2();
    efficiency2D[6][5]->Sumw2();
    efficiency2D[7][5]->Sumw2();


    if(isEleFinalState){
      nEv_muonTag[5] = t2->Draw(Form("GenPart_pt[GenPart_e2FromB_index]:GenPart_eta[GenPart_e2FromB_index] >>%s", efficiency2D[0][5]->GetName()), (llGenMassCut.at(3)).c_str());
      nEv_genAcc[5] = t2->Draw(Form("GenPart_pt[GenPart_e2FromB_index]:GenPart_eta[GenPart_e2FromB_index] >>%s", efficiency2D[1][5]->GetName()), 
			       (cut_genAcc+" && "+llGenMassCut.at(3)).c_str());
      nEv_genEff[5] = t2->Draw(Form("GenPart_pt[GenPart_e2FromB_index]:GenPart_eta[GenPart_e2FromB_index] >>%s", efficiency2D[2][5]->GetName()), 
			       (cut_genAcc+" && "+cut_genEff+" && "+llGenMassCut.at(3)).c_str());
      nEv_recoEff[5] = t2->Draw(Form("GenPart_pt[GenPart_e2FromB_index]:GenPart_eta[GenPart_e2FromB_index] >>%s", efficiency2D[3][5]->GetName()), 
				(cut_genAcc+" && "+cut_genEff+" && "+llMassCut.at(3)+" && "+llGenMassCut.at(3)).c_str());
      nEv_chargeEff[5] = t2->Draw(Form("GenPart_pt[GenPart_e2FromB_index]:GenPart_eta[GenPart_e2FromB_index] >>%s", efficiency2D[4][5]->GetName()), 
				  (cut_genAcc+" && "+cut_genEff+" && "+cut_chargeEff+" && "+llMassCut.at(3)+" && "+llGenMassCut.at(3)).c_str());
      nEv_alphaEff[5] = t2->Draw(Form("GenPart_pt[GenPart_e2FromB_index]:GenPart_eta[GenPart_e2FromB_index] >>%s", efficiency2D[5][5]->GetName()), 
				 (cut_genAcc+" && "+cut_genEff+" && "+cut_chargeEff+" && "+cut_alphaEff+" && "+llMassCut.at(3)+" && "+llGenMassCut.at(3)).c_str());
      nEv_vtxCLEff[5] = t2->Draw(Form("GenPart_pt[GenPart_e2FromB_index]:GenPart_eta[GenPart_e2FromB_index] >>%s", efficiency2D[6][5]->GetName()), 
				 (cut_genAcc+" && "+cut_genEff+" && "+cut_chargeEff+" && "+cut_alphaEff+" && "+cut_vtxCLEff+" && "+llMassCut.at(3)+" && "+llGenMassCut.at(3)).c_str());
      nEv_DCAEff[5] = t2->Draw(Form("GenPart_pt[GenPart_e2FromB_index]:GenPart_eta[GenPart_e2FromB_index] >>%s", efficiency2D[7][5]->GetName()), 
			       (cut_genAcc+" && "+cut_genEff+" && "+cut_chargeEff+" && "+cut_alphaEff+" && "+cut_vtxCLEff+" && "+cut_DCAEff+" && "+llMassCut.at(3)+" && "+llGenMassCut.at(3)).c_str());
    }
    else{
      nEv_muonTag[5] = t2->Draw(Form("GenPart_pt[GenPart_mu2FromB_index]:GenPart_eta[GenPart_mu2FromB_index] >>%s", efficiency2D[0][5]->GetName()), (llGenMassCut.at(3)).c_str());
      nEv_genAcc[5] = t2->Draw(Form("GenPart_pt[GenPart_mu2FromB_index]:GenPart_eta[GenPart_mu2FromB_index] >>%s", efficiency2D[1][5]->GetName()), 
			       (cut_genAcc+" && "+llGenMassCut.at(3)).c_str());
      nEv_genEff[5] = t2->Draw(Form("GenPart_pt[GenPart_mu2FromB_index]:GenPart_eta[GenPart_mu2FromB_index] >>%s", efficiency2D[2][5]->GetName()), 
			       (cut_genAcc+" && "+cut_genEff+" && "+llGenMassCut.at(3)).c_str());
      nEv_recoEff[5] = t2->Draw(Form("GenPart_pt[GenPart_mu2FromB_index]:GenPart_eta[GenPart_mu2FromB_index] >>%s", efficiency2D[3][5]->GetName()), 
				(cut_genAcc+" && "+cut_genEff+" && "+llMassCut.at(3)+" && "+llGenMassCut.at(3)).c_str());
      nEv_chargeEff[5] = t2->Draw(Form("GenPart_pt[GenPart_mu2FromB_index]:GenPart_eta[GenPart_mu2FromB_index] >>%s", efficiency2D[4][5]->GetName()), 
				  (cut_genAcc+" && "+cut_genEff+" && "+cut_chargeEff+" && "+llMassCut.at(3)+" && "+llGenMassCut.at(3)).c_str());
      nEv_alphaEff[5] = t2->Draw(Form("GenPart_pt[GenPart_mu2FromB_index]:GenPart_eta[GenPart_mu2FromB_index] >>%s", efficiency2D[5][5]->GetName()), 
				 (cut_genAcc+" && "+cut_genEff+" && "+cut_chargeEff+" && "+cut_alphaEff+" && "+llMassCut.at(3)+" && "+llGenMassCut.at(3)).c_str());
      nEv_vtxCLEff[5] = t2->Draw(Form("GenPart_pt[GenPart_mu2FromB_index]:GenPart_eta[GenPart_mu2FromB_index] >>%s", efficiency2D[6][5]->GetName()), 
				 (cut_genAcc+" && "+cut_genEff+" && "+cut_chargeEff+" && "+cut_alphaEff+" && "+cut_vtxCLEff+" && "+llMassCut.at(3)+" && "+llGenMassCut.at(3)).c_str());
      nEv_DCAEff[5] = t2->Draw(Form("GenPart_pt[GenPart_mu2FromB_index]:GenPart_eta[GenPart_mu2FromB_index] >>%s", efficiency2D[7][5]->GetName()), 
			       (cut_genAcc+" && "+cut_genEff+" && "+cut_chargeEff+" && "+cut_alphaEff+" && "+cut_vtxCLEff+" && "+cut_DCAEff+" && "+llMassCut.at(3)+" && "+llGenMassCut.at(3)).c_str());
    }



    //1D plotting;
    TH1F* efficiency[8][6];
    TH2F* efficiencyRatio[8][6];
    for(int ij=0; ij<6; ++ij){
      if(ij == 5){
	efficiency[0][5] = new TH1F(Form((cutFlow.at(0)+"_JPsimassBin_%.2f-%.2f").c_str(), llMassBoundary.at(3), llMassBoundary.at(3+1)), "", 20, 0., 10.);
	efficiency[1][5] = new TH1F(Form((cutFlow.at(1)+"_JPsimassBin_%.2f-%.2f").c_str(), llMassBoundary.at(3), llMassBoundary.at(3+1)), "", 20, 0., 10.);
	efficiency[2][5] = new TH1F(Form((cutFlow.at(2)+"_JPsimassBin_%.2f-%.2f").c_str(), llMassBoundary.at(3), llMassBoundary.at(3+1)), "", 20, 0., 10.);
	efficiency[3][5] = new TH1F(Form((cutFlow.at(3)+"_JPsimassBin_%.2f-%.2f").c_str(), llMassBoundary.at(3), llMassBoundary.at(3+1)), "", 20, 0., 10.);
	efficiency[4][5] = new TH1F(Form((cutFlow.at(4)+"_JPsimassBin_%.2f-%.2f").c_str(), llMassBoundary.at(3), llMassBoundary.at(3+1)), "", 20, 0., 10.);
	efficiency[5][5] = new TH1F(Form((cutFlow.at(5)+"_JPsimassBin_%.2f-%.2f").c_str(), llMassBoundary.at(3), llMassBoundary.at(3+1)), "", 20, 0., 10.);
	efficiency[6][5] = new TH1F(Form((cutFlow.at(6)+"_JPsimassBin_%.2f-%.2f").c_str(), llMassBoundary.at(3), llMassBoundary.at(3+1)), "", 20, 0., 10.);
	efficiency[7][5] = new TH1F(Form((cutFlow.at(7)+"_JPsimassBin_%.2f-%.2f").c_str(), llMassBoundary.at(3), llMassBoundary.at(3+1)), "", 20, 0., 10.);
	
	efficiencyRatio[0][5] = (TH2F*)(efficiency2D[0][5]->Clone(Form(("ratioEff_"+cutFlow.at(0)+"_2D_JPsimassBin_%.2f-%.2f").c_str(), llMassBoundary.at(3), llMassBoundary.at(3+1))) );
	efficiencyRatio[1][5] = (TH2F*)(efficiency2D[1][5]->Clone(Form(("ratioEff_"+cutFlow.at(1)+"_2D_JPsimassBin_%.2f-%.2f").c_str(), llMassBoundary.at(3), llMassBoundary.at(3+1))) );
	efficiencyRatio[2][5] = (TH2F*)(efficiency2D[2][5]->Clone(Form(("ratioEff_"+cutFlow.at(2)+"_2D_JPsimassBin_%.2f-%.2f").c_str(), llMassBoundary.at(3), llMassBoundary.at(3+1))) );
	efficiencyRatio[3][5] = (TH2F*)(efficiency2D[3][5]->Clone(Form(("ratioEff_"+cutFlow.at(3)+"_2D_JPsimassBin_%.2f-%.2f").c_str(), llMassBoundary.at(3), llMassBoundary.at(3+1))) );
	efficiencyRatio[4][5] = (TH2F*)(efficiency2D[4][5]->Clone(Form(("ratioEff_"+cutFlow.at(4)+"_2D_JPsimassBin_%.2f-%.2f").c_str(), llMassBoundary.at(3), llMassBoundary.at(3+1))) );
	efficiencyRatio[5][5] = (TH2F*)(efficiency2D[5][5]->Clone(Form(("ratioEff_"+cutFlow.at(3)+"_2D_JPsimassBin_%.2f-%.2f").c_str(), llMassBoundary.at(3), llMassBoundary.at(3+1))) );
	efficiencyRatio[6][5] = (TH2F*)(efficiency2D[6][5]->Clone(Form(("ratioEff_"+cutFlow.at(6)+"_2D_JPsimassBin_%.2f-%.2f").c_str(), llMassBoundary.at(3), llMassBoundary.at(3+1))) );
	efficiencyRatio[7][5] = (TH2F*)(efficiency2D[7][5]->Clone(Form(("ratioEff_"+cutFlow.at(7)+"_2D_JPsimassBin_%.2f-%.2f").c_str(), llMassBoundary.at(3), llMassBoundary.at(3+1))) );     
      }
      else{
	efficiency[0][ij] = new TH1F(Form((cutFlow.at(0)+"_llmassBin_%.2f-%.2f").c_str(), llMassBoundary.at(ij), llMassBoundary.at(ij+1)), "", 20, 0., 10.);
	efficiency[1][ij] = new TH1F(Form((cutFlow.at(1)+"_llmassBin_%.2f-%.2f").c_str(), llMassBoundary.at(ij), llMassBoundary.at(ij+1)), "", 20, 0., 10.);
	efficiency[2][ij] = new TH1F(Form((cutFlow.at(2)+"_llmassBin_%.2f-%.2f").c_str(), llMassBoundary.at(ij), llMassBoundary.at(ij+1)), "", 20, 0., 10.);
	efficiency[3][ij] = new TH1F(Form((cutFlow.at(3)+"_llmassBin_%.2f-%.2f").c_str(), llMassBoundary.at(ij), llMassBoundary.at(ij+1)), "", 20, 0., 10.);
	efficiency[4][ij] = new TH1F(Form((cutFlow.at(4)+"_llmassBin_%.2f-%.2f").c_str(), llMassBoundary.at(ij), llMassBoundary.at(ij+1)), "", 20, 0., 10.);
	efficiency[5][ij] = new TH1F(Form((cutFlow.at(5)+"_llmassBin_%.2f-%.2f").c_str(), llMassBoundary.at(ij), llMassBoundary.at(ij+1)), "", 20, 0., 10.);
	efficiency[6][ij] = new TH1F(Form((cutFlow.at(6)+"_llmassBin_%.2f-%.2f").c_str(), llMassBoundary.at(ij), llMassBoundary.at(ij+1)), "", 20, 0., 10.);
	efficiency[7][ij] = new TH1F(Form((cutFlow.at(7)+"_llmassBin_%.2f-%.2f").c_str(), llMassBoundary.at(ij), llMassBoundary.at(ij+1)), "", 20, 0., 10.);
	
	efficiencyRatio[0][ij] = (TH2F*)(efficiency2D[0][ij]->Clone(Form(("ratioEff_"+cutFlow.at(0)+"_2D_llmassBin_%.2f-%.2f").c_str(), llMassBoundary.at(ij), llMassBoundary.at(ij+1))) );
	efficiencyRatio[1][ij] = (TH2F*)(efficiency2D[1][ij]->Clone(Form(("ratioEff_"+cutFlow.at(1)+"_2D_llmassBin_%.2f-%.2f").c_str(), llMassBoundary.at(ij), llMassBoundary.at(ij+1))) );
	efficiencyRatio[2][ij] = (TH2F*)(efficiency2D[2][ij]->Clone(Form(("ratioEff_"+cutFlow.at(2)+"_2D_llmassBin_%.2f-%.2f").c_str(), llMassBoundary.at(ij), llMassBoundary.at(ij+1))) );
	efficiencyRatio[3][ij] = (TH2F*)(efficiency2D[3][ij]->Clone(Form(("ratioEff_"+cutFlow.at(3)+"_2D_llmassBin_%.2f-%.2f").c_str(), llMassBoundary.at(ij), llMassBoundary.at(ij+1))) );
	efficiencyRatio[4][ij] = (TH2F*)(efficiency2D[4][ij]->Clone(Form(("ratioEff_"+cutFlow.at(4)+"_2D_llmassBin_%.2f-%.2f").c_str(), llMassBoundary.at(ij), llMassBoundary.at(ij+1))) );
	efficiencyRatio[5][ij] = (TH2F*)(efficiency2D[5][ij]->Clone(Form(("ratioEff_"+cutFlow.at(5)+"_2D_llmassBin_%.2f-%.2f").c_str(), llMassBoundary.at(ij), llMassBoundary.at(ij+1))) );
	efficiencyRatio[6][ij] = (TH2F*)(efficiency2D[6][ij]->Clone(Form(("ratioEff_"+cutFlow.at(6)+"_2D_llmassBin_%.2f-%.2f").c_str(), llMassBoundary.at(ij), llMassBoundary.at(ij+1))) );
	efficiencyRatio[7][ij] = (TH2F*)(efficiency2D[7][ij]->Clone(Form(("ratioEff_"+cutFlow.at(7)+"_2D_llmassBin_%.2f-%.2f").c_str(), llMassBoundary.at(ij), llMassBoundary.at(ij+1))) );
      }

      efficiency[0][ij] = (TH1F*)(efficiency2D[0][ij]->ProjectionY())->Clone(efficiency[0][ij]->GetName());
      efficiency[1][ij] = (TH1F*)(efficiency2D[1][ij]->ProjectionY())->Clone(efficiency[1][ij]->GetName());
      efficiency[2][ij] = (TH1F*)(efficiency2D[2][ij]->ProjectionY())->Clone(efficiency[2][ij]->GetName());
      efficiency[3][ij] = (TH1F*)(efficiency2D[3][ij]->ProjectionY())->Clone(efficiency[3][ij]->GetName());
      efficiency[4][ij] = (TH1F*)(efficiency2D[4][ij]->ProjectionY())->Clone(efficiency[4][ij]->GetName());
      efficiency[5][ij] = (TH1F*)(efficiency2D[5][ij]->ProjectionY())->Clone(efficiency[5][ij]->GetName());
      efficiency[6][ij] = (TH1F*)(efficiency2D[6][ij]->ProjectionY())->Clone(efficiency[6][ij]->GetName());
      efficiency[7][ij] = (TH1F*)(efficiency2D[7][ij]->ProjectionY())->Clone(efficiency[7][ij]->GetName());

      efficiencyRatio[0][ij]->Divide(efficiency2D[0][ij]);
      efficiencyRatio[1][ij]->Divide(efficiency2D[0][ij]);
      efficiencyRatio[2][ij]->Divide(efficiency2D[0][ij]);
      efficiencyRatio[3][ij]->Divide(efficiency2D[0][ij]);
      efficiencyRatio[4][ij]->Divide(efficiency2D[0][ij]);
      efficiencyRatio[5][ij]->Divide(efficiency2D[0][ij]);
      efficiencyRatio[6][ij]->Divide(efficiency2D[0][ij]);
      efficiencyRatio[7][ij]->Divide(efficiency2D[0][ij]);
    }

    ///plotting 
    TLegend *legM = new TLegend(0.72,0.58,0.93,0.93,NULL,"brNDC");
    legM->SetTextFont(42);
    legM->SetFillColor(kWhite);
    legM->SetLineColor(kWhite);
    legM->SetShadowColor(kWhite);
    //legM->SetFillStyle(0);
    legM->SetTextSize(0.03);


    int iColors[8] = {kYellow+3, kViolet, kRed, kGreen+1, kRed+1, kBlue, kCyan+2, kMagenta};

    for(int ij=0; ij<8; ++ij){
      for(int kl=0; kl<6; ++kl){
	efficiency[ij][kl]->SetLineColor(iColors[ij]);
	efficiency[ij][kl]->SetLineWidth(2);
      }
      legM->AddEntry(efficiency[ij][0], cutFlow.at(ij).c_str(), "l");
    }


    TCanvas* c1[5];
    for(int ij=0; ij<5; ++ij){
      legM->SetHeader(Form("mass bin [%.2f-%.2f]", llMassBoundary.at(ij), llMassBoundary.at(ij+1)), "c" );

      c1[ij] = new TCanvas();
      c1[ij]->cd();
      gPad->SetLogy();
      if(ij != 3){
	if(isEleFinalState){
	  efficiency[0][ij]->GetXaxis()->SetTitle(" ele2 pT (GeV)");
	  efficiency[0][ij]->GetYaxis()->SetTitle(" nEvents (Kee)");
	}
	else{
	  efficiency[0][ij]->GetXaxis()->SetTitle(" muon2 pT (GeV)");
	  efficiency[0][ij]->GetYaxis()->SetTitle(" nEvents (Kmumu)");
	}
	efficiency[0][ij]->GetXaxis()->SetTitleOffset(1.05);
	efficiency[0][ij]->GetYaxis()->SetTitleOffset(1.05);
	efficiency[0][ij]->Draw("hist");
	for(int kl=1; kl<8;++kl) efficiency[kl][ij]->Draw("hist, same");
	legM->Draw("hist, same");
	if(isEleFinalState)
	  c1[ij]->Print(Form("plots/EffxAcc/vsPt/NonResonant_eeK_massBin_%.2f-%.2f.png", llMassBoundary.at(ij), llMassBoundary.at(ij+1)), "png");
	else
	  c1[ij]->Print(Form("plots/EffxAcc/vsPt/NonResonant_mumuK_massBin_%.2f-%.2f.png", llMassBoundary.at(ij), llMassBoundary.at(ij+1)), "png");
      }
      else{
	if(isEleFinalState){
	  efficiency[0][5]->GetXaxis()->SetTitle(" ele2 pT (GeV)");
	  efficiency[0][5]->GetYaxis()->SetTitle(" nEvents (Kee)");
	}
	else{
	  efficiency[0][5]->GetXaxis()->SetTitle(" muon2 pT (GeV)");
	  efficiency[0][5]->GetYaxis()->SetTitle(" nEvents (Kmumu)");
	}
	efficiency[0][5]->GetXaxis()->SetTitleOffset(1.05);
	efficiency[0][5]->GetYaxis()->SetTitleOffset(1.05);
	efficiency[0][5]->Draw("hist");
	for(int kl=1; kl<8;++kl) efficiency[kl][5]->Draw("hist, same");
	legM->Draw("hist, same");
	if(isEleFinalState)
	  c1[ij]->Print(Form("plots/EffxAcc/vsPt/JPsi_eeK_massBin_%.2f-%.2f.png", llMassBoundary.at(ij), llMassBoundary.at(ij+1)), "png");
	else
	  c1[ij]->Print(Form("plots/EffxAcc/vsPt/JPsi_mumuK_massBin_%.2f-%.2f.png", llMassBoundary.at(ij), llMassBoundary.at(ij+1)), "png");
      }
    }


    //print 2D eff
    //efficiencyRatio[0][ij]
    TCanvas* c2[8][5];
    for(int kl=0; kl<4; ++kl){
    for(int ij=0; ij<5; ++ij){

      c2[kl][ij] = new TCanvas();
      c2[kl][ij]->cd();
      if(ij != 3){
        if(isEleFinalState){
          efficiencyRatio[kl][ij]->GetXaxis()->SetTitle(" ele2 #eta");
          efficiencyRatio[kl][ij]->GetYaxis()->SetTitle(" ele2 pT (GeV)");
        }
        else{
          efficiencyRatio[kl][ij]->GetXaxis()->SetTitle(" muon2 #eta");
          efficiencyRatio[kl][ij]->GetYaxis()->SetTitle(" muon2 pT (GeV)");
        }
        efficiencyRatio[kl][ij]->GetYaxis()->SetTitleOffset(1.05);
        efficiencyRatio[kl][ij]->Draw("colz");
        if(isEleFinalState)
          c2[kl][ij]->Print(Form(("plots/EffxAcc/Eff_"+cutFlow.at(kl)+"_over_"+cutFlow.at(0)+"_NonResonant_eeK_massBin_%.2f-%.2f.png").c_str(), llMassBoundary.at(ij), llMassBoundary.at(ij+1)), "png");
        else
          c2[kl][ij]->Print(Form(("plots/EffxAcc/Eff_"+cutFlow.at(kl)+"_over_"+cutFlow.at(0)+"_NonResonant_mumuK_massBin_%.2f-%.2f.png").c_str(), llMassBoundary.at(ij), llMassBoundary.at(ij+1)), "png");
      }
      else{
        if(isEleFinalState){
          efficiencyRatio[kl][5]->GetXaxis()->SetTitle(" ele2 #eta");
          efficiencyRatio[kl][5]->GetYaxis()->SetTitle(" ele2 pT (GeV)");
        }
        else{
          efficiencyRatio[kl][5]->GetXaxis()->SetTitle(" muon2 #eta");
          efficiencyRatio[kl][5]->GetYaxis()->SetTitle(" muon2 pT (GeV)");
        }
        efficiencyRatio[kl][5]->GetYaxis()->SetTitleOffset(1.05);
        efficiencyRatio[kl][5]->Draw("colz");
        if(isEleFinalState)
	  c2[kl][ij]->Print(Form(("plots/EffxAcc/Eff_"+cutFlow.at(kl)+"_over_"+cutFlow.at(0)+"_JPsi_eeK_massBin_%.2f-%.2f.png").c_str(), llMassBoundary.at(ij), llMassBoundary.at(ij+1)), "png");
        else
	  c2[kl][ij]->Print(Form(("plots/EffxAcc/Eff_"+cutFlow.at(kl)+"_over_"+cutFlow.at(0)+"_JPsi_mumuK_massBin_%.2f-%.2f.png").c_str(), llMassBoundary.at(ij), llMassBoundary.at(ij+1)), "png");
      }
    }
    }



    //print table

    std::string outFileName = "plots/EffxAcc/Kmumu_finalState.txt";
    if(isEleFinalState) outFileName = "plots/EffxAcc/Kee_finalState.txt";

    std::ofstream outFileLong(outFileName.c_str(), std::ios::out);

  //upgrate to 2D plot for future or something better than printout
  for(int ij=0; ij<5; ++ij){
    // std::cout << " \n mass bin \t nonResonant = " << llGenMassCut.at(ij) << "\t        \t Resonant (J/Psi) \t   \t \t" << std::endl;
    outFileLong << " \n \t nonResonant ll mass-bin [" << llMassBoundary.at(ij) << "-"<< llMassBoundary.at(ij+1) << "]  \t\t\t Resonant (J/Psi) \t" << std::endl;
    outFileLong << " muonTag = \t " << nEv_muonTag[ij] << " \t\t\t" << nEv_muonTag[5]  << " \t " << std::endl;

    outFileLong << " genAcc = \t " << nEv_genAcc[ij] << " \t /muonTag = " << nEv_genAcc[ij]/nEv_muonTag[ij]
		<< " \t " << nEv_genAcc[5] << " \t  /muonTag = " << nEv_genAcc[5]/nEv_muonTag[5] 
		<< " \t double ratio = " << (nEv_genAcc[ij]/nEv_muonTag[ij])/(nEv_genAcc[5]/nEv_muonTag[5]) << "\n";

    outFileLong << " x genEff = \t " << nEv_genEff[ij] << " \t /muonTag = " << nEv_genEff[ij]/nEv_muonTag[ij]
		<< " \t " << nEv_genEff[5] << " \t /muonTag = " << nEv_genEff[5]/nEv_muonTag[5]
		<< " \t double ratio = " << (nEv_genEff[ij]/nEv_muonTag[ij])/(nEv_genEff[5]/nEv_muonTag[5]) << "\n";

    outFileLong << " x recoEff = \t " << nEv_recoEff[ij] << " \t /muonTag = " << nEv_recoEff[ij]/nEv_muonTag[ij]
		<< " \t " << nEv_recoEff[5] << " \t /muonTag = " << nEv_recoEff[5]/nEv_muonTag[5]
		<< " \t  double ratio = " << (nEv_recoEff[ij]/nEv_muonTag[ij])/(nEv_recoEff[5]/nEv_muonTag[5]) << "\n";

    outFileLong << " x chargeEff = \t " << nEv_chargeEff[ij] << " \t /muonTag = " << nEv_chargeEff[ij]/nEv_muonTag[ij]
		<< " \t " << nEv_chargeEff[5] << " \t /muonTag = " << nEv_chargeEff[5]/nEv_muonTag[5]
		<< " \t double ratio = " << (nEv_chargeEff[ij]/nEv_muonTag[ij])/(nEv_chargeEff[5]/nEv_muonTag[5]) << "\n";

    outFileLong << " x alphaEff = \t " << nEv_alphaEff[ij] << " \t /muonTag = " << nEv_alphaEff[ij]/nEv_muonTag[ij]
		<< " \t " << nEv_alphaEff[5] << " \t /muonTag = " << nEv_alphaEff[5]/nEv_muonTag[5]
		<< " \t double ratio = " << (nEv_alphaEff[ij]/nEv_muonTag[ij])/(nEv_alphaEff[5]/nEv_muonTag[5]) << "\n";

    outFileLong << " x vtxCLEff = \t " << nEv_vtxCLEff[ij] << " \t /muonTag = " << nEv_vtxCLEff[ij]/nEv_muonTag[ij]
		<< " \t " << nEv_vtxCLEff[5] << " \t /muonTag = " << nEv_vtxCLEff[5]/nEv_muonTag[5]
		<< " \t double ratio = " << (nEv_vtxCLEff[ij]/nEv_muonTag[ij])/(nEv_vtxCLEff[5]/nEv_muonTag[5]) << "\n";
    
    outFileLong << " x DCAEff = \t " << nEv_DCAEff[ij] << " \t /muonTag = " << nEv_DCAEff[ij]/nEv_muonTag[ij]
		<< " \t " << nEv_DCAEff[5]  << " \t /muonTag = " << nEv_DCAEff[5]/nEv_muonTag[5]
		<< " \t double ratio = " << (nEv_DCAEff[ij]/nEv_muonTag[ij])/(nEv_DCAEff[5]/nEv_muonTag[5]) << std::endl;
  }
  outFileLong.close();



  std::cout << " ciao " << std::endl;



}
