//g++ -Wall -o analyzeCharged_fastDATA `root-config --cflags --glibs` -lRooFitCore analyzeCharged_fastDATA.cpp
//
//to run on ele all dataset 
// ./analyzeCharged_fastDATA 1 -1
//to run on muon runA
// ./analyzeCharged_fastDATA 0  runA
// options are: isEleFinalState (1, 0)   dataset (-1, runA, runB)  nMaxEvents (-1, N)

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>

#include "TApplication.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TCut.h"
#include "TError.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TKey.h"
#include "TLegend.h"
#include "TMath.h"
#include "TPad.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

#include "RooAddPdf.h"
#include "RooCBShape.h"
#include "RooChebychev.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooGaussian.h"
#include "RooMsgService.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"

using namespace RooFit;

const int kBToKllMax = 100000;


int main(int argc, char *argv[]){

  if(argc < 2) {
    std::cout << " Missing arguments " << std::endl;
    return -1;
  }
  int isEleFinalState = atoi(argv[1]);
  std::string dataset = "-1";
  int nMaxEvents = -1;
  if(argc > 2) dataset = argv[2];
  if(argc > 3) nMaxEvents = atoi(argv[3]);


  gROOT->Reset();
  gROOT->Macro("./setStyle.C");
  gSystem->Load("libRooFit") ;
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  TChain* t1 = new TChain("Events");
  //new prod Kee
  if(isEleFinalState){
    if(dataset == "runA" || dataset == "-1"){
      t1->Add("/vols/cms/tstreble/BPH/BToKee_ntuple/BPHParking1_2018A_18_08_14/*root");
      t1->Add("/vols/cms/tstreble/BPH/BToKee_ntuple/BPHParking2_2018A_18_08_14/*root");
      t1->Add("/vols/cms/tstreble/BPH/BToKee_ntuple/BPHParking3_2018A_18_08_14/*root");
      t1->Add("/vols/cms/tstreble/BPH/BToKee_ntuple/BPHParking4_2018A_18_08_14/*root");
      t1->Add("/vols/cms/tstreble/BPH/BToKee_ntuple/BPHParking5_2018A_18_08_14/*root");
      t1->Add("/vols/cms/tstreble/BPH/BToKee_ntuple/BPHParking6_2018A_18_08_14/*root");
    }
    if(dataset == "runB" || dataset == "-1"){
      t1->Add("/vols/cms/tstreble/BPH/BToKee_ntuple/BPHParking1_2018B_18_08_14/*root");
      t1->Add("/vols/cms/tstreble/BPH/BToKee_ntuple/BPHParking2_2018B_18_08_14/*root");
      t1->Add("/vols/cms/tstreble/BPH/BToKee_ntuple/BPHParking3_2018B_18_08_14/*root");
      t1->Add("/vols/cms/tstreble/BPH/BToKee_ntuple/BPHParking4_2018B_18_08_14/*root");
      t1->Add("/vols/cms/tstreble/BPH/BToKee_ntuple/BPHParking5_2018B_18_08_14/*root");
    }
  }
  else{
    if(dataset == "runA" || dataset == "-1"){
      t1->Add("/vols/cms/tstreble/BPH/BToKmumu_ntuple/BPHParking1_2018A_18_08_14/*root");
      t1->Add("/vols/cms/tstreble/BPH/BToKmumu_ntuple/BPHParking2_2018A_18_08_14/*root");
      t1->Add("/vols/cms/tstreble/BPH/BToKmumu_ntuple/BPHParking3_2018A_18_08_14/*root");
      t1->Add("/vols/cms/tstreble/BPH/BToKmumu_ntuple/BPHParking4_2018A_18_08_14/*root");
      t1->Add("/vols/cms/tstreble/BPH/BToKmumu_ntuple/BPHParking5_2018A_18_08_14/*root");
      t1->Add("/vols/cms/tstreble/BPH/BToKmumu_ntuple/BPHParking6_2018A_18_08_14/*root");
    }
    if(dataset == "runB" || dataset == "-1"){
      t1->Add("/vols/cms/tstreble/BPH/BToKmumu_ntuple/BPHParking1_2018B_18_08_14/*root");
      t1->Add("/vols/cms/tstreble/BPH/BToKmumu_ntuple/BPHParking2_2018B_18_08_14/*root");
      t1->Add("/vols/cms/tstreble/BPH/BToKmumu_ntuple/BPHParking3_2018B_18_08_14/*root");
      t1->Add("/vols/cms/tstreble/BPH/BToKmumu_ntuple/BPHParking4_2018B_18_08_14/*root");
      t1->Add("/vols/cms/tstreble/BPH/BToKmumu_ntuple/BPHParking5_2018B_18_08_14/*root");
    }
  }

  int nEvts = t1->GetEntries();
  std::cout << " #initial n. events: " << nEvts << std::endl;
  
  int MuonTag_index = -1;
  int BToKll_sel_index = -1;
  float BToKll_lep1_charge[kBToKllMax];
  float BToKll_lep2_charge[kBToKllMax];
  float BToKll_cosAlpha[kBToKllMax];
  float BToKll_CL_vtx[kBToKllMax];
  float BToKll_kaon_DCASig[kBToKllMax];
  float BToKll_Lxy[kBToKllMax];
  float BToKll_llKFit_ll_mass[kBToKllMax];
  float BToKll_mass[kBToKllMax];
  float BToKll_pt[kBToKllMax];
  float BToKll_kaon_pt[kBToKllMax];
  int passed_2trk = -1;
  
  t1->SetBranchStatus("*", 0);
  if(isEleFinalState){
    t1->SetBranchStatus("Muon_sel_index", 1);            t1->SetBranchAddress("Muon_sel_index", &MuonTag_index);
    t1->SetBranchStatus("BToKee_sel_index", 1);          t1->SetBranchAddress("BToKee_sel_index", &BToKll_sel_index);
    t1->SetBranchStatus("BToKee_ele1_charge", 1);        t1->SetBranchAddress("BToKee_ele1_charge", &BToKll_lep1_charge);
    t1->SetBranchStatus("BToKee_ele2_charge", 1);        t1->SetBranchAddress("BToKee_ele2_charge", &BToKll_lep2_charge);
    t1->SetBranchStatus("BToKee_cosAlpha", 1);           t1->SetBranchAddress("BToKee_cosAlpha", &BToKll_cosAlpha);
    t1->SetBranchStatus("BToKee_CL_vtx", 1);             t1->SetBranchAddress("BToKee_CL_vtx", &BToKll_CL_vtx);
    t1->SetBranchStatus("BToKee_kaon_DCASig", 1);        t1->SetBranchAddress("BToKee_kaon_DCASig", &BToKll_kaon_DCASig);
    t1->SetBranchStatus("BToKee_Lxy", 1);                t1->SetBranchAddress("BToKee_Lxy", &BToKll_Lxy);
    t1->SetBranchStatus("BToKee_eeKFit_ee_mass", 1);     t1->SetBranchAddress("BToKee_eeKFit_ee_mass", &BToKll_llKFit_ll_mass);
    t1->SetBranchStatus("BToKee_mass", 1);               t1->SetBranchAddress("BToKee_mass", &BToKll_mass);
    t1->SetBranchStatus("BToKee_pt", 1);                 t1->SetBranchAddress("BToKee_pt", &BToKll_pt);
    t1->SetBranchStatus("BToKee_kaon_pt", 1);            t1->SetBranchAddress("BToKee_kaon_pt", &BToKll_kaon_pt);
    t1->SetBranchStatus("BToKee_eeRefit", 1);            t1->SetBranchAddress("BToKee_eeRefit", &passed_2trk);
  }
  else{
    t1->SetBranchStatus("Muon_probe_index", 1);            t1->SetBranchAddress("Muon_probe_index", &MuonTag_index);
    t1->SetBranchStatus("BToKmumu_sel_index", 1);          t1->SetBranchAddress("BToKmumu_sel_index", &BToKll_sel_index);
    t1->SetBranchStatus("BToKmumu_mu1_charge", 1);         t1->SetBranchAddress("BToKmumu_mu1_charge", &BToKll_lep1_charge);
    t1->SetBranchStatus("BToKmumu_mu2_charge", 1);         t1->SetBranchAddress("BToKmumu_mu2_charge", &BToKll_lep2_charge);
    t1->SetBranchStatus("BToKmumu_cosAlpha", 1);           t1->SetBranchAddress("BToKmumu_cosAlpha", &BToKll_cosAlpha);
    t1->SetBranchStatus("BToKmumu_CL_vtx", 1);             t1->SetBranchAddress("BToKmumu_CL_vtx", &BToKll_CL_vtx);
    t1->SetBranchStatus("BToKmumu_kaon_DCASig", 1);        t1->SetBranchAddress("BToKmumu_kaon_DCASig", &BToKll_kaon_DCASig);
    t1->SetBranchStatus("BToKmumu_Lxy", 1);                t1->SetBranchAddress("BToKmumu_Lxy", &BToKll_Lxy);
    t1->SetBranchStatus("BToKmumu_mumuKFit_mumu_mass", 1); t1->SetBranchAddress("BToKmumu_mumuKFit_mumu_mass", &BToKll_llKFit_ll_mass);
    t1->SetBranchStatus("BToKmumu_mass", 1);               t1->SetBranchAddress("BToKmumu_mass", &BToKll_mass);
    t1->SetBranchStatus("BToKmumu_pt", 1);                 t1->SetBranchAddress("BToKmumu_pt", &BToKll_pt);
    t1->SetBranchStatus("BToKmumu_kaon_pt", 1);            t1->SetBranchAddress("BToKmumu_kaon_pt", &BToKll_kaon_pt);
    t1->SetBranchStatus("BToKmumu_mumuRefit", 1);          t1->SetBranchAddress("BToKmumu_mumuRefit", &passed_2trk);
  }


  //selections 
  /*
  std::string cut_muonTag = "Muon_sel_index != -1";
  std::string cut_Recocandidate = "BToKee_sel_index != -1";
  std::string cut_chargeEff = "BToKee_ele1_charge[BToKee_sel_index]*BToKee_ele2_charge[BToKee_sel_index] < 0.";
  std::string cut_alphaEff = "BToKee_cosAlpha[BToKee_sel_index] > 0.999";
  std::string cut_vtxCLEff = "BToKee_CL_vtx[BToKee_sel_index] > 0.1";
  std::string cut_LxyEff  = "BToKee_Lxy[BToKee_sel_index] > 6";

  if(!isEleFinalState){
    cut_muonTag = "Muon_probe_index != -1";
    cut_Recocandidate = "BToKmumu_sel_index != -1";
    cut_chargeEff = "BToKmumu_mu1_charge[BToKmumu_sel_index]*BToKmumu_mu2_charge[BToKmumu_sel_index] < 0.";
    cut_alphaEff = "BToKmumu_cosAlpha[BToKmumu_sel_index] > 0.999";
    cut_vtxCLEff = "BToKmumu_CL_vtx[BToKmumu_sel_index] > 0.1";
    cut_LxyEff  = "BToKmumu_Lxy[BToKmumu_sel_index] > 6";
  }
  */


  std::vector<float> llMassBoundary;
  llMassBoundary.push_back(0.);
  llMassBoundary.push_back(1.);
  llMassBoundary.push_back(2.5);
  llMassBoundary.push_back(2.9);
  llMassBoundary.push_back(3.2);
  llMassBoundary.push_back(3.58);

  /*
  std::vector<std::string> llMassCut;
  for(int ij=0; ij<5; ++ij){
    std::string cut = Form("BToKee_eeKFit_ee_mass[BToKee_sel_index] > %.2f && BToKee_eeKFit_ee_mass[BToKee_sel_index] < %.2f",
                           llMassBoundary.at(ij), llMassBoundary.at(ij+1));
    if(!isEleFinalState){
      cut = Form("BToKmumu_mumuKFit_mumu_mass[BToKmumu_sel_index] > %.2f && BToKmumu_mumuKFit_mumu_mass[BToKmumu_sel_index] < %.2f",
                 llMassBoundary.at(ij), llMassBoundary.at(ij+1));
    }
    llMassCut.push_back(cut);
  }
  */

  ///histos: 1 per bin plus inclusive
  TH1F* hAlpha[6];
  TH1F* hCLVtx[6];
  TH1F* hDCASig[6];
  TH1F* hLxy[6];
  TH1F* hllMass[6];
  TH2F* hllMass_vs_Bmass[6];
  TH1F* hBmass[6];

  for(int ij=0; ij<6; ++ij){
    hAlpha[ij] = new TH1F(Form("hAlpha_bin%d", ij), "", 500, 0, 1.1);
    hAlpha[ij]->Sumw2();
    hAlpha[ij]->SetLineColor(kRed);
    hAlpha[ij]->SetLineWidth(2);

    hCLVtx[ij] = new TH1F(Form("hCLVtx_%d", ij), "", 100, 0., 1.);
    hCLVtx[ij]->Sumw2();
    hCLVtx[ij]->SetLineColor(kRed);
    hCLVtx[ij]->SetLineWidth(2);
 
    hDCASig[ij] = new TH1F(Form("hDCASig_%d", ij), "", 100, -50., 50.);
    hDCASig[ij]->Sumw2();
    hDCASig[ij]->SetLineColor(kRed);
    hDCASig[ij]->SetLineWidth(2);

    hLxy[ij] = new TH1F(Form("hLxy_%d", ij), "", 100, 0., 100.);
    hLxy[ij]->Sumw2();
    hLxy[ij]->SetLineColor(kRed);
    hLxy[ij]->SetLineWidth(2);

    hllMass[ij] = new TH1F(Form("hllMass_%d", ij), "", 750, 0., 15.);
    hllMass[ij]->Sumw2();
    hllMass[ij]->SetLineColor(kRed);
    hllMass[ij]->SetLineWidth(2);

    hllMass_vs_Bmass[ij] = new TH2F(Form("hllMass_vs_Bmass_%d", ij), "", 500, 0., 15., 500, 0., 15.);
    hllMass_vs_Bmass[ij]->Sumw2();
    hllMass_vs_Bmass[ij]->SetMarkerColor(kRed);
    hllMass_vs_Bmass[ij]->SetMarkerStyle(20);

    hBmass[ij] = new TH1F(Form("Bmass_%d", ij), "", 750, 0., 15.); // 75, 4.5, 6.);
    hBmass[ij]->Sumw2();
    hBmass[ij]->SetLineColor(kRed);
    hBmass[ij]->SetLineWidth(2);
  }


  float nEv_muonTag[5] = {0.};
  float nEv_recoCand[5] = {0.};
  float nEv_chargeEff[5] = {0.};
  float nEv_alphaEff[5] = {0.};
  float nEv_vtxCLEff[5] = {0.};
  float nEv_LxyEff[5] = {0.};

  if(nMaxEvents == -1) nMaxEvents = nEvts;
  for(int iEvt = 0; iEvt<nMaxEvents; ++iEvt){

    if(iEvt%500000 == 0) std::cout << " >>> processing event " << iEvt << " " << 1.*iEvt/nEvts*100. << std::endl;

    t1->GetEntry(iEvt);
    if(MuonTag_index == -1) continue;
    ++nEv_muonTag[0];

    if(BToKll_sel_index == -1) continue; 
    ++nEv_recoCand[0];

    if(BToKll_lep1_charge[BToKll_sel_index]*BToKll_lep2_charge[BToKll_sel_index] > 0.) continue;

    /*
    if(BToKll_kaon_pt[BToKll_sel_index] < 1.5 || BToKll_pt[BToKll_sel_index] < 10.) continue;
    */

    float llInvMass = BToKll_llKFit_ll_mass[BToKll_sel_index];
    int massBin = -1;
    for(unsigned int kl=0; kl<llMassBoundary.size()-1; ++kl){
      if(llInvMass >= llMassBoundary[kl] && llInvMass < llMassBoundary[kl+1]){
	massBin = kl;
	break;
      }
    }

    if(massBin == -1) continue;
    ++nEv_chargeEff[massBin];

    hAlpha[massBin]->Fill(BToKll_cosAlpha[BToKll_sel_index]);
    hAlpha[5]->Fill(BToKll_cosAlpha[BToKll_sel_index]);
    //if(BToKll_cosAlpha[BToKll_sel_index] < 0.999) continue;
    if(BToKll_cosAlpha[BToKll_sel_index] < 0.99) continue;
    ++nEv_alphaEff[massBin];

    hCLVtx[massBin]->Fill(BToKll_CL_vtx[BToKll_sel_index]);
    hCLVtx[5]->Fill(BToKll_CL_vtx[BToKll_sel_index]);
    if(BToKll_CL_vtx[BToKll_sel_index] < 0.1) continue;
    ++nEv_vtxCLEff[massBin];    

    hDCASig[massBin]->Fill(BToKll_kaon_DCASig[BToKll_sel_index]);
    hLxy[massBin]->Fill(BToKll_Lxy[BToKll_sel_index]);
    hDCASig[5]->Fill(BToKll_kaon_DCASig[BToKll_sel_index]);
    hLxy[5]->Fill(BToKll_Lxy[BToKll_sel_index]);
    if(BToKll_Lxy[BToKll_sel_index] < 6.) continue;
    ++nEv_LxyEff[massBin];


    hllMass[massBin]->Fill(llInvMass);
    hllMass_vs_Bmass[massBin]->Fill(BToKll_mass[BToKll_sel_index], llInvMass);
    hBmass[massBin]->Fill(BToKll_mass[BToKll_sel_index]);

    hllMass[5]->Fill(llInvMass);
    hllMass_vs_Bmass[5]->Fill(BToKll_mass[BToKll_sel_index], llInvMass);
    hBmass[5]->Fill(BToKll_mass[BToKll_sel_index]);

  }//loop over events


  std::string outName = "outMassHistos_DATA_Kee.root";
  if(!isEleFinalState) outName = "outMassHistos_DATA_Kmumu.root";
  TFile outMassHistos(outName.c_str(), "recreate");
  outMassHistos.cd();

  std::cout << " ***** summary ***** "<< std::endl;
  for(int ij=0; ij<6; ++ij){
    hAlpha[ij]->Write(hAlpha[ij]->GetName());
    hCLVtx[ij]->Write(hCLVtx[ij]->GetName());
    hDCASig[ij]->Write(hDCASig[ij]->GetName());
    hLxy[ij]->Write(hLxy[ij]->GetName());
    hllMass[ij]->Write(hllMass[ij]->GetName());
    hllMass_vs_Bmass[ij]->Write(hllMass_vs_Bmass[ij]->GetName());
    hBmass[ij]->Write(hBmass[ij]->GetName());

    if(ij > 4) continue;
    std::cout << "\n massBin: " << llMassBoundary[ij] << " - " << llMassBoundary[ij+1]
              << " \n \t recoEvts = " << nEv_chargeEff[ij] << " \t alphaCut = " << nEv_alphaEff[ij]
	      << " \t vtxCLCut = " << nEv_vtxCLEff[ij] << " \t LxyCut = " << nEv_LxyEff[ij] << std::endl;

  }
  outMassHistos.Close();

  TLegend *legM = new TLegend(0.70,0.70,0.98,0.95,NULL,"brNDC");
  legM->SetTextFont(42);
  legM->SetFillColor(kWhite);
  legM->SetLineColor(kWhite);
  legM->SetShadowColor(kWhite);
  legM->SetFillStyle(0);
  legM->SetTextSize(0.05);
  legM->AddEntry(hAlpha[3], "DATA", "l");
  

  TCanvas* c1 = new TCanvas();
  c1->cd();
  gPad->SetLogy();
  hAlpha[3]->GetXaxis()->SetTitle("cosAlpha");
  hAlpha[3]->Draw("hist");
  legM->Draw("same");
  if(isEleFinalState){
    c1->Print("plots/Kee_cosAlpha_DA.png", "png");
    c1->Print("plots/Kee_cosAlpha_DA.pdf", "pdf");
    c1->Print("plots/Kee_cosAlpha_DA.root", "root");
  }
  else{
    c1->Print("plots/Kmumu_cosAlpha_DA.png", "png");
    c1->Print("plots/Kmumu_cosAlpha_DA.pdf", "pdf");
    c1->Print("plots/Kmumu_cosAlpha_DA.root", "root");
  }

  TCanvas* c1b = new TCanvas();
  c1b->cd();
  gPad->SetLogy();
  hCLVtx[3]->GetXaxis()->SetTitle("CL_vtx");
  hCLVtx[3]->Draw("hist");
  legM->Draw("same");
  if(isEleFinalState){
    c1b->Print("plots/Kee_CL_vtx_DA.png", "png");
    c1b->Print("plots/Kee_CL_vtx_DA.pdf", "pdf");
    c1b->Print("plots/Kee_CL_vtx_DA.root", "root");
  }
  else{
    c1b->Print("plots/Kmumu_CL_vtx_DA.png", "png");
    c1b->Print("plots/Kmumu_CL_vtx_DA.pdf", "pdf");
    c1b->Print("plots/Kmumu_CL_vtx_DA.root", "root");
  }

  TCanvas* c1c = new TCanvas();
  hDCASig[3]->GetXaxis()->SetTitle("kaon DCA SIP");
  hDCASig[3]->Draw("hist");
  legM->Draw("same");
  if(isEleFinalState){
    c1c->Print("plots/Kee_kaonDCA_DA.png", "png");
    c1c->Print("plots/Kee_kaonDCA_DA.pdf", "pdf");
    c1c->Print("plots/Kee_kaonDCA_DA.root", "root");
  }
  else{
    c1c->Print("plots/Kmumu_kaonDCA_DA.png", "png");
    c1c->Print("plots/Kmumu_kaonDCA_DA.pdf", "pdf");
    c1c->Print("plots/Kmumu_kaonDCA_DA.root", "root");
  }
  
  TCanvas* c1d = new TCanvas();
  gPad->SetLogy();
  hLxy[3]->GetXaxis()->SetTitle("Lxy IP(cm)");
  hLxy[3]->Draw("hist");
  legM->Draw("same");
  if(isEleFinalState){
    c1d->Print("plots/Kee_Lxy_DA.png", "png");
    c1d->Print("plots/Kee_Lxy_DA.pdf", "pdf");
    c1d->Print("plots/Kee_Lxy_DA.root", "root");
  }
  else{
    c1d->Print("plots/Kmumu_Lxy_DA.png", "png");
    c1d->Print("plots/Kmumu_Lxy_DA.pdf", "pdf");
    c1d->Print("plots/Kmumu_Lxy_DA.root", "root");
  }

  TCanvas* c1e = new TCanvas();
  gPad->SetLogy();
  if(isEleFinalState)
    hllMass[3]->GetXaxis()->SetTitle("Mee mass");
  else hllMass[3]->GetXaxis()->SetTitle("Mmumu mass");
  hllMass[3]->Draw("hist");
  legM->Draw("same");
  if(isEleFinalState){
    c1e->Print("plots/Kee_Mee_DA.png", "png");
    c1e->Print("plots/Kee_Mee_DA.pdf", "pdf");
    c1e->Print("plots/Kee_Mee_DA.root", "root");
  }
  else{
    c1e->Print("plots/Kmumu_Mee_DA.png", "png");
    c1e->Print("plots/Kmumu_Mee_DA.pdf", "pdf");
    c1e->Print("plots/Kmumu_Mee_DA.root", "root");
  }

  TCanvas* c1f = new TCanvas();
  if(isEleFinalState){
    hllMass_vs_Bmass[3]->GetXaxis()->SetTitle("Kee mass");
    hllMass_vs_Bmass[3]->GetYaxis()->SetTitle("ee mass");
  }
  else{
    hllMass_vs_Bmass[3]->GetXaxis()->SetTitle("Kmumu mass");
    hllMass_vs_Bmass[3]->GetYaxis()->SetTitle("mumu mass");
  }
  hllMass_vs_Bmass[3]->Draw("");
  legM->Draw("same");
  if(isEleFinalState){
    c1f->Print("plots/Kee_Mee_vsBmass_DA.png", "png");
    c1f->Print("plots/Kee_Mee_vsBmass_DA.pdf", "pdf");
    c1f->Print("plots/Kee_Mee_vsBmass_DA.root", "root");
  }
  else{
    c1f->Print("plots/Kmumu_Mmumu_vsBmass_DA.png", "png");
    c1f->Print("plots/Kmumu_Mmumu_vsBmass_DA.pdf", "pdf");
    c1f->Print("plots/Kmumu_Mmumu_vsBmass_DA.root", "root");
  }

  TCanvas* c1g = new TCanvas();
  gPad->SetLogy();
  if(isEleFinalState)  hBmass[3]->GetXaxis()->SetTitle("B(Kee) mass");
  else hBmass[3]->GetXaxis()->SetTitle("B(Kmumu) mass");
  hBmass[3]->Draw("hist");
  legM->Draw("same");
  if(isEleFinalState){
    c1g->Print("plots/Kee_mass_DA.png", "png");
    c1g->Print("plots/Kee_mass_DA.pdf", "pdf");
    c1g->Print("plots/Kee_mass_DA.root", "root");
  }
  else{
    c1g->Print("plots/Kmumu_mass_DA.png", "png");
    c1g->Print("plots/Kmumu_mass_DA.pdf", "pdf");
    c1g->Print("plots/Kmumu_mass_DA.root", "root");
  }



  //comment out for the moment
  ///////////////////////////////////////////////////////////
  //now fitting

  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

  float nEv_postFit[5] = {0.};
  float nEvError_postFit[5] = {0.};

  RooWorkspace w("w");
  w.factory("x[0, 10]");
  //w.factory("x[4.5, 6.]");                                                                                                                                           

  w.factory("nbackground[10000, 0, 10000]");
  w.factory("nsignal[100, 0.0, 10000.0]");

  for(int ij=0; ij<5; ++ij){
    //w.factory("Gaussian::smodel(x,mu[5.3,4.5,6],sigma[0.05,0,0.2])");                                                                                                
    //w.factory("RooCBShape::smodel(x,m[5.3,4.5,6],s[0.1,0.,1.],a[1.2,0.,3.],n[1,0.1,6.])");
    w.factory("RooCBShape::smodel(x,m[5.3,0.,10.],s[0.1,0.,1.],a[1.2,0.,3.],n[1,0.1,6.])");
    //w.factory("RooCBShape::CBall(x[0,15], mean[11000,13000], sigma[5000,200000], alpha[0,10000],n[0,100000])");                                                      
    RooAbsPdf * smodel = w.pdf("smodel");

    w.factory("Exponential::bmodel(x,tau[-2,-3,0])");
    RooAbsPdf * bmodel = w.pdf("bmodel");

    w.factory("SUM::model(nbackground*bmodel, nsignal*smodel)");
    RooAbsPdf * model = w.pdf("model");

    RooDataHist hBMass("hBMass", "hBMass", *w.var("x"), Import(*(hBmass[ij])));

    RooFitResult * r = model->fitTo(hBMass, Minimizer("Minuit2"),Save(true));

    RooPlot * plot = w.var("x")->frame();
    if(isEleFinalState){
      plot->SetXTitle("Kee mass (GeV)");
    }
    else{
      plot->SetXTitle("K#mu#mu mass (GeV)");
    }
    plot->SetTitle("");
    plot->SetAxisRange(4.,6);
    hBMass.plotOn(plot);
    model->plotOn(plot);
    model->plotOn(plot, Components("bmodel"),LineStyle(kDashed));
    model->plotOn(plot, Components("smodel"),LineColor(kRed));

    TCanvas * cc = new TCanvas();
    cc->SetLogy(0);
    plot->Draw();
    if(isEleFinalState) cc->Print(Form("plots/Bmass_DATA/Kee_%s.png",hBmass[ij]->GetName()), "png");
    else cc->Print(Form("plots/Bmass_DATA/Kmumu_%s.png",hBmass[ij]->GetName()), "png");

    RooRealVar* parS = (RooRealVar*) r->floatParsFinal().find("nsignal");
    RooRealVar* parB = (RooRealVar*) r->floatParsFinal().find("nbackground");
    nEv_postFit[ij] = parS->getValV();
    nEvError_postFit[ij] = parS->getError();

    std::cout << " selection signal events = \t " << parS->getValV() << " error = " << parS->getError()
              << " bkg events = " << parB->getValV() << " error = " << parB->getError() << std::endl;
  }

  std::cout << " ***** summary ***** "<< std::endl;
  for(int ij=0; ij<5; ++ij){

    std::cout << "\n category = " << hBmass[ij]->GetName()
              << " \n \t recoEvts = " << nEv_LxyEff[ij] << " postFit " << nEv_postFit[ij] << "+/-" << nEvError_postFit[ij] << std::endl;
  }

}
