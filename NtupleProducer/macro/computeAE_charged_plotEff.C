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
#include "TLegend.h"

using namespace std;

void computeAE_charged_plotEff(TString nonResonantFile, TString ResonantFile, int isEleFinalState){

  gROOT->Reset();
  gROOT->Macro("setStyle.C");
  gROOT->SetBatch(kTRUE);

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

  vector<TH2F*> eff2D;
  eff2D.resize(cut_name.size());
  vector<vector<TH2F*>> efficiency2D;
  efficiency2D.resize(llMassCut.size()+1, eff2D);

  vector<TH1F*> th1;
  th1.resize(cut_name.size());
  vector<vector<TH1F*>> th1_pt;
  th1_pt.resize(llMassCut.size()+1, th1);

  TString gen_index = "GenPart_e2FromB_index";
  if(!isEleFinalState) gen_index = "GenPart_mu2FromB_index";

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


      if(i<llMassCut.size()){ //Non-resonant

	nEv_cuts[i][j] = t1->GetEntries(cuts[j]);

	TString name2D = "ratioEff_" + cut_name[j]+ Form("_2D_llmassBin_%.2f-%.2f", llMassBoundary[i], llMassBoundary[i+1]);
	TString name1D = cut_name[j]+ Form("_llmassBin_%.2f-%.2f", llMassBoundary[i], llMassBoundary[i+1]);
	efficiency2D[i][j] = new TH2F(name2D, "", 15, -3., 3., 20, 0., 10.);
	efficiency2D[i][j]->Sumw2();
	t1->Draw("GenPart_pt["+gen_index+"]:GenPart_eta["+gen_index+"] >>"+efficiency2D[i][j]->GetName(),cuts[j],"goff");
	th1_pt[i][j] = (TH1F*)(efficiency2D[i][j]->ProjectionY())->Clone(name1D);

      }

      else{ ; //Resonant
	nEv_cuts[i][j] = t2->GetEntries(cuts[j]);

	TString name2D = "ratioEff_" + cut_name[j]+ Form("_2D_JPsimassBin_%.2f-%.2f", llMassBoundary[JPsi_bin], llMassBoundary[JPsi_bin+1]);
	TString name1D = cut_name[j]+ Form("_JPsimassBin_%.2f-%.2f", llMassBoundary[JPsi_bin], llMassBoundary[JPsi_bin+1]);

	efficiency2D[i][j] = new TH2F(name2D, "", 15, -3., 3., 20, 0., 10.);
	efficiency2D[i][j]->Sumw2();	
	t2->Draw("GenPart_pt["+gen_index+"]:GenPart_eta["+gen_index+"] >>"+efficiency2D[i][j]->GetName(),cuts[j],"goff");
	th1_pt[i][j] = (TH1F*)(efficiency2D[i][j]->ProjectionY())->Clone(name1D);

      }

      if(j>0) efficiency2D[i][j]->Divide(efficiency2D[i][0]);

    }
  }


  ///plotting
  TLegend *legM = new TLegend(0.72,0.58,0.93,0.93,NULL,"brNDC");
  legM->SetTextFont(42);
  legM->SetFillColor(kWhite);
  legM->SetLineColor(kWhite);
  legM->SetShadowColor(kWhite);
  //legM->SetFillStyle(0);
  legM->SetTextSize(0.03);


  vector<int> iColors = {kYellow+3, kViolet, kRed, kGreen+1, kRed+1, kBlue, kCyan+2, kMagenta};

  for(unsigned j=0; j<cut_name.size(); j++){
    for(unsigned i=0; i<llMassCut.size()+1; i++){
      th1_pt[i][j]->SetLineColor(iColors[j]);
      th1_pt[i][j]->SetLineWidth(2);
    }
    legM->AddEntry(th1_pt[0][j], cut_name[j], "l");
  }


  TCanvas* c;
  for(unsigned i=0; i<llMassCut.size()+1; i++){
    if(i==JPsi_bin) continue;

    if(i<llMassCut.size()) legM->SetHeader(Form("mass bin [%.2f-%.2f]", llMassBoundary[i], llMassBoundary[i+1]), "c" );
    else legM->SetHeader(Form("mass bin [%.2f-%.2f]", llMassBoundary[JPsi_bin], llMassBoundary[JPsi_bin+1]), "c" );
    c = new TCanvas();
    c->cd();
    gPad->SetLogy();

    if(isEleFinalState){
      th1_pt[i][0]->GetXaxis()->SetTitle(" ele2 pT (GeV)");
      th1_pt[i][0]->GetYaxis()->SetTitle(" nEvents (Kee)");
    }
    else{
      th1_pt[i][0]->GetXaxis()->SetTitle(" muon2 pT (GeV)");
      th1_pt[i][0]->GetYaxis()->SetTitle(" nEvents (Kmumu)");
    }

    th1_pt[i][0]->GetXaxis()->SetTitleOffset(1.05);
    th1_pt[i][0]->GetYaxis()->SetTitleOffset(1.05);
    th1_pt[i][0]->Draw("histgoff");

    for(unsigned int j=1;j<cut_name.size();j++) th1_pt[i][j]->Draw("hist, same goff");
    legM->Draw("hist, same goff");

    TString filename = "plots/EffxAcc/vsPt/";
    if(i<llMassCut.size()) filename+="NonResonant_";
    else filename+="JPsi_";
    if(isEleFinalState) filename+="eeK_";
    else filename+="mumuK_";
    if(i<llMassCut.size()) filename+=Form("massBin_%.2f-%.2f.png", llMassBoundary[i], llMassBoundary[i+1]);
    else filename+=Form("massBin_%.2f-%.2f.png", llMassBoundary[JPsi_bin], llMassBoundary[JPsi_bin+1]);
    c->Print(filename,"png");
  }



  //print 2D eff

  for(unsigned i=0; i<llMassCut.size()+1; i++){
    for(unsigned int j=0;j<cut_name.size();j++){

      c = new TCanvas();
      c->cd();

      if(isEleFinalState){
	efficiency2D[i][j]->GetXaxis()->SetTitle(" ele2 #eta");
	efficiency2D[i][j]->GetYaxis()->SetTitle(" ele2 pT (GeV)");
      }
      else{
	efficiency2D[i][j]->GetXaxis()->SetTitle(" muon2 #eta");
	efficiency2D[i][j]->GetYaxis()->SetTitle(" muon2 pT (GeV)");
      }

      efficiency2D[i][j]->GetYaxis()->SetTitleOffset(1.05);
      efficiency2D[i][j]->Draw("colz goff");

      TString filename = "plots/EffxAcc/";
      if(j==0) filename+=cut_name[0]+"_";
      else filename+="Eff_"+cut_name[j]+"_over_"+cut_name[0]+"_";
      if(i<llMassCut.size()) filename+="NonResonant_";
      else filename+="JPsi_";
      if(isEleFinalState) filename+="eeK_";
      else filename+="mumuK_";
      if(i<llMassCut.size()) filename+=Form("massBin_%.2f-%.2f.png", llMassBoundary[i], llMassBoundary[i+1]);
      else filename+=Form("massBin_%.2f-%.2f.png", llMassBoundary[JPsi_bin], llMassBoundary[JPsi_bin+1]);
      c->Print(filename,"png");

    }
  }


  //print table

  TString outFileName = "plots/EffxAcc/Kmumu_finalState.txt";
  if(isEleFinalState) outFileName = "plots/EffxAcc/Kee_finalState.txt";

  std::ofstream outFileLong(outFileName, std::ios::out);

  for(unsigned i=0; i<llMassCut.size(); i++){
    outFileLong << " \n   \t nonResonant ll mass-bin [" << llMassBoundary.at(i) << "-"<< llMassBoundary.at(i+1) << "] \t       \t Resonant (J/Psi) \t   \t \t" << std::endl;

    for(unsigned j=0; j<cut_name.size(); j++){
      if(j==0){
	outFileLong << " " << cut_name[j] << " = \t " << nEv_cuts[i][j] << " \t              \t \t        \t " << nEv_cuts[llMassCut.size()][j]  << " \t " << std::endl;
      }
      else{
	double eff = nEv_cuts[i][j]/nEv_cuts[i][0];
	double eff_JPsi = nEv_cuts[llMassCut.size()][j]/nEv_cuts[llMassCut.size()][0];
	outFileLong << " " << cut_name[j] << " = \t " << nEv_cuts[i][j] << " \t /muonTag = " << eff
		    << "               \t " << nEv_cuts[llMassCut.size()][j] << " \t  /muonTag = " << eff_JPsi
		    << " \t double ratio = \t" << eff/eff_JPsi << "\n";
      }

    }

  }

  outFileLong.close();

  std::cout << " ciao " << std::endl;



}
