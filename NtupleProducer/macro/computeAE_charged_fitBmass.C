//compute double ratio nnReso/JPsi  event counts following all levels of selections
//muonTag + 
//gen Acceptance +
//gen Efficiency (reco final state matched to gen level) +
//selections for the analysis in data  

//should be ok for electron final state 
//example to run
//electron final state
//root -l computeAE_charged_fitBmass.C'("/vols/cms/amartell/BParking/ntuPROD/newNANO_20Aug/ntu_BToKee_18_08_14.root" , "/vols/cms/amartell/BParking/ntuPROD/newNANO_20Aug/ntu_BToKJPsiee_18_08_14.root", 1, true)' 
//muon final state
//root -l computeAE_charged_fitBmass.C'("/vols/cms/amartell/BParking/ntuPROD/newNANO_20Aug/ntu_BToKmumu_18_09_3.root" , "/vols/cms/amartell/BParking/ntuPROD/newNANO_20Aug/ntu_BToKJPsimumu_18_09_3.root", 0, true)' 


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
#include "TGraphErrors.h"
#include "TPad.h"
#include "TMultiGraph.h"
#include "TStyle.h"
#include "TChain.h"

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooCBShape.h"
#include "RooAddPdf.h"
#include "RooWorkspace.h"
#include "RooFitResult.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"

using namespace RooFit;


void computeAE_charged_fitBmass(std::string nonResonantFile, std::string ResonantFile, int isEleFinalState, bool foldGenMassBin=true){

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
  std::string cut_Recocandidate = "BToKee_sel_index != -1";
  std::string cut_chargeEff = "BToKee_ele1_charge[BToKee_sel_index]*BToKee_ele2_charge[BToKee_sel_index] < 0.";
  std::string cut_alphaEff = "BToKee_cosAlpha[BToKee_sel_index] > 0.99";
  std::string cut_vtxCLEff = "BToKee_CL_vtx[BToKee_sel_index] > 0.1";  
  std::string cut_LxyEff  = "BToKee_Lxy[BToKee_sel_index] > 6";

  if(!isEleFinalState){
    cut_muonTag = "Muon_probe_index != -1";
    cut_Recocandidate = "BToKmumu_sel_index != -1";
    cut_chargeEff = "BToKmumu_mu1_charge[BToKmumu_sel_index]*BToKmumu_mu2_charge[BToKmumu_sel_index] < 0.";
    cut_alphaEff = "BToKmumu_cosAlpha[BToKmumu_sel_index] > 0.99";
    cut_vtxCLEff = "BToKmumu_CL_vtx[BToKmumu_sel_index] > 0.1";
    cut_LxyEff  = "BToKmumu_Lxy[BToKmumu_sel_index] > 6";
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
    std::string cut = Form("BToKee_eeKFit_ee_mass[BToKee_sel_index] > %.2f && BToKee_eeKFit_ee_mass[BToKee_sel_index] < %.2f",
                           llMassBoundary.at(ij), llMassBoundary.at(ij+1));
    std::string gencut = Form("BToKee_gen_eeMass > %.2f && BToKee_gen_eeMass < %.2f",
                              llMassBoundary.at(ij), llMassBoundary.at(ij+1));
    if(!isEleFinalState){
      cut = Form("BToKmumu_mumuKFit_mumu_mass[BToKmumu_sel_index] > %.2f && BToKmumu_mumuKFit_mumu_mass[BToKmumu_sel_index] < %.2f",
		 llMassBoundary.at(ij), llMassBoundary.at(ij+1));
      gencut = Form("BToKmumu_gen_mumuMass > %.2f && BToKmumu_gen_mumuMass < %.2f",
		    llMassBoundary.at(ij), llMassBoundary.at(ij+1));

    }

    llMassCut.push_back(cut);
    if(foldGenMassBin) llGenMassCut.push_back(gencut);
    else llGenMassCut.push_back(" 1 == 1");
  }

  //non resonant first - JPsi last
  float nEv_muonTag[6] = {0.};
  float nEv_recoCand[6] = {0.};
  float nEv_chargeEff[6] = {0.};
  float nEv_alphaEff[6] = {0.};
  float nEv_vtxCLEff[6] = {0.};
  float nEv_LxyEff[6] = {0.};

  TH1F* h_Bmass_llbin[6];
  h_Bmass_llbin[0] = new TH1F(Form("h_Bmass_JPsi_%.2f-%.2f", llMassBoundary.at(3), llMassBoundary.at(3+1)), "", 1000, 0., 10.);
  if(isEleFinalState){
    nEv_muonTag[0] = t2->Draw("Muon_sel_index", (cut_muonTag+" && "+llGenMassCut.at(3)).c_str(), "goff");
    nEv_LxyEff[0] = t2->Draw(Form("BToKee_mass[BToKee_sel_index] >> %s", h_Bmass_llbin[0]->GetName()),
    (cut_muonTag+" && "+cut_Recocandidate+" && "+cut_chargeEff+" && "+cut_alphaEff+" && "+cut_vtxCLEff+" && "+cut_LxyEff+" && "+llMassCut.at(3)+" && "+llGenMassCut.at(3)).c_str());
  }
  else{
    nEv_muonTag[0] = t2->Draw("Muon_sel_index", (cut_muonTag+" && "+llGenMassCut.at(3)).c_str(), "goff");
    nEv_LxyEff[0] = t2->Draw(Form("BToKmumu_mass[BToKmumu_sel_index] >> %s", h_Bmass_llbin[0]->GetName()),
    (cut_muonTag+" && "+cut_Recocandidate+" && "+cut_chargeEff+" && "+cut_alphaEff+" && "+cut_vtxCLEff+" && "+cut_LxyEff+" && "+llMassCut.at(3)+" && "+llGenMassCut.at(3)).c_str());
  }
  
  std::cout << " tot JPsi_MC muonTag events = " << nEv_muonTag[0] << std::endl;
  std::cout << " tot JPsi_MC endSelection Events = " << nEv_LxyEff[0] << std::endl;

  for(int ij=0; ij<5; ++ij){
    h_Bmass_llbin[ij+1] = new TH1F(Form("h_Bmass_bin_%.2f-%.2f", llMassBoundary.at(ij), llMassBoundary.at(ij+1)), "", 1000, 0., 10.);

    if(isEleFinalState){
      nEv_muonTag[ij+1] = t1->Draw("Muon_sel_index", (cut_muonTag+" && "+llGenMassCut.at(ij)).c_str(), "goff");
      nEv_LxyEff[ij+1] = t1->Draw(Form("BToKee_mass[BToKee_sel_index] >> %s", h_Bmass_llbin[ij+1]->GetName()),
      (cut_muonTag+" && "+cut_Recocandidate+" && "+cut_chargeEff+" && "+cut_alphaEff+" && "+cut_vtxCLEff+" && "+cut_LxyEff+" && "+llMassCut.at(ij)+" && "+llGenMassCut.at(ij)).c_str());
    }
    else{
      nEv_muonTag[ij+1] = t1->Draw("Muon_sel_index", (cut_muonTag+" && "+llGenMassCut.at(ij)).c_str(), "goff");
      nEv_LxyEff[ij+1] = t1->Draw(Form("BToKmumu_mass[BToKmumu_sel_index] >> %s", h_Bmass_llbin[ij+1]->GetName()),
      (cut_muonTag+" && "+cut_Recocandidate+" && "+cut_chargeEff+" && "+cut_alphaEff+" && "+cut_vtxCLEff+" && "+cut_LxyEff+" && "+llMassCut.at(ij)+" && "+llGenMassCut.at(ij)).c_str());
    }
    std::cout << " tot nonResonant_MC muonTag events = " << nEv_muonTag[ij+1] << " in mass bin " << llMassCut.at(ij) << std::endl;
    std::cout << " tot nonResonant_MC endSelection Events = " << nEv_LxyEff[ij+1] << std::endl;
  }//loop


  std::string outName = "outMassHistos_MC_ee.root";
  if(!isEleFinalState) outName ="outMassHistos_MC_mumu.root";
  TFile outMassHistos(outName.c_str(), "recreate");
  outMassHistos.cd();
  for(int ij=0; ij<6; ++ij)
    h_Bmass_llbin[ij]->Write(h_Bmass_llbin[ij]->GetName());
  outMassHistos.Close();



  //now fitting
  float nEv_postFit[6] = {0.};
  float nEvError_postFit[6] = {0.};

  RooWorkspace w("w");    
  w.factory("x[0, 10]");  
  //w.factory("x[4.5, 6.]");  

  w.factory("nbackground[10000, 0, 10000]");   
  w.factory("nsignal[100, 0.0, 10000.0]");

  for(int ij=0; ij<6; ++ij){
    //w.factory("Gaussian::smodel(x,mu[5.3,4.5,6],sigma[0.05,0,0.2])");
    w.factory("RooCBShape::smodel(x,m[5.3,4.5,6],s[0.1,0.,1.],a[1.2,0.,3.],n[1,0.1,6.])");
    //w.factory("RooCBShape::CBall(x[0,15], mean[11000,13000], sigma[5000,200000], alpha[0,10000],n[0,100000])");
    RooAbsPdf * smodel = w.pdf("smodel");
    
    w.factory("Exponential::bmodel(x,tau[-2,-3,0])");
    RooAbsPdf * bmodel = w.pdf("bmodel");
    
    w.factory("SUM::model(nbackground*bmodel, nsignal*smodel)");
    RooAbsPdf * model = w.pdf("model");
    
    RooDataHist hBMass("hBMass", "hBMass", *w.var("x"), Import(*(h_Bmass_llbin[ij])));

    RooFitResult * r = model->fitTo(hBMass, Minimizer("Minuit2"),Save(true));

    RooPlot * plot = w.var("x")->frame();
    if(isEleFinalState){
      if(ij == 0) plot->SetXTitle("K(JPsi)ee mass (GeV)");
      else plot->SetXTitle("Kee mass (GeV)");
    }
    else{
      if(ij == 0) plot->SetXTitle("K(JPsi)#mu#mu mass (GeV)");
      else plot->SetXTitle("K#mu#mu mass (GeV)");
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
    if(isEleFinalState) cc->Print(Form("plots/Bmass_MC/Kee_%s.png",h_Bmass_llbin[ij]->GetName()), "png");
    else cc->Print(Form("plots/Bmass_MC/Kmumu_%s.png",h_Bmass_llbin[ij]->GetName()), "png");

    RooRealVar* parS = (RooRealVar*) r->floatParsFinal().find("nsignal");
    RooRealVar* parB = (RooRealVar*) r->floatParsFinal().find("nbackground");
    nEv_postFit[ij] = parS->getValV();
    nEvError_postFit[ij] = parS->getError();

    std::cout << " JPsi selection signal events = \t " << parS->getValV() << " error = " << parS->getError() 
	      << " bkg events = " << parB->getValV() << " error = " << parB->getError() << std::endl;
  }
  
  float errb = pow(nEvError_postFit[0]/nEv_postFit[0], 2);

  std::cout << " ***** summary ***** "<< std::endl;
  for(int ij=0; ij<6; ++ij){

    float errRatio = pow(nEvError_postFit[ij]/nEv_postFit[ij], 2) + errb;
    errRatio = sqrt(errRatio) * nEvError_postFit[ij]/nEvError_postFit[0];

    std::cout << "\n \n  category: " << h_Bmass_llbin[ij]->GetName()
              << " \t integral muonTag = " << nEv_muonTag[ij]
              << " evtCount " << nEv_LxyEff[ij] << " postFit " << nEv_postFit[ij] <<"+/-"<<nEvError_postFit[ij]
              << "\n \t\t\t eff(/muonTag)    evtCount " << nEv_LxyEff[ij]/nEv_muonTag[ij]
              << " postFit " << nEv_postFit[ij]/nEv_muonTag[ij] << "+/-" << nEvError_postFit[ij]/nEv_muonTag[ij]
              << "\n  \t\t\t double ratio wrt JPsi    evtCount " << (nEv_LxyEff[ij]/nEv_muonTag[ij])/(nEv_LxyEff[0]/nEv_muonTag[0])
              << " postFit " << (nEv_postFit[ij]/nEv_muonTag[ij])/(nEv_postFit[0]/nEv_muonTag[0]) << "+/-"<<errRatio<< std::endl;
  }

}

