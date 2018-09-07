//plot and fit BToKll 

//example to run
//root -l fitBmass_fromHistos.C'(1, "/vols/cms/amartell/BParking/data/outMassHistos_Kee_runAB_allStatTightestSelections.root")' 

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


void fitBmass_fromHistos(int isEleFinalState, std::string inFile){

  gROOT->Reset();
  gROOT->Macro("setStyle.C");

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  std::cout << " isEleFinalState = " << isEleFinalState << std::endl;


  TFile* inF = TFile::Open(inFile.c_str());

  std::vector<float> llMassBoundary;
  llMassBoundary.push_back(0.);
  llMassBoundary.push_back(1.);
  llMassBoundary.push_back(2.5);
  llMassBoundary.push_back(2.9);
  llMassBoundary.push_back(3.2);
  llMassBoundary.push_back(3.58);


  TH1F* h_Bmass[6];
  for(int ij=0; ij<6; ++ij){
    h_Bmass[ij] = (TH1F*)inF->Get(Form("Bmass_%d", ij))->Clone(Form("h_Bmass_%d", ij));
    //    h_Bmass[ij]->GetXaxis()->SetRangeUser(4.5, 6.);
  }//loop


  //now fitting
  float nEv_postFit[6] = {0.};
  float nEvError_postFit[6] = {0.};
  float nBkg_postFit[6] = {0.};
  float nBkgError_postFit[6] = {0.};
  float chi2[6] = {0.};

  for(int ij=0; ij<6; ++ij){
    RooWorkspace w("w");    
    //w.factory("x[0, 10]");  
    w.factory("x[4.5, 6.]");  
    
    w.factory("nbackground[10000, 0, 100000]");   
    w.factory("nsignal[100, 0.0, 10000]");
    
    w.factory("Gaussian::smodel(x,mu[5.3,4.5,6],sigma[0.05,0,0.2])");
    //w.factory("RooCBShape::smodel(x,m[5.3,4.5,6],s[0.1,0.,1.],a[1.2,0.,3.],n[1,0.1,6.])");

    //w.factory("RooCBShape::CBall(x[0,15], mean[11000,13000], sigma[5000,200000], alpha[0,10000],n[0,100000])");
    RooAbsPdf * smodel = w.pdf("smodel");
    
    w.factory("Exponential::bmodel(x,tau[-2,-3,0])");
    RooAbsPdf * bmodel = w.pdf("bmodel");
    
    w.factory("SUM::model(nbackground * bmodel, nsignal * smodel)");
    RooAbsPdf * model = w.pdf("model");
    
    RooDataHist hBMass("hBMass", "hBMass", *w.var("x"), Import(*(h_Bmass[ij])));
    w.Print();

    RooFitResult * r = model->fitTo(hBMass, Minimizer("Minuit2"),Save(true));
    std::cout << " fit status = " << r->status() << std::endl;

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
    chi2[ij] = plot->chiSquare();


    TCanvas * cc = new TCanvas();
    cc->SetLogy(0);
    plot->Draw();
    if(isEleFinalState) cc->Print(Form("plots/Bmass_DA/Kee_%s.png",h_Bmass[ij]->GetName()), "png");
    else cc->Print(Form("plots/Bmass_DA/Kmumu_%s.png",h_Bmass[ij]->GetName()), "png");

    RooRealVar* parS = (RooRealVar*) r->floatParsFinal().find("nsignal");
    RooRealVar* parB = (RooRealVar*) r->floatParsFinal().find("nbackground");
    nEv_postFit[ij] = parS->getValV();
    nEvError_postFit[ij] = parS->getError();
    nBkg_postFit[ij] = parB->getValV();
    nBkgError_postFit[ij] = parB->getError();

    std::cout << " JPsi selection signal events = \t " << parS->getValV() << " error = " << parS->getError() 
	      << " bkg events = " << parB->getValV() << " error = " << parB->getError() << std::endl;
  }
  

  std::cout << " ***** summary ***** "<< std::endl;
  for(int ij=0; ij<6; ++ij){

    std::cout << "\n \n  category: " << h_Bmass[ij]->GetName()
	      << "\n \t signal = " << nEv_postFit[ij] << "+/-" << nEvError_postFit[ij]
              << "\n \t bkg = " << nBkg_postFit[ij] << "+/-" << nBkgError_postFit[ij] << " chi2 = " << chi2[ij] << std::endl;
  }

}

