// Author: T. Strebler (IC)
// Date:   31 May 2018
//
// Add new variables from NanoAOD tree for BToKpipi analysis
// Can be either included in a interpreted macro or compiled in c++
// (use `root-config --glibs --cflags`)
//



#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"
#include "TLorentzVector.h"

#include "NanoAODTree.h"


float D0Mass_ = 1.864;
float BuMass_ = 5.279;
float KaonMass_ = 0.493677;
float PionMass_ = 0.139570;


int main(int argc, char** argv) {

  string status_sample = *(argv + 1);
  bool isMC = false;
  bool isData = false;
  if (status_sample.compare("mc") == 0) isMC = true;
  if (status_sample.compare("data") == 0) isData = true;


  string output;
  for (int i = 1; i < argc; ++i) {
    if(std::string(argv[i]) == "--output") {
      if (i + 1 < argc) {
	output = argv[i+1];
	break;
      } else {
	std::cerr << "--output option requires one argument." << std::endl;
	return 1;
      }      
    }  
  }
  if(output==""){
    std::cerr << "--output argument required" << std::endl;
    return 1;
  }
    
  

  string input;
  for (int i = 1; i < argc; ++i) {
    if(std::string(argv[i]) == "--input") {
      if (i + 1 < argc) {
	input = argv[i+1];
	break;
      } else {
	std::cerr << "--intput option requires one argument." << std::endl;
	return 1;
      }      
    }  
  }
  if(input==""){
    std::cerr << "--input argument required" << std::endl;
    return 1;
  }


  bool overwrite;
  for (int i = 1; i < argc; ++i) {
    if(std::string(argv[i]) == "--overwrite") {
      overwrite = true;
      break;
    }
  }


  bool saveFullNanoAOD;
  for (int i = 1; i < argc; ++i) {
    if(std::string(argv[i]) == "--saveFullNanoAOD") {
      saveFullNanoAOD = true;
      break;
    }
  }
 

  TFile* f_new = TFile::Open(output.c_str());
  if(f_new!=0 && !overwrite){
    cout<<output<<" already exists, please delete it before converting again"<<endl;
    return 0;
  }
  f_new = TFile::Open(output.c_str(),"RECREATE");


  TChain* oldtree = new TChain("Events");
  oldtree->Add(input.c_str());
  NanoAODTree* tree = new NanoAODTree(oldtree);

  TTree* tree_new=new TTree("BToKpipiTree","BToKpipiTree");
  if(saveFullNanoAOD)
    tree_new=tree->GetTree()->CloneTree(0);


  //New branches

  int _Muon_sel_index = -1;
  int _BToKpipi_sel_index = -1;

  tree_new->Branch("Muon_sel_index",&_Muon_sel_index,"Muon_sel_index/I");
  tree_new->Branch("BToKpipi_sel_index",&_BToKpipi_sel_index,"BToKpipi_sel_index/I");

  int _GenPart_BToKpipi_index = -1;
  int _GenPart_D0FromB_index = -1;
  int _GenPart_piFromB_index = -1;
  int _GenPart_KFromD0_index = -1;
  int _GenPart_piFromD0_index = -1;
  int _BToKpipi_gen_index = -1;

  if(isMC){

    tree_new->Branch("GenPart_BToKpipi_index",&_GenPart_BToKpipi_index,"GenPart_BToKpipi_index/I");
    tree_new->Branch("GenPart_D0FromB_index",&_GenPart_D0FromB_index,"GenPart_D0FromB_index/I");
    tree_new->Branch("GenPart_piFromB_index",&_GenPart_piFromB_index,"GenPart_piFromB_index/I");
    tree_new->Branch("GenPart_KFromD0_index",&_GenPart_KFromD0_index,"GenPart_KFromD0_index/I");
    tree_new->Branch("GenPart_piFromD0_index",&_GenPart_piFromD0_index,"GenPart_piFromD0_index/I");
    tree_new->Branch("BToKpipi_gen_index",&_BToKpipi_gen_index,"BToKpipi_gen_index/I");

  }

  int nentries = tree->GetEntries();
  cout<<"Nentries="<<nentries<<endl;
  cout<<"isMC="<<isMC<<endl;

  for (int iEntry = 0; iEntry < nentries ; iEntry++){

    tree->GetEntry(iEntry);

    //if(iEntry%10000==0) cout<<"Entry #"<<iEntry<<" "<< int(100*float(iEntry)/nentries)<<"%"<<endl;
    cout<<"Entry #"<<iEntry<<" "<< int(100*float(iEntry)/nentries)<<"%"<<endl;

    _Muon_sel_index = -1;
    _BToKpipi_sel_index = -1;

    _GenPart_BToKpipi_index = -1;
    _GenPart_D0FromB_index = -1;
    _GenPart_piFromB_index = -1;
    _GenPart_KFromD0_index = -1;
    _GenPart_piFromD0_index = -1;
    _BToKpipi_gen_index = -1;

    //Select the muon

    int nMuon = tree->nMuon;

    for(int i_mu=0; i_mu<nMuon; i_mu++){

      //Selections on the muon to refine (in particular trigger matching)
      if( 1 ){
	_Muon_sel_index = i_mu;
	break; //Take leading muon passing the selections (muons are pt-ordered)
      }

    }

    if(_Muon_sel_index <0){
      tree_new->Fill();
      continue;
    }

    //Select the BToKpipi candidate with reco criteria

    int nBToKpipi = tree->nBToKpipi;
    float best_D0_mass = -1.;
    float best_Bu_mass = -1.;

    for(int i_BToKpipi=0; i_BToKpipi<nBToKpipi; nBToKpipi++){            

      if(tree->BToKpipi_piBu_charge[i_BToKpipi]*tree->Muon_charge[_Muon_sel_index]>0) continue; //Only consider BToKpipi with opposite charge to muon
      
      float Kpi_mass = tree->BToKpipi_Kpi_mass[i_BToKpipi];
      float Kpi_CL_vtx = tree->BToKpipi_Kpi_CL_vtx[i_BToKpipi];
	
      //D0 selection
      if( !(best_D0_mass < 0. 
	    || abs(best_D0_mass-Kpi_mass)<1e-3 //Several BToKpipi can share the same D0->Kpi
	    || abs(Kpi_mass-D0Mass_) < abs(best_D0_mass-D0Mass_)) )       
	continue;

      //if( Kpi_CL_vtx < min_CL) continue; //cut on Kpi vtx refitting

      float B_mass = tree->BToKpipi_mass[i_BToKpipi];
      float B_CL_vtx = tree->BToKpipi_CL_vtx[i_BToKpipi];
      
      if( !(best_Bu_mass < 0. 
	    || abs(B_mass-BuMass_) < abs(best_Bu_mass-BuMass_)) )       
	continue;
      
      best_D0_mass = Kpi_mass;
      best_Bu_mass = B_mass;
      _BToKpipi_sel_index = i_BToKpipi;

    }


    
    //Select the BToKpipi candidate based on gen matching

    if(isMC){
      
      int nGenPart = tree->nGenPart;
      
      for(int i_Bu=0; i_Bu<nGenPart; i_Bu++){

	if(abs(tree->GenPart_pdgId[i_Bu])==521){

	  for(int i_gen=0; i_gen<nGenPart; i_gen++){
	    int pdgId = tree->GenPart_pdgId[i_gen];
	    int mother_index = tree->GenPart_genPartIdxMother[i_gen];
	    if(abs(pdgId)==421 && mother_index == i_Bu)
	      _GenPart_D0FromB_index = i_gen;
	    else if(abs(pdgId)==211 && mother_index == i_Bu)
	      _GenPart_piFromB_index = i_gen;
	    if(_GenPart_D0FromB_index>=0 && _GenPart_piFromB_index>=0){
	      _GenPart_BToKpipi_index = i_Bu;
	      break;
	    }
	  }	  

	}
	
	if(_GenPart_BToKpipi_index>=0) break;

      }


      
      for(int i_gen=0; i_gen<nGenPart; i_gen++){

	int pdgId = tree->GenPart_pdgId[i_gen];
	int mother_index = tree->GenPart_genPartIdxMother[i_gen];
	if(abs(pdgId)==321 && mother_index == _GenPart_D0FromB_index)
	  _GenPart_KFromD0_index = i_gen;
	if(abs(pdgId)==211 && mother_index == _GenPart_D0FromB_index)
	  _GenPart_piFromD0_index = i_gen;
	if(_GenPart_KFromD0_index>=0 && _GenPart_piFromD0_index>=0) break;

      }


      TLorentzVector gen_piFromB_tlv;
      TLorentzVector gen_KFromD0_tlv;
      TLorentzVector gen_piFromD0_tlv;
      
      gen_piFromB_tlv.SetPtEtaPhiM(tree->GenPart_pt[_GenPart_piFromB_index],
				   tree->GenPart_eta[_GenPart_piFromB_index],
				   tree->GenPart_phi[_GenPart_piFromB_index],
				   PionMass_);
      gen_KFromD0_tlv.SetPtEtaPhiM(tree->GenPart_pt[_GenPart_KFromD0_index],
				   tree->GenPart_eta[_GenPart_KFromD0_index],
				   tree->GenPart_phi[_GenPart_KFromD0_index],
				   KaonMass_);
      gen_piFromD0_tlv.SetPtEtaPhiM(tree->GenPart_pt[_GenPart_piFromD0_index],
				    tree->GenPart_eta[_GenPart_piFromD0_index],
				    tree->GenPart_phi[_GenPart_piFromD0_index],
				    PionMass_);

      float best_dR = -1.;

      for(int i_BToKpipi=0; i_BToKpipi<nBToKpipi; nBToKpipi++){

	TLorentzVector piFromB_tlv;
	TLorentzVector KFromD0_tlv;
	TLorentzVector piFromD0_tlv;
	
	piFromB_tlv.SetPtEtaPhiM(tree->BToKpipi_piBu_pt[i_BToKpipi],
				 tree->BToKpipi_piBu_eta[i_BToKpipi],
				 tree->BToKpipi_piBu_phi[i_BToKpipi],
				 PionMass_);
	KFromD0_tlv.SetPtEtaPhiM(tree->BToKpipi_kaon_pt[i_BToKpipi],
				 tree->BToKpipi_kaon_eta[i_BToKpipi],
				 tree->BToKpipi_kaon_phi[i_BToKpipi],
				 KaonMass_);
	piFromD0_tlv.SetPtEtaPhiM(tree->BToKpipi_piD0_pt[i_BToKpipi],
				  tree->BToKpipi_piD0_eta[i_BToKpipi],
				  tree->BToKpipi_piD0_phi[i_BToKpipi],
				  PionMass_);

	float dR_piFromB = piFromB_tlv.DeltaR(gen_piFromB_tlv);
	float dR_KFromD0 = KFromD0_tlv.DeltaR(gen_KFromD0_tlv);
	float dR_piFromD0 = piFromD0_tlv.DeltaR(gen_piFromD0_tlv);

	float dR_tot = dR_piFromB + dR_KFromD0 + dR_piFromD0; //In case several BToKpipi matches, take the closest one in dR_tot

	if( dR_piFromB<0.1 && dR_KFromD0<0.1 && dR_piFromD0
	    && (best_dR<0. || dR_tot<best_dR) ){
	  best_dR = dR_tot;
	  _BToKpipi_gen_index = i_BToKpipi;	  
	}

      }

    }

    tree_new->Fill();

  }


  f_new->cd();
  if(!saveFullNanoAOD) tree_new->AddFriend("Events",input.c_str());

  tree_new->Write();
  f_new->Close();
  return 0;

}



