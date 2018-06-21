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

float PFCand_pt_Min_ = 2.;
float piBu_pt_Min_ = 2.;
float Kpi_CL_vtx_Min_ = 0.05;
float Kpi_abs_dz_Max_ = 1;
float Kpi_pt_Min_ = 5.;
float B_CL_vtx_Min_ = 0.;
float B_Lxy_Min_ = 5.;
float B_cosAlpha_Min_ = 0.999;


int main(int argc, char** argv) {

  string status_sample = *(argv + 1);
  bool isMC = false;
  bool isData = false;
  if (status_sample.compare("mc") == 0) isMC = true;
  if (status_sample.compare("data") == 0) isData = true;

  bool isBPHParking = false;
  for (int i = 1; i < argc; ++i) {
    if(std::string(argv[i]) == "BPHParking") {
      isBPHParking = true;
      break;
    }
  }

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


  string sign;
  for (int i = 1; i < argc; ++i) {
    if(std::string(argv[i]) == "--sign") {
      if (i + 1 < argc) {
	sign = argv[i+1];
	break;
      } else {
	std::cerr << "--sign option requires one argument." << std::endl;
	return 1;
      }
    }
  }


  bool saveFullNanoAOD = true;
  /*for (int i = 1; i < argc; ++i) {
    if(std::string(argv[i]) == "--saveFullNanoAOD") {
      saveFullNanoAOD = true;
      break;
    }
    }*/
 

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
  int nentries = oldtree->GetEntries();
  if(saveFullNanoAOD)
    tree_new=oldtree->GetTree()->CloneTree(0);


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
  int _BToKpipi_D0_gen_index = -1; //In case we can find at least the D0 candidates matched

  if(isMC){

    tree_new->Branch("GenPart_BToKpipi_index",&_GenPart_BToKpipi_index,"GenPart_BToKpipi_index/I");
    tree_new->Branch("GenPart_D0FromB_index",&_GenPart_D0FromB_index,"GenPart_D0FromB_index/I");
    tree_new->Branch("GenPart_piFromB_index",&_GenPart_piFromB_index,"GenPart_piFromB_index/I");
    tree_new->Branch("GenPart_KFromD0_index",&_GenPart_KFromD0_index,"GenPart_KFromD0_index/I");
    tree_new->Branch("GenPart_piFromD0_index",&_GenPart_piFromD0_index,"GenPart_piFromD0_index/I");
    tree_new->Branch("BToKpipi_gen_index",&_BToKpipi_gen_index,"BToKpipi_gen_index/I");
    tree_new->Branch("BToKpipi_D0_gen_index",&_BToKpipi_D0_gen_index,"BToKpipi_D0_gen_index/I");

  }


  bool _HLT_Mu8p5_IP3p5 = false;
  bool _HLT_Mu10p5_IP3p5 = false;
  bool _HLT_Mu9_IP6 = false;
  bool _HLT_Mu8_IP3 = false;
  bool _HLT_BPHParking = false;

  if(isBPHParking){

    tree_new->Branch("HLT_Mu8p5_IP3p5",&_HLT_Mu8p5_IP3p5,"HLT_Mu8p5_IP3p5/O");
    tree_new->Branch("HLT_Mu10p5_IP3p5",&_HLT_Mu10p5_IP3p5,"HLT_Mu10p5_IP3p5/O");
    tree_new->Branch("HLT_Mu9_IP6",&_HLT_Mu9_IP6,"HLT_Mu9_IP6/O");
    tree_new->Branch("HLT_Mu8_IP3",&_HLT_Mu8_IP3,"HLT_Mu8_IP3/O");
    tree_new->Branch("HLT_BPHParking",&_HLT_BPHParking,"HLT_BPHParking/O");

  }

  int _PFCand_genPiFromB_index = -1;
  int _PFCand_genKFromD0_index = -1;
  int _PFCand_genPiFromD0_index = -1;

  tree_new->Branch("PFCand_genPiFromB_index",&_PFCand_genPiFromB_index,"PFCand_genPiFromB_index/I");
  tree_new->Branch("PFCand_genKFromD0_index",&_PFCand_genKFromD0_index,"PFCand_genKFromD0_index/I");
  tree_new->Branch("PFCand_genPiFromD0_index",&_PFCand_genPiFromD0_index,"PFCand_genPiFromD0_index/I");


  cout<<"Nentries="<<nentries<<endl;
  cout<<"isMC="<<isMC<<endl;

  for (int iEntry = 0; iEntry < nentries ; iEntry++){

    tree->GetEntry(iEntry);

    if(iEntry%10000==0) cout<<"Entry #"<<iEntry<<" "<< int(100*float(iEntry)/nentries)<<"%"<<endl;


    _Muon_sel_index = -1;
    _BToKpipi_sel_index = -1;

    _GenPart_BToKpipi_index = -1;
    _GenPart_D0FromB_index = -1;
    _GenPart_piFromB_index = -1;
    _GenPart_KFromD0_index = -1;
    _GenPart_piFromD0_index = -1;
    _BToKpipi_gen_index = -1;
    _BToKpipi_D0_gen_index = -1;

    _HLT_Mu8p5_IP3p5 = false;
    _HLT_Mu10p5_IP3p5 = false;
    _HLT_Mu9_IP6 = false;
    _HLT_Mu8_IP3 = false;
    _HLT_BPHParking = false;

    _PFCand_genPiFromB_index = -1;
    _PFCand_genKFromD0_index = -1;
    _PFCand_genPiFromD0_index = -1;

    //Select the muon

    int nMuon = tree->nMuon;

    for(int i_mu=0; i_mu<nMuon; i_mu++){

      bool isMuSel = tree->Muon_softId[i_mu];

      //Trigger selection + matching

      //Only for isBPHParking for now
      if(isMC){
	isMuSel &= true;
      }


      if(isBPHParking){

	_HLT_Mu8p5_IP3p5 = tree->HLT_Mu8p5_IP3p5_part0
	  || tree->HLT_Mu8p5_IP3p5_part1
	  || tree->HLT_Mu8p5_IP3p5_part2
	  || tree->HLT_Mu8p5_IP3p5_part3
	  || tree->HLT_Mu8p5_IP3p5_part4
	  || tree->HLT_Mu8p5_IP3p5_part5;
	_HLT_Mu10p5_IP3p5 = tree->HLT_Mu10p5_IP3p5_part0
	  || tree->HLT_Mu10p5_IP3p5_part1
	  || tree->HLT_Mu10p5_IP3p5_part2
	  || tree->HLT_Mu10p5_IP3p5_part3
	  || tree->HLT_Mu10p5_IP3p5_part4
	  || tree->HLT_Mu10p5_IP3p5_part5;
	_HLT_Mu9_IP6 = tree->HLT_Mu9_IP6_part0
	  || tree->HLT_Mu9_IP6_part1
	  || tree->HLT_Mu9_IP6_part2
	  || tree->HLT_Mu9_IP6_part3
	  || tree->HLT_Mu9_IP6_part4
	  || tree->HLT_Mu9_IP6_part5;
	_HLT_Mu8_IP3 = tree->HLT_Mu8_IP3_part0
	  || tree->HLT_Mu8_IP3_part1
	  || tree->HLT_Mu8_IP3_part2
	  || tree->HLT_Mu8_IP3_part3
	  || tree->HLT_Mu8_IP3_part4
	  || tree->HLT_Mu8_IP3_part5;
	_HLT_BPHParking = _HLT_Mu8p5_IP3p5 || _HLT_Mu10p5_IP3p5 || _HLT_Mu9_IP6 || _HLT_Mu8_IP3;

	TLorentzVector mu;
	mu.SetPtEtaPhiM(tree->Muon_pt[i_mu],tree->Muon_eta[i_mu],tree->Muon_phi[i_mu],tree->Muon_mass[i_mu]);

	bool isTrigMatched = false;
	int nTrigObj = tree->nTrigObj;
	for(unsigned int i_trig = 0; i_trig<nTrigObj; i_trig++){

	  if( tree->TrigObj_id[i_trig]==13 && ((tree->TrigObj_filterBits[i_trig])>>3)&1 ){
	    TLorentzVector trig;
	    trig.SetPtEtaPhiM(tree->TrigObj_pt[i_trig],tree->TrigObj_eta[i_trig],tree->TrigObj_phi[i_trig],0);
	    float dR = mu.DeltaR(trig);
	    if(dR<0.1){
	      isTrigMatched = true;
	      break;
	    }
	  }

	}

	isMuSel &= (_HLT_BPHParking && isTrigMatched);

      }


      if( isMuSel ){
	_Muon_sel_index = i_mu;
	break; //Take leading muon passing the selections (muons are pt-ordered)
      }

    }

    if(_Muon_sel_index <0){
      //Let's skim events which do not have a tag muon (including trigger)
      //tree_new->Fill();
      continue;
    }

    //Select the BToKpipi candidate with reco criteria

    int nBToKpipi = tree->nBToKpipi;
    float best_D0_mass = -1.;
    float best_Bu_mass = -1.;

    for(int i_BToKpipi=0; i_BToKpipi<nBToKpipi; i_BToKpipi++){            

      if(sign=="isOS" && tree->BToKpipi_piBu_charge[i_BToKpipi]*tree->Muon_charge[_Muon_sel_index]>0) continue; //Only consider BToKpipi with opposite charge to muon
      if(sign=="isSS" && tree->BToKpipi_piBu_charge[i_BToKpipi]*tree->Muon_charge[_Muon_sel_index]<0) continue; //Only keep SS for control region
      if(tree->BToKpipi_piBu_charge[i_BToKpipi]*tree->BToKpipi_kaon_charge[i_BToKpipi]<0) continue; //Remove D0/D0bar suppressed decay

      float Kpi_CL_vtx = tree->BToKpipi_Kpi_CL_vtx[i_BToKpipi];
      float delta_dz = fabs(tree->BToKpipi_kaon_dz[i_BToKpipi]-tree->BToKpipi_piD0_dz[i_BToKpipi]);
      float Kpi_pt = tree->BToKpipi_Kpi_pt[i_BToKpipi];
      
      float B_CL_vtx = tree->BToKpipi_CL_vtx[i_BToKpipi];
      float B_Lxy = tree->BToKpipi_Lxy[i_BToKpipi];
      float B_cosAlpha = tree->BToKpipi_cosAlpha[i_BToKpipi];

      float piBu_pt = tree->BToKpipi_piBu_pt[i_BToKpipi];
      float piD0_pt = tree->BToKpipi_piD0_pt[i_BToKpipi];
      float kaon_pt = tree->BToKpipi_piD0_pt[i_BToKpipi];

      if( piBu_pt < piBu_pt_Min_ || piD0_pt < PFCand_pt_Min_ || kaon_pt < PFCand_pt_Min_ ) continue;

      if( Kpi_CL_vtx < Kpi_CL_vtx_Min_) continue; //cut on Kpi vtx refitting
      if(delta_dz > Kpi_abs_dz_Max_) continue;
      if(Kpi_pt < Kpi_pt_Min_) continue;

      if(B_CL_vtx <= B_CL_vtx_Min_) continue;
      if(B_Lxy < B_Lxy_Min_) continue;
      if(B_cosAlpha < B_cosAlpha_Min_) continue;

      float Kpi_mass = tree->BToKpipi_Kpi_mass[i_BToKpipi];

      //D0 selection
      if ( !(best_CL_D0 < 0.
	     || abs(best_CL_D0-Kpi_CL_vtx)<1e-3 //Several BToKpipi can share the same D0->Kpi
	     || Kpi_CL_vtx > best_CL_D0) )
	continue;

      float B_mass = tree->BToKpipi_mass[i_BToKpipi];
      
      if( !(best_Bu_mass < 0. 
	    || abs(B_mass-BuMass_) < abs(best_Bu_mass-BuMass_)) )       
	continue;
      
      best_CL_D0 = Kpi_CL_vtx;
      best_Bu_mass = B_mass;
      _BToKpipi_sel_index = i_BToKpipi;

    }




    if(isMC){
      
      int nGenPart = tree->nGenPart;

      for(int i_Bu=0; i_Bu<nGenPart; i_Bu++){

	_GenPart_D0FromB_index = -1;
	_GenPart_piFromB_index = -1;

	if(abs(tree->GenPart_pdgId[i_Bu])==521){

	  bool additional_daughters = false;

	  for(int i_gen=0; i_gen<nGenPart; i_gen++){
	    int pdgId = tree->GenPart_pdgId[i_gen];
	    int mother_index = tree->GenPart_genPartIdxMother[i_gen];
	    if(mother_index == i_Bu){
	      if(abs(pdgId)==421 && _GenPart_D0FromB_index<0) _GenPart_D0FromB_index = i_gen;
	      else if(abs(pdgId)==211 && _GenPart_piFromB_index<0) _GenPart_piFromB_index = i_gen;
	      else if(abs(pdgId)!=22) additional_daughters = true;
	    }
	  }

	  if(_GenPart_D0FromB_index>=0 && _GenPart_piFromB_index>=0 && !additional_daughters){
	    _GenPart_BToKpipi_index = i_Bu;
	    break;
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


      //Select the BToKpipi candidate based on gen matching

      for(int i_BToKpipi=0; i_BToKpipi<nBToKpipi; i_BToKpipi++){

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

	if( dR_piFromB<0.1 && dR_KFromD0<0.1 && dR_piFromD0<0.1
	    && (best_dR<0. || dR_tot<best_dR) ){
	  best_dR = dR_tot;
	  _BToKpipi_gen_index = i_BToKpipi;
	  _BToKpipi_D0_gen_index = i_BToKpipi;
	}

	else if(dR_KFromD0<0.1 && dR_piFromD0<0.1){
	  _BToKpipi_D0_gen_index = i_BToKpipi;
	}

      }


      //Select the PFCand based on gen matching

      float best_dR_piFromB = -1.;
      float best_dR_KFromD0 = -1.;
      float best_dR_piFromD0 = -1.;

      int nPFCand = tree->nPFCand;
      for(int i_pf=0;i_pf<nPFCand;i_pf++){

	if(abs(tree->PFCand_pdgId[i_pf])!=211) continue;

	TLorentzVector PFCand_tlv;
	PFCand_tlv.SetPtEtaPhiM(tree->PFCand_pt[i_pf],tree->PFCand_eta[i_pf],tree->PFCand_phi[i_pf],tree->PFCand_mass[i_pf]);

	float dR_piFromB = PFCand_tlv.DeltaR(gen_piFromB_tlv);

	if(dR_piFromB<0.1 && (best_dR_piFromB<0. || dR_piFromB<best_dR_piFromB)){
	  _PFCand_genPiFromB_index = i_pf;
	  best_dR_piFromB = dR_piFromB;
	}

	float dR_KFromD0 = PFCand_tlv.DeltaR(gen_KFromD0_tlv);
	if(dR_KFromD0<0.1 && (best_dR_KFromD0<0. || dR_KFromD0<best_dR_KFromD0)){
	  _PFCand_genKFromD0_index = i_pf;
	  best_dR_KFromD0 = dR_KFromD0;
	}

	float dR_piFromD0 = PFCand_tlv.DeltaR(gen_piFromD0_tlv);
	if(dR_piFromD0<0.1 && (best_dR_piFromD0<0. || dR_piFromD0<best_dR_piFromD0)){
	  _PFCand_genPiFromD0_index = i_pf;
	  best_dR_piFromD0 = dR_piFromD0;
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



