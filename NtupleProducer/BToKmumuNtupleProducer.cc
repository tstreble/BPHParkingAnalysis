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


float JPsiMass_ = 3.0969;
float BuMass_ = 5.279;
float KaonMass_ = 0.493677;
float MuonMass_ = 0.10565837;


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

  //Always saveFullNanoAOD info because we are skimming
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

  TTree* tree_new=new TTree("BToKmumuTree","BToKmumuTree");
  if(saveFullNanoAOD)
    tree_new=tree->GetTree()->CloneTree(0);


  //New branches

  int _Muon_sel_index = -1;
  int _BToKmumu_sel_index = -1;

  tree_new->Branch("Muon_sel_index",&_Muon_sel_index,"Muon_sel_index/I");
  tree_new->Branch("BToKmumu_sel_index",&_BToKmumu_sel_index,"BToKmumu_sel_index/I");

  int _GenPart_BToKmumu_index = -1;
  int _GenPart_JPsiFromB_index = -1;
  int _GenPart_KFromB_index = -1;
  int _GenPart_mu1FromJPsi_index = -1;
  int _GenPart_mu2FromJPsi_index = -1;
  int _BToKmumu_gen_index = -1;

  if(isMC){

    tree_new->Branch("GenPart_BToKmumu_index",&_GenPart_BToKmumu_index,"GenPart_BToKmumu_index/I");
    tree_new->Branch("GenPart_JPsiFromB_index",&_GenPart_JPsiFromB_index,"GenPart_JPsiFromB_index/I");
    tree_new->Branch("GenPart_KFromB_index",&_GenPart_KFromB_index,"GenPart_KFromB_index/I");
    tree_new->Branch("GenPart_mu1FromJPsi_index",&_GenPart_mu1FromJPsi_index,"GenPart_mu1FromJPsi_index/I");
    tree_new->Branch("GenPart_mu2FromJPsi_index",&_GenPart_mu2FromJPsi_index,"GenPart_mu2FromJPsi_index/I");
    tree_new->Branch("BToKmumu_gen_index",&_BToKmumu_gen_index,"BToKmumu_gen_index/I");

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


  int nentries = tree->GetEntries();
  cout<<"Nentries="<<nentries<<endl;
  cout<<"isMC="<<isMC<<endl;

  for (int iEntry = 0; iEntry < nentries ; iEntry++){

    tree->GetEntry(iEntry);

    if(iEntry%10000==0) cout<<"Entry #"<<iEntry<<" "<< int(100*float(iEntry)/nentries)<<"%"<<endl;

    _Muon_sel_index = -1;
    _BToKmumu_sel_index = -1;

    _GenPart_BToKmumu_index = -1;
    _GenPart_JPsiFromB_index = -1;
    _GenPart_KFromB_index = -1;
    _GenPart_mu1FromJPsi_index = -1;
    _GenPart_mu2FromJPsi_index = -1;
    _BToKmumu_gen_index = -1;

    _HLT_Mu8p5_IP3p5 = false;
    _HLT_Mu10p5_IP3p5 = false;
    _HLT_Mu9_IP6 = false;
    _HLT_Mu8_IP3 = false;
    _HLT_BPHParking = false;

    //Select the muon

    int nMuon = tree->nMuon;

    for(int i_mu=0; i_mu<nMuon; i_mu++){

      bool isMuSel = tree->Muon_softId[i_mu];

      //Trigger selection + matching

      //Only for isBPHParking for now
      if(isMC){
	isMuSel = true;
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

      //Should implement something to avoid overlap with BToKmm final state (nMuon<2 won't work because only muons with pT>3 GeV are stored in default NanoAOD)
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

    //Select the BToKmumu candidate with reco criteria

    int nBToKmumu = tree->nBToKmumu;
    float best_JPsi_mass = -1.;
    float best_Bu_mass = -1.;

    for(int i_BToKmumu=0; i_BToKmumu<nBToKmumu; i_BToKmumu++){            

      if(tree->BToKmumu_kaon_charge[i_BToKmumu]*tree->Muon_charge[_Muon_sel_index]>0) continue; //Only consider BToKmumu with opposite charge to muon
      
      TLorentzVector mu_sel;
      TLorentzVector mu1;
      TLorentzVector mu2;

      mu_sel.SetPtEtaPhiM(tree->Muon_pt[_Muon_sel_index],
			  tree->Muon_eta[_Muon_sel_index],
			  tree->Muon_phi[_Muon_sel_index],
			  tree->Muon_mass[_Muon_sel_index]);

      mu1.SetPtEtaPhiM(tree->BToKmumu_mu1_pt[i_BToKmumu],
		       tree->BToKmumu_mu1_eta[i_BToKmumu],
		       tree->BToKmumu_mu1_phi[i_BToKmumu],
		       MuonMass_);

      mu2.SetPtEtaPhiM(tree->BToKmumu_mu2_pt[i_BToKmumu],
		       tree->BToKmumu_mu2_eta[i_BToKmumu],
		       tree->BToKmumu_mu2_phi[i_BToKmumu],
		       MuonMass_);

      if(mu1.DeltaR(mu_sel)<0.1 || mu2.DeltaR(mu_sel)<0.1) continue; //Avoid any ambibuity between tag muon and probe muons

      float mumu_mass = tree->BToKmumu_mumu_mass[i_BToKmumu];
      float mumu_CL_vtx = tree->BToKmumu_mumu_CL_vtx[i_BToKmumu];
	
      //JPsi selection
      if( !(best_JPsi_mass < 0. 
	    || abs(best_JPsi_mass-mumu_mass)<1e-3 //Several BToKmumu can share the same JPsi->mumu
	    || abs(mumu_mass-JPsiMass_) < abs(best_JPsi_mass-JPsiMass_)) )       
	continue;

      //if( mumu_CL_vtx < min_CL) continue; //cut on mumu vtx refitting

      float B_mass = tree->BToKmumu_mass[i_BToKmumu];
      float B_CL_vtx = tree->BToKmumu_CL_vtx[i_BToKmumu];
      
      if( !(best_Bu_mass < 0. 
	    || abs(B_mass-BuMass_) < abs(best_Bu_mass-BuMass_)) )       
	continue;
      
      best_JPsi_mass = mumu_mass;
      best_Bu_mass = B_mass;
      _BToKmumu_sel_index = i_BToKmumu;

    }


    
    //Select the BToKmumu candidate based on gen matching

    if(isMC){
      
      int nGenPart = tree->nGenPart;
      
      for(int i_Bu=0; i_Bu<nGenPart; i_Bu++){

	if(abs(tree->GenPart_pdgId[i_Bu])==521){

	  for(int i_gen=0; i_gen<nGenPart; i_gen++){
	    int pdgId = tree->GenPart_pdgId[i_gen];
	    int mother_index = tree->GenPart_genPartIdxMother[i_gen];
	    if(abs(pdgId)==443 && mother_index == i_Bu)
	      _GenPart_JPsiFromB_index = i_gen;
	    else if(abs(pdgId)==321 && mother_index == i_Bu)
	      _GenPart_KFromB_index = i_gen;
	    if(_GenPart_JPsiFromB_index>=0 && _GenPart_KFromB_index>=0){
	      _GenPart_BToKmumu_index = i_Bu;
	      break;
	    }
	  }	  

	}
	
	if(_GenPart_BToKmumu_index>=0) break;

      }


      
      for(int i_gen=0; i_gen<nGenPart; i_gen++){

	int pdgId = tree->GenPart_pdgId[i_gen];
	int mother_index = tree->GenPart_genPartIdxMother[i_gen];
	if(abs(pdgId)==13 && mother_index == _GenPart_JPsiFromB_index && _GenPart_mu1FromJPsi_index<0)
	  _GenPart_mu1FromJPsi_index = i_gen;
	else if(abs(pdgId)==13 && mother_index == _GenPart_JPsiFromB_index)
	  _GenPart_mu2FromJPsi_index = i_gen;
	if(_GenPart_mu1FromJPsi_index>=0 && _GenPart_mu2FromJPsi_index>=0) break;

      }

      //mu1FromJPsi stored a leading daughter
      if(tree->GenPart_pt[_GenPart_mu2FromJPsi_index]>tree->GenPart_pt[_GenPart_mu1FromJPsi_index]){
	int i_temp = _GenPart_mu1FromJPsi_index;
	_GenPart_mu1FromJPsi_index = _GenPart_mu2FromJPsi_index;
	_GenPart_mu2FromJPsi_index = i_temp;
      }


      TLorentzVector gen_KFromB_tlv;
      TLorentzVector gen_mu1FromJPsi_tlv;
      TLorentzVector gen_mu2FromJPsi_tlv;
      
      gen_KFromB_tlv.SetPtEtaPhiM(tree->GenPart_pt[_GenPart_KFromB_index],
				   tree->GenPart_eta[_GenPart_KFromB_index],
				   tree->GenPart_phi[_GenPart_KFromB_index],
				   KaonMass_);
      gen_mu1FromJPsi_tlv.SetPtEtaPhiM(tree->GenPart_pt[_GenPart_mu1FromJPsi_index],
				       tree->GenPart_eta[_GenPart_mu1FromJPsi_index],
				       tree->GenPart_phi[_GenPart_mu1FromJPsi_index],
				       MuonMass_);
      gen_mu2FromJPsi_tlv.SetPtEtaPhiM(tree->GenPart_pt[_GenPart_mu2FromJPsi_index],
				       tree->GenPart_eta[_GenPart_mu2FromJPsi_index],
				       tree->GenPart_phi[_GenPart_mu2FromJPsi_index],
				       MuonMass_);

      float best_dR = -1.;

      for(int i_BToKmumu=0; i_BToKmumu<nBToKmumu; i_BToKmumu++){

	TLorentzVector kaon_tlv;
	TLorentzVector mu1_tlv;
	TLorentzVector mu2_tlv;
	
	kaon_tlv.SetPtEtaPhiM(tree->BToKmumu_kaon_pt[i_BToKmumu],
			      tree->BToKmumu_kaon_eta[i_BToKmumu],
			      tree->BToKmumu_kaon_phi[i_BToKmumu],
			      KaonMass_);
	mu1_tlv.SetPtEtaPhiM(tree->BToKmumu_mu1_pt[i_BToKmumu],
			      tree->BToKmumu_mu1_eta[i_BToKmumu],
			      tree->BToKmumu_mu1_phi[i_BToKmumu],
			      MuonMass_);
	mu2_tlv.SetPtEtaPhiM(tree->BToKmumu_mu2_pt[i_BToKmumu],
			      tree->BToKmumu_mu2_eta[i_BToKmumu],
			      tree->BToKmumu_mu2_phi[i_BToKmumu],
			      MuonMass_);

	float dR_KFromB = kaon_tlv.DeltaR(gen_KFromB_tlv);
	float dR_mu1FromJPsi = min(mu1_tlv.DeltaR(gen_mu1FromJPsi_tlv),mu2_tlv.DeltaR(gen_mu1FromJPsi_tlv));
	float dR_mu2FromJPsi = min(mu1_tlv.DeltaR(gen_mu2FromJPsi_tlv),mu2_tlv.DeltaR(gen_mu2FromJPsi_tlv));
	//Should check that same objects not selected twice

	float dR_tot = dR_KFromB + dR_mu1FromJPsi + dR_mu2FromJPsi; //In case several BToKmumu matches, take the closest one in dR_tot

	if( dR_KFromB<0.1 && dR_mu1FromJPsi<0.1 && dR_mu2FromJPsi<0.1
	    && (best_dR<0. || dR_tot<best_dR) ){
	  best_dR = dR_tot;
	  _BToKmumu_gen_index = i_BToKmumu;	  
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



