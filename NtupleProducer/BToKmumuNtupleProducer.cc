// Author: T. Strebler (IC)
// Date:   31 May 2018
//
// Add new variables from NanoAOD tree for BToKpipi analysis
// Can be either included in a interpreted macro or compiled in c++
// (use `root-config --glibs --cflags`)
//



#include <iostream>
#include <fstream>
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


  float genMuPtCut = -1.;
  for (int i = 1; i < argc; ++i) {
    if(std::string(argv[i]) == "--genMuPtCut") {
      if (i + 1 < argc) {
	genMuPtCut = std::stof(argv[i+1]);
	break;
      } else {
	std::cerr << "--genMuPtCut option requires one argument." << std::endl;
	return 1;
      }
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
  string listFilesTXT;
  bool inputTXT = false;
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
    else if(std::string(argv[i]) == "--inputTXT") {
      if (i + 1 < argc) {
        inputTXT = true;
        listFilesTXT = argv[i+1];
        break;
      } else {
	std::cerr << "--intputTXT option requires one argument." << std::endl;
        return 1;
      }
    }
  }

  if(listFilesTXT == "" && input == ""){
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

  bool isResonant = false;
  for (int i = 1; i < argc; ++i) {
    if(std::string(argv[i]) == "Resonant") {
      isResonant = true;
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
  if(inputTXT){
    string reader;
    std::ifstream inFileLong;
    inFileLong.open(listFilesTXT.c_str(), std::ios::in);

    while(!inFileLong.eof()){
      inFileLong >> reader;
      std::cout << " Adding " << reader << std::endl;
      oldtree->Add(reader.c_str());
    }
  }
  else   oldtree->Add(input.c_str());
  NanoAODTree* tree = new NanoAODTree(oldtree);

  TTree* tree_new=new TTree("BToKmumuTree","BToKmumuTree");
  if(saveFullNanoAOD)
    tree_new=tree->GetTree()->CloneTree(0);


  //New branches

  int _BToKmumu_sel_index = -1;
  int _Muon_sel_index = -1; //Probe muon with selection algo.
  int _Muon_probe_index = -1; //Probe muon for Acc.xEff. = _Muon_sel_index in data

  tree_new->Branch("BToKmumu_sel_index",&_BToKmumu_sel_index,"BToKmumu_sel_index/I");
  tree_new->Branch("Muon_sel_index",&_Muon_sel_index,"Muon_sel_index/I");
  tree_new->Branch("Muon_probe_index",&_Muon_probe_index,"Muon_probe_index/I");


  int _GenPart_BToKmumu_index = -1;
  int _GenPart_JPsiFromB_index = -1;
  int _GenPart_KFromB_index = -1;
  int _GenPart_mu1FromB_index = -1;
  int _GenPart_mu2FromB_index = -1;
  int _BToKmumu_gen_index = -1;
  float _BToKmumu_gen_mumuMass = -1;
  float _BToKmumu_gen_mass = -1;

  int _Muon_mu1FromB_index = -1;
  int _Muon_mu2FromB_index = -1;
  int _PFCand_genKFromB_index = -1;

  if(isMC){

    tree_new->Branch("GenPart_BToKmumu_index",&_GenPart_BToKmumu_index,"GenPart_BToKmumu_index/I");
    tree_new->Branch("GenPart_JPsiFromB_index",&_GenPart_JPsiFromB_index,"GenPart_JPsiFromB_index/I");
    tree_new->Branch("GenPart_KFromB_index",&_GenPart_KFromB_index,"GenPart_KFromB_index/I");
    tree_new->Branch("GenPart_mu1FromB_index",&_GenPart_mu1FromB_index,"GenPart_mu1FromB_index/I");
    tree_new->Branch("GenPart_mu2FromB_index",&_GenPart_mu2FromB_index,"GenPart_mu2FromB_index/I");
    tree_new->Branch("BToKmumu_gen_index",&_BToKmumu_gen_index,"BToKmumu_gen_index/I");
    tree_new->Branch("BToKmumu_gen_mumuMass",&_BToKmumu_gen_mumuMass,"BToKmumu_gen_mumuMass/F");
    tree_new->Branch("BToKmumu_gen_mass",&_BToKmumu_gen_mass,"BToKmumu_gen_mass/F");
    tree_new->Branch("Muon_mu1FromB_index",&_Muon_mu1FromB_index,"Muon_mu1FromB_index/I");
    tree_new->Branch("Muon_mu2FromB_index",&_Muon_mu2FromB_index,"Muon_mu2FromB_index/I");
    tree_new->Branch("PFCand_genKFromB_index",&_PFCand_genKFromB_index,"PFCand_genKFromB_index/I");

  }


  bool _HLT_Mu8p5_IP3p5 = false;
  bool _HLT_Mu10p5_IP3p5 = false;
  bool _HLT_Mu9_IP6 = false;
  bool _HLT_Mu8_IP3 = false;
  bool _HLT_BPHParking = false;

  bool _Muon_isHLT_BPHParking[kMuonMax];

  if(isBPHParking){

    tree_new->Branch("HLT_Mu8p5_IP3p5",&_HLT_Mu8p5_IP3p5,"HLT_Mu8p5_IP3p5/O");
    tree_new->Branch("HLT_Mu10p5_IP3p5",&_HLT_Mu10p5_IP3p5,"HLT_Mu10p5_IP3p5/O");
    tree_new->Branch("HLT_Mu9_IP6",&_HLT_Mu9_IP6,"HLT_Mu9_IP6/O");
    tree_new->Branch("HLT_Mu8_IP3",&_HLT_Mu8_IP3,"HLT_Mu8_IP3/O");
    tree_new->Branch("HLT_BPHParking",&_HLT_BPHParking,"HLT_BPHParking/O");

    tree_new->Branch("Muon_isHLT_BPHParking",_Muon_isHLT_BPHParking,"Muon_isHLT_BPHParking[nMuon]/O");

  }


  int nentries = tree->GetEntries();
  cout<<"Nentries="<<nentries<<endl;
  cout<<"isMC="<<isMC<<endl;

  for (int iEntry = 0; iEntry < nentries ; iEntry++){

    int out = tree->GetEntry(iEntry);
    if(out<0){
      cout<<"Error retrievieng entry #"<<iEntry<<endl;
      return -1;
    }

    if(iEntry%10000==0) cout<<"Entry #"<<iEntry<<" "<< int(100*float(iEntry)/nentries)<<"%"<<endl;

    _BToKmumu_sel_index = -1;
    _Muon_sel_index = -1;
    _Muon_probe_index = -1;

    _GenPart_BToKmumu_index = -1;
    _GenPart_JPsiFromB_index = -1;
    _GenPart_KFromB_index = -1;
    _GenPart_mu1FromB_index = -1;
    _GenPart_mu2FromB_index = -1;
    _BToKmumu_gen_index = -1;
    _BToKmumu_gen_mumuMass = -1;
    _BToKmumu_gen_mass = -1;
    _Muon_mu1FromB_index = -1;
    _Muon_mu2FromB_index = -1;
    _PFCand_genKFromB_index = -1;

    _HLT_Mu8p5_IP3p5 = false;
    _HLT_Mu10p5_IP3p5 = false;
    _HLT_Mu9_IP6 = false;
    _HLT_Mu8_IP3 = false;
    _HLT_BPHParking = false;


    int nMuon = tree->nMuon;

    for(int i_mu=0; i_mu<nMuon; i_mu++){

      //Trigger selection + matching

      //Only for isBPHParking for now
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

	_Muon_isHLT_BPHParking[i_mu] = isTrigMatched;

      }

    }


    //Select the BToKmumu candidate with reco criteria

    int nBToKmumu = tree->nBToKmumu;
    float best_B_CL_vtx = -1.;

    for(int i_BToKmumu=0; i_BToKmumu<nBToKmumu; i_BToKmumu++){            

      //Disabled for now
      //if(tree->BToKmumu_kaon_charge[i_BToKmumu]*tree->Muon_charge[_Muon_sel_index]>0) continue; //Only consider BToKmumu with opposite charge to muon
      
      if(tree->BToKmumu_mu1_charge[i_BToKmumu] * tree->BToKmumu_mu2_charge[i_BToKmumu] > 0.) continue;

      float B_CL_vtx = tree->BToKmumu_CL_vtx[i_BToKmumu];
      
      if( best_B_CL_vtx < 0. || B_CL_vtx>best_B_CL_vtx ){
	best_B_CL_vtx = B_CL_vtx;
	_BToKmumu_sel_index = i_BToKmumu;
      }

    }

    //Take as probe muon leading soft ID muon + trigger-matched if BPHParking data

    if(isData && _BToKmumu_sel_index<0) continue;


    if(_BToKmumu_sel_index>=0){
      
      for(int i_mu=0; i_mu<nMuon; i_mu++){

	if(i_mu==tree->BToKmumu_mu1_index[_BToKmumu_sel_index] || i_mu==tree->BToKmumu_mu2_index[_BToKmumu_sel_index]) continue;

	if(tree->Muon_softId[i_mu] && tree->Muon_pt[i_mu] > 8.
	   && (!isBPHParking || _Muon_isHLT_BPHParking[i_mu])){
	  _Muon_sel_index = i_mu;
	  break;
	}

      }

    }


    //!!! Can have selected B->Kmumu even without additional probe muons
    //Needs to ask both BToKmumu_sel_index>=0 && Muon_sel_index>=0 for analysis on probe side
    //BToKmumu_sel_index not reset to -1 for potential analysis of the probe side


    
    //Select the BToKmumu candidate based on gen matching

    if(isMC){
      
      int nGenPart = tree->nGenPart;
      
      if(isResonant){

	for(int i_Bu=0; i_Bu<nGenPart; i_Bu++){

	  _GenPart_JPsiFromB_index = -1;
	  _GenPart_mu1FromB_index = -1;
	  _GenPart_mu2FromB_index = -1;
	  _GenPart_KFromB_index = -1;

	  if(abs(tree->GenPart_pdgId[i_Bu])==521){

	    for(int i_gen=0; i_gen<nGenPart; i_gen++){

	      int pdgId = tree->GenPart_pdgId[i_gen];
	      int mother_index = tree->GenPart_genPartIdxMother[i_gen];

	      int mu1_index = -1;
	      int mu2_index = -1;

	      if(abs(pdgId)==443 && mother_index == i_Bu){

		for(int j_gen=0; j_gen<nGenPart; j_gen++){

		  int pdgId = tree->GenPart_pdgId[j_gen];
		  int mother_index = tree->GenPart_genPartIdxMother[j_gen];
		  if(abs(pdgId)==13 && mother_index == i_gen && mu1_index<0)
		    mu1_index = j_gen;
		  else if(abs(pdgId)==13 && mother_index == i_gen)
		    mu2_index = j_gen;
		  if(mu1_index>=0 && mu2_index>=0) break;

		}

		if(mu1_index>=0 && mu2_index>=0){
		  _GenPart_JPsiFromB_index = i_gen;
		  _GenPart_mu1FromB_index = mu1_index;
		  _GenPart_mu2FromB_index = mu2_index;
		}
		else break;

	      }

	      else if(abs(pdgId)==321 && mother_index == i_Bu) _GenPart_KFromB_index = i_gen;

	      else if(mother_index == i_Bu) break; //Additional B decay products

	    }

	  }

	  if(_GenPart_JPsiFromB_index>=0 && _GenPart_KFromB_index>=0){
	    _GenPart_BToKmumu_index = i_Bu;
	    break;
	  }

	}	

      }


      else{ //Non-resonant

	for(int i_Bu=0; i_Bu<nGenPart; i_Bu++){

	  _GenPart_JPsiFromB_index = -1;
	  _GenPart_mu1FromB_index = -1;
	  _GenPart_mu2FromB_index = -1;
	  _GenPart_KFromB_index = -1;

	  if(abs(tree->GenPart_pdgId[i_Bu])==521){

	    for(int i_gen=0; i_gen<nGenPart; i_gen++){

	      int pdgId = tree->GenPart_pdgId[i_gen];
	      int mother_index = tree->GenPart_genPartIdxMother[i_gen];
	      if(abs(pdgId)==13 && mother_index == i_Bu && _GenPart_mu1FromB_index<0)
		_GenPart_mu1FromB_index = i_gen;
	      else if(abs(pdgId)==13 && mother_index == i_Bu)
		_GenPart_mu2FromB_index = i_gen;
	      else if(abs(pdgId)==321 && mother_index == i_Bu)
		_GenPart_KFromB_index = i_gen;
	      else if(mother_index == i_Bu) break;

	    }

	  }//if B

	  if(_GenPart_mu1FromB_index>=0 && _GenPart_mu2FromB_index>=0 && _GenPart_KFromB_index>=0){
	    _GenPart_BToKmumu_index = i_Bu;
	    break;
	  }
	  
	} //loop over B

      }


      if(_GenPart_BToKmumu_index>=0){

	//mu1FromB stored a leading daughter
	if(tree->GenPart_pt[_GenPart_mu2FromB_index]>tree->GenPart_pt[_GenPart_mu1FromB_index]){
	  int i_temp = _GenPart_mu1FromB_index;
	  _GenPart_mu1FromB_index = _GenPart_mu2FromB_index;
	  _GenPart_mu2FromB_index = i_temp;
	}

	TLorentzVector gen_KFromB_tlv;
	TLorentzVector gen_mu1FromB_tlv;
	TLorentzVector gen_mu2FromB_tlv;

	gen_KFromB_tlv.SetPtEtaPhiM(tree->GenPart_pt[_GenPart_KFromB_index],
				    tree->GenPart_eta[_GenPart_KFromB_index],
				    tree->GenPart_phi[_GenPart_KFromB_index],
				    KaonMass_);
	gen_mu1FromB_tlv.SetPtEtaPhiM(tree->GenPart_pt[_GenPart_mu1FromB_index],
				      tree->GenPart_eta[_GenPart_mu1FromB_index],
				      tree->GenPart_phi[_GenPart_mu1FromB_index],
				      MuonMass_);
	gen_mu2FromB_tlv.SetPtEtaPhiM(tree->GenPart_pt[_GenPart_mu2FromB_index],
				      tree->GenPart_eta[_GenPart_mu2FromB_index],
				      tree->GenPart_phi[_GenPart_mu2FromB_index],
				      MuonMass_);

	_BToKmumu_gen_mumuMass = (gen_mu1FromB_tlv+gen_mu2FromB_tlv).Mag();
	_BToKmumu_gen_mass = (gen_mu1FromB_tlv+gen_mu2FromB_tlv+gen_KFromB_tlv).Mag();

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
	  float dR_mu1FromB = min(mu1_tlv.DeltaR(gen_mu1FromB_tlv),mu2_tlv.DeltaR(gen_mu1FromB_tlv));
	  float dR_mu2FromB = min(mu1_tlv.DeltaR(gen_mu2FromB_tlv),mu2_tlv.DeltaR(gen_mu2FromB_tlv));
	  //Should check that same objects not selected twice

	  float dR_tot = dR_KFromB + dR_mu1FromB + dR_mu2FromB; //In case several BToKmumu matches, take the closest one in dR_tot

	  if( dR_KFromB<0.1 && dR_mu1FromB<0.1 && dR_mu2FromB<0.1
	      && (best_dR<0. || dR_tot<best_dR) ){
	    best_dR = dR_tot;
	    _BToKmumu_gen_index = i_BToKmumu;
	  }

	}

	float best_dR_mu1FromB = -1.;
	float best_dR_mu2FromB = -1.;
	float best_dR_KFromB = -1.;

	for(int i_mu=0; i_mu<nMuon; i_mu++){

	  TLorentzVector mu_tlv;
	  mu_tlv.SetPtEtaPhiM(tree->Muon_pt[i_mu],
			      tree->Muon_eta[i_mu],
			      tree->Muon_phi[i_mu],
			      MuonMass_);

	  float dR_mu1FromB = mu_tlv.DeltaR(gen_mu1FromB_tlv);
	  float dR_mu2FromB = mu_tlv.DeltaR(gen_mu2FromB_tlv);

	  if(dR_mu1FromB<0.1 && (best_dR_mu1FromB<0. || dR_mu1FromB<best_dR_mu1FromB)){
	    _Muon_mu1FromB_index = i_mu;
	    best_dR_mu1FromB = dR_mu1FromB;
	  }
	  if(dR_mu2FromB<0.1 && (best_dR_mu2FromB<0. || dR_mu2FromB<best_dR_mu2FromB)){
	    _Muon_mu2FromB_index = i_mu;
	    best_dR_mu2FromB = dR_mu2FromB;
	  }

	}

	int nPFCand = tree->nPFCand;
	for(int i_pf=0;i_pf<nPFCand;i_pf++){

	  if(abs(tree->PFCand_pdgId[i_pf])!=211) continue;

	  TLorentzVector PFCand_tlv;
	  PFCand_tlv.SetPtEtaPhiM(tree->PFCand_pt[i_pf],tree->PFCand_eta[i_pf],tree->PFCand_phi[i_pf],tree->PFCand_mass[i_pf]);

	  float dR_KFromB = PFCand_tlv.DeltaR(gen_KFromB_tlv);

	  if(dR_KFromB<0.1 && (best_dR_KFromB<0. || dR_KFromB<best_dR_KFromB)){
	    _PFCand_genKFromB_index = i_pf;
	    best_dR_KFromB = dR_KFromB;
	  }

	}


	//Gen muon pt filter on probe
	
	bool isTagMuonHighPt = false;
	
	for(int i_gen=0; i_gen<nGenPart; i_gen++){
	  
	  if(i_gen==_GenPart_mu1FromB_index || i_gen==_GenPart_mu2FromB_index) continue;
	  
	  int pdgId = tree->GenPart_pdgId[i_gen];
	  if(abs(pdgId)!=13) continue;
	  
	  TLorentzVector gen_tagMu_tlv;
	  gen_tagMu_tlv.SetPtEtaPhiM(tree->GenPart_pt[i_gen],
				     tree->GenPart_eta[i_gen],
				     tree->GenPart_phi[i_gen],
				     MuonMass_);

	  if(gen_tagMu_tlv.DeltaR(gen_mu1FromB_tlv)>0.1 && gen_tagMu_tlv.DeltaR(gen_mu2FromB_tlv)>0.1 //In case there are several copies of same muon (with FSR for instance)
	     && gen_tagMu_tlv.Pt()>genMuPtCut){
	    isTagMuonHighPt = true;
	    break;
	  }

	}

	if(!isTagMuonHighPt) continue; //Skip events where the gen filter in MC has been applied to the probe muon

      }



      //Require additional probe muon passing soft ID for Acc.xEff. symmetric with electrons

      bool isProbeMuonSoftID = false;
      for(int i_mu=0; i_mu<nMuon; i_mu++){

	if(i_mu==_Muon_mu1FromB_index || i_mu==_Muon_mu2FromB_index) continue;
	if(tree->Muon_softId[i_mu] && tree->Muon_pt[i_mu] > 8.){
	  isProbeMuonSoftID = true;
	  _Muon_probe_index = i_mu;
	  break;
	}

      }

      if(!isProbeMuonSoftID) continue; //Skip events where there is no probe muon passing the soft ID

    }


    if(isBPHParking){
      _Muon_probe_index = _Muon_sel_index;
    }


    tree_new->Fill();

  }


  f_new->cd();
  if(!saveFullNanoAOD) tree_new->AddFriend("Events",input.c_str());

  tree_new->Write();
  f_new->Close();
  return 0;

}



