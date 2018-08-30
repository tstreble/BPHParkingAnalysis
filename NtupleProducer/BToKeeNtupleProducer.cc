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
float ElectronMass_ = 0.5109989e-3;

/* possible selctions cut for signal
float min_CL = 0.1;
float min_BcosAlpha = 0.9;
float min_KDCA = 6;
float min_BIP = 1;
*/

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

  bool isResonant = false;
  for (int i = 1; i < argc; ++i) {
    if(std::string(argv[i]) == "Resonant") {
      isResonant = true;
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

  TTree* tree_new=new TTree("BToKeeTree","BToKeeTree");
  if(saveFullNanoAOD)
    tree_new=tree->GetTree()->CloneTree(0);


  //New branches

  int _Muon_sel_index = -1;
  int _BToKee_sel_index = -1;

  tree_new->Branch("Muon_sel_index",&_Muon_sel_index,"Muon_sel_index/I");
  tree_new->Branch("BToKee_sel_index",&_BToKee_sel_index,"BToKee_sel_index/I");

  int _GenPart_BToKee_index = -1;
  int _GenPart_JPsiFromB_index = -1;
  int _GenPart_KFromB_index = -1;
  int _GenPart_e1FromB_index = -1;
  int _GenPart_e2FromB_index = -1;
  int _BToKee_gen_index = -1;
  float _BToKee_gen_eeMass = -1;
  float _BToKee_gen_mass = -1;

  int _Electron_e1FromB_index = -1;
  int _Electron_e2FromB_index = -1;
  int _PFCand_genKFromB_index = -1;

  if(isMC){

    tree_new->Branch("GenPart_BToKee_index",&_GenPart_BToKee_index,"GenPart_BToKee_index/I");
    tree_new->Branch("GenPart_JPsiFromB_index",&_GenPart_JPsiFromB_index,"GenPart_JPsiFromB_index/I");
    tree_new->Branch("GenPart_KFromB_index",&_GenPart_KFromB_index,"GenPart_KFromB_index/I");
    tree_new->Branch("GenPart_e1FromB_index",&_GenPart_e1FromB_index,"GenPart_e1FromB_index/I");
    tree_new->Branch("GenPart_e2FromB_index",&_GenPart_e2FromB_index,"GenPart_e2FromB_index/I");
    tree_new->Branch("BToKee_gen_index",&_BToKee_gen_index,"BToKee_gen_index/I");
    tree_new->Branch("Electron_e1FromB_index",&_Electron_e1FromB_index,"Electron_e1FromB_index/I");
    tree_new->Branch("Electron_e2FromB_index",&_Electron_e2FromB_index,"Electron_e2FromB_index/I");
    tree_new->Branch("PFCand_genKFromB_index",&_PFCand_genKFromB_index,"PFCand_genKFromB_index/I");
    tree_new->Branch("BToKee_gen_eeMass",&_BToKee_gen_eeMass,"BToKee_gen_eeMass/F");
    tree_new->Branch("BToKee_gen_mass",&_BToKee_gen_mass,"BToKee_gen_mass/F");
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

    _Muon_sel_index = -1;
    _BToKee_sel_index = -1;

    _GenPart_BToKee_index = -1;
    _GenPart_JPsiFromB_index = -1;
    _GenPart_KFromB_index = -1;
    _GenPart_e1FromB_index = -1;
    _GenPart_e2FromB_index = -1;
    _BToKee_gen_index = -1;
    _Electron_e1FromB_index = -1;
    _Electron_e2FromB_index = -1;
    _PFCand_genKFromB_index = -1;
    _BToKee_gen_eeMass = -1;
    _BToKee_gen_mass = -1;

    _HLT_Mu8p5_IP3p5 = false;
    _HLT_Mu10p5_IP3p5 = false;
    _HLT_Mu9_IP6 = false;
    _HLT_Mu8_IP3 = false;
    _HLT_BPHParking = false;

    //Trigger selection + matching

    int nMuon = tree->nMuon;

    for(int i_mu=0; i_mu<nMuon; i_mu++){

      //Only trigger for isBPHParking for now
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


    //Select the BToKee candidate with reco criteria

    int nBToKee = tree->nBToKee;
    float best_B_CL_vtx = -1.;

    for(int i_BToKee=0; i_BToKee<nBToKee; i_BToKee++){            

      //Disabled for now
      //if(tree->BToKee_kaon_charge[i_BToKee]*tree->Muon_charge[_Muon_sel_index]>0) continue; //Only consider BToKee with opposite charge to muon
      
      //can be too restrictive for the moment should study effect of charge flip...
      /*
	if(BToKee_ele1_charge[i_BToKee] != Electron_charge[BToKee_ele1_charge[i_BToKee]] ||
	BToKee_ele1_charge[i_BToKee] == -1 || BToKee_ele1_charge[i_BToKee] > nElectron) continue;
	if(BToKee_ele2_charge[i_BToKee] != Electron_charge[BToKee_ele2_charge[i_BToKee]] ||
	BToKee_ele2_charge[i_BToKee] == -1 || BToKee_ele2_charge[i_BToKee] > nElectron) continue;   
      */

      if(tree->BToKee_ele1_charge[i_BToKee] * tree->BToKee_ele2_charge[i_BToKee] > 0.) continue;

      float B_CL_vtx = tree->BToKee_CL_vtx[i_BToKee];
      
      if( best_B_CL_vtx < 0. || B_CL_vtx>best_B_CL_vtx ){      
	best_B_CL_vtx = B_CL_vtx;
	_BToKee_sel_index = i_BToKee;
      }

    }


    if(isData && _BToKee_sel_index<0) continue;

    
    //Select the BToKee candidate based on gen matching

    if(isMC){
      
      int nGenPart = tree->nGenPart;
      
      if(isResonant){

	for(int i_Bu=0; i_Bu<nGenPart; i_Bu++){

	  _GenPart_JPsiFromB_index = -1;
	  _GenPart_e1FromB_index = -1;
	  _GenPart_e2FromB_index = -1;
	  _GenPart_KFromB_index = -1;

	  if(abs(tree->GenPart_pdgId[i_Bu])==521){

	    for(int i_gen=0; i_gen<nGenPart; i_gen++){

	      int pdgId = tree->GenPart_pdgId[i_gen];
	      int mother_index = tree->GenPart_genPartIdxMother[i_gen];
	      int e1_index = -1;
	      int e2_index = -1;

	      if(abs(pdgId)==443 && mother_index == i_Bu){

		for(int j_gen=0; j_gen<nGenPart; j_gen++){
		  int pdgId = tree->GenPart_pdgId[j_gen];
		  int mother_index = tree->GenPart_genPartIdxMother[j_gen];
		  if(abs(pdgId)==11 && mother_index == i_gen && e1_index<0)
		    e1_index = j_gen;
		  else if(abs(pdgId)==11 && mother_index == i_gen)
		    e2_index = j_gen;
		  if(e1_index>=0 && e2_index>=0) break;
		}

		if(e1_index>=0 && e2_index>=0){
		  _GenPart_JPsiFromB_index = i_gen;
		  _GenPart_e1FromB_index = e1_index;
		  _GenPart_e2FromB_index = e2_index;
		}
		else break;

	      }

	      else if(abs(pdgId)==321 && mother_index == i_Bu) _GenPart_KFromB_index = i_gen;

	      else if(mother_index == i_Bu) break; //Additional B decay products

	    }

	  }

	  if(_GenPart_JPsiFromB_index>=0 && _GenPart_KFromB_index>=0){
	    _GenPart_BToKee_index = i_Bu;
	    break;
	  }

	}

      }//resonant

      else{

        for(int i_Bu=0; i_Bu<nGenPart; i_Bu++){

	  _GenPart_JPsiFromB_index = -1;
	  _GenPart_e1FromB_index = -1;
	  _GenPart_e2FromB_index = -1;
	  _GenPart_KFromB_index = -1;

          if(abs(tree->GenPart_pdgId[i_Bu])==521){

            for(int i_gen=0; i_gen<nGenPart; i_gen++){

              int pdgId = tree->GenPart_pdgId[i_gen];
              int mother_index = tree->GenPart_genPartIdxMother[i_gen];
              if(abs(pdgId)==11 && mother_index == i_Bu && _GenPart_e1FromB_index<0)
                _GenPart_e1FromB_index = i_gen;
              else if(abs(pdgId)==11 && mother_index == i_Bu)
                _GenPart_e2FromB_index = i_gen;
              else if(abs(pdgId)==321 && mother_index == i_Bu)
                _GenPart_KFromB_index = i_gen;
	      else if(mother_index == i_Bu) break; //Additional B decay products

	    }
          }//if B

	  if(_GenPart_e1FromB_index>=0 && _GenPart_e2FromB_index>=0 && _GenPart_KFromB_index>=0){
	    _GenPart_BToKee_index = i_Bu;
	    break;
	  }

        }//loop over B
      }


      if(_GenPart_BToKee_index>=0){

	//e1FromB stored a leading daughter
	if(tree->GenPart_pt[_GenPart_e2FromB_index]>tree->GenPart_pt[_GenPart_e1FromB_index]){
	  int i_temp = _GenPart_e1FromB_index;
	  _GenPart_e1FromB_index = _GenPart_e2FromB_index;
	  _GenPart_e2FromB_index = i_temp;
	}

	TLorentzVector gen_KFromB_tlv;
	TLorentzVector gen_e1FromB_tlv;
	TLorentzVector gen_e2FromB_tlv;

	gen_KFromB_tlv.SetPtEtaPhiM(tree->GenPart_pt[_GenPart_KFromB_index],
				    tree->GenPart_eta[_GenPart_KFromB_index],
				    tree->GenPart_phi[_GenPart_KFromB_index],
				    KaonMass_);
	gen_e1FromB_tlv.SetPtEtaPhiM(tree->GenPart_pt[_GenPart_e1FromB_index],
				     tree->GenPart_eta[_GenPart_e1FromB_index],
				     tree->GenPart_phi[_GenPart_e1FromB_index],
				     ElectronMass_);
	gen_e2FromB_tlv.SetPtEtaPhiM(tree->GenPart_pt[_GenPart_e2FromB_index],
				     tree->GenPart_eta[_GenPart_e2FromB_index],
				     tree->GenPart_phi[_GenPart_e2FromB_index],
				     ElectronMass_);

	_BToKee_gen_eeMass = (gen_e1FromB_tlv+gen_e2FromB_tlv).Mag();
	_BToKee_gen_mass = (gen_e1FromB_tlv+gen_e2FromB_tlv+gen_KFromB_tlv).Mag();

	float best_dR = -1.;

	for(int i_BToKee=0; i_BToKee<nBToKee; i_BToKee++){

	  TLorentzVector kaon_tlv;
	  TLorentzVector ele1_tlv;
	  TLorentzVector ele2_tlv;

	  kaon_tlv.SetPtEtaPhiM(tree->BToKee_kaon_pt[i_BToKee],
				tree->BToKee_kaon_eta[i_BToKee],
				tree->BToKee_kaon_phi[i_BToKee],
				KaonMass_);
	  ele1_tlv.SetPtEtaPhiM(tree->BToKee_ele1_pt[i_BToKee],
				tree->BToKee_ele1_eta[i_BToKee],
				tree->BToKee_ele1_phi[i_BToKee],
				ElectronMass_);
	  ele2_tlv.SetPtEtaPhiM(tree->BToKee_ele2_pt[i_BToKee],
				tree->BToKee_ele2_eta[i_BToKee],
				tree->BToKee_ele2_phi[i_BToKee],
				ElectronMass_);

	  float dR_KFromB = kaon_tlv.DeltaR(gen_KFromB_tlv);
	  float dR_e1FromB = min(ele1_tlv.DeltaR(gen_e1FromB_tlv),ele2_tlv.DeltaR(gen_e1FromB_tlv));
	  float dR_e2FromB = min(ele1_tlv.DeltaR(gen_e2FromB_tlv),ele2_tlv.DeltaR(gen_e2FromB_tlv));
	  //Should check that same objects not selected twice

	  float dR_tot = dR_KFromB + dR_e1FromB + dR_e2FromB; //In case several BToKee matches, take the closest one in dR_tot

	  if( dR_KFromB<0.1 && dR_e1FromB<0.1 && dR_e2FromB<0.1
	      && (best_dR<0. || dR_tot<best_dR) ){
	    best_dR = dR_tot;
	    _BToKee_gen_index = i_BToKee;
	  }

	}//matched to reco


	float best_dR_e1FromB = -1.;
	float best_dR_e2FromB = -1.;
	float best_dR_KFromB = -1.;

	int nElectron = tree->nElectron;
	for(int i_ele=0; i_ele<nElectron; i_ele++){

	  TLorentzVector ele_tlv;
	  ele_tlv.SetPtEtaPhiM(tree->Electron_pt[i_ele],
			       tree->Electron_eta[i_ele],
			       tree->Electron_phi[i_ele],
			       ElectronMass_);

	  float dR_e1FromB = ele_tlv.DeltaR(gen_e1FromB_tlv);
	  float dR_e2FromB = ele_tlv.DeltaR(gen_e2FromB_tlv);

	  if(dR_e1FromB<0.1 && (best_dR_e1FromB<0. || dR_e1FromB<best_dR_e1FromB)){
	    _Electron_e1FromB_index = i_ele;
	    best_dR_e1FromB = dR_e1FromB;
	  }
	  if(dR_e2FromB<0.1 && (best_dR_e2FromB<0. || dR_e2FromB<best_dR_e2FromB)){
	    _Electron_e2FromB_index = i_ele;
	    best_dR_e2FromB = dR_e2FromB;
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

      }

    }//gen info


    //Require probe muon passing soft ID for Acc.xEff.
    //+ trigger-matched in BPHParking dataset

    bool isProbeMuonSoftID = false;

    TLorentzVector ele1_tlv;
    ele1_tlv.SetPtEtaPhiM(tree->BToKee_ele1_pt[_BToKee_sel_index],
			  tree->BToKee_ele1_eta[_BToKee_sel_index],
			  tree->BToKee_ele1_phi[_BToKee_sel_index],
			  ElectronMass_);
    TLorentzVector ele2_tlv;
    ele2_tlv.SetPtEtaPhiM(tree->BToKee_ele2_pt[_BToKee_sel_index],
			  tree->BToKee_ele2_eta[_BToKee_sel_index],
			  tree->BToKee_ele2_phi[_BToKee_sel_index],
			  ElectronMass_);


    for(int i_mu=0; i_mu<nMuon; i_mu++){

      TLorentzVector mu;
      mu.SetPtEtaPhiM(tree->Muon_pt[i_mu],tree->Muon_eta[i_mu],tree->Muon_phi[i_mu],tree->Muon_mass[i_mu]);

      //Anti-matching with electrons to be safe
      if(mu.DeltaR(ele1_tlv)<0.1 || mu.DeltaR(ele2_tlv)<0.1) continue;

      if(tree->Muon_softId[i_mu] && tree->Muon_pt[i_mu] > 8. && (!isBPHParking || _Muon_isHLT_BPHParking[i_mu])){
	isProbeMuonSoftID = true;
	_Muon_sel_index = i_mu;
	break;
      }
    }

    if(!isProbeMuonSoftID) continue; //Skip events where there is no probe muon passing the soft ID


    tree_new->Fill();

  }


  f_new->cd();
  if(!saveFullNanoAOD) tree_new->AddFriend("Events",input.c_str());

  tree_new->Write();
  f_new->Close();
  return 0;

}



