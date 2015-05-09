#define DetectorEffects_cxx
#include "../interface/DetectorEffects.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
// -----------------------------------------------------------------------------------------
// Fro relval TTBar 
// Removed           //if( (DeltaPt_*GENJET_p4.Pt()) > 200.0 && pfMetRawPt > 200.0) {                                                                                    
// Removed genlevel taus cuts 
// Removed genlevel muon cuts
// Add them back when running on DY -> mumu samples. 
// -----------------------------------------------------------------------------------------


void DetectorEffects::Loop(TString outfilename)
{
  using namespace std;
  const double m_pi = 3.1415926535 ;
  Setetaphi();
  OutputFileName = outfilename;
  f = new TFile(OutputFileName,"RECREATE");
  outTree_ = new TTree("outTree_","outTree_");
  
  MakeBranches();
  MakeHistos();
  
  std::cout<<" begin loop"<<std::endl;
  bool debug = false;
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Clear();
    
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    if(debug) std::cout<< "n muons = "<<genjetngenMuons<<std::endl;
    if( genjetngenMuons ==2 ) {
      //if(true){
      
      int ntaus=0;
      int nmuons=0;
      int neles=0;
      for (auto ipart=0; ipart < (*genParId).size(); ipart++){
	if ( (*genParId)[ipart] == 15 || (*genParId)[ipart] == -15 ) 	ntaus++;
	if ( (*genParId)[ipart] == 13 || (*genParId)[ipart] == -13 ) 	nmuons++;
	if ( (*genParId)[ipart] == 11 || (*genParId)[ipart] == -11 ) 	neles++;
      }
      
      if ( ntaus  == 0 ){
	//if(!(ntaus==0 && nmuons==0 && neles==0)){
	//	if(true){
	
	//compute dPt before filling the histograms 
	Float_t tmpsumdPt=0;
	for (auto i = 0;i < AK5nJet; i++){
	  if((*AK5genjetEn)[i]<-900) continue; // don't use jet when there is no matched genjet for a given PFJet
	  RECOJET_p4.SetPtEtaPhiE( (*AK5jetPt)[i],
				   (*AK5jetEta)[i],
				   (*AK5jetPhi)[i],
				   (*AK5jetEn)[i] );
	  GENJET_p4.SetPxPyPzE( (*AK5genjetPx)[i],
				(*AK5genjetPy)[i],
				(*AK5genjetPz)[i],
				(*AK5genjetEn)[i]);
	  
	  Float_t tmpDeltaPt_ = (GENJET_p4.Pt()-RECOJET_p4.Pt());
	  tmpsumdPt =+ (tmpDeltaPt_) ;
	}//compute dPt before filling the histograms ends here 
	
	if(debug) std::cout<<" value of sum dpt = "<<tmpsumdPt<<std::endl;
	
	Float_t HT=0;
	Float_t sumdPt =0.;
	for (auto j = 0;j < AK5nJet; j++){
	  if((*AK5genjetEn)[j]<-900) continue; // don't use jet when there is no matched genjet for a given PFJet
	  TLorentzVector RECOJET_p4;
	  RECOJET_p4.SetPtEtaPhiE( (*AK5jetPt)[j]*(*AK5jecFactor)[j],
				   (*AK5jetEta)[j],
				   (*AK5jetPhi)[j],
				   (*AK5jetEn)[j] );
	  
	  
	  TLorentzVector GENJET_p4;
	  GENJET_p4.SetPxPyPzE( (*AK5genjetPx)[j],
				(*AK5genjetPy)[j],
				(*AK5genjetPz)[j],
				(*AK5genjetEn)[j]);
	  
	  bool isEventFound = FindEvent(GENJET_p4.Eta(), GENJET_p4.Phi());
	  if(!isEventFound) continue; // this will fill only those events which are in vicinity of the ECAL Holes. 
	  
	  //if( RECOJET_p4.Eta() > 2.4) continue;
	  //if( RECOJET_p4.Eta() < -2.4) continue;
	  
	  DeltaPt_ = (GENJET_p4.Pt()-RECOJET_p4.Pt())/GENJET_p4.Pt();
	  
	  HT =+ RECOJET_p4.Pt(); // compute sum pT of all the jets
	  sumdPt =+ (DeltaPt_*GENJET_p4.Pt()) ;
	  
	  //Fill branches to New Tree;
	  dpT_.push_back((double)DeltaPt_*GENJET_p4.Pt()) ;
	  dpT_Over_pT_.push_back((double)DeltaPt_);
	  GenJetpT_.push_back((double)GENJET_p4.Pt());
	  GenJeteta_.push_back((double)GENJET_p4.Eta());
	  GenJetphi_.push_back((double)GENJET_p4.Phi());
	  
	  
	  // This will make sure that only outliers will be filled in the histograms
	  //if( tmpsumdPt > 250.0 && pfMetRawPt > 250.0) {
	  //if( (DeltaPt_*GENJET_p4.Pt()) > 200.0 && pfMetRawPt > 200.0) {
	  delta->Fill(DeltaPt_);
	  dpT->Fill((DeltaPt_*GENJET_p4.Pt()));
	  deltapt_vs_eta->Fill( GENJET_p4.Eta(),DeltaPt_);
	  deltapt_vs_phi->Fill( GENJET_p4.Phi(),DeltaPt_);
	  dpTOverpT_vs_MET->Fill(DeltaPt_,pfMetRawPt);
	  //eta_vs_phi_profile_dPt->Fill(GENJET_p4.Eta(), GENJET_p4.Phi(), (GENJET_p4.Energy() - RECOJET_p4.Energy()));
	  eta_vs_phi_profile_dPt->Fill(GENJET_p4.Eta(), GENJET_p4.Phi(), DeltaPt_);
	  MET_vs_dPt->Fill(pfMetRawPt,(DeltaPt_*GENJET_p4.Pt()));
	  dpT_vs_genJetpT->Fill(DeltaPt_*GENJET_p4.Pt(),GENJET_p4.Pt());
	  dpT_vs_recoJetpT->Fill(DeltaPt_*GENJET_p4.Pt(),RECOJET_p4.Pt());
	  
	  GenHAD_vs_RecoChHAD->Fill((*AK5genjetHAD)[j]/GENJET_p4.Energy(), (*AK5jetCHadEF)[j]);
	  GenEM_vs_RecoEM->Fill((*AK5genjetEM)[j]/GENJET_p4.Energy(), (*AK5jetPhoEF)[j]);
	  GenHAD_vs_RecoHAD->Fill((*AK5genjetHAD)[j]/GENJET_p4.Energy(), (*AK5jetCHadEF)[j] + (*AK5jetNHadEF)[j]);
	  
	  
	  //Dilute with dphi (MET, jet pT) < 0.5
	  double phi1 = GENJET_p4.Phi();
	  double phi2 = pfMetRawPhi;;
	  double result = phi1 - phi2;
	  if(result>m_pi) result -= 2*m_pi;
	  else if (result <= -m_pi) result += 2*m_pi ;
	  else result = result;
	  dPhi_vs_MET->Fill(result, pfMetRawPt);
	  
	  std::cout<<" dphi = "<<result
		   <<" MET = "<<pfMetRawPt
		   <<" Jet pt,eta,phi = ("<<RECOJET_p4.Pt()<<", "<<RECOJET_p4.Eta()<<", "<<RECOJET_p4.Phi()<<") "
		   <<" genJet pt,eta,phi ("<<GENJET_p4.Pt()<<", "<<GENJET_p4.Eta()<<", "<<GENJET_p4.Phi()<<") "
		   <<" dpT "<<DeltaPt_*GENJET_p4.Pt()
		   <<std::endl;
	  
	  if(result <0.5 && result > -0.5){
	    eta_vs_phi_profile_dPt_Dilute_dphi->Fill(GENJET_p4.Eta(), GENJET_p4.Phi(), DeltaPt_);
	    MET_vs_dPt_Dilute_dphi->Fill(pfMetRawPt,(DeltaPt_*GENJET_p4.Pt()));
	    GenEM_vs_RecoEM_Dilute_dPhi->Fill((*AK5genjetEM)[j]/GENJET_p4.Energy(), (*AK5jetPhoEF)[j]);
	    GenHAD_vs_RecoHAD_Dilute_dPhi->Fill((*AK5genjetHAD)[j]/GENJET_p4.Energy(), (*AK5jetCHadEF)[j] + (*AK5jetNHadEF)[j]);
	    GenHAD_vs_RecoChHAD_Dilute_dPhi->Fill((*AK5genjetHAD)[j]/GENJET_p4.Energy(), (*AK5jetCHadEF)[j]);
	  }
	  
	  // Dilute with dpT > 200 GeV && MET > 200 GeV
	  if( (DeltaPt_*GENJET_p4.Pt()) > 200.0 && pfMetRawPt > 200.0) { 
	    //  if( (DeltaPt_*GENJET_p4.Pt()) > 50.0) { 
	    eta_vs_phi_profile_dPt_DilutedpT->Fill(GENJET_p4.Eta(), GENJET_p4.Phi(), DeltaPt_);
	    dPhi_MET_Jet->Fill(result);
	    dPhi_vs_MET_DilutedpT->Fill(result, pfMetRawPt);
	    double dpt = DeltaPt_*GENJET_p4.Pt();
	    dPhi_vs_METOverdpT_DilutedpT->Fill(result,pfMetRawPt/dpt);
	    GenEM_vs_RecoEM_Dilute_dpT->Fill((*AK5genjetEM)[j]/GENJET_p4.Energy(), (*AK5jetPhoEF)[j]);
	    GenHAD_vs_RecoHAD_Dilute_dpT->Fill((*AK5genjetHAD)[j]/GENJET_p4.Energy(), (*AK5jetCHadEF)[j] + (*AK5jetNHadEF)[j]);
	    GenHAD_vs_RecoChHAD_Dilute_dpT->Fill((*AK5genjetHAD)[j]/GENJET_p4.Energy(), (*AK5jetCHadEF)[j]);
	    if(true)	      std::cout<<info_runId
				       <<":"<<info_lumiSection
				       <<":"<<info_eventId
				       <<std::endl;
	    
	  }
	  
	  
	  
	}//for (auto j = 0;j < AK5nJet; j++){
	
       	MET_vs_HT->Fill(pfMetRawPt,HT);
	MET_vs_SumdPt->Fill(pfMetRawPt,sumdPt);
	SumdpT->Fill(sumdPt);
	MET->Fill(pfMetRawPt);
		
	// Fill New Branches;
	MET_ = (float) pfMetRawPt;
	HT_ = (float) HT;
	SumdpT_ = (float) sumdPt;
	
	Run_ = info_runId;
	Lumi_ = info_lumiSection;
	Event_ = info_eventId;
	outTree_->Fill(); 
      }// ntaus > 0
      }    // exact two gen muons
      
    if(false){
      for(auto imu=0;imu<nMu;imu++){
	std::cout<<" pT, eta,phi = ("<<(*muPt)[imu]<<", "<<(*muEta)[imu]<<", "<<(*muPhi)[imu]<<")"
		 <<" GM,TM = ("<<(*isGlobalMuon)[imu]<<", "<<(*isTrackerMuon)[imu]<<")"
		 <<"muTrkLayers = "<<(*muTrkLayers)[imu]
		 <<"muPixelHits = "<<(*muPixelHits)[imu]
		 <<"muHits = "<<(*muHits)[imu]
		 <<"muMatches = "<<(*muMatches)[imu]
		 <<"muSegID = "<<(*muSegID)[imu]
		 <<"muNSegs = "<<(*muNSegs)[imu]
		 <<std::endl;
      }
    }
      //electrons.pVector = (TLorentzVector*) RK_Electron_4Momentum->At(i);
      //if(jentry%100000==1) 
    std::cout<<"---------------- event number ---------------- = "<<jentry<<std::endl;
    }
    
  
  
  f->cd();
  //outTree_->Write();
  
  delta->Write();
  MET->Write();
  SumdpT->Write();
  
  dpT->Write();
  deltapt_vs_eta->Write();
  deltapt_vs_phi->Write();
  dpTOverpT_vs_MET->Write();
  eta_vs_phi_profile_dPt->Write();
  MET_vs_HT->Write();
  MET_vs_SumdPt->Write();
  MET_vs_dPt->Write();
  eta_vs_phi_profile_dPt_DilutedpT->Write();
  eta_vs_phi_profile_dPt_Dilute_dphi->Write();
  MET_vs_dPt_Dilute_dphi->Write();
  dpT_vs_genJetpT->Write();
  dpT_vs_recoJetpT->Write();
  GenEM_vs_RecoEM->Write();
  GenHAD_vs_RecoChHAD->Write();
  GenHAD_vs_RecoHAD->Write();
  dPhi_MET_Jet->Write();
  
  dPhi_vs_MET_DilutedpT->Write();
  dPhi_vs_METOverdpT_DilutedpT->Write();
  

  dPhi_vs_dpT->Write();
  dPhi_vs_MET->Write();
  GenEM_vs_RecoEM_Dilute_dPhi->Write();
  GenHAD_vs_RecoHAD_Dilute_dPhi->Write();
  GenHAD_vs_RecoChHAD_Dilute_dPhi->Write();
  GenEM_vs_RecoEM_Dilute_dpT->Write();
  GenHAD_vs_RecoHAD_Dilute_dpT->Write();
  GenHAD_vs_RecoChHAD_Dilute_dpT->Write();

  f->Close();
  std::cout<<" file closed and job finished"<<std::endl;
}




void DetectorEffects::MakeBranches(){
  
  outTree_->Branch("MET_",&MET_,"MET_/F");
  outTree_->Branch("HT_",&HT_,"HT_/F");
  outTree_->Branch("SumdpT_",&SumdpT_,"SumdpT_/F");
  
  outTree_->Branch("dpT_","std::vector<double>",&dpT_);
  outTree_->Branch("dpT_Over_pT_","std::vector<double>",&dpT_Over_pT_);
  outTree_->Branch("GenJetpT_","std::vector<double>",&GenJetpT_);
  outTree_->Branch("GenJeteta_","std::vector<double>",&GenJeteta_);
  outTree_->Branch("GenJetphi_","std::vector<double>",&GenJetphi_);
  outTree_->Branch("Run_",&Run_,"Run_/I");
  outTree_->Branch("Lumi_",&Lumi_,"Lumi_/I");
  outTree_->Branch("Event_",&Event_,"Event_/I");
}


void DetectorEffects::MakeHistos(){
  delta = new TH1F("delta","delta;(p_{T}^{genJet}-p_{T}^{recoJet})/p_{T}^{genJet};# of entries",120,-3,3.);
  MET   = new TH1F("MET","MET;MET;# of events", 750, 0,1500);
  SumdpT = new TH1F("SumdpT","SumdpT;#Sigma #Delta p_{T};# of Events",500,-1000,1000);
  dpT = new TH1F("dpT","dpT; #Delta p_{T};# of Events",500,-1000,1000);
  dPhi_MET_Jet = new TH1F("dPhi_MET_Jet","dPhi_MET_Jet;#Delta #phi (MET,jet);# of Events", 100,-20,20);
  
  dpTOverpT_vs_MET = new TH2F("dpTOverpT_vs_MET","dpTOverpT_vs_MET;#Delta p_{T}/p_{T};MET",90,-3.0,3.0,500,0,1000);
  
  deltapt_vs_eta = new TH2F("deltapt_vs_eta","deltapt_vs_eta; #eta_{genJet} ; #Delta p_{T}/p_{T}^{genJet}",100,-3.14,3.14, 400,-10.,10.);
  deltapt_vs_phi = new TH2F("deltapt_vs_phi","deltapt_vs_phi; #phi_{genJet}; #Delta p_{T}/p_{T}^{genJet} ",100,-3.14,3.14, 400,-10.,10.);
  eta_vs_phi_profile_dPt = new TProfile2D("eta_vs_phi_profile_dPt","eta_vs_phi_profile_dPt;#eta^{genJet};#phi_{genJet}",100,-3.14,3.14, 100,-3.14,3.14); 
  MET_vs_HT = new TH2F("MET_vs_HT","MET_vs_HT;MET; H_{T}",200,0,1000,200,0,1000);
  MET_vs_SumdPt = new TH2F("MET_vs_SumdPt","MET_vs_SumdPt;MET; #Sigma #Delta p_{T}",400,-1000,1000,400,-1000,1000);  
  MET_vs_dPt    = new TH2F("MET_vs_dPt","MET_vs_dPt;MET;#Delta p_{T}",200,0,1000,200,0,1000);
  
  dpT_vs_genJetpT = new TH2F("dpT_vs_genJetpT","dpT_vs_genJetpT;#Delta p_{T}; p_{T}^{genJet}",1000,-1000,1000, 1000,-1000,1000);
  dpT_vs_recoJetpT = new TH2F("dpT_vs_recoJetpT","dpT_vs_recoJetpT;#Delta p_{T}; p_{T}^{recoJet}",1000,-1000,1000, 1000,-1000,1000);
  
  eta_vs_phi_profile_dPt_DilutedpT = new TProfile2D("eta_vs_phi_profile_dPt_DilutedpT","eta_vs_phi_profile_dPt_DilutedpT;#eta^{genJet};#phi_{genJet}",100,-3.14,3.14, 100,-3.14,3.14); 
  eta_vs_phi_profile_dPt_Dilute_dphi = new TProfile2D("eta_vs_phi_profile_dPt_Dilute_dphi","eta_vs_phi_profile_dPt_Dilute_dphi;#eta^{genJet};#phi_{genJet}",100,-3.14,3.14, 100,-3.14,3.14); 
  MET_vs_dPt_Dilute_dphi    = new TH2F("MET_vs_dPt_Dilute_dphi","MET_vs_dPt_Dilute_dphi;MET;#Delta p_{T}",200,0,1000,200,0,1000);
  
  GenEM_vs_RecoEM = new TH2F("GenEM_vs_RecoEM","GenEM_vs_RecoEM;EM_{gen};EM_{reco}",100,0,1,100,0,1);
  GenHAD_vs_RecoChHAD = new TH2F("GenHAD_vs_RecoChHAD","GenHAD_vs_RecoChHAD;Had_{gen};ChHad_{reco}",100,0,1,100,0,1);
  GenHAD_vs_RecoHAD = new TH2F("GenHAD_vs_RecoHAD","GenHAD_vs_RecoHAD;Had_{gen};Had_{reco}",100,0,1,100,0,1);
  
  dPhi_vs_MET = new TH2F("dPhi_vs_MET","dPhi_vs_MET;#Delta #phi; MET",100,-20,20,500,0,1000);
  dPhi_vs_dpT = new TH2F("dPhi_vs_dpT","dPhi_vs_dpT;#Delta #phi;#Delta p_{T}",100,-20,20,200,0,1000);
  

  dPhi_vs_MET_DilutedpT = new TH2F("dPhi_vs_MET_DilutedpT","dPhi_vs_MET_DilutedpT;#Delta #phi; MET",100,-20,20,500,0,1000);

  dPhi_vs_METOverdpT_DilutedpT = new TH2F("dPhi_vs_METOverdpT_DilutedpT","dPhi_vs_METOverdpT_DilutedpT;#Delta #phi; MET/#Delta p_{T}",100,-20,20,100,-1,1);
  
  GenEM_vs_RecoEM_Dilute_dPhi = new TH2F("GenEM_vs_RecoEM_Dilute_dPhi","GenEM_vs_RecoEM_Dilute_dPhi;EM_{gen};EM_{reco}",100,0,1,100,0,1);
  GenHAD_vs_RecoHAD_Dilute_dPhi = new TH2F("GenHAD_vs_RecoHAD_Dilute_dPhi","GenHAD_vs_RecoHAD_Dilute_dPhi;Had_{gen};Had_{reco}",100,0,1,100,0,1);
  GenHAD_vs_RecoChHAD_Dilute_dPhi = new TH2F("GenHAD_vs_RecoChHAD_Dilute_dPhi","GenHAD_vs_RecoChHAD_Dilute_dPhi;Had_{gen};Had_{reco}",100,0,1,100,0,1);

  GenEM_vs_RecoEM_Dilute_dpT = new TH2F("GenEM_vs_RecoEM_Dilute_dpT","GenEM_vs_RecoEM_Dilute_dpT;EM_{gen};EM_{reco}",100,0,1,100,0,1);
  GenHAD_vs_RecoHAD_Dilute_dpT = new TH2F("GenHAD_vs_RecoHAD_Dilute_dpT","GenHAD_vs_RecoHAD_Dilute_dpT;Had_{gen};Had_{reco}",100,0,1,100,0,1);
  GenHAD_vs_RecoChHAD_Dilute_dpT = new TH2F("GenHAD_vs_RecoChHAD_Dilute_dpT","GenHAD_vs_RecoChHAD_Dilute_dpT;Had_{gen};Had_{reco}",100,0,1,100,0,1);

}



void DetectorEffects::Clear(){
  MET_ = 0.;
  HT_ = 0.;
  SumdpT_ = 0.;
  dpT_.clear();
  dpT_Over_pT_.clear();
  GenJetpT_.clear();
  GenJeteta_.clear();
  GenJetphi_.clear();
  Run_ = 0;
  Lumi_ = 0;
  Event_ = 0;
}


/*
bool DetectorEffects::FindEvent(Float_t eta1, Float_t phi1){
  ifstream badchannels;
  std::string line;
  float eta;
  float phi;
  float deta;
  float dphi;
  bool result = false;
  badchannels.open ("../../src/data/ECALHoles.txt",ios::in);
  
  std::cout<<" file is open "<<badchannels.is_open()<<std::endl;
  
  if(badchannels.is_open()){
    while (getline(badchannels,line)){
      
      std::istringstream  issbadchannels(line);
      issbadchannels >> eta;
      issbadchannels >> phi;
      deta = eta-eta1;
      dphi = phi-phi1;
      float dr = TMath::Sqrt(deta*deta + dphi*dphi);
      if(dr<0.15) {
	std::cout<<" eta, phi,dr = "<<eta1<<", "<<phi1
		 <<", "<<eta<<", "<<phi
		 <<", "<<dr<<std::endl;
	result = true;
	break;
      }
    }  
  }
  return result;
}
*/


bool DetectorEffects::FindEvent(Float_t eta1, Float_t phi1){
  float eta;
  float phi;
  float deta;
  float dphi;
  bool result = false;
  
  for (int i =0; i<(int)etaVec.size();i++){
    eta = etaVec[i];  phi = phiVec[i];
    deta = eta-eta1;
    dphi = phi-phi1;
    float dr = TMath::Sqrt(deta*deta + dphi*dphi);
    if(dr<0.15) {
      std::cout<<" eta, phi,dr = "<<eta1<<", "<<phi1
	       <<", "<<eta<<", "<<phi
	       <<", "<<dr<<std::endl;
      result = true;
      break;
    }
  }  
  return result;
}


void DetectorEffects::Setetaphi(){
  etaVec.clear(); phiVec.clear();
etaVec.push_back(-1.470);  phiVec.push_back( +1.841);
etaVec.push_back(-1.470);  phiVec.push_back( +2.190);
etaVec.push_back(-1.470);  phiVec.push_back( +2.295);
etaVec.push_back(-1.470);  phiVec.push_back( +2.679);
etaVec.push_back(-1.470);  phiVec.push_back( +3.081);
etaVec.push_back(-1.470);  phiVec.push_back( +3.133);
etaVec.push_back(-1.470);  phiVec.push_back( -3.133);
etaVec.push_back(-1.470);  phiVec.push_back( -3.115);
etaVec.push_back(-1.470);  phiVec.push_back( -3.098);
etaVec.push_back(-1.470);  phiVec.push_back( -3.081);
etaVec.push_back(-1.470);  phiVec.push_back( -3.063);
etaVec.push_back(-1.470);  phiVec.push_back( -1.649);
etaVec.push_back(-1.470);  phiVec.push_back( -1.632);
etaVec.push_back(-1.470);  phiVec.push_back( -1.614);
etaVec.push_back(-1.453);  phiVec.push_back( +0.777);
etaVec.push_back(-1.453);  phiVec.push_back( +1.300);
etaVec.push_back(-1.453);  phiVec.push_back( -3.133);
etaVec.push_back(-1.453);  phiVec.push_back( -3.115);
etaVec.push_back(-1.453);  phiVec.push_back( -3.098);
etaVec.push_back(-1.453);  phiVec.push_back( -3.081);
etaVec.push_back(-1.453);  phiVec.push_back( -3.063);
etaVec.push_back(-1.436);  phiVec.push_back( -3.133);
etaVec.push_back(-1.436);  phiVec.push_back( -3.115);
etaVec.push_back(-1.436);  phiVec.push_back( -3.098);
etaVec.push_back(-1.436);  phiVec.push_back( -3.081);
etaVec.push_back(-1.436);  phiVec.push_back( -3.063);
etaVec.push_back(-1.436);  phiVec.push_back( -1.580);
etaVec.push_back(-1.418);  phiVec.push_back( -3.133);
etaVec.push_back(-1.418);  phiVec.push_back( -3.115);
etaVec.push_back(-1.418);  phiVec.push_back( -3.098);
etaVec.push_back(-1.418);  phiVec.push_back( -3.081);
etaVec.push_back(-1.418);  phiVec.push_back( -3.063);
etaVec.push_back(-1.418);  phiVec.push_back( -1.597);
etaVec.push_back(-1.401);  phiVec.push_back( -3.133);
etaVec.push_back(-1.401);  phiVec.push_back( -3.115);
etaVec.push_back(-1.401);  phiVec.push_back( -3.098);
etaVec.push_back(-1.401);  phiVec.push_back( -3.081);
etaVec.push_back(-1.401);  phiVec.push_back( -3.063);
etaVec.push_back(-1.401);  phiVec.push_back( -1.580);
etaVec.push_back(-1.383);  phiVec.push_back( +1.510);
etaVec.push_back(-1.314);  phiVec.push_back( -0.166);
etaVec.push_back(-1.314);  phiVec.push_back( -0.148);
etaVec.push_back(-1.314);  phiVec.push_back( -0.131);
etaVec.push_back(-1.314);  phiVec.push_back( -0.113);
etaVec.push_back(-1.314);  phiVec.push_back( -0.096);
etaVec.push_back(-1.296);  phiVec.push_back( -2.697);
etaVec.push_back(-1.296);  phiVec.push_back( -2.679);
etaVec.push_back(-1.296);  phiVec.push_back( -2.662);
etaVec.push_back(-1.296);  phiVec.push_back( -2.644);
etaVec.push_back(-1.296);  phiVec.push_back( -2.627);
etaVec.push_back(-1.279);  phiVec.push_back( -2.697);
etaVec.push_back(-1.279);  phiVec.push_back( -2.679);
etaVec.push_back(-1.279);  phiVec.push_back( -2.662);
etaVec.push_back(-1.279);  phiVec.push_back( -2.644);
etaVec.push_back(-1.279);  phiVec.push_back( -2.627);
etaVec.push_back(-1.279);  phiVec.push_back( -1.161);
etaVec.push_back(-1.262);  phiVec.push_back( -2.697);
etaVec.push_back(-1.262);  phiVec.push_back( -2.679);
etaVec.push_back(-1.262);  phiVec.push_back( -2.662);
etaVec.push_back(-1.262);  phiVec.push_back( -2.644);
etaVec.push_back(-1.262);  phiVec.push_back( -2.627);
etaVec.push_back(-1.209);  phiVec.push_back( +2.923);
etaVec.push_back(-1.209);  phiVec.push_back( -1.213);
etaVec.push_back(-1.209);  phiVec.push_back( -1.196);
etaVec.push_back(-1.209);  phiVec.push_back( -1.178);
etaVec.push_back(-1.209);  phiVec.push_back( -1.161);
etaVec.push_back(-1.209);  phiVec.push_back( -1.143);
etaVec.push_back(-1.192);  phiVec.push_back( -1.213);
etaVec.push_back(-1.192);  phiVec.push_back( -1.196);
etaVec.push_back(-1.192);  phiVec.push_back( -1.178);
etaVec.push_back(-1.192);  phiVec.push_back( -1.161);
etaVec.push_back(-1.192);  phiVec.push_back( -1.143);
etaVec.push_back(-1.174);  phiVec.push_back( -1.213);
etaVec.push_back(-1.174);  phiVec.push_back( -1.196);
etaVec.push_back(-1.174);  phiVec.push_back( -1.178);
etaVec.push_back(-1.174);  phiVec.push_back( -1.161);
etaVec.push_back(-1.174);  phiVec.push_back( -1.143);
etaVec.push_back(-1.157);  phiVec.push_back( -1.213);
etaVec.push_back(-1.157);  phiVec.push_back( -1.196);
etaVec.push_back(-1.157);  phiVec.push_back( -1.178);
etaVec.push_back(-1.157);  phiVec.push_back( -1.161);
etaVec.push_back(-1.157);  phiVec.push_back( -1.143);
etaVec.push_back(-1.140);  phiVec.push_back( +2.138);
etaVec.push_back(-1.140);  phiVec.push_back( -1.213);
etaVec.push_back(-1.140);  phiVec.push_back( -1.196);
etaVec.push_back(-1.140);  phiVec.push_back( -1.178);
etaVec.push_back(-1.140);  phiVec.push_back( -1.161);
etaVec.push_back(-1.140);  phiVec.push_back( -1.143);
etaVec.push_back(-1.087);  phiVec.push_back( +3.115);
etaVec.push_back(-1.070);  phiVec.push_back( -1.161);
etaVec.push_back(-1.000);  phiVec.push_back( -0.340);
etaVec.push_back(-0.983);  phiVec.push_back( +2.295);
etaVec.push_back(-0.966);  phiVec.push_back( +1.562);
etaVec.push_back(-0.948);  phiVec.push_back( +1.614);
etaVec.push_back(-0.948);  phiVec.push_back( +2.016);
etaVec.push_back(-0.948);  phiVec.push_back( +2.033);
etaVec.push_back(-0.948);  phiVec.push_back( +2.051);
etaVec.push_back(-0.948);  phiVec.push_back( +2.068);
etaVec.push_back(-0.948);  phiVec.push_back( +2.086);
etaVec.push_back(-0.948);  phiVec.push_back( +2.976);
etaVec.push_back(-0.948);  phiVec.push_back( +2.993);
etaVec.push_back(-0.948);  phiVec.push_back( +3.011);
etaVec.push_back(-0.948);  phiVec.push_back( +3.028);
etaVec.push_back(-0.948);  phiVec.push_back( +3.046);
etaVec.push_back(-0.931);  phiVec.push_back( +2.016);
etaVec.push_back(-0.931);  phiVec.push_back( +2.033);
etaVec.push_back(-0.931);  phiVec.push_back( +2.051);
etaVec.push_back(-0.931);  phiVec.push_back( +2.068);
etaVec.push_back(-0.931);  phiVec.push_back( +2.086);
etaVec.push_back(-0.931);  phiVec.push_back( +2.976);
etaVec.push_back(-0.931);  phiVec.push_back( +2.993);
etaVec.push_back(-0.931);  phiVec.push_back( +3.011);
etaVec.push_back(-0.931);  phiVec.push_back( +3.028);
etaVec.push_back(-0.931);  phiVec.push_back( +3.046);
etaVec.push_back(-0.914);  phiVec.push_back( +2.016);
etaVec.push_back(-0.914);  phiVec.push_back( +2.033);
etaVec.push_back(-0.914);  phiVec.push_back( +2.051);
etaVec.push_back(-0.914);  phiVec.push_back( +2.068);
etaVec.push_back(-0.914);  phiVec.push_back( +2.086);
etaVec.push_back(-0.914);  phiVec.push_back( +2.976);
etaVec.push_back(-0.914);  phiVec.push_back( +2.993);
etaVec.push_back(-0.914);  phiVec.push_back( +3.011);
etaVec.push_back(-0.914);  phiVec.push_back( +3.028);
etaVec.push_back(-0.914);  phiVec.push_back( +3.046);
etaVec.push_back(-0.896);  phiVec.push_back( +2.016);
etaVec.push_back(-0.896);  phiVec.push_back( +2.033);
etaVec.push_back(-0.896);  phiVec.push_back( +2.051);
etaVec.push_back(-0.896);  phiVec.push_back( +2.068);
etaVec.push_back(-0.896);  phiVec.push_back( +2.086);
etaVec.push_back(-0.896);  phiVec.push_back( +2.976);
etaVec.push_back(-0.896);  phiVec.push_back( +2.993);
etaVec.push_back(-0.896);  phiVec.push_back( +3.011);
etaVec.push_back(-0.896);  phiVec.push_back( +3.028);
etaVec.push_back(-0.896);  phiVec.push_back( +3.046);
etaVec.push_back(-0.879);  phiVec.push_back( +0.113);
etaVec.push_back(-0.879);  phiVec.push_back( +2.016);
etaVec.push_back(-0.879);  phiVec.push_back( +2.033);
etaVec.push_back(-0.879);  phiVec.push_back( +2.051);
etaVec.push_back(-0.879);  phiVec.push_back( +2.068);
etaVec.push_back(-0.879);  phiVec.push_back( +2.086);
etaVec.push_back(-0.879);  phiVec.push_back( +2.976);
etaVec.push_back(-0.879);  phiVec.push_back( +2.993);
etaVec.push_back(-0.879);  phiVec.push_back( +3.011);
etaVec.push_back(-0.879);  phiVec.push_back( +3.028);
etaVec.push_back(-0.879);  phiVec.push_back( +3.046);
etaVec.push_back(-0.844);  phiVec.push_back( +0.951);
etaVec.push_back(-0.826);  phiVec.push_back( +2.714);
etaVec.push_back(-0.774);  phiVec.push_back( +2.068);
etaVec.push_back(-0.774);  phiVec.push_back( -2.278);
etaVec.push_back(-0.774);  phiVec.push_back( -0.777);
etaVec.push_back(-0.774);  phiVec.push_back( -0.759);
etaVec.push_back(-0.774);  phiVec.push_back( -0.742);
etaVec.push_back(-0.774);  phiVec.push_back( -0.724);
etaVec.push_back(-0.774);  phiVec.push_back( -0.707);
etaVec.push_back(-0.757);  phiVec.push_back( -0.777);
etaVec.push_back(-0.757);  phiVec.push_back( -0.759);
etaVec.push_back(-0.757);  phiVec.push_back( -0.742);
etaVec.push_back(-0.757);  phiVec.push_back( -0.724);
etaVec.push_back(-0.757);  phiVec.push_back( -0.707);
etaVec.push_back(-0.739);  phiVec.push_back( +1.108);
etaVec.push_back(-0.739);  phiVec.push_back( -0.777);
etaVec.push_back(-0.739);  phiVec.push_back( -0.759);
etaVec.push_back(-0.739);  phiVec.push_back( -0.742);
etaVec.push_back(-0.739);  phiVec.push_back( -0.724);
etaVec.push_back(-0.739);  phiVec.push_back( -0.707);
etaVec.push_back(-0.722);  phiVec.push_back( -0.777);
etaVec.push_back(-0.722);  phiVec.push_back( -0.759);
etaVec.push_back(-0.722);  phiVec.push_back( -0.742);
etaVec.push_back(-0.722);  phiVec.push_back( -0.724);
etaVec.push_back(-0.722);  phiVec.push_back( -0.707);
etaVec.push_back(-0.705);  phiVec.push_back( +2.662);
etaVec.push_back(-0.705);  phiVec.push_back( -0.777);
etaVec.push_back(-0.705);  phiVec.push_back( -0.759);
etaVec.push_back(-0.705);  phiVec.push_back( -0.742);
etaVec.push_back(-0.705);  phiVec.push_back( -0.724);
etaVec.push_back(-0.705);  phiVec.push_back( -0.707);
etaVec.push_back(-0.652);  phiVec.push_back( +3.098);
etaVec.push_back(-0.652);  phiVec.push_back( -0.602);
etaVec.push_back(-0.583);  phiVec.push_back( +0.951);
etaVec.push_back(-0.583);  phiVec.push_back( +1.929);
etaVec.push_back(-0.513);  phiVec.push_back( +2.714);
etaVec.push_back(-0.513);  phiVec.push_back( +2.731);
etaVec.push_back(-0.513);  phiVec.push_back( +2.749);
etaVec.push_back(-0.513);  phiVec.push_back( +2.766);
etaVec.push_back(-0.513);  phiVec.push_back( +2.784);
etaVec.push_back(-0.513);  phiVec.push_back( -2.260);
etaVec.push_back(-0.513);  phiVec.push_back( -2.190);
etaVec.push_back(-0.513);  phiVec.push_back( -1.038);
etaVec.push_back(-0.496);  phiVec.push_back( +2.714);
etaVec.push_back(-0.496);  phiVec.push_back( +2.731);
etaVec.push_back(-0.496);  phiVec.push_back( +2.749);
etaVec.push_back(-0.496);  phiVec.push_back( +2.766);
etaVec.push_back(-0.496);  phiVec.push_back( +2.784);
etaVec.push_back(-0.479);  phiVec.push_back( +1.318);
etaVec.push_back(-0.479);  phiVec.push_back( +2.714);
etaVec.push_back(-0.479);  phiVec.push_back( +2.731);
etaVec.push_back(-0.479);  phiVec.push_back( +2.749);
etaVec.push_back(-0.479);  phiVec.push_back( +2.766);
etaVec.push_back(-0.479);  phiVec.push_back( +2.784);
etaVec.push_back(-0.461);  phiVec.push_back( +2.714);
etaVec.push_back(-0.461);  phiVec.push_back( +2.731);
etaVec.push_back(-0.461);  phiVec.push_back( +2.749);
etaVec.push_back(-0.461);  phiVec.push_back( +2.766);
etaVec.push_back(-0.461);  phiVec.push_back( +2.784);
etaVec.push_back(-0.444);  phiVec.push_back( +2.714);
etaVec.push_back(-0.444);  phiVec.push_back( +2.731);
etaVec.push_back(-0.444);  phiVec.push_back( +2.749);
etaVec.push_back(-0.444);  phiVec.push_back( +2.766);
etaVec.push_back(-0.444);  phiVec.push_back( +2.784);
etaVec.push_back(-0.426);  phiVec.push_back( +2.278);
etaVec.push_back(-0.426);  phiVec.push_back( -1.946);
etaVec.push_back(-0.409);  phiVec.push_back( +2.033);
etaVec.push_back(-0.409);  phiVec.push_back( +2.958);
etaVec.push_back(-0.391);  phiVec.push_back( +1.580);
etaVec.push_back(-0.391);  phiVec.push_back( +2.278);
etaVec.push_back(-0.357);  phiVec.push_back( +3.081);
etaVec.push_back(-0.357);  phiVec.push_back( -2.155);
etaVec.push_back(-0.339);  phiVec.push_back( +0.183);
etaVec.push_back(-0.339);  phiVec.push_back( +0.201);
etaVec.push_back(-0.339);  phiVec.push_back( +0.218);
etaVec.push_back(-0.339);  phiVec.push_back( +0.236);
etaVec.push_back(-0.339);  phiVec.push_back( +0.253);
etaVec.push_back(-0.322);  phiVec.push_back( +0.183);
etaVec.push_back(-0.322);  phiVec.push_back( +0.201);
etaVec.push_back(-0.322);  phiVec.push_back( +0.218);
etaVec.push_back(-0.322);  phiVec.push_back( +0.236);
etaVec.push_back(-0.322);  phiVec.push_back( +0.253);
etaVec.push_back(-0.322);  phiVec.push_back( -1.213);
etaVec.push_back(-0.305);  phiVec.push_back( +0.183);
etaVec.push_back(-0.305);  phiVec.push_back( +0.201);
etaVec.push_back(-0.305);  phiVec.push_back( +0.218);
etaVec.push_back(-0.305);  phiVec.push_back( +0.236);
etaVec.push_back(-0.305);  phiVec.push_back( +0.253);
etaVec.push_back(-0.287);  phiVec.push_back( +0.183);
etaVec.push_back(-0.287);  phiVec.push_back( +0.201);
etaVec.push_back(-0.287);  phiVec.push_back( +0.218);
etaVec.push_back(-0.287);  phiVec.push_back( +0.236);
etaVec.push_back(-0.287);  phiVec.push_back( +0.253);
etaVec.push_back(-0.287);  phiVec.push_back( -1.422);
etaVec.push_back(-0.270);  phiVec.push_back( +0.183);
etaVec.push_back(-0.270);  phiVec.push_back( +0.201);
etaVec.push_back(-0.270);  phiVec.push_back( +0.218);
etaVec.push_back(-0.270);  phiVec.push_back( +0.236);
etaVec.push_back(-0.270);  phiVec.push_back( +0.253);
etaVec.push_back(-0.252);  phiVec.push_back( +0.096);
etaVec.push_back(-0.252);  phiVec.push_back( +0.113);
etaVec.push_back(-0.252);  phiVec.push_back( +0.131);
etaVec.push_back(-0.252);  phiVec.push_back( +0.148);
etaVec.push_back(-0.252);  phiVec.push_back( +0.166);
etaVec.push_back(-0.235);  phiVec.push_back( +0.096);
etaVec.push_back(-0.235);  phiVec.push_back( +0.113);
etaVec.push_back(-0.235);  phiVec.push_back( +0.131);
etaVec.push_back(-0.235);  phiVec.push_back( +0.148);
etaVec.push_back(-0.235);  phiVec.push_back( +0.166);
etaVec.push_back(-0.235);  phiVec.push_back( -1.876);
etaVec.push_back(-0.218);  phiVec.push_back( +0.096);
etaVec.push_back(-0.218);  phiVec.push_back( +0.113);
etaVec.push_back(-0.218);  phiVec.push_back( +0.131);
etaVec.push_back(-0.218);  phiVec.push_back( +0.148);
etaVec.push_back(-0.218);  phiVec.push_back( +0.166);
etaVec.push_back(-0.218);  phiVec.push_back( +0.864);
etaVec.push_back(-0.200);  phiVec.push_back( +0.096);
etaVec.push_back(-0.200);  phiVec.push_back( +0.113);
etaVec.push_back(-0.200);  phiVec.push_back( +0.131);
etaVec.push_back(-0.200);  phiVec.push_back( +0.148);
etaVec.push_back(-0.200);  phiVec.push_back( +0.166);
etaVec.push_back(-0.183);  phiVec.push_back( +0.096);
etaVec.push_back(-0.183);  phiVec.push_back( +0.113);
etaVec.push_back(-0.183);  phiVec.push_back( +0.131);
etaVec.push_back(-0.183);  phiVec.push_back( +0.148);
etaVec.push_back(-0.183);  phiVec.push_back( +0.166);
etaVec.push_back(-0.183);  phiVec.push_back( -3.098);
etaVec.push_back(-0.183);  phiVec.push_back( -1.806);
etaVec.push_back(-0.165);  phiVec.push_back( +0.497);
etaVec.push_back(-0.165);  phiVec.push_back( +0.707);
etaVec.push_back(-0.165);  phiVec.push_back( +0.724);
etaVec.push_back(-0.165);  phiVec.push_back( +0.742);
etaVec.push_back(-0.165);  phiVec.push_back( +0.759);
etaVec.push_back(-0.165);  phiVec.push_back( +0.777);
etaVec.push_back(-0.165);  phiVec.push_back( -2.609);
etaVec.push_back(-0.165);  phiVec.push_back( -2.592);
etaVec.push_back(-0.165);  phiVec.push_back( -2.574);
etaVec.push_back(-0.165);  phiVec.push_back( -2.557);
etaVec.push_back(-0.165);  phiVec.push_back( -2.539);
etaVec.push_back(-0.148);  phiVec.push_back( +0.707);
etaVec.push_back(-0.148);  phiVec.push_back( +0.724);
etaVec.push_back(-0.148);  phiVec.push_back( +0.742);
etaVec.push_back(-0.148);  phiVec.push_back( +0.759);
etaVec.push_back(-0.148);  phiVec.push_back( +0.777);
etaVec.push_back(-0.148);  phiVec.push_back( -2.609);
etaVec.push_back(-0.148);  phiVec.push_back( -2.592);
etaVec.push_back(-0.148);  phiVec.push_back( -2.574);
etaVec.push_back(-0.148);  phiVec.push_back( -2.557);
etaVec.push_back(-0.148);  phiVec.push_back( -2.539);
etaVec.push_back(-0.131);  phiVec.push_back( +0.707);
etaVec.push_back(-0.131);  phiVec.push_back( +0.724);
etaVec.push_back(-0.131);  phiVec.push_back( +0.742);
etaVec.push_back(-0.131);  phiVec.push_back( +0.759);
etaVec.push_back(-0.131);  phiVec.push_back( +0.777);
etaVec.push_back(-0.131);  phiVec.push_back( -2.609);
etaVec.push_back(-0.131);  phiVec.push_back( -2.592);
etaVec.push_back(-0.131);  phiVec.push_back( -2.574);
etaVec.push_back(-0.131);  phiVec.push_back( -2.557);
etaVec.push_back(-0.131);  phiVec.push_back( -2.539);
etaVec.push_back(-0.113);  phiVec.push_back( +0.707);
etaVec.push_back(-0.113);  phiVec.push_back( +0.724);
etaVec.push_back(-0.113);  phiVec.push_back( +0.742);
etaVec.push_back(-0.113);  phiVec.push_back( +0.759);
etaVec.push_back(-0.113);  phiVec.push_back( +0.777);
etaVec.push_back(-0.113);  phiVec.push_back( -2.609);
etaVec.push_back(-0.113);  phiVec.push_back( -2.592);
etaVec.push_back(-0.113);  phiVec.push_back( -2.574);
etaVec.push_back(-0.113);  phiVec.push_back( -2.557);
etaVec.push_back(-0.113);  phiVec.push_back( -2.539);
etaVec.push_back(-0.113);  phiVec.push_back( -0.689);
etaVec.push_back(-0.096);  phiVec.push_back( +0.393);
etaVec.push_back(-0.096);  phiVec.push_back( +0.707);
etaVec.push_back(-0.096);  phiVec.push_back( +0.724);
etaVec.push_back(-0.096);  phiVec.push_back( +0.742);
etaVec.push_back(-0.096);  phiVec.push_back( +0.759);
etaVec.push_back(-0.096);  phiVec.push_back( +0.777);
etaVec.push_back(-0.096);  phiVec.push_back( +3.063);
etaVec.push_back(-0.096);  phiVec.push_back( -2.609);
etaVec.push_back(-0.096);  phiVec.push_back( -2.592);
etaVec.push_back(-0.096);  phiVec.push_back( -2.574);
etaVec.push_back(-0.096);  phiVec.push_back( -2.557);
etaVec.push_back(-0.096);  phiVec.push_back( -2.539);
etaVec.push_back(-0.078);  phiVec.push_back( +0.881);
etaVec.push_back(-0.078);  phiVec.push_back( +0.899);
etaVec.push_back(-0.078);  phiVec.push_back( +0.916);
etaVec.push_back(-0.078);  phiVec.push_back( +0.934);
etaVec.push_back(-0.078);  phiVec.push_back( +0.951);
etaVec.push_back(-0.061);  phiVec.push_back( +0.881);
etaVec.push_back(-0.061);  phiVec.push_back( +0.899);
etaVec.push_back(-0.061);  phiVec.push_back( +0.916);
etaVec.push_back(-0.061);  phiVec.push_back( +0.934);
etaVec.push_back(-0.061);  phiVec.push_back( +0.951);
etaVec.push_back(-0.043);  phiVec.push_back( +0.881);
etaVec.push_back(-0.043);  phiVec.push_back( +0.899);
etaVec.push_back(-0.043);  phiVec.push_back( +0.916);
etaVec.push_back(-0.043);  phiVec.push_back( +0.934);
etaVec.push_back(-0.043);  phiVec.push_back( +0.951);
etaVec.push_back(-0.026);  phiVec.push_back( +0.881);
etaVec.push_back(-0.026);  phiVec.push_back( +0.899);
etaVec.push_back(-0.026);  phiVec.push_back( +0.916);
etaVec.push_back(-0.026);  phiVec.push_back( +0.934);
etaVec.push_back(-0.026);  phiVec.push_back( +0.951);
etaVec.push_back(-0.026);  phiVec.push_back( -2.976);
etaVec.push_back(-0.009);  phiVec.push_back( +0.881);
etaVec.push_back(-0.009);  phiVec.push_back( +0.899);
etaVec.push_back(-0.009);  phiVec.push_back( +0.916);
etaVec.push_back(-0.009);  phiVec.push_back( +0.934);
etaVec.push_back(-0.009);  phiVec.push_back( +0.951);
etaVec.push_back(-0.009);  phiVec.push_back( +3.063);
etaVec.push_back(-0.009);  phiVec.push_back( -1.894);
etaVec.push_back(+0.061);  phiVec.push_back( -0.759);
etaVec.push_back(+0.096);  phiVec.push_back( +1.754);
etaVec.push_back(+0.096);  phiVec.push_back( +1.772);
etaVec.push_back(+0.096);  phiVec.push_back( +1.789);
etaVec.push_back(+0.096);  phiVec.push_back( +1.806);
etaVec.push_back(+0.096);  phiVec.push_back( +1.824);
etaVec.push_back(+0.096);  phiVec.push_back( +2.697);
etaVec.push_back(+0.096);  phiVec.push_back( -1.388);
etaVec.push_back(+0.096);  phiVec.push_back( -1.126);
etaVec.push_back(+0.113);  phiVec.push_back( +1.754);
etaVec.push_back(+0.113);  phiVec.push_back( +1.772);
etaVec.push_back(+0.113);  phiVec.push_back( +1.789);
etaVec.push_back(+0.113);  phiVec.push_back( +1.806);
etaVec.push_back(+0.113);  phiVec.push_back( +1.824);
etaVec.push_back(+0.113);  phiVec.push_back( -0.497);
etaVec.push_back(+0.131);  phiVec.push_back( +1.754);
etaVec.push_back(+0.131);  phiVec.push_back( +1.772);
etaVec.push_back(+0.131);  phiVec.push_back( +1.789);
etaVec.push_back(+0.131);  phiVec.push_back( +1.806);
etaVec.push_back(+0.131);  phiVec.push_back( +1.824);
etaVec.push_back(+0.131);  phiVec.push_back( +2.627);
etaVec.push_back(+0.131);  phiVec.push_back( -0.916);
etaVec.push_back(+0.148);  phiVec.push_back( +1.754);
etaVec.push_back(+0.148);  phiVec.push_back( +1.772);
etaVec.push_back(+0.148);  phiVec.push_back( +1.789);
etaVec.push_back(+0.148);  phiVec.push_back( +1.806);
etaVec.push_back(+0.148);  phiVec.push_back( +1.824);
etaVec.push_back(+0.165);  phiVec.push_back( +1.754);
etaVec.push_back(+0.165);  phiVec.push_back( +1.772);
etaVec.push_back(+0.165);  phiVec.push_back( +1.789);
etaVec.push_back(+0.165);  phiVec.push_back( +1.806);
etaVec.push_back(+0.165);  phiVec.push_back( +1.824);
etaVec.push_back(+0.287);  phiVec.push_back( -1.056);
etaVec.push_back(+0.305);  phiVec.push_back( -0.759);
etaVec.push_back(+0.339);  phiVec.push_back( +2.801);
etaVec.push_back(+0.444);  phiVec.push_back( +1.876);
etaVec.push_back(+0.444);  phiVec.push_back( -0.393);
etaVec.push_back(+0.531);  phiVec.push_back( +0.742);
etaVec.push_back(+0.566);  phiVec.push_back( -0.079);
etaVec.push_back(+0.583);  phiVec.push_back( -0.812);
etaVec.push_back(+0.583);  phiVec.push_back( -0.777);
etaVec.push_back(+0.618);  phiVec.push_back( +0.655);
etaVec.push_back(+0.635);  phiVec.push_back( +0.079);
etaVec.push_back(+0.652);  phiVec.push_back( +0.881);
etaVec.push_back(+0.687);  phiVec.push_back( -2.173);
etaVec.push_back(+0.722);  phiVec.push_back( -1.422);
etaVec.push_back(+0.739);  phiVec.push_back( -2.208);
etaVec.push_back(+0.774);  phiVec.push_back( +0.113);
etaVec.push_back(+0.774);  phiVec.push_back( +0.497);
etaVec.push_back(+0.774);  phiVec.push_back( -2.697);
etaVec.push_back(+0.774);  phiVec.push_back( -1.702);
etaVec.push_back(+0.774);  phiVec.push_back( -0.480);
etaVec.push_back(+0.792);  phiVec.push_back( -0.567);
etaVec.push_back(+0.792);  phiVec.push_back( -0.183);
etaVec.push_back(+0.879);  phiVec.push_back( +0.358);
etaVec.push_back(+0.879);  phiVec.push_back( +0.375);
etaVec.push_back(+0.879);  phiVec.push_back( +0.393);
etaVec.push_back(+0.879);  phiVec.push_back( +0.410);
etaVec.push_back(+0.879);  phiVec.push_back( +0.428);
etaVec.push_back(+0.879);  phiVec.push_back( +1.667);
etaVec.push_back(+0.879);  phiVec.push_back( +1.684);
etaVec.push_back(+0.879);  phiVec.push_back( +1.702);
etaVec.push_back(+0.879);  phiVec.push_back( +1.719);
etaVec.push_back(+0.879);  phiVec.push_back( +1.737);
etaVec.push_back(+0.879);  phiVec.push_back( +2.801);
etaVec.push_back(+0.879);  phiVec.push_back( +2.819);
etaVec.push_back(+0.879);  phiVec.push_back( +2.836);
etaVec.push_back(+0.879);  phiVec.push_back( +2.854);
etaVec.push_back(+0.879);  phiVec.push_back( +2.871);
etaVec.push_back(+0.896);  phiVec.push_back( +1.667);
etaVec.push_back(+0.896);  phiVec.push_back( +1.684);
etaVec.push_back(+0.896);  phiVec.push_back( +1.702);
etaVec.push_back(+0.896);  phiVec.push_back( +1.719);
etaVec.push_back(+0.896);  phiVec.push_back( +1.737);
etaVec.push_back(+0.896);  phiVec.push_back( +2.801);
etaVec.push_back(+0.896);  phiVec.push_back( +2.819);
etaVec.push_back(+0.896);  phiVec.push_back( +2.836);
etaVec.push_back(+0.896);  phiVec.push_back( +2.854);
etaVec.push_back(+0.896);  phiVec.push_back( +2.871);
etaVec.push_back(+0.896);  phiVec.push_back( -1.963);
etaVec.push_back(+0.914);  phiVec.push_back( +1.091);
etaVec.push_back(+0.914);  phiVec.push_back( +1.667);
etaVec.push_back(+0.914);  phiVec.push_back( +1.684);
etaVec.push_back(+0.914);  phiVec.push_back( +1.702);
etaVec.push_back(+0.914);  phiVec.push_back( +1.719);
etaVec.push_back(+0.914);  phiVec.push_back( +1.737);
etaVec.push_back(+0.914);  phiVec.push_back( +2.801);
etaVec.push_back(+0.914);  phiVec.push_back( +2.819);
etaVec.push_back(+0.914);  phiVec.push_back( +2.836);
etaVec.push_back(+0.914);  phiVec.push_back( +2.854);
etaVec.push_back(+0.914);  phiVec.push_back( +2.871);
etaVec.push_back(+0.931);  phiVec.push_back( +1.667);
etaVec.push_back(+0.931);  phiVec.push_back( +1.684);
etaVec.push_back(+0.931);  phiVec.push_back( +1.702);
etaVec.push_back(+0.931);  phiVec.push_back( +1.719);
etaVec.push_back(+0.931);  phiVec.push_back( +1.737);
etaVec.push_back(+0.931);  phiVec.push_back( +2.801);
etaVec.push_back(+0.931);  phiVec.push_back( +2.819);
etaVec.push_back(+0.931);  phiVec.push_back( +2.836);
etaVec.push_back(+0.931);  phiVec.push_back( +2.854);
etaVec.push_back(+0.931);  phiVec.push_back( +2.871);
etaVec.push_back(+0.931);  phiVec.push_back( +3.098);
etaVec.push_back(+0.948);  phiVec.push_back( +1.667);
etaVec.push_back(+0.948);  phiVec.push_back( +1.684);
etaVec.push_back(+0.948);  phiVec.push_back( +1.702);
etaVec.push_back(+0.948);  phiVec.push_back( +1.719);
etaVec.push_back(+0.948);  phiVec.push_back( +1.737);
etaVec.push_back(+0.948);  phiVec.push_back( +2.801);
etaVec.push_back(+0.948);  phiVec.push_back( +2.819);
etaVec.push_back(+0.948);  phiVec.push_back( +2.836);
etaVec.push_back(+0.948);  phiVec.push_back( +2.854);
etaVec.push_back(+0.948);  phiVec.push_back( +2.871);
etaVec.push_back(+0.966);  phiVec.push_back( +0.620);
etaVec.push_back(+0.966);  phiVec.push_back( +0.637);
etaVec.push_back(+0.966);  phiVec.push_back( +0.655);
etaVec.push_back(+0.966);  phiVec.push_back( +0.672);
etaVec.push_back(+0.966);  phiVec.push_back( +0.689);
etaVec.push_back(+0.983);  phiVec.push_back( +0.620);
etaVec.push_back(+0.983);  phiVec.push_back( +0.637);
etaVec.push_back(+0.983);  phiVec.push_back( +0.655);
etaVec.push_back(+0.983);  phiVec.push_back( +0.672);
etaVec.push_back(+0.983);  phiVec.push_back( +0.689);
etaVec.push_back(+1.000);  phiVec.push_back( +0.620);
etaVec.push_back(+1.000);  phiVec.push_back( +0.637);
etaVec.push_back(+1.000);  phiVec.push_back( +0.655);
etaVec.push_back(+1.000);  phiVec.push_back( +0.672);
etaVec.push_back(+1.000);  phiVec.push_back( +0.689);
etaVec.push_back(+1.018);  phiVec.push_back( +0.620);
etaVec.push_back(+1.018);  phiVec.push_back( +0.637);
etaVec.push_back(+1.018);  phiVec.push_back( +0.655);
etaVec.push_back(+1.018);  phiVec.push_back( +0.672);
etaVec.push_back(+1.018);  phiVec.push_back( +0.689);
etaVec.push_back(+1.018);  phiVec.push_back( +2.417);
etaVec.push_back(+1.035);  phiVec.push_back( +0.620);
etaVec.push_back(+1.035);  phiVec.push_back( +0.637);
etaVec.push_back(+1.035);  phiVec.push_back( +0.655);
etaVec.push_back(+1.035);  phiVec.push_back( +0.672);
etaVec.push_back(+1.035);  phiVec.push_back( +0.689);
etaVec.push_back(+1.053);  phiVec.push_back( -3.133);
etaVec.push_back(+1.053);  phiVec.push_back( -3.115);
etaVec.push_back(+1.053);  phiVec.push_back( -3.098);
etaVec.push_back(+1.053);  phiVec.push_back( -3.081);
etaVec.push_back(+1.053);  phiVec.push_back( -3.063);
etaVec.push_back(+1.070);  phiVec.push_back( -3.133);
etaVec.push_back(+1.070);  phiVec.push_back( -3.115);
etaVec.push_back(+1.070);  phiVec.push_back( -3.098);
etaVec.push_back(+1.070);  phiVec.push_back( -3.081);
etaVec.push_back(+1.070);  phiVec.push_back( -3.063);
etaVec.push_back(+1.087);  phiVec.push_back( -3.133);
etaVec.push_back(+1.087);  phiVec.push_back( -3.115);
etaVec.push_back(+1.087);  phiVec.push_back( -3.098);
etaVec.push_back(+1.087);  phiVec.push_back( -3.081);
etaVec.push_back(+1.087);  phiVec.push_back( -3.063);
etaVec.push_back(+1.105);  phiVec.push_back( -3.133);
etaVec.push_back(+1.105);  phiVec.push_back( -3.115);
etaVec.push_back(+1.105);  phiVec.push_back( -3.098);
etaVec.push_back(+1.105);  phiVec.push_back( -3.081);
etaVec.push_back(+1.105);  phiVec.push_back( -3.063);
etaVec.push_back(+1.122);  phiVec.push_back( +3.063);
etaVec.push_back(+1.122);  phiVec.push_back( -3.133);
etaVec.push_back(+1.122);  phiVec.push_back( -3.115);
etaVec.push_back(+1.122);  phiVec.push_back( -3.098);
etaVec.push_back(+1.122);  phiVec.push_back( -3.081);
etaVec.push_back(+1.122);  phiVec.push_back( -3.063);
etaVec.push_back(+1.122);  phiVec.push_back( -2.871);
etaVec.push_back(+1.140);  phiVec.push_back( +0.026);
etaVec.push_back(+1.140);  phiVec.push_back( +2.365);
etaVec.push_back(+1.140);  phiVec.push_back( -0.846);
etaVec.push_back(+1.157);  phiVec.push_back( -1.230);
etaVec.push_back(+1.174);  phiVec.push_back( +1.667);
etaVec.push_back(+1.174);  phiVec.push_back( -1.894);
etaVec.push_back(+1.192);  phiVec.push_back( +1.091);
etaVec.push_back(+1.227);  phiVec.push_back( +2.435);
etaVec.push_back(+1.227);  phiVec.push_back( -0.515);
etaVec.push_back(+1.227);  phiVec.push_back( -0.497);
etaVec.push_back(+1.227);  phiVec.push_back( -0.480);
etaVec.push_back(+1.227);  phiVec.push_back( -0.463);
etaVec.push_back(+1.227);  phiVec.push_back( -0.445);
etaVec.push_back(+1.244);  phiVec.push_back( +1.649);
etaVec.push_back(+1.244);  phiVec.push_back( -0.986);
etaVec.push_back(+1.244);  phiVec.push_back( -0.515);
etaVec.push_back(+1.244);  phiVec.push_back( -0.497);
etaVec.push_back(+1.244);  phiVec.push_back( -0.480);
etaVec.push_back(+1.244);  phiVec.push_back( -0.463);
etaVec.push_back(+1.244);  phiVec.push_back( -0.445);
etaVec.push_back(+1.262);  phiVec.push_back( -0.515);
etaVec.push_back(+1.262);  phiVec.push_back( -0.497);
etaVec.push_back(+1.262);  phiVec.push_back( -0.480);
etaVec.push_back(+1.262);  phiVec.push_back( -0.463);
etaVec.push_back(+1.262);  phiVec.push_back( -0.445);
etaVec.push_back(+1.279);  phiVec.push_back( -0.515);
etaVec.push_back(+1.279);  phiVec.push_back( -0.497);
etaVec.push_back(+1.279);  phiVec.push_back( -0.480);
etaVec.push_back(+1.279);  phiVec.push_back( -0.463);
etaVec.push_back(+1.279);  phiVec.push_back( -0.445);
etaVec.push_back(+1.296);  phiVec.push_back( -0.515);
etaVec.push_back(+1.296);  phiVec.push_back( -0.497);
etaVec.push_back(+1.296);  phiVec.push_back( -0.480);
etaVec.push_back(+1.296);  phiVec.push_back( -0.463);
etaVec.push_back(+1.296);  phiVec.push_back( -0.445);
etaVec.push_back(+1.296);  phiVec.push_back( -0.288);
etaVec.push_back(+1.314);  phiVec.push_back( -0.201);
etaVec.push_back(+1.401);  phiVec.push_back( +2.714);
etaVec.push_back(+1.401);  phiVec.push_back( +2.731);
etaVec.push_back(+1.401);  phiVec.push_back( +2.749);
etaVec.push_back(+1.401);  phiVec.push_back( +2.766);
etaVec.push_back(+1.401);  phiVec.push_back( +2.784);
etaVec.push_back(+1.401);  phiVec.push_back( -2.609);
etaVec.push_back(+1.401);  phiVec.push_back( -2.592);
etaVec.push_back(+1.401);  phiVec.push_back( -2.574);
etaVec.push_back(+1.401);  phiVec.push_back( -2.557);
etaVec.push_back(+1.401);  phiVec.push_back( -2.539);
etaVec.push_back(+1.401);  phiVec.push_back( -0.253);
etaVec.push_back(+1.401);  phiVec.push_back( -0.236);
etaVec.push_back(+1.401);  phiVec.push_back( -0.218);
etaVec.push_back(+1.401);  phiVec.push_back( -0.201);
etaVec.push_back(+1.401);  phiVec.push_back( -0.183);
etaVec.push_back(+1.418);  phiVec.push_back( +0.375);
etaVec.push_back(+1.418);  phiVec.push_back( +2.714);
etaVec.push_back(+1.418);  phiVec.push_back( +2.731);
etaVec.push_back(+1.418);  phiVec.push_back( +2.749);
etaVec.push_back(+1.418);  phiVec.push_back( +2.766);
etaVec.push_back(+1.418);  phiVec.push_back( +2.784);
etaVec.push_back(+1.418);  phiVec.push_back( -2.609);
etaVec.push_back(+1.418);  phiVec.push_back( -2.592);
etaVec.push_back(+1.418);  phiVec.push_back( -2.574);
etaVec.push_back(+1.418);  phiVec.push_back( -2.557);
etaVec.push_back(+1.418);  phiVec.push_back( -2.539);
etaVec.push_back(+1.418);  phiVec.push_back( -0.253);
etaVec.push_back(+1.418);  phiVec.push_back( -0.236);
etaVec.push_back(+1.418);  phiVec.push_back( -0.218);
etaVec.push_back(+1.418);  phiVec.push_back( -0.201);
etaVec.push_back(+1.418);  phiVec.push_back( -0.183);
etaVec.push_back(+1.436);  phiVec.push_back( +2.714);
etaVec.push_back(+1.436);  phiVec.push_back( +2.731);
etaVec.push_back(+1.436);  phiVec.push_back( +2.749);
etaVec.push_back(+1.436);  phiVec.push_back( +2.766);
etaVec.push_back(+1.436);  phiVec.push_back( +2.784);
etaVec.push_back(+1.436);  phiVec.push_back( -2.609);
etaVec.push_back(+1.436);  phiVec.push_back( -2.592);
etaVec.push_back(+1.436);  phiVec.push_back( -2.574);
etaVec.push_back(+1.436);  phiVec.push_back( -2.557);
etaVec.push_back(+1.436);  phiVec.push_back( -2.539);
etaVec.push_back(+1.436);  phiVec.push_back( -0.253);
etaVec.push_back(+1.436);  phiVec.push_back( -0.236);
etaVec.push_back(+1.436);  phiVec.push_back( -0.218);
etaVec.push_back(+1.436);  phiVec.push_back( -0.201);
etaVec.push_back(+1.436);  phiVec.push_back( -0.183);
etaVec.push_back(+1.453);  phiVec.push_back( +2.714);
etaVec.push_back(+1.453);  phiVec.push_back( +2.731);
etaVec.push_back(+1.453);  phiVec.push_back( +2.749);
etaVec.push_back(+1.453);  phiVec.push_back( +2.766);
etaVec.push_back(+1.453);  phiVec.push_back( +2.784);
etaVec.push_back(+1.453);  phiVec.push_back( -0.253);
etaVec.push_back(+1.453);  phiVec.push_back( -0.236);
etaVec.push_back(+1.453);  phiVec.push_back( -0.218);
etaVec.push_back(+1.453);  phiVec.push_back( -0.201);
etaVec.push_back(+1.453);  phiVec.push_back( -0.183);
etaVec.push_back(+1.470);  phiVec.push_back( +2.714);
etaVec.push_back(+1.470);  phiVec.push_back( +2.731);
etaVec.push_back(+1.470);  phiVec.push_back( +2.749);
etaVec.push_back(+1.470);  phiVec.push_back( +2.766);
etaVec.push_back(+1.470);  phiVec.push_back( +2.784);
etaVec.push_back(+1.470);  phiVec.push_back( +2.801);
etaVec.push_back(+1.470);  phiVec.push_back( -0.253);
etaVec.push_back(+1.470);  phiVec.push_back( -0.236);
etaVec.push_back(+1.470);  phiVec.push_back( -0.218);
etaVec.push_back(+1.470);  phiVec.push_back( -0.201);
etaVec.push_back(+1.470);  phiVec.push_back( -0.183);
etaVec.push_back(+1.507);  phiVec.push_back( -3.092);
etaVec.push_back(+1.526);  phiVec.push_back( -3.091);
etaVec.push_back(+1.546);  phiVec.push_back( -3.090);
etaVec.push_back(+1.566);  phiVec.push_back( -3.089);
etaVec.push_back(+1.587);  phiVec.push_back( -3.088);
etaVec.push_back(+1.574);  phiVec.push_back( -2.871);
etaVec.push_back(+1.595);  phiVec.push_back( -2.865);
etaVec.push_back(+1.615);  phiVec.push_back( -2.860);
etaVec.push_back(+1.637);  phiVec.push_back( -2.854);
etaVec.push_back(+1.658);  phiVec.push_back( -2.847);
etaVec.push_back(+1.597);  phiVec.push_back( +2.679);
etaVec.push_back(+1.588);  phiVec.push_back( +2.660);
etaVec.push_back(+1.899);  phiVec.push_back( +2.757);
etaVec.push_back(+1.996);  phiVec.push_back( +2.748);
etaVec.push_back(+1.555);  phiVec.push_back( -2.162);
etaVec.push_back(+1.674);  phiVec.push_back( +2.256);
etaVec.push_back(+1.584);  phiVec.push_back( -2.157);
etaVec.push_back(+1.854);  phiVec.push_back( +2.401);
etaVec.push_back(+1.689);  phiVec.push_back( +2.238);
etaVec.push_back(+2.307);  phiVec.push_back( -3.027);
etaVec.push_back(+1.563);  phiVec.push_back( +2.020);
etaVec.push_back(+2.113);  phiVec.push_back( +2.198);
etaVec.push_back(+2.107);  phiVec.push_back( +1.978);
etaVec.push_back(+2.835);  phiVec.push_back( +2.413);
etaVec.push_back(+2.782);  phiVec.push_back( +2.357);
etaVec.push_back(+1.522);  phiVec.push_back( +1.703);
etaVec.push_back(+1.525);  phiVec.push_back( +1.662);
etaVec.push_back(+1.609);  phiVec.push_back( -1.605);
etaVec.push_back(+1.695);  phiVec.push_back( +1.462);
etaVec.push_back(+1.505);  phiVec.push_back( -1.461);
etaVec.push_back(+1.843);  phiVec.push_back( -1.412);
etaVec.push_back(+2.164);  phiVec.push_back( -1.348);
etaVec.push_back(+2.884);  phiVec.push_back( -0.900);
etaVec.push_back(+2.835);  phiVec.push_back( -0.838);
etaVec.push_back(+2.891);  phiVec.push_back( -0.784);
etaVec.push_back(+2.782);  phiVec.push_back( -0.785);
etaVec.push_back(+2.835);  phiVec.push_back( -0.729);
etaVec.push_back(+2.052);  phiVec.push_back( +1.229);
etaVec.push_back(+1.820);  phiVec.push_back( +0.961);
etaVec.push_back(+1.567);  phiVec.push_back( -0.996);
etaVec.push_back(+1.817);  phiVec.push_back( +0.782);
etaVec.push_back(+1.798);  phiVec.push_back( +0.801);
etaVec.push_back(+1.780);  phiVec.push_back( +0.820);
etaVec.push_back(+1.761);  phiVec.push_back( +0.837);
etaVec.push_back(+1.743);  phiVec.push_back( +0.855);
etaVec.push_back(+1.551);  phiVec.push_back( +1.009);
etaVec.push_back(+1.798);  phiVec.push_back( +0.764);
etaVec.push_back(+1.780);  phiVec.push_back( +0.783);
etaVec.push_back(+1.762);  phiVec.push_back( +0.801);
etaVec.push_back(+1.744);  phiVec.push_back( +0.819);
etaVec.push_back(+1.726);  phiVec.push_back( +0.836);
etaVec.push_back(+1.779);  phiVec.push_back( +0.746);
etaVec.push_back(+1.761);  phiVec.push_back( +0.765);
etaVec.push_back(+1.744);  phiVec.push_back( +0.784);
etaVec.push_back(+1.727);  phiVec.push_back( +0.801);
etaVec.push_back(+1.709);  phiVec.push_back( +0.818);
etaVec.push_back(+1.760);  phiVec.push_back( +0.729);
etaVec.push_back(+1.743);  phiVec.push_back( +0.748);
etaVec.push_back(+1.726);  phiVec.push_back( +0.766);
etaVec.push_back(+1.709);  phiVec.push_back( +0.784);
etaVec.push_back(+1.693);  phiVec.push_back( +0.801);
etaVec.push_back(+1.740);  phiVec.push_back( +0.713);
etaVec.push_back(+1.724);  phiVec.push_back( +0.731);
etaVec.push_back(+1.708);  phiVec.push_back( +0.750);
etaVec.push_back(+1.692);  phiVec.push_back( +0.767);
etaVec.push_back(+1.676);  phiVec.push_back( +0.785);
etaVec.push_back(+1.663);  phiVec.push_back( -0.764);
etaVec.push_back(+1.678);  phiVec.push_back( -0.747);
etaVec.push_back(+1.694);  phiVec.push_back( -0.729);
etaVec.push_back(+1.709);  phiVec.push_back( -0.711);
etaVec.push_back(+1.724);  phiVec.push_back( -0.692);
etaVec.push_back(+1.739);  phiVec.push_back( -0.671);
etaVec.push_back(+1.754);  phiVec.push_back( -0.651);
etaVec.push_back(+1.769);  phiVec.push_back( -0.631);
etaVec.push_back(+1.783);  phiVec.push_back( -0.609);
etaVec.push_back(+1.798);  phiVec.push_back( -0.587);
etaVec.push_back(+1.646);  phiVec.push_back( -0.749);
etaVec.push_back(+1.661);  phiVec.push_back( -0.732);
etaVec.push_back(+1.676);  phiVec.push_back( -0.714);
etaVec.push_back(+1.691);  phiVec.push_back( -0.696);
etaVec.push_back(+1.706);  phiVec.push_back( -0.677);
etaVec.push_back(+1.719);  phiVec.push_back( -0.656);
etaVec.push_back(+1.734);  phiVec.push_back( -0.636);
etaVec.push_back(+1.748);  phiVec.push_back( -0.616);
etaVec.push_back(+1.762);  phiVec.push_back( -0.595);
etaVec.push_back(+1.776);  phiVec.push_back( -0.573);
etaVec.push_back(+1.630);  phiVec.push_back( -0.734);
etaVec.push_back(+1.644);  phiVec.push_back( -0.717);
etaVec.push_back(+1.658);  phiVec.push_back( -0.699);
etaVec.push_back(+1.673);  phiVec.push_back( -0.681);
etaVec.push_back(+1.687);  phiVec.push_back( -0.663);
etaVec.push_back(+1.700);  phiVec.push_back( -0.642);
etaVec.push_back(+1.714);  phiVec.push_back( -0.622);
etaVec.push_back(+1.727);  phiVec.push_back( -0.602);
etaVec.push_back(+1.741);  phiVec.push_back( -0.581);
etaVec.push_back(+1.754);  phiVec.push_back( -0.560);
etaVec.push_back(+1.613);  phiVec.push_back( -0.719);
etaVec.push_back(+1.627);  phiVec.push_back( -0.702);
etaVec.push_back(+1.641);  phiVec.push_back( -0.685);
etaVec.push_back(+1.655);  phiVec.push_back( -0.667);
etaVec.push_back(+1.668);  phiVec.push_back( -0.649);
etaVec.push_back(+1.681);  phiVec.push_back( -0.628);
etaVec.push_back(+1.694);  phiVec.push_back( -0.609);
etaVec.push_back(+1.707);  phiVec.push_back( -0.589);
etaVec.push_back(+1.720);  phiVec.push_back( -0.568);
etaVec.push_back(+1.733);  phiVec.push_back( -0.547);
etaVec.push_back(+1.596);  phiVec.push_back( -0.705);
etaVec.push_back(+1.610);  phiVec.push_back( -0.689);
etaVec.push_back(+1.623);  phiVec.push_back( -0.671);
etaVec.push_back(+1.637);  phiVec.push_back( -0.654);
etaVec.push_back(+1.650);  phiVec.push_back( -0.635);
etaVec.push_back(+1.662);  phiVec.push_back( -0.615);
etaVec.push_back(+1.674);  phiVec.push_back( -0.595);
etaVec.push_back(+1.687);  phiVec.push_back( -0.576);
etaVec.push_back(+1.700);  phiVec.push_back( -0.555);
etaVec.push_back(+1.712);  phiVec.push_back( -0.535);
etaVec.push_back(+1.849);  phiVec.push_back( -0.072);
etaVec.push_back(+1.520);  phiVec.push_back( -0.766);
etaVec.push_back(+1.533);  phiVec.push_back( -0.751);
etaVec.push_back(+1.546);  phiVec.push_back( -0.736);
etaVec.push_back(+1.559);  phiVec.push_back( -0.721);
etaVec.push_back(+1.572);  phiVec.push_back( -0.705);
etaVec.push_back(+1.584);  phiVec.push_back( -0.687);
etaVec.push_back(+1.597);  phiVec.push_back( -0.671);
etaVec.push_back(+1.610);  phiVec.push_back( -0.653);
etaVec.push_back(+1.622);  phiVec.push_back( -0.636);
etaVec.push_back(+1.635);  phiVec.push_back( -0.618);
etaVec.push_back(+1.506);  phiVec.push_back( -0.753);
etaVec.push_back(+1.519);  phiVec.push_back( -0.738);
etaVec.push_back(+1.531);  phiVec.push_back( -0.723);
etaVec.push_back(+1.544);  phiVec.push_back( -0.708);
etaVec.push_back(+1.557);  phiVec.push_back( -0.692);
etaVec.push_back(+1.568);  phiVec.push_back( -0.675);
etaVec.push_back(+1.580);  phiVec.push_back( -0.658);
etaVec.push_back(+1.592);  phiVec.push_back( -0.641);
etaVec.push_back(+1.605);  phiVec.push_back( -0.623);
etaVec.push_back(+1.617);  phiVec.push_back( -0.605);
etaVec.push_back(+1.551);  phiVec.push_back( -0.662);
etaVec.push_back(+1.563);  phiVec.push_back( -0.646);
etaVec.push_back(+1.575);  phiVec.push_back( -0.629);
etaVec.push_back(+1.587);  phiVec.push_back( -0.611);
etaVec.push_back(+1.599);  phiVec.push_back( -0.593);
etaVec.push_back(+1.535);  phiVec.push_back( -0.650);
etaVec.push_back(+1.547);  phiVec.push_back( -0.634);
etaVec.push_back(+1.558);  phiVec.push_back( -0.617);
etaVec.push_back(+1.570);  phiVec.push_back( -0.600);
etaVec.push_back(+1.581);  phiVec.push_back( -0.582);
etaVec.push_back(+1.519);  phiVec.push_back( -0.638);
etaVec.push_back(+1.530);  phiVec.push_back( -0.622);
etaVec.push_back(+1.542);  phiVec.push_back( -0.605);
etaVec.push_back(+1.553);  phiVec.push_back( -0.588);
etaVec.push_back(+1.564);  phiVec.push_back( -0.571);
etaVec.push_back(+1.508);  phiVec.push_back( -0.622);
etaVec.push_back(+1.518);  phiVec.push_back( -0.606);
etaVec.push_back(+1.529);  phiVec.push_back( -0.590);
etaVec.push_back(+1.539);  phiVec.push_back( -0.573);
etaVec.push_back(+1.550);  phiVec.push_back( -0.555);
etaVec.push_back(+1.605);  phiVec.push_back( -0.441);
etaVec.push_back(+1.492);  phiVec.push_back( -0.611);
etaVec.push_back(+1.502);  phiVec.push_back( -0.595);
etaVec.push_back(+1.512);  phiVec.push_back( -0.579);
etaVec.push_back(+1.523);  phiVec.push_back( -0.562);
etaVec.push_back(+1.533);  phiVec.push_back( -0.545);
etaVec.push_back(-1.495);  phiVec.push_back( -2.975);
etaVec.push_back(-1.508);  phiVec.push_back( -3.131);
etaVec.push_back(-1.530);  phiVec.push_back( -2.948);
etaVec.push_back(-1.563);  phiVec.push_back( +2.764);
etaVec.push_back(-1.825);  phiVec.push_back( -3.127);
etaVec.push_back(-1.957);  phiVec.push_back( -2.965);
etaVec.push_back(-1.744);  phiVec.push_back( -2.323);
etaVec.push_back(-1.743);  phiVec.push_back( -2.287);
etaVec.push_back(-1.584);  phiVec.push_back( +2.157);
etaVec.push_back(-1.567);  phiVec.push_back( +2.145);
etaVec.push_back(-1.950);  phiVec.push_back( -2.471);
etaVec.push_back(-1.969);  phiVec.push_back( -2.496);
etaVec.push_back(-1.987);  phiVec.push_back( -2.522);
etaVec.push_back(-2.006);  phiVec.push_back( -2.549);
etaVec.push_back(-2.024);  phiVec.push_back( -2.578);
etaVec.push_back(-1.974);  phiVec.push_back( -2.451);
etaVec.push_back(-1.994);  phiVec.push_back( -2.477);
etaVec.push_back(-2.014);  phiVec.push_back( -2.503);
etaVec.push_back(-2.034);  phiVec.push_back( -2.531);
etaVec.push_back(-2.053);  phiVec.push_back( -2.559);
etaVec.push_back(-1.999);  phiVec.push_back( -2.430);
etaVec.push_back(-2.020);  phiVec.push_back( -2.456);
etaVec.push_back(-2.041);  phiVec.push_back( -2.483);
etaVec.push_back(-2.062);  phiVec.push_back( -2.511);
etaVec.push_back(-2.082);  phiVec.push_back( -2.540);
etaVec.push_back(-2.024);  phiVec.push_back( -2.408);
etaVec.push_back(-2.047);  phiVec.push_back( -2.434);
etaVec.push_back(-2.069);  phiVec.push_back( -2.461);
etaVec.push_back(-2.090);  phiVec.push_back( -2.489);
etaVec.push_back(-2.112);  phiVec.push_back( -2.519);
etaVec.push_back(-2.050);  phiVec.push_back( -2.385);
etaVec.push_back(-2.073);  phiVec.push_back( -2.411);
etaVec.push_back(-2.096);  phiVec.push_back( -2.439);
etaVec.push_back(-2.119);  phiVec.push_back( -2.467);
etaVec.push_back(-2.142);  phiVec.push_back( -2.497);
etaVec.push_back(-1.820);  phiVec.push_back( +2.180);
etaVec.push_back(-2.063);  phiVec.push_back( -2.205);
etaVec.push_back(-1.894);  phiVec.push_back( -1.893);
etaVec.push_back(-1.513);  phiVec.push_back( -1.761);
etaVec.push_back(-1.505);  phiVec.push_back( -1.681);
etaVec.push_back(-2.061);  phiVec.push_back( -1.732);
etaVec.push_back(-2.108);  phiVec.push_back( -1.628);
etaVec.push_back(-1.853);  phiVec.push_back( +1.587);
etaVec.push_back(-2.187);  phiVec.push_back( -1.549);
etaVec.push_back(-2.228);  phiVec.push_back( -1.548);
etaVec.push_back(-2.270);  phiVec.push_back( -1.547);
etaVec.push_back(-2.315);  phiVec.push_back( -1.545);
etaVec.push_back(-2.361);  phiVec.push_back( -1.544);
etaVec.push_back(-1.722);  phiVec.push_back( -1.532);
etaVec.push_back(-2.185);  phiVec.push_back( -1.508);
etaVec.push_back(-2.226);  phiVec.push_back( -1.505);
etaVec.push_back(-2.268);  phiVec.push_back( -1.503);
etaVec.push_back(-2.312);  phiVec.push_back( -1.500);
etaVec.push_back(-2.359);  phiVec.push_back( -1.497);
etaVec.push_back(-2.182);  phiVec.push_back( -1.469);
etaVec.push_back(-2.222);  phiVec.push_back( -1.464);
etaVec.push_back(-2.264);  phiVec.push_back( -1.460);
etaVec.push_back(-2.308);  phiVec.push_back( -1.455);
etaVec.push_back(-2.354);  phiVec.push_back( -1.450);
etaVec.push_back(-2.177);  phiVec.push_back( -1.429);
etaVec.push_back(-2.217);  phiVec.push_back( -1.423);
etaVec.push_back(-2.259);  phiVec.push_back( -1.417);
etaVec.push_back(-2.302);  phiVec.push_back( -1.410);
etaVec.push_back(-2.348);  phiVec.push_back( -1.403);
etaVec.push_back(-2.171);  phiVec.push_back( -1.390);
etaVec.push_back(-2.210);  phiVec.push_back( -1.383);
etaVec.push_back(-2.251);  phiVec.push_back( -1.375);
etaVec.push_back(-2.294);  phiVec.push_back( -1.366);
etaVec.push_back(-2.339);  phiVec.push_back( -1.357);
etaVec.push_back(-2.542);  phiVec.push_back( -0.906);
etaVec.push_back(-2.506);  phiVec.push_back( +0.706);
etaVec.push_back(-1.743);  phiVec.push_back( -0.855);
etaVec.push_back(-1.761);  phiVec.push_back( -0.837);
etaVec.push_back(-1.551);  phiVec.push_back( +1.009);
etaVec.push_back(-1.572);  phiVec.push_back( +0.967);
etaVec.push_back(-1.555);  phiVec.push_back( +0.979);
etaVec.push_back(-2.072);  phiVec.push_back( +0.019);
etaVec.push_back(-1.962);  phiVec.push_back( +0.144);
etaVec.push_back(-1.635);  phiVec.push_back( +0.618);
etaVec.push_back(-1.622);  phiVec.push_back( +0.636);
etaVec.push_back(-1.610);  phiVec.push_back( +0.653);
etaVec.push_back(-1.597);  phiVec.push_back( +0.671);
etaVec.push_back(-1.584);  phiVec.push_back( +0.687);
etaVec.push_back(-1.617);  phiVec.push_back( +0.605);
etaVec.push_back(-1.605);  phiVec.push_back( +0.623);
etaVec.push_back(-1.592);  phiVec.push_back( +0.641);
etaVec.push_back(-1.580);  phiVec.push_back( +0.658);
etaVec.push_back(-1.568);  phiVec.push_back( +0.675);
etaVec.push_back(-1.599);  phiVec.push_back( +0.593);
etaVec.push_back(-1.587);  phiVec.push_back( +0.611);
etaVec.push_back(-1.575);  phiVec.push_back( +0.629);
etaVec.push_back(-1.563);  phiVec.push_back( +0.646);
etaVec.push_back(-1.551);  phiVec.push_back( +0.662);
etaVec.push_back(-1.581);  phiVec.push_back( +0.582);
etaVec.push_back(-1.570);  phiVec.push_back( +0.600);
etaVec.push_back(-1.558);  phiVec.push_back( +0.617);
etaVec.push_back(-1.547);  phiVec.push_back( +0.634);
etaVec.push_back(-1.535);  phiVec.push_back( +0.650);
etaVec.push_back(-1.564);  phiVec.push_back( +0.571);
etaVec.push_back(-1.553);  phiVec.push_back( +0.588);
etaVec.push_back(-1.542);  phiVec.push_back( +0.605);
etaVec.push_back(-1.530);  phiVec.push_back( +0.622);
etaVec.push_back(-1.519);  phiVec.push_back( +0.638);
etaVec.push_back(-1.522);  phiVec.push_back( +0.112);
etaVec.push_back(-1.504);  phiVec.push_back( -0.089);
}
