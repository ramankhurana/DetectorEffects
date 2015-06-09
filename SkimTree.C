void SkimTree(TString inputfilename, TString outputfilename){
  TFile* inputfile = new TFile(inputfilename,"READ");
  inputfile->cd();
  tree->cd();
  TTree* newtree = tree->CopyTree("pfpatgenMetPt_<30","",1000000000,0);
  TFile* outputfile = new TFile(outputfilename,"RECREATE");
  outputfile->cd();
  newtree->Write();
  outputfile->Write();

}
