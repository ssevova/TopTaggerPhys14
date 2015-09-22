//================================================================================================
//
// Perform top MVA ID training
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                          // access to gROOT, entry point to ROOT system
#include <TSystem.h>                        // interface to OS
#include <TStyle.h>                         // class to handle ROOT plotting styles
#include <TFile.h>                          // file handle class
#include <TTree.h>                          // class to access ntuples
#include <TLorentzVector.h>                 // 4-vector class
#include <TVector2.h>                       // 2-vector class
#include <TMVA/Factory.h>                   // ROOT MVA library
#include <vector>                           // STL vector class
#include <iostream>                         // standard I/O
#include <iomanip>                          // functions to format standard I/O
#include <fstream>                          // functions for file I/O
#include <string>                           // C++ string class
#include <cmath>                            // C++ math library
#include <cassert>

#endif

using namespace std;

void mvaTrain()
{
  //--------------------------------------------------------------------------------------------------------------
  // Settings
  //==============================================================================================================
  const bool   doHighPt  = true;
  const string infilename("trainingbits_semilep.root");

  string outfilename;
  // if(doHighPt){
  //   outfilename = "toptrainingbits_noMdrop_noDR_highpt.root";
  // }
  // else{
  //   outfilename = "toptrainingbits_noMdrop_noDR_lowpt.root";
  // }

  outfilename = "toptrainingbits_baseMVA.root";
  const double TOPMASSLOW  = 100;
  const double TOPMASSHIGH = 300;
  
  //--------------------------------------------------------------------------------------------------------------
  // Set up output file
  //==============================================================================================================
  
  unsigned int isSig, b_mis, w_mis, wb_mis, evtType, eventNum;
  TLorentzVector *vjet1=0, *vjet2=0, *vjet3=0;
  float qgid1, qgid2, detaj1b, detaj2b, dphij1b, dphij2b, drj1b, drj2b, bjcsv, jet1csv, jet2csv, mdrop, mtop, prob, cost, cosTS1, cosTS2; 
   
  //
  // trees for TMVA training
  //
  TTree *t_sig = new TTree("signal", "signal");            t_sig->SetDirectory(0);
  TTree *t_bkg1 = new TTree("background1", "background1"); t_bkg1->SetDirectory(0);
  TTree *t_bkg2 = new TTree("background2", "background2"); t_bkg2->SetDirectory(0);
  TTree *t_bkg3 = new TTree("background3", "background3"); t_bkg3->SetDirectory(0);

  t_sig->Branch("isSig",  &isSig);
  t_sig->Branch("b_mis",  &b_mis);
  t_sig->Branch("w_mis",  &w_mis);
  t_sig->Branch("wb_mis", &wb_mis);
  t_sig->Branch("mtop", &mtop);
  // t_sig->Branch("mdrop", &mdrop);
  t_sig->Branch("qgid1",  &qgid1);
  t_sig->Branch("qgid2",  &qgid2);
  // t_sig->Branch("detaj1b",&detaj1b);
  // t_sig->Branch("detaj2b",&detaj2b);
  t_sig->Branch("dphij1b",&dphij1b);
  t_sig->Branch("dphij2b",&dphij2b);
  t_sig->Branch("drj1b",  &drj1b);
  t_sig->Branch("drj2b",  &drj2b);
  t_sig->Branch("bjcsv",  &bjcsv);
  t_sig->Branch("jet1csv",  &jet1csv);
  t_sig->Branch("jet2csv",  &jet2csv);
  // t_sig->Branch("prob",   &prob);
  // t_sig->Branch("cost",   &cost);
  // t_sig->Branch("cosTS1",  &cosTS1);
  // t_sig->Branch("cosTS2",  &cosTS2);   
 
  t_bkg1->Branch("isSig",  &isSig);
  t_bkg1->Branch("b_mis",  &b_mis);
  t_bkg1->Branch("w_mis",  &w_mis);
  t_bkg1->Branch("wb_mis", &wb_mis);
  t_bkg1->Branch("mtop",  &mtop);
  // t_bkg1->Branch("mdrop",  &mdrop);
  t_bkg1->Branch("qgid1",  &qgid1);
  t_bkg1->Branch("qgid2",  &qgid2);
  // t_bkg1->Branch("detaj1b",&detaj1b);
  // t_bkg1->Branch("detaj2b",&detaj2b);
  t_bkg1->Branch("dphij1b",&dphij1b);
  t_bkg1->Branch("dphij2b",&dphij2b);
  t_bkg1->Branch("drj1b",  &drj1b);
  t_bkg1->Branch("drj2b",  &drj2b);
  t_bkg1->Branch("bjcsv",  &bjcsv);
  t_bkg1->Branch("jet1csv",  &jet1csv);
  t_bkg1->Branch("jet2csv",  &jet2csv);
  // t_bkg1->Branch("prob",  &prob);
  // t_bkg1->Branch("cost",  &cost);
  // t_bkg1->Branch("cosTS1",  &cosTS1);
  // t_bkg1->Branch("cosTS2",  &cosTS2);     

  
  t_bkg2->Branch("isSig",  &isSig);
  t_bkg2->Branch("b_mis",  &b_mis);
  t_bkg2->Branch("w_mis",  &w_mis);
  t_bkg2->Branch("wb_mis", &wb_mis);
  t_bkg2->Branch("mtop",  &mtop);
  // t_bkg2->Branch("mdrop",  &mdrop);
  t_bkg2->Branch("qgid1",  &qgid1);
  t_bkg2->Branch("qgid2",  &qgid2);
  // t_bkg2->Branch("detaj1b",&detaj1b);
  // t_bkg2->Branch("detaj2b",&detaj2b);
  t_bkg2->Branch("dphij1b",&dphij1b);
  t_bkg2->Branch("dphij2b",&dphij2b);
  t_bkg2->Branch("drj1b",  &drj1b);
  t_bkg2->Branch("drj2b",  &drj2b);
  t_bkg2->Branch("bjcsv",  &bjcsv);
  t_bkg2->Branch("jet1csv",  &jet1csv);
  t_bkg2->Branch("jet2csv",  &jet2csv);
  // t_bkg2->Branch("prob",     &prob);
  // t_bkg2->Branch("cost",     &cost);
  // t_bkg2->Branch("cosTS1",  &cosTS1);
  // t_bkg2->Branch("cosTS2",  &cosTS2);     
  
  t_bkg3->Branch("isSig",  &isSig);
  t_bkg3->Branch("b_mis",  &b_mis);
  t_bkg3->Branch("w_mis",  &w_mis);
  t_bkg3->Branch("wb_mis", &wb_mis);
  t_bkg3->Branch("mtop",  &mtop);
  // t_bkg3->Branch("mdrop",  &mdrop);
  t_bkg3->Branch("qgid1",  &qgid1);
  t_bkg3->Branch("qgid2",  &qgid2);
  // t_bkg3->Branch("detaj1b",&detaj1b);
  // t_bkg3->Branch("detaj2b",&detaj2b);
  t_bkg3->Branch("dphij1b",&dphij1b);
  t_bkg3->Branch("dphij2b",&dphij2b);
  t_bkg3->Branch("drj1b",  &drj1b);
  t_bkg3->Branch("drj2b",  &drj2b);
  t_bkg3->Branch("bjcsv",  &bjcsv);
  t_bkg3->Branch("jet1csv",  &jet1csv);
  t_bkg3->Branch("jet2csv",  &jet2csv);
  // t_bkg3->Branch("prob",     &prob);
  // t_bkg3->Branch("cost",     &cost);
  // t_bkg3->Branch("cosTS1",  &cosTS1);
  // t_bkg3->Branch("cosTS2",  &cosTS2);
  
  //--------------------------------------------------------------------------------------------------------------
  // Process input file
  //==============================================================================================================
  std::cout << "Processing " << infilename << "..." << std::endl;
    
  TFile *infile = TFile::Open(infilename.c_str()); assert(infile);
  TTree *intree = (TTree*)infile->Get("Events");   assert(intree);
  
  intree->SetBranchAddress("isSig",   &isSig);
  intree->SetBranchAddress("b_mis",   &b_mis);
  intree->SetBranchAddress("w_mis",   &w_mis);
  intree->SetBranchAddress("wb_mis",  &wb_mis);
  intree->SetBranchAddress("eventNum",&eventNum);   
  intree->SetBranchAddress("evtType", &evtType);
  intree->SetBranchAddress("mtop",   &mtop);
  intree->SetBranchAddress("qgid1",   &qgid1);
  intree->SetBranchAddress("qgid2",   &qgid2);
  // intree->SetBranchAddress("mdrop",   &mdrop);
  // intree->SetBranchAddress("detaj1b", &detaj1b);
  // intree->SetBranchAddress("detaj2b", &detaj2b);
  intree->SetBranchAddress("dphij1b", &dphij1b);
  intree->SetBranchAddress("dphij2b", &dphij2b);
  intree->SetBranchAddress("drj1b",   &drj1b);
  intree->SetBranchAddress("drj2b",   &drj2b);
  intree->SetBranchAddress("bjcsv",   &bjcsv);
  intree->SetBranchAddress("jet1csv",   &jet1csv);
  intree->SetBranchAddress("jet2csv",   &jet2csv);
  intree->SetBranchAddress("vjet1",   &vjet1);
  intree->SetBranchAddress("vjet2",   &vjet2);
  intree->SetBranchAddress("vjet3",   &vjet3);
  // intree->SetBranchAddress("prob",   &prob);
  // intree->SetBranchAddress("cost",   &cost);
  // intree->SetBranchAddress("cosTS1",  &cosTS1);
  // intree->SetBranchAddress("cosTS2",  &cosTS2);
  
  for(unsigned int ientry=0; ientry<intree->GetEntries(); ientry++) {  
    intree->GetEntry(ientry);
    if(evtType == 0) continue;  //skip events reserved for testing
    
    TLorentzVector dijet;
    dijet = (*vjet1 + *vjet2);

    if(mtop<TOPMASSLOW || mtop>TOPMASSHIGH) continue;
    // if(doHighPt) {
    //   if(dijet.Pt()<160) continue;
    // }
    // else{
    //   if(dijet.Pt()>160) continue;
    // }
  
    if     (isSig==1)   t_sig ->Fill();  
    else if(b_mis==1)   t_bkg1->Fill();  
    else if(w_mis==1)   t_bkg2->Fill();  
    else if(wb_mis==1)  t_bkg3->Fill();  
  }
  
  delete infile;
  infile=0, intree=0;


  //--------------------------------------------------------------------------------------------------------------
  // Perform TMVA training
  //==============================================================================================================
    
  std::cout << std::endl;
  std::cout << "============= time for TMVA ==================" << std::endl;
  std::cout << std::endl;

  TFile *tf = new TFile(outfilename.c_str(), "RECREATE");
  TMVA::Factory * factory = new TMVA::Factory(outfilename.c_str(), tf, 
                                              "!V:!Silent:Color:DrawProgressBar:Transformations=I");
  factory->AddSignalTree(t_sig);
  factory->AddBackgroundTree(t_bkg1);
  factory->AddBackgroundTree(t_bkg2);
  factory->AddBackgroundTree(t_bkg3);
  factory->AddSpectator("isSig",  'I');
  factory->AddSpectator("b_mis",  'I');
  factory->AddSpectator("w_mis",  'I');
  factory->AddSpectator("wb_mis",  'I');
  factory->AddSpectator("mtop"    ,'F');
  factory->AddVariable("qgid1"   ,'F');
  factory->AddVariable("qgid2"   ,'F');
  // factory->AddVariable("mdrop",   'F');
  // factory->AddVariable("detaj1b", 'F');
  // factory->AddVariable("detaj2b", 'F');
  factory->AddVariable("dphij1b", 'F');
  factory->AddVariable("dphij2b", 'F');
  factory->AddVariable("drj1b",   'F');
  factory->AddVariable("drj2b",   'F');
  factory->AddVariable("bjcsv",   'F');
  factory->AddVariable("jet1csv",   'F');
  factory->AddVariable("jet2csv",   'F');
  // factory->AddVariable("prob",    'F');
  // factory->AddVariable("cost",    'F');
  // factory->AddVariable("cosTS1",   'F');
  // factory->AddVariable("cosTS2",   'F');
 
  factory->BookMethod(TMVA::Types::kBDT, "BDTG", 
		      "!H:!V:NTrees=1000:BoostType=Grad:UseBaggedBoost=F:nCuts=20000:Shrinkage=0.05:MaxDepth=3:UseYesNoLeaf=F:MinNodeSize=5");
  //factory->BookMethod(TMVA::Types::kBDT, "BDTG", 
  //                    "!H:!V:NTrees=1000:BoostType=Grad:UseBaggedGrad=F:nCuts=20000:Shrinkage=0.05:MaxDepth=3:UseYesNoLeaf=F:nEventsMin=200");
  TCut mycut("");
  factory->PrepareTrainingAndTestTree(mycut,"nTrain_Signal=100000:nTrain_Background=200000:nTest_Signal=5000:nTest_Background=5000");
  factory->TrainAllMethods();
  factory->TestAllMethods();
  factory->EvaluateAllMethods();
  tf->Close();  
}


