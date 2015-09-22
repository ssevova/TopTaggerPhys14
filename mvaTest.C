//================================================================================================
//
// Create bacon bits for testing top MVA ID
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                          // access to gROOT, entry point to ROOT system
#include <TSystem.h>                        // interface to OS
#include <TFile.h>                          // file handle class
#include <TTree.h>                          // class to access ntuples
#include <TLorentzVector.h>
#include <vector>                           // STL vector class
#include <iostream>                         // standard I/O
#include <iomanip>                          // functions to format standard I/O
#include <fstream>                          // functions for file I/O
#include <string>                           // C++ string class
#include <cmath>                            // C++ math library
#include <cassert>

#include <TMVA/Reader.h>

#endif

using namespace std;


//=== MAIN MACRO =================================================================================================

void mvaTest()

{
  //--------------------------------------------------------------------------------------------------------------
  // Settings
  //==============================================================================================================
   
  const string infilename("trainingbits_semilep.root");
  const string outfilename("testingbits_baseMVA_tt1l.root");
 
  const double TOPMASSLOW  = 100;
  const double TOPMASSHIGH = 300;

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code
  //==============================================================================================================
    
  //
  vector<string> weightfilesv;
  weightfilesv.push_back("weights/toptrainingbits_baseMVA.root_BDTG.weights.xml");

  unsigned int bkgType;
  unsigned int isSig, b_mis, w_mis, wb_mis, evtType, eventNum;
  int dummy, dummy2, dummy3, dummy4;
  float qgid1, qgid2, detaj1b, detaj2b, dphij1b, dphij2b, drj1b, drj2b, bjcsv, jet1csv, jet2csv, mdrop, cost, prob, cosTS1, cosTS2;
  float weight;
  float topmva;
  float mtop;
  TLorentzVector *q1vec=0, *q2vec=0, *q3vec=0;
  TLorentzVector *vjet1=0, *vjet2=0, *vjet3=0;

  TFile *outFile = new TFile(outfilename.c_str(),"RECREATE");
  TTree *outTree = new TTree("Events","Events");  

  outTree->Branch("eventNum",  &eventNum , "eventNum/i");
  outTree->Branch("bkgType",   &bkgType,   "bkgType/i");
  outTree->Branch("isSig",     &isSig,     "isSig/i");
  outTree->Branch("b_mis",     &b_mis,     "b_mis/i");
  outTree->Branch("w_mis",     &w_mis,     "w_mis/i");
  outTree->Branch("wb_mis",    &wb_mis,    "wb_mis/i");
  outTree->Branch("weight",    &weight,    "weight/F");
  outTree->Branch("mtop",      &mtop,      "mtop/F");
  outTree->Branch("qgid1",     &qgid1,     "qgid1/F");
  outTree->Branch("qgid2",     &qgid2,     "qgid2/F");
  //outTree->Branch("mdrop",     &mdrop,     "mdrop/F");
  //outTree->Branch("detaj1b",   &detaj1b,   "detaj1b/F");
  //outTree->Branch("detaj2b",   &detaj2b,   "detaj2b/F");
  outTree->Branch("dphij1b",   &dphij1b,   "dphij1b/F");
  outTree->Branch("dphij2b",   &dphij2b,   "dphij2b/F");
  outTree->Branch("drj1b",     &drj1b,     "drj1b/F");
  outTree->Branch("drj2b",     &drj2b,     "drj2b/F");
  outTree->Branch("topmva",    &topmva,    "topmva/F");
  outTree->Branch("bjcsv",     &bjcsv,     "bjcsv/F");
  outTree->Branch("jet1csv",   &jet1csv,   "jet1csv/F");
  outTree->Branch("jet2csv",   &jet2csv,   "jet2csv/F");
  // outTree->Branch("prob",      &prob,      "prob/F");
  // outTree->Branch("cost",      &cost,      "cost/F");
  // outTree->Branch("cosTS1",     &cosTS1,     "cosTS1/F");
  // outTree->Branch("cosTS2",     &cosTS2,     "cosTS2/F");
  outTree->Branch("q1vec", "TLorentzVector", &q1vec);
  outTree->Branch("q2vec", "TLorentzVector", &q2vec);
  outTree->Branch("q3vec", "TLorentzVector", &q3vec);
  outTree->Branch("vjet1", "TLorentzVector", &vjet1);
  outTree->Branch("vjet2", "TLorentzVector", &vjet2);
  outTree->Branch("vjet3", "TLorentzVector", &vjet3);
  

  
  TMVA::Reader *readers[weightfilesv.size()];
  for(unsigned int iw=0; iw<weightfilesv.size(); iw++) {
    readers[iw] = new TMVA::Reader("");
    readers[iw]->AddSpectator("isSig",&dummy);
    readers[iw]->AddSpectator("b_mis", &dummy2);
    readers[iw]->AddSpectator("w_mis", &dummy3);
    readers[iw]->AddSpectator("wb_mis",&dummy4);
    readers[iw]->AddSpectator("mtop",  &mtop);
    readers[iw]->AddVariable("qgid1",  &qgid1);
    readers[iw]->AddVariable("qgid2",  &qgid2);
    //readers[iw]->AddVariable("mdrop",  &mdrop);
    //readers[iw]->AddVariable("detaj1b",&detaj1b);
    //readers[iw]->AddVariable("detaj2b",&detaj2b);
    readers[iw]->AddVariable("dphij1b",&dphij1b);
    readers[iw]->AddVariable("dphij2b",&dphij2b);
    readers[iw]->AddVariable("drj1b",  &drj1b);
    readers[iw]->AddVariable("drj2b",  &drj2b);
    readers[iw]->AddVariable("bjcsv",  &bjcsv);
    readers[iw]->AddVariable("jet1csv", &jet1csv);
    readers[iw]->AddVariable("jet2csv", &jet2csv);
    // readers[iw]->AddVariable("prob",    &prob);
    // readers[iw]->AddVariable("cost",    &cost);
    // readers[iw]->AddVariable("cosTS1",    &cosTS1);
    // readers[iw]->AddVariable("cosTS2",    &cosTS2);
    readers[iw]->BookMVA("BDTG", weightfilesv[iw].c_str());
  }
 
  std::cout << "Processing " << infilename << "..." << std::endl;
    
  TFile *infile = TFile::Open(infilename.c_str()); assert(infile);
  TTree *intree = (TTree*)infile->Get("Events");   assert(intree);

  intree->SetBranchAddress("bkgType", &bkgType);
  intree->SetBranchAddress("isSig",  &isSig);
  intree->SetBranchAddress("b_mis",  &b_mis);
  intree->SetBranchAddress("w_mis",  &w_mis);
  intree->SetBranchAddress("wb_mis", &wb_mis);
  intree->SetBranchAddress("weight", &weight);
  intree->SetBranchAddress("eventNum", &eventNum);
  intree->SetBranchAddress("evtType",  &evtType);
  intree->SetBranchAddress("mtop",     &mtop);
  //intree->SetBranchAddress("mdrop",  &mdrop);
  intree->SetBranchAddress("qgid1",    &qgid1);
  intree->SetBranchAddress("qgid2",    &qgid2);
  //intree->SetBranchAddress("detaj1b",&detaj1b);
  //intree->SetBranchAddress("detaj2b",&detaj2b);
  intree->SetBranchAddress("dphij1b",  &dphij1b);
  intree->SetBranchAddress("dphij2b",  &dphij2b);
  intree->SetBranchAddress("drj1b",    &drj1b);
  intree->SetBranchAddress("drj2b",    &drj2b);
  intree->SetBranchAddress("bjcsv",    &bjcsv);
  intree->SetBranchAddress("jet1csv",  &jet1csv);
  intree->SetBranchAddress("jet2csv",  &jet2csv);
  // intree->SetBranchAddress("prob",     &prob);
  // intree->SetBranchAddress("cost",     &cost);
  // intree->SetBranchAddress("cosTS1",    &cosTS1);
  // intree->SetBranchAddress("cosTS2",    &cosTS2);
  intree->SetBranchAddress("q1vec",    &q1vec);
  intree->SetBranchAddress("q2vec",    &q2vec);
  intree->SetBranchAddress("q3vec",    &q3vec);
  intree->SetBranchAddress("vjet1",    &vjet1);
  intree->SetBranchAddress("vjet2",    &vjet2);
  intree->SetBranchAddress("vjet3",    &vjet3);
      
  
  
  for(unsigned int ientry=1; ientry<intree->GetEntries(); ientry++) {
    intree->GetEntry(ientry);

    if(evtType == 1){ continue; } //skips events reserved for training 
    
    if(mtop<TOPMASSLOW || mtop>TOPMASSHIGH) continue;
       
    dummy =(int)isSig;
    dummy2 =(int)b_mis;
    dummy3 =(int)w_mis;
    dummy4 =(int)wb_mis;

    topmva = readers[0]->EvaluateMVA("BDTG");
    
    outTree->Fill();

  }
  
  delete infile;
  infile=0, intree=0;
  
  outFile->Write();
  outFile->Close();
}
