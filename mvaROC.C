//================================================================================================
//
// Make ROC curves
//
//________________________________________________________________________________________________

using namespace std;

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                          // access to gROOT, entry point to ROOT system
#include <TSystem.h>                        // interface to OS
#include <TStyle.h>                         // class to handle ROOT plotting styles
#include <TFile.h>                          // file handle class
#include <TTree.h>                          // class to access ntuples
#include <TH1D.h>                           // 1D histogram class
#include <TH2D.h>                           // 2D histogram class
#include <TGraphErrors.h>                   // graph class
#include <vector>                           // STL vector class
#include <iostream>                         // standard I/O
#include <iomanip>                          // functions to format standard I/O
#include <fstream>                          // functions for file I/O
#include <string>                           // C++ string class
#include <cmath>                            // C++ math library
#include <utility>
#include <cassert>
#include <TLorentzVector.h>

#include "CPlot.hh"                         // helper class for plots
#include "KStyle.hh"                        // style settings for drawing
#endif


//=== MAIN MACRO =================================================================================================
// if  compiling through root ...
void mvaROC()
// if compiling with g++
//int  main()
{
  //--------------------------------------------------------------------------------------------------------------
  // Settings
  //==============================================================================================================
  gStyle->SetTitleOffset(1.700,"Y");
  
  
  vector<string> infilename; 
  string outputDir;
  infilename.push_back("TestingBits/testingbits_baseMVA_tt1l.root");
  infilename.push_back("TestingBits/testingbits_baseMVA+prob_tt1l.root");
  infilename.push_back("TestingBits/testingbits_baseMVA+prob+cost_tt1l.root");
  
  // infilename.push_back("TestingBits/testingbits_inclZjets.root");
  // infilename.push_back("TestingBits/testingbits_noMdrop_inclZjets.root");
  // infilename.push_back("TestingBits/testingbits_inclZjets_btag.root");
  // infilename.push_back("TestingBits/testingbits_noMdrop_inclZjets_btag.root");

  // // N-1 Plotting
  
  // infilename.push_back("TestingBits/testingbits_tt1l.root");
  // infilename.push_back("TestingBits/testingbits_noCosTS_tt1l.root");
  // infilename.push_back("TestingBits/testingbits_nodR_tt1l.root");
  // infilename.push_back("TestingBits/testingbits_nodPhi_tt1l.root");
  
  // infilename.push_back("TestingBits/testingbits_kinfit_tt1l.root");
  // infilename.push_back("TestingBits/testingbits_noQGL_tt1l.root");
  // infilename.push_back("TestingBits/testingbits_cosTS_tt1l.root");


  
  // infilename.push_back("TestingBits/testingbits_nocost_tt1l.root");
  // infilename.push_back("TestingBits/testingbits_noprob_tt1l.root");
  
  // infilename.push_back("TestingBits/testingbits_tt2lbkg.root");
  // infilename.push_back("TestingBits/testingbits_noCosTS_tt2lbkg.root");
  // infilename.push_back("TestingBits/testingbits_nodR_tt2lbkg.root");
  // infilename.push_back("TestingBits/testingbits_nodPhi_tt2lbkg.root");
  // infilename.push_back("TestingBits/testingbits_cosTS_tt2lbkg.root");
  // infilename.push_back("TestingBits/testingbits_nocost_tt2lbkg.root");
  // infilename.push_back("TestingBits/testingbits_noprob_tt2lbkg.root");

  // // N=1 Plotting
  // infilename.push_back("TestingBits/testingbits_tt1l.root");
  // infilename.push_back("TestingBits/testingbits_cosTS_tt1l.root");
  // infilename.push_back("TestingBits/testingbits_cost_tt1l.root");
  // infilename.push_back("TestingBits/testingbits_prob_tt1l.root");
  
  // infilename.push_back("TestingBits/testingbits_kinfit_cosTS_tt2lbkg.root");  
  // infilename.push_back("TestingBits/testingbits_cosTS_tt2lbkg.root");
  // infilename.push_back("TestingBits/testingbits_cost_tt2lbkg.root");
  // infilename.push_back("TestingBits/testingbits_prob_tt2lbkg.root");


  // // tt1l vs. tt2l vs. Zjets
  // infilename.push_back("TestingBits/testingbits_tt1l.root");
  // infilename.push_back("TestingBits/testingbits_kinfit_cosTS_tt2lbkg.root");
  // infilename.push_back("TestingBits/testingbits_zjetsbkg.root");
  // infilename.push_back("TestingBits/testingbits_zjetsbkg_btag.root");

  //
  
  //infilename.push_back("TestingBits/testingbits_noMdrop.root");
  // infilename.push_back("TestingBits/testingbits_noMdrop_lowpt.root");
  // infilename.push_back("testingbits_noMdrop_noDeta_lowpt.root");
  // infilename.push_back("testingbits_noMdrop_noDphi_lowpt.root");
  // infilename.push_back("testingbits_noMdrop_noDR_lowpt.root");

  // infilename.push_back("TestingBits/testingbits_noMdrop_highpt.root");
  // infilename.push_back("testingbits_noMdrop_noDeta_highpt.root");
  // infilename.push_back("testingbits_noMdrop_noDphi_highpt.root");
  // infilename.push_back("testingbits_noMdrop_noDR_highpt.root");
  
  
  // infilename.push_back("testingbits_nobmis_bkg.root");
  // infilename.push_back("TestingBits/testingbits_noQGL.root");
  // infilename.push_back("TestingBits/testingbits_noMdrop.root");
  // infilename.push_back("TestingBits/testingbits_noDeta.root");
  // infilename.push_back("TestingBits/testingbits_noDphi.root");
  // infilename.push_back("TestingBits/testingbits_noDR.root");

  //infilename.push_back("TestingBits/testingbits_onlyQGL.root");
  //infilename.push_back("TestingBits/testingbits_onlyMdrop.root");
  //infilename.push_back("TestingBits/testingbits_onlyDeta.root");
  //infilename.push_back("TestingBits/testingbits_onlyDphi.root");
  //infilename.push_back("TestingBits/testingbits_onlyDR.root");

  outputDir = "ROC";

  gSystem->mkdir(outputDir.c_str(),true);
  CPlot::sOutDir = outputDir;
   

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code
  //==============================================================================================================
  //
  // Set up output signal and background trees
  //
  
  float sig_detaj1b, sig_detaj2b, sig_dphij1b, sig_dphij2b, sig_drj1b, sig_drj2b;
  float sig_qgid1, sig_qgid2;
  float sig_bjcsv, sig_jet1csv, sig_jet2csv;
  float sig_weight;
  float sig_prob, sig_cost;
  float sig_cosTS, sig_cosTS1, sig_cosTS2;
  float sig_mva, sig_mtop, sig_mw;
  float sig_genTopPt;
  unsigned int sig_sample;
  
  unsigned int bkg_type;
  unsigned int bkg_bmis, bkg_wmis, bkg_wbmis;
  float bkg_detaj1b, bkg_detaj2b, bkg_dphij1b, bkg_dphij2b, bkg_drj1b, bkg_drj2b;
  float bkg_qgid1, bkg_qgid2;
  float bkg_bjcsv, bkg_jet1csv, bkg_jet2csv;
  float bkg_weight;
  float bkg_prob, bkg_cost;
  float bkg_cosTS, bkg_cosTS1, bkg_cosTS2;
  float bkg_mva, bkg_mtop, bkg_mw;

  unsigned int bkg_sample;
  
  TFile *outFile = new TFile("SigBkgROC.root","RECREATE");
  TTree *sigTree = new TTree("sigEvents","sigEvents");
  TTree *bkgTree = new TTree("bkgEvents","bkgEvents");
  
  sigTree->Branch("sig_qgid1",  &sig_qgid1,  "sig_qgid1/F");
  sigTree->Branch("sig_qgid2",  &sig_qgid2,  "sig_qgid2/F");
  sigTree->Branch("sig_detaj1b",&sig_detaj1b,"sig_detaj1b/F");
  sigTree->Branch("sig_detaj2b",&sig_detaj2b,"sig_detaj2b/F");
  sigTree->Branch("sig_dphij1b",&sig_dphij1b,"sig_dphij1b/F");
  sigTree->Branch("sig_dphij2b",&sig_dphij2b,"sig_dphij2b/F");
  sigTree->Branch("sig_drj1b",  &sig_drj1b,  "sig_drj1b/F");
  sigTree->Branch("sig_drj2b",  &sig_drj2b,  "sig_drj2b/F");
  sigTree->Branch("sig_bjcsv",    &sig_bjcsv,    "sig_bjcsv/F");
  sigTree->Branch("sig_jet1csv",  &sig_jet1csv,  "sig_jet1csv/F");
  sigTree->Branch("sig_jet2csv",  &sig_jet2csv,  "sig_jet2csv/F");
  sigTree->Branch("sig_weight",   &sig_weight,   "sig_weight/F");
  sigTree->Branch("sig_prob",     &sig_prob,     "sig_prob/F");
  sigTree->Branch("sig_cost",     &sig_cost,     "sig_cost/F");
  // sigTree->Branch("sig_cosTS",     &sig_cosTS,     "sig_cosTS/F");
  // sigTree->Branch("sig_cosTS1",     &sig_cosTS1,     "sig_cosTS1/F");
  // sigTree->Branch("sig_cosTS2",     &sig_cosTS2,     "sig_cosTS2/F");
  sigTree->Branch("sig_mva",      &sig_mva,      "sig_mva/F");
  sigTree->Branch("sig_mtop",      &sig_mtop,      "sig_mtop/F");
  sigTree->Branch("sig_mw",      &sig_mw,      "sig_mw/F");
  sigTree->Branch("sig_genTopPt", &sig_genTopPt, "sig_genTopPt/F");
  sigTree->Branch("sig_sample",   &sig_sample,   "sig_sample/i");
  
  
  bkgTree->Branch("bkg_bmis",   &bkg_bmis,   "bkg_bmis/i");
  bkgTree->Branch("bkg_wmis",   &bkg_wmis,   "bkg_wmis/i");
  bkgTree->Branch("bkg_wbmis",  &bkg_wbmis,  "bkg_wbmis/i");
  bkgTree->Branch("bkg_qgid1",  &bkg_qgid1,  "bkg_qgid1/F");
  bkgTree->Branch("bkg_qgid2",  &bkg_qgid2,  "bkg_qgid2/F");
  bkgTree->Branch("bkg_detaj1b",&bkg_detaj1b,"bkg_detaj1b/F");
  bkgTree->Branch("bkg_detaj2b",&bkg_detaj2b,"bkg_detaj2b/F");
  bkgTree->Branch("bkg_dphij1b",&bkg_dphij1b,"bkg_dphij1b/F");
  bkgTree->Branch("bkg_dphij2b",&bkg_dphij2b,"bkg_dphij2b/F");
  bkgTree->Branch("bkg_drj1b",  &bkg_drj1b,  "bkg_drj1b/F");
  bkgTree->Branch("bkg_drj2b",  &bkg_drj2b,  "bkg_drj2b/F");
  bkgTree->Branch("bkg_bjcsv",    &bkg_bjcsv,    "bkg_bjcsv/F");
  bkgTree->Branch("bkg_jet1csv",  &bkg_jet1csv,  "bkg_jet1csv/F");
  bkgTree->Branch("bkg_jet2csv",  &bkg_jet2csv,  "bkg_jet2csv/F");
  bkgTree->Branch("bkg_weight",   &bkg_weight,   "bkg_weight/F");
  bkgTree->Branch("bkg_prob",     &bkg_prob,     "bkg_prob/F");
  bkgTree->Branch("bkg_cost",     &bkg_cost,     "bkg_cost/F");
  // bkgTree->Branch("bkg_cosTS",     &bkg_cosTS,     "bkg_cosTS/F");
  // bkgTree->Branch("bkg_cosTS1",   &bkg_cosTS1,   "bkg_cosTS1/F");
  // bkgTree->Branch("bkg_cosTS2",   &bkg_cosTS2,   "bkg_cosTS2/F");
  bkgTree->Branch("bkg_type",     &bkg_type,     "bkg_type/i");
  bkgTree->Branch("bkg_mva",      &bkg_mva,      "bkg_mva/F");
  bkgTree->Branch("bkg_mtop",     &bkg_mtop,     "bkg_mtop/F");
  bkgTree->Branch("bkg_mw",       &bkg_mw,       "bkg_mw/F");
  bkgTree->Branch("bkg_sample",   &bkg_sample,   "bkg_sample/i");  
  
  //Declare histograms
  char hname[100];

 			   
  vector <TH1D*> htopmvav, hmdropv, hdetaj1bv, hdetaj2bv,hdphij1bv,hdphij2bv, hdrj1bv,hdrj2bv, hqgid1v, hqgid2v, hbjcsv, hj1csv, hj2csv, htopmvaBkg;
  for(unsigned int icombo=0; icombo<9; icombo++){
    sprintf(hname,"hbjcsv_%i",icombo);   hbjcsv.push_back(new TH1D(hname,"",50,0,1));      hbjcsv[icombo]->Sumw2();
    sprintf(hname,"hj1csv_%i",icombo);   hj1csv.push_back(new TH1D(hname,"",50,0,1));      hj1csv[icombo]->Sumw2();
    sprintf(hname,"hj2csv_%i",icombo);   hj2csv.push_back(new TH1D(hname,"",50,0,1));      hj2csv[icombo]->Sumw2();
    sprintf(hname,"hqgid1_%i",icombo);   hqgid1v.push_back(new TH1D(hname,"",50,0,1));     hqgid1v[icombo]->Sumw2();
    sprintf(hname,"hqgid2_%i",icombo);   hqgid2v.push_back(new TH1D(hname,"",50,0,1));     hqgid2v[icombo]->Sumw2();
    sprintf(hname,"htopmva_%i",icombo);  htopmvav.push_back(new TH1D(hname,"",50,-1,1));   htopmvav[icombo]->Sumw2();
    sprintf(hname,"hmdrop_%i",icombo);   hmdropv.push_back(new TH1D(hname,"",50,0,1));     hmdropv[icombo]->Sumw2();
    sprintf(hname,"hdetaj1b_%i",icombo); hdetaj1bv.push_back(new TH1D(hname,"",50,0,4));   hdetaj1bv[icombo]->Sumw2();
    sprintf(hname,"hdetaj2b_%i",icombo); hdetaj2bv.push_back(new TH1D(hname,"",50,0,4));   hdetaj2bv[icombo]->Sumw2();
    sprintf(hname,"hdphij1b_%i",icombo); hdphij1bv.push_back(new TH1D(hname,"",40,0,3.2));   hdphij1bv[icombo]->Sumw2();
    sprintf(hname,"hdphij2b_%i",icombo); hdphij2bv.push_back(new TH1D(hname,"",40,0,3.2));   hdphij2bv[icombo]->Sumw2();
    sprintf(hname,"hdrj1b_%i",icombo);   hdrj1bv.push_back(new TH1D(hname,"",50,0.8,4));   hdrj1bv[icombo]->Sumw2();
    sprintf(hname,"hdrj2b_%i",icombo);   hdrj2bv.push_back(new TH1D(hname,"",50,0.8,4));   hdrj2bv[icombo]->Sumw2();
    sprintf(hname,"htopmvaBkg_%i",icombo); htopmvaBkg.push_back(new TH1D(hname,"",50,-1,1)); htopmvaBkg[icombo]->Sumw2();
  }

  vector<string> labelsv;
  vector<int> colorsv, linesv;
  // labelsv.push_back("Z(#nu#nu)+jets bkg (with Mdrop)");         colorsv.push_back(kAzure-8);    linesv.push_back(7);
  // labelsv.push_back("Z(#nu#nu)+jets bkg (no Mdrop)");           colorsv.push_back(kViolet+9);   linesv.push_back(7);
  // labelsv.push_back("Z(#nu#nu)+jets & b-tag bkg (with Mdrop)"); colorsv.push_back(kOrange+7);        linesv.push_back(7);
  // labelsv.push_back("Z(#nu#nu)+jets & b-tag bkg (no Mdrop)");   colorsv.push_back(kRed);        linesv.push_back(1);

  // // N-1 Plotting (including kinematic fitter variables and cos theta star)
  labelsv.push_back("t#bar{t}(1l) base MVA");                  colorsv.push_back(kGreen+2);   linesv.push_back(1);
  labelsv.push_back("t#bar{t}(1l) base MVA+prob");             colorsv.push_back(kPink+4);   linesv.push_back(1);
  labelsv.push_back("t#bar{t}(1l) base MVA+prob+cost");        colorsv.push_back(kBlue-4);   linesv.push_back(1);
  // labelsv.push_back("t#bar{t} 1L [no cos(#theta*)]");        colorsv.push_back(kGreen-3);   linesv.push_back(7);
  // labelsv.push_back("t#bar{t} 1L [no QG])");        colorsv.push_back(kGreen-3);   linesv.push_back(4);
  // labelsv.push_back("t#bar{t} 1L [no kin fit]");        colorsv.push_back(kGreen-3);   linesv.push_back(9);

   // labelsv.push_back("t#bar{t} 1L (no probability)"); colorsv.push_back(kGreen-3);   linesv.push_back(9);
  
  
  // labelsv.push_back("t#bar{t} 2L");                  colorsv.push_back(kPink+4);    linesv.push_back(1);
  // labelsv.push_back("t#bar{t} 2L [no cos(#theta*)]");        colorsv.push_back(kPink+9);    linesv.push_back(7);
  // labelsv.push_back("t#bar{t} 2L [no #DeltaR]");        colorsv.push_back(kPink+9);    linesv.push_back(7);
  // labelsv.push_back("t#bar{t} 2L [no #Delta#phi]");        colorsv.push_back(kPink+9);    linesv.push_back(7);
  
  // labelsv.push_back("t#bar{t} 2L [no QGL]");        colorsv.push_back(kPink+9);    linesv.push_back(4);
  // labelsv.push_back("t#bar{t} 2L [no kin fit]");        colorsv.push_back(kPink+9);   linesv.push_back(9);
   // labelsv.push_back("t#bar{t} 2L (no cos#theta*)");  colorsv.push_back(kPink+9);    linesv.push_back(4);
  
  // labelsv.push_back("t#bar{t} 2L (no probability)"); colorsv.push_back(kPink+9);    linesv.push_back(9);
  
  // // N=1 Plotting (kin. fitter vars & cos theta star)
  // labelsv.push_back("t#bar{t} 1L");                  colorsv.push_back(kGreen+2);   linesv.push_back(1);
  // labelsv.push_back("t#bar{t} 1L [cos(#theta*)]");          colorsv.push_back(kGreen-3);   linesv.push_back(4);
  // labelsv.push_back("t#bar{t} 1L [cost]");                  colorsv.push_back(kGreen-3);   linesv.push_back(7);
  // labelsv.push_back("t#bar{t} 1L [probability]");           colorsv.push_back(kGreen-3);   linesv.push_back(9);

  // labelsv.push_back("t#bar{t} 2L");                  colorsv.push_back(kPink+4);    linesv.push_back(1);
  // labelsv.push_back("t#bar{t} 2L [cos(#theta*)]");          colorsv.push_back(kPink+9);   linesv.push_back(4);
  // labelsv.push_back("t#bar{t} 2L [cost]");                  colorsv.push_back(kPink+9);   linesv.push_back(7);
  // labelsv.push_back("t#bar{t} 2L [probability]");           colorsv.push_back(kPink+9);   linesv.push_back(9);

  // // tt1l vs tt2l vs Zjets
  // labelsv.push_back("t#bar{t} semi-lepton");  colorsv.push_back(kViolet+2); linesv.push_back(7);
  // labelsv.push_back("t#bar{t} di-lepton");  colorsv.push_back(kPink+9);   linesv.push_back(7);
  // labelsv.push_back("Z(#nu#nu) + jets"); colorsv.push_back(kAzure+8);  linesv.push_back(7);
  // labelsv.push_back("Z(#nu#nu) + jets [b-tag requirement]"); colorsv.push_back(kBlue-4);  linesv.push_back(7);
  
  // labelsv.push_back("no qgl");                colorsv.push_back(kViolet+2);       linesv.push_back(7);
  // labelsv.push_back("no #zeta");              colorsv.push_back(kGreen+2);        linesv.push_back(7);
 
  //labelsv.push_back("QGL+#Delta#eta+#Delta#phi+#DeltaR (Low p_{T})");  colorsv.push_back(kBlue-6);     linesv.push_back(7);
  // labelsv.push_back("QGL+#Delta#phi+#DeltaR (Low p_{T})");        colorsv.push_back(kOrange+7);       linesv.push_back(1);
  // labelsv.push_back("QGL+#Delta#eta+#DeltaR (Low p_{T})");        colorsv.push_back(kPink+9);        linesv.push_back(1);
  // labelsv.push_back("QGL+#Delta#eta+#Delta#phi (Low p_{T})");     colorsv.push_back(kAzure+8);         linesv.push_back(1);
  // labelsv.push_back("QGL+#Delta#eta+#Delta#phi+#DeltaR (High p_{T})");  colorsv.push_back(kCyan-6);     linesv.push_back(7);
  // labelsv.push_back("QGL+#Delta#phi+#DeltaR (High p_{T})");        colorsv.push_back(kOrange+7);       linesv.push_back(7);
  // labelsv.push_back("QGL+#Delta#eta+#DeltaR (High p_{T})");        colorsv.push_back(kPink+9);        linesv.push_back(7);
  // labelsv.push_back("QGL+#Delta#eta+#Delta#phi (High p_{T})");     colorsv.push_back(kAzure+8);         linesv.push_back(7);
  // labelsv.push_back("QGL+#Delta #eta");        colorsv.push_back(kOrange+7);       linesv.push_back(7);
  // labelsv.push_back("QGL+#Delta #phi");        colorsv.push_back(kAzure+1);        linesv.push_back(7);
  // labelsv.push_back("QGL+#Delta R");           colorsv.push_back(kBlue-4);         linesv.push_back(1);

  
 
  // labelsv.push_back("all vars (w&b mis bkg)"); colorsv.push_back(kViolet-5);     linesv.push_back(1);
  // labelsv.push_back("- qgl (w&b mis bkg)");    colorsv.push_back(kGray+3);       linesv.push_back(1);
  // labelsv.push_back("only qgl");                 colorsv.push_back(kAzure-8);      linesv.push_back(7);
  // labelsv.push_back("only #zeta");               colorsv.push_back(kRed);         linesv.push_back(7);
  // labelsv.push_back("only #Delta#eta");         colorsv.push_back(kOrange-3);    linesv.push_back(7);
  // labelsv.push_back("only #Delta#phi");         colorsv.push_back(kGreen);      linesv.push_back(7);
  // labelsv.push_back("only #Delta R");            colorsv.push_back(kPink+9);      linesv.push_back(7);
  //labelsv.push_back("all variables");                 colorsv.push_back(kViolet+2);     linesv.push_back(1);

  vector<float> mvaWP;
  mvaWP.push_back(-0.45);//0.42); // 60% sig eff for baseMVA
  mvaWP.push_back(-0.40);//0.40); // 60% sig eff for base+prob
  mvaWP.push_back(-0.42);//0.42); // 60% sig eff for base+prob+cost
 
  vector<double> cutsv;
  for(unsigned int i=0; i<21; i++) {
    cutsv.push_back(-1. + 2.0*i/20.);
  }
  vector<double> sigpass[labelsv.size()];
  vector<double> sigfail[labelsv.size()];
  vector<double> bkgpass[labelsv.size()];
  vector<double> bkgfail[labelsv.size()];
  
  vector <double> Npass[labelsv.size()];
  vector <double> Nfail[labelsv.size()];
  
  vector <double> ptlow, pthigh;
  ptlow.push_back(0);   pthigh.push_back(20);
  ptlow.push_back(20);  pthigh.push_back(40);
  ptlow.push_back(40);  pthigh.push_back(60);
  ptlow.push_back(60);  pthigh.push_back(80);
  ptlow.push_back(80);  pthigh.push_back(100);
  ptlow.push_back(100); pthigh.push_back(120);
  ptlow.push_back(120); pthigh.push_back(150);
  ptlow.push_back(150); pthigh.push_back(180);
  ptlow.push_back(180); pthigh.push_back(220);
  ptlow.push_back(220); pthigh.push_back(260);
  ptlow.push_back(260); pthigh.push_back(300);
  ptlow.push_back(300); pthigh.push_back(340);
  ptlow.push_back(340); pthigh.push_back(400);
    
  for(unsigned int iw=0; iw<labelsv.size(); iw++) {
    for(unsigned int ic=0; ic<cutsv.size(); ic++) {
      sigpass[iw].push_back(0);
      sigfail[iw].push_back(0);
      bkgpass[iw].push_back(0);
      bkgfail[iw].push_back(0);
    }

    for(unsigned int ic=0; ic<ptlow.size(); ic++) {
      Npass[iw].push_back(0);
      Nfail[iw].push_back(0);
    }
  }

  
  vector<TGraphErrors*> TopEff;
  vector<TGraphErrors*> rocsv;

  float sigmva[labelsv.size()];
  float bkgmva[labelsv.size()];
  
  for(unsigned int ifile=0; ifile < infilename.size(); ifile++){
    
    cout << "file: " << ifile << endl;
    //
    //input file vars
    //
    bool isQGL = false, isMdrop = false, isDeta = false, isDphi = false, isDR = false;
    unsigned int isSig, b_mis, w_mis, wb_mis;
    unsigned int bkgType, eventNum;
    float weight;
    float mtop, detaj1b, detaj2b, dphij1b, dphij2b, drj1b, drj2b;
    float bjcsv, jet1csv, jet2csv;
    float qgid1, qgid2, mdrop;
    float cost, prob, cosTS, cosTS1, cosTS2;
    float topmva;
    TLorentzVector *vjet1=0, *vjet2=0, *vjet3=0;
    TLorentzVector *q1vec=0, *q2vec=0, *q3vec=0;
 
        
    std::cout << "Processing " << infilename[ifile] << "..." << std::endl;
    
    TFile *infile = TFile::Open(infilename[ifile].c_str()); assert(infile);
    TTree *intree = (TTree*)infile->Get("Events");   assert(intree);

   
    intree->SetBranchAddress("eventNum", &eventNum);
   
    if(intree->GetBranchStatus("bkgType")==1) {
      intree->SetBranchAddress("bkgType", &bkgType);
    }
    intree->SetBranchAddress("isSig",   &isSig);
    intree->SetBranchAddress("b_mis",   &b_mis);
    intree->SetBranchAddress("w_mis",   &w_mis);
    intree->SetBranchAddress("wb_mis",  &wb_mis);
    intree->SetBranchAddress("topmva",  &topmva); 
    intree->SetBranchAddress("bjcsv",   &bjcsv);
    intree->SetBranchAddress("jet1csv", &jet1csv);
    intree->SetBranchAddress("jet2csv", &jet2csv);
    intree->SetBranchAddress("weight",  &weight);

    
    if(intree->GetBranchStatus("qgid1")==1 && intree->GetBranchStatus("qgid2")==1) {
      isQGL=true;
      intree->SetBranchAddress("qgid1",   &qgid1);
      intree->SetBranchAddress("qgid2",   &qgid2);
    }
    if(intree->GetBranchStatus("mdrop")==1){
      isMdrop=true;
      intree->SetBranchAddress("mdrop",   &mdrop);	
    }
    if(intree->GetBranchStatus("detaj1b")==1 && intree->GetBranchStatus("detaj2b")==1){
      isDeta=true;
      intree->SetBranchAddress("detaj1b", &detaj1b);
      intree->SetBranchAddress("detaj2b", &detaj2b);
    }
    if(intree->GetBranchStatus("dphij1b")==1 && intree->GetBranchStatus("dphij2b")==1){
      isDphi=true;
      intree->SetBranchAddress("dphij1b", &dphij1b);
      intree->SetBranchAddress("dphij2b", &dphij2b);
    }
    if(intree->GetBranchStatus("drj1b")==1 && intree->GetBranchStatus("drj2b")==1){
      isDR=true;
      intree->SetBranchAddress("drj1b",   &drj1b);
      intree->SetBranchAddress("drj2b",   &drj2b);
    }
    if(intree->GetBranchStatus("prob")==1){
      intree->SetBranchAddress("prob",  &prob);
    }
    if(intree->GetBranchStatus("cost")==1){
      intree->SetBranchAddress("cost",  &cost);
    }
    if(intree->GetBranchStatus("cosTS")==1){
      intree->SetBranchAddress("cosTS", &cosTS);
    }
    if(intree->GetBranchStatus("cosTS1")==1){
      intree->SetBranchAddress("cosTS1", &cosTS1);
      intree->SetBranchAddress("cosTS2", &cosTS2);
    }
    if(intree->GetBranchStatus("vjet1")==1){
      intree->SetBranchAddress("vjet1", &vjet1);
      intree->SetBranchAddress("vjet2", &vjet2);
      intree->SetBranchAddress("vjet3", &vjet3);
    }
    if(intree->GetBranchStatus("q1vec")==1){
      intree->SetBranchAddress("q1vec", &q1vec);
      intree->SetBranchAddress("q2vec", &q2vec);
      intree->SetBranchAddress("q3vec", &q3vec);
    }
    vector<float> comboMVA;
    vector < std::pair <int,int> > entryNum;
    
    for(unsigned int ientry=0; ientry<intree->GetEntries(); ientry++) {

      
      unsigned int evt2Num=0;
      unsigned int evt1Num=0;

      intree->GetEntry(ientry+1);
      evt2Num = eventNum;

      intree->GetEntry(ientry);
      evt1Num = eventNum;
          
      comboMVA.push_back(topmva);
      entryNum.push_back(make_pair(ientry,isSig));
           
      if(evt2Num != evt1Num){
	
	// cout << "-------------- new event -------------" << endl;
	float tmpMVA         = -999;
	int   tmpSigEntryNum = -999; 
	int   tmpBkgEntryNum = -999;
	float m_jjj=-999, m_jj=-999;
	float pt_jjj=0;
	
	for(unsigned int i=0; i<comboMVA.size(); ++i){
	  
	  if(entryNum[i].second==1){// signal 

	    sigmva[ifile]  = comboMVA[i];
	    tmpSigEntryNum = entryNum[i].first; 

	  }else{                    // background
	  	    
	    bool next = false;
	    if(ifile>3){
	      if(bkgType==3) next=true;
	    }
	    if (next) continue;
	    
	    if(comboMVA[i]>tmpMVA){
	      tmpMVA = comboMVA[i];
	      tmpBkgEntryNum = entryNum[i].first;
	    }
	  }
	}
	
	if(tmpSigEntryNum != -999){
	  intree->GetEntry(tmpSigEntryNum);
	  TLorentzVector  j1 = *vjet1;
	  TLorentzVector  j2 = *vjet2;
	  TLorentzVector  j3 = *vjet3;

	  TLorentzVector q1 = *q1vec; TLorentzVector q2 = *q2vec; TLorentzVector q3 = *q3vec;

	  pt_jjj = (q1+q2+q3).Pt();
	  
	  m_jj  = (j1+j2).M();
	  m_jjj = (j1+j2+j3).M();
	  	  
	  sig_detaj1b = detaj1b;
	  sig_detaj2b = detaj2b;
	  sig_dphij1b = dphij1b;
	  sig_dphij2b = dphij2b;
	  sig_drj1b   = drj1b;
	  sig_drj2b   = drj2b;
	  sig_qgid1   = qgid1;
	  sig_qgid2   = qgid2;
	  sig_bjcsv   = bjcsv;
	  sig_jet1csv = jet1csv;
	  sig_jet2csv = jet2csv;
	  sig_weight  = weight;
	  sig_prob    = prob;
	  sig_cost    = cost;
	  // sig_cosTS   = cosTS;
	  // sig_cosTS1   = cosTS1;
	  // sig_cosTS2   = cosTS2;
	  sig_mva     = sigmva[ifile];
	  sig_mtop    = m_jjj;
	  sig_mw      = m_jj;
	  sig_genTopPt = pt_jjj;
	  sig_sample  = ifile;

	  sigTree     ->Fill();

	  // cout << "signal event mva score: " << sig_mva << endl;
	}
	
	if (tmpBkgEntryNum != -999) {
	  intree->GetEntry(tmpBkgEntryNum);
	  TLorentzVector  j1 = *vjet1;
	  TLorentzVector  j2 = *vjet2;
	  TLorentzVector  j3 = *vjet3;
	  
	  m_jj  = (j1+j2).M();
	  m_jjj = (j1+j2+j3).M();
	 	  
	  bkgmva[ifile] = tmpMVA;
	  bkg_bmis      = b_mis;
	  bkg_wmis      = w_mis;
	  bkg_wbmis     = wb_mis;
	  bkg_detaj1b = detaj1b;
	  bkg_detaj2b = detaj2b;
	  bkg_dphij1b = dphij1b;
	  bkg_dphij2b = dphij2b;
	  bkg_drj1b   = drj1b;
	  bkg_drj2b   = drj2b;
	  bkg_qgid1   = qgid1;
	  bkg_qgid2   = qgid2;
	  bkg_bjcsv   = bjcsv;
	  bkg_jet1csv = jet1csv;
	  bkg_jet2csv = jet2csv;
	  bkg_weight  = weight;
	  bkg_prob    = prob;
	  bkg_cost    = cost;
	  // bkg_cosTS   = cosTS;
	  // bkg_cosTS1  = cosTS1;
	  // bkg_cosTS2  = cosTS2;
	  bkg_type    = bkgType;
	  bkg_mva     = bkgmva[ifile];
	  bkg_mtop    = m_jjj;
	  bkg_mw      = m_jj;
	  bkg_sample  = ifile;

	  bkgTree     ->Fill();
	  
	  // cout << "background event mva score: " << bkg_mva << endl;
      }
	
	comboMVA.clear();
	entryNum.clear();
	
	
    }
	
      
      if(isSig == 1){
      	htopmvav[0]  ->Fill(topmva);
      	hbjcsv[0]    ->Fill(bjcsv);
      	hj1csv[0]    ->Fill(jet1csv);
      	hj2csv[0]    ->Fill(jet2csv);
      	if(isQGL){
      	  hqgid1v[0]   ->Fill(qgid1);
      	  hqgid2v[0]   ->Fill(qgid2);
      	}
      	if(isMdrop){
      	  hmdropv[0]   ->Fill(mdrop);
      	}
      	if(isDeta){
      	  hdetaj1bv[0] ->Fill(detaj1b);
      	  hdetaj2bv[0] ->Fill(detaj2b);
      	}
      	if(isDphi){
      	  hdphij1bv[0] ->Fill(dphij1b);
      	  hdphij2bv[0] ->Fill(dphij2b);
      	}
      	if(isDR){
      	  hdrj1bv[0]   ->Fill(drj1b);
      	  hdrj2bv[0]   ->Fill(drj2b);
      	}
      }
      else if(b_mis == 1){
      	htopmvav[1]  ->Fill(topmva);
      	hbjcsv[1]    ->Fill(bjcsv);
      	hj1csv[1]    ->Fill(jet1csv);
      	hj2csv[1]    ->Fill(jet2csv);

      	if(isQGL){
      	  hqgid1v[1]   ->Fill(qgid1);
      	  hqgid2v[1]   ->Fill(qgid2);
      	}
      	if(isMdrop){
      	  hmdropv[1]   ->Fill(mdrop);
      	}
      	if(isDeta){
      	  hdetaj1bv[1] ->Fill(detaj1b);
      	  hdetaj2bv[1] ->Fill(detaj2b);
      	}
      	if(isDphi){
      	  hdphij1bv[1] ->Fill(dphij1b);
      	  hdphij2bv[1] ->Fill(dphij2b);
      	}
      	if(isDR){
      	  hdrj1bv[1]   ->Fill(drj1b);
      	  hdrj2bv[1]   ->Fill(drj2b);
      	}
    	
      }
      else if(w_mis == 1){
      	htopmvav[2]  ->Fill(topmva);
      	hbjcsv[2]    ->Fill(bjcsv);
      	hj1csv[2]    ->Fill(jet1csv);
      	hj2csv[2]    ->Fill(jet2csv);

      	if(isQGL){
      	  hqgid1v[2]   ->Fill(qgid1);
      	  hqgid2v[2]   ->Fill(qgid2);
      	}
      	if(isMdrop){
      	  hmdropv[2]   ->Fill(mdrop);
      	}
      	if(isDeta){
      	  hdetaj1bv[2] ->Fill(detaj1b);
      	  hdetaj2bv[2] ->Fill(detaj2b);
      	}
      	if(isDphi){
      	  hdphij1bv[2] ->Fill(dphij1b);
      	  hdphij2bv[2] ->Fill(dphij2b);
      	}
      	if(isDR){
      	  hdrj1bv[2]   ->Fill(drj1b);
      	  hdrj2bv[2]   ->Fill(drj2b);
      	}	
      }
      else if(wb_mis == 1){
      	htopmvav[3]  ->Fill(topmva);
      	hbjcsv[3]    ->Fill(bjcsv);
      	hj1csv[3]    ->Fill(jet1csv);
      	hj2csv[3]    ->Fill(jet2csv);

      	if(isQGL){
      	  hqgid1v[3]   ->Fill(qgid1);
      	  hqgid2v[3]   ->Fill(qgid2);
      	}
      	if(isMdrop){
      	  hmdropv[3]   ->Fill(mdrop);
      	}
      	if(isDeta){
      	  hdetaj1bv[3] ->Fill(detaj1b);
      	  hdetaj2bv[3] ->Fill(detaj2b);
      	}
      	if(isDphi){
      	  hdphij1bv[3] ->Fill(dphij1b);
      	  hdphij2bv[3] ->Fill(dphij2b);
      	}
      	if(isDR){
      	  hdrj1bv[3]   ->Fill(drj1b);
      	  hdrj2bv[3]   ->Fill(drj2b);
      	}		
      }

      htopmvaBkg[0] ->Add(htopmvav[1],htopmvav[2]);
      htopmvaBkg[0] ->Add(htopmvav[3]);
    }

    gStyle->SetPalette(1);
    TCanvas *c = MakeCanvas("c","c",800,800);
    c->SetTickx(1);
    c->SetTicky(1);
    
    char ylabel[100];
    float norm_topmva_sig;
    norm_topmva_sig = 1.0/htopmvav[0]->Integral();
    htopmvav[0]->Scale(norm_topmva_sig);
    float norm_topmva_bkg;
    norm_topmva_bkg = 1.0/htopmvaBkg[0]->Integral();//(htopmvav[1]->Integral() + htopmvav[2]->Integral() + htopmvav[3]->Integral()); 
    vector <float> norm_j1csv;
    for(unsigned int k=1; k < 4; ++k){
      norm_j1csv.push_back(1.0/hj1csv[k]->Integral());
      hj1csv[k]->Scale(norm_j1csv[k]);
      htopmvav[k]->Scale(norm_topmva_bkg);
    }
      
      
      sprintf(ylabel,"Fraction / %.2f",htopmvav[ifile]->GetBinWidth(1));
      sprintf(hname,"topmva_%i",ifile);
      CPlot plottopmva(hname,"Top MVA","MVA score",ylabel);
      plottopmva.AddHist1D(htopmvav[0],"Signal Top","hist",kRed+2,1,3004);
      plottopmva.AddToStack(htopmvav[1],"B mismatch",kMagenta-6,kMagenta-6);
      plottopmva.AddToStack(htopmvav[2],"W mismatch",kCyan-6,kCyan-6);
      plottopmva.AddToStack(htopmvav[3],"B & W mismatch",kBlue-6,kBlue-6);
      plottopmva.TransLegend(-0.35,+0.02);
      plottopmva.Draw(c,true,"png");

      sprintf(ylabel,"Fraction / %.2f", hj1csv[ifile]->GetBinWidth(1));
      sprintf(hname,"j1csv_%i",ifile);
      CPlot plotj1csv(hname,"jet1 csv","CSVv2+IVF",ylabel);
      plotj1csv.AddHist1D(hj1csv[0],"Signal Top","hist",kRed+2,1,3004);
      plotj1csv.AddToStack(hj1csv[1],"B mismatch",kMagenta-6,kMagenta-6);
      plotj1csv.AddToStack(hj1csv[2],"W mismatch",kCyan-6,kCyan-6);
      plotj1csv.AddToStack(hj1csv[3],"B & W mismatch",kBlue-6,kBlue-6);
      plotj1csv.TransLegend(-0.35,+0.02);
      plotj1csv.Draw(c,true,"png");

      delete infile;
      infile=0, intree=0;
  }
  outFile->Write();
  outFile->Close();
  //
  // Read in signal and background trees for ROC curve
  //
  vector<TH1D*> hSigTopPt, hSigMtop, hSigMw, hBkgMtop, hBkgMw;
  for(unsigned int it=0; it<labelsv.size(); it++){
     sprintf(hname,"gentoppt_%i",it);  hSigTopPt.push_back(new TH1D(hname,"",50,0,400)); hSigTopPt[it]->Sumw2();
     sprintf(hname,"SigMtop_%i",it);hSigMtop.push_back(new TH1D(hname,"",50,100,400)); hSigMtop[it]->Sumw2();
     sprintf(hname,"SigMw_%i",it); hSigMw.push_back(new TH1D(hname,"",50,0,250));       hSigMw[it]->Sumw2();
     sprintf(hname,"BkgMtop_%i",it); hBkgMtop.push_back(new TH1D(hname,"",50,0,400));   hBkgMtop[it]->Sumw2();
     sprintf(hname,"BkgMw_%i",it);hBkgMw.push_back(new TH1D(hname,"",50,0,250));       hBkgMw[it]->Sumw2();
  }

  TFile *inSigBkgfile = TFile::Open("SigBkgROC.root");         assert(inSigBkgfile);
  TTree *inSigtree = (TTree*)inSigBkgfile->Get("sigEvents");   assert(inSigtree);
  TTree *inBkgtree = (TTree*)inSigBkgfile->Get("bkgEvents");   assert(inBkgtree);
  
  inSigtree->SetBranchAddress("sig_qgid1",  &sig_qgid1);
  inSigtree->SetBranchAddress("sig_qgid2",  &sig_qgid2);
  inSigtree->SetBranchAddress("sig_detaj1b",&sig_detaj1b);
  inSigtree->SetBranchAddress("sig_detaj2b",&sig_detaj2b);
  inSigtree->SetBranchAddress("sig_dphij1b",&sig_dphij1b);
  inSigtree->SetBranchAddress("sig_dphij2b",&sig_dphij2b);
  inSigtree->SetBranchAddress("sig_drj1b",  &sig_drj1b);
  inSigtree->SetBranchAddress("sig_drj2b",  &sig_drj2b);
  inSigtree->SetBranchAddress("sig_bjcsv",    &sig_bjcsv);
  inSigtree->SetBranchAddress("sig_jet1csv",  &sig_jet1csv);
  inSigtree->SetBranchAddress("sig_jet2csv",  &sig_jet2csv);
  inSigtree->SetBranchAddress("sig_weight",   &sig_weight);
  inSigtree->SetBranchAddress("sig_prob",     &sig_prob);
  inSigtree->SetBranchAddress("sig_cost",     &sig_cost);
  // inSigtree->SetBranchAddress("sig_cosTS",    &sig_cosTS);
  // inSigtree->SetBranchAddress("sig_cosTS1",    &sig_cosTS1);
  // inSigtree->SetBranchAddress("sig_cosTS2",    &sig_cosTS2);
  inSigtree->SetBranchAddress("sig_mtop",     &sig_mtop);
  inSigtree->SetBranchAddress("sig_mw",       &sig_mw);
  inSigtree->SetBranchAddress("sig_mva",      &sig_mva);
  inSigtree->SetBranchAddress("sig_genTopPt", &sig_genTopPt);
  inSigtree->SetBranchAddress("sig_sample",   &sig_sample);
  
  inBkgtree->SetBranchAddress("bkg_bmis",  &bkg_bmis);
  inBkgtree->SetBranchAddress("bkg_wmis",  &bkg_wmis);
  inBkgtree->SetBranchAddress("bkg_wbmis",&bkg_wbmis);
  inBkgtree->SetBranchAddress("bkg_qgid1",  &bkg_qgid1);
  inBkgtree->SetBranchAddress("bkg_qgid2",  &bkg_qgid2);
  inBkgtree->SetBranchAddress("bkg_detaj1b",&bkg_detaj1b);
  inBkgtree->SetBranchAddress("bkg_detaj2b",&bkg_detaj2b);
  inBkgtree->SetBranchAddress("bkg_dphij1b",&bkg_dphij1b);
  inBkgtree->SetBranchAddress("bkg_dphij2b",&bkg_dphij2b);
  inBkgtree->SetBranchAddress("bkg_drj1b",  &bkg_drj1b);
  inBkgtree->SetBranchAddress("bkg_drj2b",  &bkg_drj2b);
  inBkgtree->SetBranchAddress("bkg_bjcsv",   &bkg_bjcsv);
  inBkgtree->SetBranchAddress("bkg_jet1csv", &bkg_jet1csv);
  inBkgtree->SetBranchAddress("bkg_jet2csv", &bkg_jet2csv);
  inBkgtree->SetBranchAddress("bkg_weight",  &bkg_weight);
  inBkgtree->SetBranchAddress("bkg_type",    &bkg_type);
  inBkgtree->SetBranchAddress("bkg_prob",     &bkg_prob);
  inBkgtree->SetBranchAddress("bkg_cost",     &bkg_cost);
  // inBkgtree->SetBranchAddress("bkg_cosTS",    &bkg_cosTS);
  // inBkgtree->SetBranchAddress("bkg_cosTS1",    &bkg_cosTS1);
  // inBkgtree->SetBranchAddress("bkg_cosTS2",    &bkg_cosTS2);
  inBkgtree->SetBranchAddress("bkg_mtop",     &bkg_mtop);
  inBkgtree->SetBranchAddress("bkg_mw",     &bkg_mw);
  inBkgtree->SetBranchAddress("bkg_mva",     &bkg_mva);
  inBkgtree->SetBranchAddress("bkg_sample",  &bkg_sample);

  
  vector < pair < float,float > > ptpaireff[labelsv.size()];
  
  for(unsigned int isig=0; isig<inSigtree->GetEntries(); isig++) {
    inSigtree->GetEntry(isig);

    if(sig_mva > mvaWP[sig_sample]){ 

      hSigTopPt[sig_sample]->Fill(sig_genTopPt);
      if(sig_mtop != -999){ hSigMtop[sig_sample]->Fill(sig_mtop); }
      if(sig_mw   != -999){ hSigMw[sig_sample]  ->Fill(sig_mw);   }
      ptpaireff[sig_sample].push_back(make_pair(1,sig_genTopPt));
    } else {
      ptpaireff[sig_sample].push_back(make_pair(0,sig_genTopPt));
    }
    
    for(unsigned int ic=0; ic<cutsv.size(); ic++) {
	bool sig_pass = sig_mva > cutsv[ic];
	if(sig_pass){ sigpass[sig_sample][ic] += 1; }
	else        { sigfail[sig_sample][ic] += 1; }
      }
  }
  
  for(unsigned int ibkg=0; ibkg<inBkgtree->GetEntries(); ibkg++){
    inBkgtree->GetEntry(ibkg);
    if(bkg_mva > mvaWP[bkg_sample]){
      if(bkg_mtop != -999){ hBkgMtop[bkg_sample]->Fill(bkg_mtop); }
      if(bkg_mw   != -999){ hBkgMw[bkg_sample]  ->Fill(bkg_mw);   }
    }
    
    for(unsigned int ic=0; ic<cutsv.size();ic++){
      bool bkg_pass = bkg_mva > cutsv[ic];
      if(bkg_pass){ bkgpass[bkg_sample][ic] += bkg_weight; }
      else        { bkgfail[bkg_sample][ic] += bkg_weight; }
    }
  }
  
  delete inSigBkgfile;
  inSigBkgfile=0, inSigtree=0, inBkgtree=0;

  for(unsigned int it=0; it<labelsv.size(); it++){
    hSigTopPt[it]->Scale(1.0/hSigTopPt[it]->Integral());
    hSigMtop[it] ->Scale(1.0/hSigMtop[it] ->Integral());
    hSigMw[it]   ->Scale(1.0/hSigMw[it]   ->Integral());
    hBkgMtop[it] ->Scale(1.0/hBkgMtop[it] ->Integral());
    hBkgMw[it]   ->Scale(1.0/hBkgMw[it]   ->Integral());

  }
  
  double sigeff[labelsv.size()][cutsv.size()];
  double bkgeff[labelsv.size()][cutsv.size()];
  for(unsigned int iw=0; iw<labelsv.size(); iw++) {
    for(unsigned int icut=0; icut<cutsv.size(); icut++) {    
      sigeff[iw][icut] = sigpass[iw][icut]/(sigpass[iw][icut]+sigfail[iw][icut]);
      bkgeff[iw][icut] = bkgpass[iw][icut]/(bkgpass[iw][icut]+bkgfail[iw][icut]);  
    }
  }

  for(unsigned int is=0; is < labelsv.size(); ++is){
    for(unsigned int ilw=0; ilw < ptlow.size(); ++ilw){
      for(unsigned int ip=0; ip < ptpaireff[is].size(); ++ip){
	if(ptpaireff[is][ip].second > ptlow[ilw] && ptpaireff[is][ip].second < pthigh[ilw]){
	  if     (ptpaireff[is][ip].first==1)     { Npass[is][ilw] += 1; }
	  else if(ptpaireff[is][ip].first==0)     { Nfail[is][ilw] += 1; }
	}
      }
    }
  }

  double topeff[labelsv.size()][ptlow.size()];
  for(unsigned int iw=0; iw < labelsv.size(); iw++){
    for(unsigned int ip=0; ip < ptlow.size(); ip++){
      topeff[iw][ip] = Npass[iw][ip]/(Npass[iw][ip]+Nfail[iw][ip]);
    }
  }

  double ptmid[ptlow.size()];
  for(unsigned int im=0; im<ptlow.size(); im++){
    ptmid[im] = (ptlow[im]+pthigh[im])/2.0;
  }

  for(unsigned int iw=0; iw<labelsv.size(); iw++) {
    rocsv.push_back(new TGraphErrors(cutsv.size(), sigeff[iw], bkgeff[iw], 0, 0));
    TopEff.push_back(new TGraphErrors(ptlow.size(), ptmid, topeff[iw], 0 ,0));    
  }
  
 

  
  //--------------------------------------------------------------------------------------------------------------
  // Make plots
  //==============================================================================================================
  
  // vector <float> norm_mdrop,norm_topmva, norm_qgid1, norm_qgid2, norm_detaj1b, norm_detaj2b, norm_dphij1b, norm_dphij2b, norm_drj1b, norm_drj2b, norm_bjcsv, norm_j1csv, norm_j2csv;
   
  // for(unsigned int j = 0; j<9; j++){
    
  //   // norm_topmva.push_back(1.0/htopmvav[j]->Integral());    htopmvav[j]->Scale(norm_topmva[j]);
  //   norm_bjcsv.push_back(1.0/hbjcsv[j]->Integral());       hbjcsv[j]->Scale(norm_bjcsv[j]);
  //   norm_j1csv.push_back(1.0/hj1csv[j]->Integral());       hj1csv[j]->Scale(norm_j1csv[j]);
  //   norm_j2csv.push_back(1.0/hj2csv[j]->Integral());       hj2csv[j]->Scale(norm_j2csv[j]);
  //   norm_mdrop.push_back(1.0/hmdropv[j]->Integral());      hmdropv[j]->Scale(norm_mdrop[j]);
  //   norm_detaj1b.push_back(1.0/hdetaj1bv[j]->Integral());  hdetaj1bv[j]->Scale(norm_detaj1b[j]);
  //   norm_detaj2b.push_back(1.0/hdetaj2bv[j]->Integral());  hdetaj2bv[j]->Scale(norm_detaj2b[j]);
  //   norm_dphij1b.push_back(1.0/hdphij1bv[j]->Integral());  hdphij1bv[j]->Scale(norm_dphij1b[j]);
  //   norm_dphij2b.push_back(1.0/hdphij2bv[j]->Integral());  hdphij2bv[j]->Scale(norm_dphij2b[j]);
  //   norm_drj1b.push_back(1.0/hdrj1bv[j]->Integral());      hdrj1bv[j]->Scale(norm_drj1b[j]);
  //   norm_drj2b.push_back(1.0/hdrj2bv[j]->Integral());      hdrj2bv[j]->Scale(norm_drj2b[j]);
  //   norm_qgid1.push_back(1.0/hqgid1v[j]->Integral());      hqgid1v[j]->Scale(norm_qgid1[j]);
  //   norm_qgid2.push_back(1.0/hqgid2v[j]->Integral());      hqgid2v[j]->Scale(norm_qgid2[j]);
  // }
				 
  //
  // ROC curve
  //  
  gStyle->SetPalette(1);
  TCanvas *c = MakeCanvas("c","c",800,800);
  c->SetTickx(1);
  c->SetTicky(1);
  
  char ylabel[100];

  CPlot plot("roc_combined","","Signal Efficiency","Background Efficiency");
  for(unsigned int iw=0; iw<labelsv.size(); iw++) {
    plot.AddGraph(rocsv[iw],labelsv[iw].c_str(),"CP",colorsv[iw],kFullDotMedium,linesv[iw]);
  }
  plot.AddLine(0,0,1,1,kBlack,7);
  plot.SetYRange(0,1);
  plot.SetXRange(0,1);
  plot.AddTextBox("CMS",0.20,0.90,0.32,0.84,0,kBlack,62);
  plot.AddTextBox("Simulation",0.20,0.84,0.38,0.79,0,kBlack,52);
  plot.AddTextBox("(13 TeV)",0.80,0.99,0.95,0.92,0,kBlack);
  plot.TransLegend(-0.4,-0.12);
  plot.Draw(c,true,"png");


  CPlot ploteff("topeffvspt","","top p_{T}","#epsilon");
  for(unsigned int iw=0; iw<labelsv.size(); iw++) {
    ploteff.AddGraph(TopEff[iw],labelsv[iw].c_str(),"CP",colorsv[iw],kFullDotLarge,7);
  }
  ploteff.SetYRange(0,1);
  ploteff.SetXRange(0,400);
  ploteff.TransLegend(-0.05,-0.02);
  ploteff.Draw(c,true,"png");
 
  for(unsigned int it=0; it<labelsv.size(); it++){
    sprintf(ylabel, "Fraction / %.2f", hSigMtop[it]->GetBinWidth(1));
    sprintf(hname,"mtop_%i",it);
    CPlot plotMtop(hname,"","m_{top} [GeV]",ylabel);
    plotMtop.AddHist1D(hBkgMtop[it], "Background","hist", kBlue-6, 1,1001);
    plotMtop.AddHist1D(hSigMtop[it], "Signal","hist",kRed+2,1,3004);
    plotMtop.TransLegend(-0.05,-0.05);
    plotMtop.Draw(c,true,"png");
    
    sprintf(ylabel, "Fraction / %.2f", hSigMw[it]->GetBinWidth(1));
    sprintf(hname,"m_W_%i",it);
    CPlot plotMw(hname,"","m_{W} [GeV]",ylabel);
    plotMw.AddHist1D(hBkgMw[it], "Background","hist", kBlue-6, 1,1001);
    plotMw.AddHist1D(hSigMw[it], "Signal","hist",kRed+2,1,3004);
    plotMw.TransLegend(-0.05,-0.05);
    plotMw.Draw(c,true,"png");

    sprintf(ylabel,"Fraction / %.2f",hSigTopPt[it]->GetBinWidth(1));
    sprintf(hname,"genTopPt_%i",it);
    CPlot plotGenTopPt(hname,"","gen top p_{T}",ylabel);
    plotGenTopPt.AddHist1D(hSigTopPt[it],"Signal Top","hist", kRed+2, 1, 3004);
    plotGenTopPt.Draw(c,true,"png");
  }

  // sprintf(ylabel,"Fraction / %.2f",hmdropv[0]->GetBinWidth(1));
  // CPlot plotmdropZjets("mdrop_zjets","","#zeta",ylabel);
  // plotmdropZjets.AddHist1D(hmdropv[0],"Z+jets","hist",kBlue,1,0);
  // plotmdropZjets.AddHist1D(hmdropv[1],"Z+jets + btag","hist",kGreen+2,1,0);
  // plotmdropZjets.AddHist1D(hmdropv[2],"Signal","hist",kRed,1,0);
  // plotmdropZjets.TransLegend(-0.05,-0.05);
  // plotmdropZjets.Draw(c,true,"png");
  
  // sprintf(ylabel,"Fraction / %.2f",hbjcsv[0]->GetBinWidth(1));
  // CPlot plotbjcsv("bjcsv_out","","CSVv2+IVF",ylabel);
  // plotbjcsv.AddHist1D(hbjcsv[0],"Signal Top","hist",kRed,1,0);
  // plotbjcsv.AddHist1D(hbjcsv[1],"B mismatch","hist",kGreen+2,1,0);
  // plotbjcsv.AddHist1D(hbjcsv[2],"W mismatch","hist",kBlue,1,0);
  // plotbjcsv.AddHist1D(hbjcsv[3],"B & W mismatch","hist",kOrange+7,1,0);
  // plotbjcsv.TransLegend(-0.05,-0.05);
  // plotbjcsv.Draw(c,true,"png");

  // sprintf(ylabel,"Fraction / %.2f",hj1csv[0]->GetBinWidth(1));
  // CPlot plotj1csv("j1csv_out","","CSVv2+IVF",ylabel);
  // plotj1csv.AddHist1D(hj1csv[0],"Signal Top","hist",kRed,1,0);
  // plotj1csv.AddHist1D(hj1csv[1],"B mismatch","hist",kGreen+2,1,0);
  // plotj1csv.AddHist1D(hj1csv[2],"W mismatch","hist",kBlue,1,0);
  // plotj1csv.AddHist1D(hj1csv[3],"B & W mismatch","hist",kOrange+7,1,0);
  // plotj1csv.TransLegend(-0.05,-0.05);
  // plotj1csv.Draw(c,true,"png");

  // sprintf(ylabel,"Fraction / %.2f",hj2csv[0]->GetBinWidth(1));
  // CPlot plotj2csv("j2csv_out","","CSVv2+IVF",ylabel);
  // plotj2csv.AddHist1D(hj2csv[0],"Signal Top","hist",kRed,1,0);
  // plotj2csv.AddHist1D(hj2csv[1],"B mismatch","hist",kGreen+2,1,0);
  // plotj2csv.AddHist1D(hj2csv[2],"W mismatch","hist",kBlue,1,0);
  // plotj2csv.AddHist1D(hj2csv[3],"B & W mismatch","hist",kOrange+7,1,0);
  // plotj2csv.TransLegend(-0.05,-0.05);
  // plotj2csv.Draw(c,true,"png");

  // sprintf(ylabel,"Fraction / %.2f",hqgid1v[0]->GetBinWidth(1));
  // CPlot plotqgid1("qgid1_out", "","Q/G",ylabel);
  // plotqgid1.AddHist1D(hqgid1v[0],"Signal Top","hist",kRed,1,0);
  // plotqgid1.AddHist1D(hqgid1v[1],"B mismatch","hist",kGreen+2,1,0);
  // plotqgid1.AddHist1D(hqgid1v[2],"W mismatch","hist",kBlue,1,0);
  // plotqgid1.AddHist1D(hqgid1v[3],"B & W mismatch","hist",kOrange+7,1,0);
  // plotqgid1.TransLegend(-0.05,-0.05);
  // plotqgid1.Draw(c,true,"png");

  // sprintf(ylabel,"Fraction / %.2f",hqgid2v[0]->GetBinWidth(1));
  // CPlot plotqgid2("qgid2_out", "","Q/G",ylabel);
  // plotqgid2.AddHist1D(hqgid2v[0],"Signal Top","hist",kRed,1,0);
  // plotqgid2.AddHist1D(hqgid2v[1],"B mismatch","hist",kGreen+2,1,0);
  // plotqgid2.AddHist1D(hqgid2v[2],"W mismatch","hist",kBlue,1,0);
  // plotqgid2.AddHist1D(hqgid2v[3],"B & W mismatch","hist",kOrange+7,1,0);
  // plotqgid2.TransLegend(-0.05,-0.05);
  // plotqgid2.Draw(c,true,"png");

  // // sprintf(ylabel,"Fraction / %.2f",hmdropv[0]->GetBinWidth(1));
  // // CPlot plotmdrop("mdrop_out", "","#zeta",ylabel);
  // // plotmdrop.AddHist1D(hmdropv[0],"Signal Top","hist",kRed,1,0);
  // // plotmdrop.AddHist1D(hmdropv[1],"B mismatch","hist",kGreen+2,1,0);
  // // plotmdrop.AddHist1D(hmdropv[2],"W mismatch","hist",kBlue,1,0);
  // // plotmdrop.AddHist1D(hmdropv[3],"B & W mismatch","hist",kOrange+7,1,0);
  // // plotmdrop.TransLegend(-0.05,-0.05);
  // // plotmdrop.Draw(c,true,"png");
  
  
  

  // sprintf(ylabel,"Fraction / %.2f",hdetaj1bv[0]->GetBinWidth(1));
  // CPlot plotdetaj1b("detaj1b_out","#Delta #eta (j_{1},B)","#Delta #eta",ylabel);
  // plotdetaj1b.AddHist1D(hdetaj1bv[0],"Signal Top","hist",kRed,1,0);
  // plotdetaj1b.AddHist1D(hdetaj1bv[1],"B mismatch","hist",kGreen+2,1,0);
  // plotdetaj1b.AddHist1D(hdetaj1bv[2],"W mismatch","hist",kBlue,1,0);
  // plotdetaj1b.AddHist1D(hdetaj1bv[3],"B & W mismatch","hist",kOrange+7,1,0);
  // plotdetaj1b.TransLegend(-0.05,-0.05);
  // plotdetaj1b.Draw(c,true,"png");

  // sprintf(ylabel,"Fraction / %.2f",hdetaj2bv[0]->GetBinWidth(1));
  // CPlot plotdetaj2b("detaj2b_out","#Delta #eta (j_{2},B)","#Delta #eta",ylabel);
  // plotdetaj2b.AddHist1D(hdetaj2bv[0],"Signal Top","hist",kRed,1,0);
  // plotdetaj2b.AddHist1D(hdetaj2bv[1],"B mismatch","hist",kGreen+2,1,0);
  // plotdetaj2b.AddHist1D(hdetaj2bv[2],"W mismatch","hist",kBlue,1,0);
  // plotdetaj2b.AddHist1D(hdetaj2bv[3],"B & W mismatch","hist",kOrange+7,1,0);
  // plotdetaj2b.TransLegend(-0.05,-0.05);
  // plotdetaj2b.Draw(c,true,"png");

  // sprintf(ylabel,"Fraction / %.2f",hdphij1bv[0]->GetBinWidth(1));
  // CPlot plotdphij1b("dphij1b_out","#Delta #phi (j_{1},B)","#Delta #phi",ylabel);
  // plotdphij1b.AddHist1D(hdphij1bv[0],"Signal Top","hist",kRed,1,0);
  // plotdphij1b.AddHist1D(hdphij1bv[1],"B mismatch","hist",kGreen+2,1,0);
  // plotdphij1b.AddHist1D(hdphij1bv[2],"W mismatch","hist",kBlue,1,0);
  // plotdphij1b.AddHist1D(hdphij1bv[3],"B & W mismatch","hist",kOrange+7,1,0);
  // plotdphij1b.SetYRange(0,0.04);
  // plotdphij1b.TransLegend(-0.05,-0.02);
  // plotdphij1b.Draw(c,true,"png");

  // sprintf(ylabel,"Fraction / %.2f",hdphij2bv[0]->GetBinWidth(1));
  // CPlot plotdphij2b("dphij2b_out","#Delta #phi (j_{2},B)","#Delta #phi",ylabel);
  // plotdphij2b.AddHist1D(hdphij2bv[0],"Signal Top","hist",kRed,1,0);
  // plotdphij2b.AddHist1D(hdphij2bv[1],"B mismatch","hist",kGreen+2,1,0);
  // plotdphij2b.AddHist1D(hdphij2bv[2],"W mismatch","hist",kBlue,1,0);
  // plotdphij2b.AddHist1D(hdphij2bv[3],"B & W mismatch","hist",kOrange+7,1,0);
  // plotdphij2b.SetYRange(0,0.04);
  // plotdphij2b.TransLegend(-0.05,-0.02);
  // plotdphij2b.Draw(c,true,"png");

  // sprintf(ylabel,"Fraction / %.2f",hdrj1bv[0]->GetBinWidth(1));
  // CPlot plotdrj1b("drj1b_out","#Delta R (j_{1},B)","#Delta R",ylabel);
  // plotdrj1b.AddHist1D(hdrj1bv[0],"Signal Top","hist",kRed,1,0);
  // plotdrj1b.AddHist1D(hdrj1bv[1],"B mismatch","hist",kGreen+2,1,0);
  // plotdrj1b.AddHist1D(hdrj1bv[2],"W mismatch","hist",kBlue,1,0);
  // plotdrj1b.AddHist1D(hdrj1bv[3],"B & W mismatch","hist",kOrange+7,1,0);
  // plotdrj1b.SetYRange(0,0.055);
  // plotdrj1b.TransLegend(-0.08,-0.05);
  // plotdrj1b.Draw(c,true,"png");

  // sprintf(ylabel,"Fraction / %.2f",hdrj2bv[0]->GetBinWidth(1));
  // CPlot plotdrj2b("drj2b_out","#Delta R (j_{2},B)","#Delta R",ylabel);
  // plotdrj2b.AddHist1D(hdrj2bv[0],"Signal Top","hist",kRed,1,0);
  // plotdrj2b.AddHist1D(hdrj2bv[1],"B mismatch","hist",kGreen+2,1,0);
  // plotdrj2b.AddHist1D(hdrj2bv[2],"W mismatch","hist",kBlue,1,0);
  // plotdrj2b.AddHist1D(hdrj2bv[3],"B & W mismatch","hist",kOrange+7,1,0);
  // plotdrj2b.SetYRange(0,0.055);
  // plotdrj2b.TransLegend(-0.08,-0.05);
  // plotdrj2b.Draw(c,true,"png");

}
