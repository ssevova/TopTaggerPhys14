//================================================================================================
//
// Performs preselection to produce bacon bits for training and testing top MVA ID
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                          // access to gROOT, entry point to ROOT system
#include <TSystem.h>                        // interface to OS
#include <TStyle.h>                   // class to handle ROOT plotting styles
#include <TFile.h>                          // file handle class
#include <TTree.h>                          // class to access ntuples
#include <TH1D.h>                           // 1D histogram class
#include <TLorentzVector.h>                 // 4-vector class
#include <TVector2.h>                       // 2-vector class
#include <vector>                           // STL vector class
#include <iostream>                         // standard I/O
#include <iomanip>                          // functions to format standard I/O
#include <fstream>                          // functions for file I/O
#include <string>                           // C++ string class
#include <cmath>                            // C++ math library
#include <cassert>
#include <utility>                          //std::pair

#include <TMVA/Reader.h>

#include "CPlot.hh"
#include "KStyle.hh"
#include "CSample.hh"

#endif

using namespace std;

//=== FUNCTION DECLARATIONS ======================================================================================

double deltaPhi(const double phi1, const double phi2) {
  double result = phi1 - phi2;
  if     (result >  TMath::Pi()) { result = result - 2*TMath::Pi(); }
  else if(result < -TMath::Pi()) { result = result + 2*TMath::Pi(); }
  return result;
}

float maximum( float a, float b, float c )
  {
    float max = ( a < b ) ? b : a;
    return ( ( max < c ) ? c : max );
  }

//=== MAIN MACRO =================================================================================================

void mvaSelect()
{
  //--------------------------------------------------------------------------------------------------------------
  // Settings
  //==============================================================================================================
  
  //
  // Preselection cuts
  //
  const double       WMASSLOW    = 60;
  const double       WMASSHIGH   = 120;
  const double       TOPMASSLOW  = 100;
  const double       TOPMASSHIGH = 300;

  //
  // input/output file
  //  
  const string infilename("/tthome/ssevova/Analysis/phys14/CMSSW_7_2_2_patch1/src/DMSAna/ttDM/baconbits/Phys14-PU20bx25_TTJets_MSDecaysCKM_central_Tune4C_AOD_ptsortbits.root");
  // const string outfilename("dummy.root");
  const string outfilename("trainingbits_tthad_ptsort.root");
    
  // plot output directory
  const string outputDir("TopPlots");

  gSystem->mkdir(outputDir.c_str(), true);
  CPlot::sOutDir = outputDir;
  const string format("png");
  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code
  //==============================================================================================================
  
  //
  // Declare variables to read in ntuple
  //
  unsigned int    runNum, lumiSec, evtNum;                   // event ID
  unsigned int    trainingEvt;
  unsigned int    triggerBits;
  unsigned int    metfilter;                                 // MET filter bits
  unsigned int    npv;                                       // number of PV / PU

  float           scale1fb;                                  // cross section scale factor per 1/fb
  float           pfmetraw,     pfmetphiraw,  dphijetmetraw; // raw PF MET
  float           pfmet,        pfmetphi,     dphijetmet;    // PF MET


  float           j1flavGen,  j2flavGen,  j3flavGen; 
  float           j1qgid,     j2qgid,     j3qgid;            // jet q/g discriminant
  float           j1csv,      j2csv,      j3csv;             // jet CSV b-tagger
  TLorentzVector *vjet1=0,      *vjet2=0,      *vjet3=0;     // jet 4-vector
  TLorentzVector *q1vec=0,      *q2vec=0,      *q3vec=0;
  TLorentzVector *vpar1=0,      *vpar2=0,      *vpar3=0;
  float Prob, Cost;
  float thadpt;

  //
  // Set up output file
  //
  
  unsigned int isSig, b_mis, w_mis, wb_mis;
  unsigned int bkgType;
  unsigned int evtType;
  unsigned int eventNum;
  
  float pttop, mtop, detaj1b, detaj2b, dphij1b, dphij2b, drj1b, drj2b;
  float qgid1, qgid2, mdrop;
  float bjcsv, jet1csv, jet2csv;
  float weight;
  float cosTS1, cosTS2;
  float prob, cost;
  
  TFile *outFile = new TFile(outfilename.c_str(),"RECREATE");
  TTree *outTree = new TTree("Events","Events");
  
  outTree->Branch("bkgType",   &bkgType,   "bkgType/i");
  outTree->Branch("isSig",  &isSig,  "isSig/i");
  outTree->Branch("b_mis",  &b_mis,  "b_mis/i");
  outTree->Branch("w_mis",  &w_mis,  "w_mis/i");
  outTree->Branch("wb_mis", &wb_mis, "wb_mis/i");
  outTree->Branch("weight",  &weight, "weight/F");
  
  outTree->Branch("mdrop",  &mdrop,  "mdrop/F");
  outTree->Branch("mtop",   &mtop,   "mtop/F");
  outTree->Branch("pttop",  &pttop,  "pttop/F");
  outTree->Branch("qgid1",  &qgid1,  "qgid1/F");
  outTree->Branch("qgid2",  &qgid2,  "qgid2/F");
  
  outTree->Branch("detaj1b",&detaj1b,"detaj1b/F");
  outTree->Branch("detaj2b",&detaj2b,"detaj2b/F");
  outTree->Branch("dphij1b",&dphij1b,"dphij1b/F");
  outTree->Branch("dphij2b",&dphij2b,"dphij2b/F");
  outTree->Branch("drj1b",  &drj1b,  "drj1b/F");
  outTree->Branch("drj2b",  &drj2b,  "drj2b/F");
  
  outTree->Branch("bjcsv",    &bjcsv,    "bjcsv/F");
  outTree->Branch("jet1csv",  &jet1csv,  "jet1csv/F");
  outTree->Branch("jet2csv",  &jet2csv,  "jet2csv/F");
  
  outTree->Branch("vjet1", "TLorentzVector", &vjet1);
  outTree->Branch("vjet2", "TLorentzVector", &vjet2);
  outTree->Branch("vjet3", "TLorentzVector", &vjet3);

  outTree->Branch("q1vec", "TLorentzVector", &q1vec);
  outTree->Branch("q2vec", "TLorentzVector", &q2vec);
  outTree->Branch("q3vec", "TLorentzVector", &q3vec);

  outTree->Branch("vpar1", "TLorentzVector", &vpar1);
  outTree->Branch("vpar2", "TLorentzVector", &vpar2);
  outTree->Branch("vpar3", "TLorentzVector", &vpar3);

  outTree->Branch("thadpt", &thadpt, "thadpt/F");
  
  outTree->Branch("evtType",  &evtType,  "evtType/i");
  outTree->Branch("eventNum", &eventNum, "eventNum/i");

  outTree->Branch("prob", &prob, "prob/F");
  outTree->Branch("cost", &cost, "cost/F");

  outTree->Branch("cosTS1", &cosTS1, "cosTS1/F");
  outTree->Branch("cosTS2", &cosTS2, "cosTS2/F");

  // 
  // Declare histograms
  //
  vector<TH1D*> hTopMassCombov, hDijetv, hTopPtv, hDijetPtv, hTriPartonPtv, hDiPartonPtv, hTriPartonMv, hDiPartonMv, hdRj1j2v;
  vector<TH2D*> hDelPt_Mtopv;
  
  char hname[100];

  
  for(unsigned int i=0; i<4; ++i){
    sprintf(hname,"hdRj1j2_%i",i);       hdRj1j2v.push_back(new TH1D(hname,"",50,0,4));                 hdRj1j2v[i]      ->Sumw2();
    sprintf(hname,"hTopMassCombo_%i",i); hTopMassCombov.push_back(new TH1D(hname,"",50,100,500));         hTopMassCombov[i]->Sumw2();
    sprintf(hname,"hTopPt_%i",i);        hTopPtv.push_back(new TH1D(hname,"",50,0,500));                  hTopPtv[i]       ->Sumw2();
    sprintf(hname,"hDijet_%i",i);        hDijetv.push_back(new TH1D(hname,"",50,0,300));                  hDijetv[i]       ->Sumw2();
    sprintf(hname,"hDijetPt_%i",i);      hDijetPtv.push_back(new TH1D(hname,"",50,0,500));                hDijetPtv[i]     ->Sumw2();
    sprintf(hname,"hTriPartonPt_%i",i);  hTriPartonPtv.push_back(new TH1D(hname,"",50,0,500));            hTriPartonPtv[i] ->Sumw2();
    sprintf(hname,"hDiPartonPt_%i",i);   hDiPartonPtv.push_back(new TH1D(hname,"",50,0,500));             hDiPartonPtv[i]  ->Sumw2();
    sprintf(hname,"hDelPt_Mtop_%i",i); hDelPt_Mtopv.push_back(new TH2D(hname,"",50,200,500,50,-100,100)); hDelPt_Mtopv[i]  ->Sumw2();
    sprintf(hname,"hTriPartonM_%i",i); hTriPartonMv.push_back(new TH1D(hname,"",50,100,300));             hTriPartonMv[i]  ->Sumw2();
    sprintf(hname,"hDiPartonM_%i",i);  hDiPartonMv.push_back(new TH1D(hname,"",50,0,300));                hDiPartonMv[i]   ->Sumw2();
  }
  
  TFile *infile=0;
  TTree *intree=0;
  
  cout << " ==> Processing " << infilename << "..." << endl;
  infile = new TFile(infilename.c_str()); assert(infile);
  intree = (TTree*)infile->Get("Events"); assert(intree);

  intree->SetBranchAddress("evtNum",        &evtNum);
  intree->SetBranchAddress("runNum",        &runNum);
  intree->SetBranchAddress("lumiSec",       &lumiSec);
  intree->SetBranchAddress("evtNum",        &evtNum);
  intree->SetBranchAddress("triggerBits",   &triggerBits);
  intree->SetBranchAddress("metfilter",     &metfilter);
  intree->SetBranchAddress("npv",           &npv);

  intree->SetBranchAddress("scale1fb",      &scale1fb);
  intree->SetBranchAddress("pfmetraw",      &pfmetraw);
  intree->SetBranchAddress("pfmetphiraw",   &pfmetphiraw);
  intree->SetBranchAddress("dphijetmetraw", &dphijetmetraw);
  intree->SetBranchAddress("pfmet",         &pfmet);
  intree->SetBranchAddress("pfmetphi",      &pfmetphi);
  intree->SetBranchAddress("dphijetmet",    &dphijetmet);
  
  intree->SetBranchAddress("j1flavGen",     &j1flavGen);
  intree->SetBranchAddress("j1qgid",        &j1qgid);
  intree->SetBranchAddress("j1csv",         &j1csv);
  intree->SetBranchAddress("vjet1",         &vjet1);

  intree->SetBranchAddress("j2flavGen",     &j2flavGen);
  intree->SetBranchAddress("j2qgid",        &j2qgid);
  intree->SetBranchAddress("j2csv",         &j2csv);
  intree->SetBranchAddress("vjet2",         &vjet2);

  intree->SetBranchAddress("j3flavGen",     &j3flavGen);
  intree->SetBranchAddress("j3qgid",        &j3qgid);
  intree->SetBranchAddress("j3csv",         &j3csv);
  intree->SetBranchAddress("vjet3",         &vjet3);

  intree->SetBranchAddress("q1vec",         &q1vec);
  intree->SetBranchAddress("q2vec",         &q2vec);
  intree->SetBranchAddress("q3vec",         &q3vec);

  intree->SetBranchAddress("vpar1",         &vpar1);
  intree->SetBranchAddress("vpar2",         &vpar2);
  intree->SetBranchAddress("vpar3",         &vpar3);

  intree->SetBranchAddress("thadpt",       &thadpt);

  intree->SetBranchAddress("Prob",          &Prob);
  intree->SetBranchAddress("Cost",          &Cost);

  intree->SetBranchAddress("trainingEvt", &trainingEvt);
  
  for(unsigned int ientry=0; ientry<intree->GetEntries(); ientry++) {
    intree->GetEntry(ientry);
        
    unsigned int eventType;
    if(trainingEvt == 1)     { eventType = 1; }
    else if(trainingEvt == 0){ eventType = 0; }
    
    if(metfilter!=0)        continue;
    
    TLorentzVector  j1 = *vjet1;
    TLorentzVector  j2 = *vjet2;
    TLorentzVector  j3 = *vjet3;

    TLorentzVector q1 = *q1vec;
    TLorentzVector q2 = *q2vec;
    TLorentzVector q3 = *q3vec;
    
    TLorentzVector  dijet = (j1+j2); 
    TLorentzVector  vtop = (dijet + j3);

    // polarization angle calculation
    j1.Boost(-1.*dijet.BoostVector());
    dijet.Boost(-1.*vtop.BoostVector());
    float costhetastar1 = TMath::Cos(j1.Angle(dijet.Vect()));

    j2.Boost(-1.*dijet.BoostVector());
    dijet.Boost(-1.*vtop.BoostVector());
    float costhetastar2 = TMath::Cos(j2.Angle(dijet.Vect()));
    
    if(fabs(j1.Eta())>2.4) continue;
    if(fabs(j2.Eta())>2.4) continue;
    if(fabs(j3.Eta())>2.4) continue;

    if(j1.Pt() < 30) continue;
    if(j2.Pt() < 30) continue;
    if(j3.Pt() < 30) continue;

    if(j1csv > 1 || j2csv > 1 || j3csv > 1) continue; //removes O(100) events with inf csv values
    double wgt = 1;
    
    bool wmatch = 
      (j1flavGen==1 && j2flavGen==-2) || (j1flavGen==-2 && j2flavGen==1) || (j1flavGen==-1 && j2flavGen==2) ||  (j1flavGen==2 && j2flavGen==-1) 
      ||  (j1flavGen==3 && j2flavGen==-4) || (j1flavGen==-4 && j2flavGen==3) || (j1flavGen==-3 && j2flavGen==4) || (j1flavGen==4 && j2flavGen==-3);

    // bool single_wmatch = fabs(j1flavGen)==1 || fabs(j1flavGen)==2 || fabs(j1flavGen)==3 || fabs(j1flavGen)==4 ||  
    //                     fabs(j2flavGen)==1 || fabs(j2flavGen)==2 || fabs(j2flavGen)==3 || fabs(j2flavGen)==4;
    bool bmatch = false;
    if( ((j1flavGen== 1 || j1flavGen== 3) && j3flavGen==-5) ||
	((j1flavGen==-1 || j1flavGen==-3) && j3flavGen== 5) ||
	((j1flavGen== 2 || j1flavGen== 4) && j3flavGen== 5) ||
	((j1flavGen==-2 || j1flavGen==-4) && j3flavGen==-5) 
	) 
      bmatch = true;

    float delPt=0;
    if(wmatch && !bmatch){
      hdRj1j2v[0]      ->Fill(j1.DeltaR(j2));
      hTopMassCombov[0]->Fill(vtop.M());           hDijetv[0]     ->Fill(dijet.M());
      hTopPtv[0]       ->Fill(vtop.Pt());          hDijetPtv[0]   ->Fill(dijet.Pt());
      hTriPartonPtv[0] ->Fill((q1+q2+q3).Pt());    hDiPartonPtv[0]->Fill((q1+q2).Pt());
    }
    else if(!wmatch && bmatch ){
      hdRj1j2v[1]      ->Fill(j1.DeltaR(j2));
      hTopMassCombov[1]->Fill(vtop.M());           hDijetv[1]     ->Fill(dijet.M());
      hTopPtv[1]       ->Fill(vtop.Pt());          hDijetPtv[1]   ->Fill(dijet.Pt());
      hTriPartonPtv[1] ->Fill((q1+q2+q3).Pt());    hDiPartonPtv[1]->Fill((q1+q2).Pt());
    }
    else if(!wmatch && !bmatch){
      hdRj1j2v[2]      ->Fill(j1.DeltaR(j2));
      hTopMassCombov[2]->Fill(vtop.M());           hDijetv[2]     ->Fill(dijet.M());
      hTopPtv[2]       ->Fill(vtop.Pt());          hDijetPtv[2]   ->Fill(dijet.Pt());
      hTriPartonPtv[2] ->Fill((q1+q2+q3).Pt());    hDiPartonPtv[2]->Fill((q1+q2).Pt());
    }
    else if( wmatch && bmatch ) {
      hdRj1j2v[3]      ->Fill(j1.DeltaR(j2));
      hTopMassCombov[3]->Fill(vtop.M());           hDijetv[3]     ->Fill(dijet.M());
      hTriPartonMv[3]  ->Fill((q1+q2+q3).M());     hDiPartonMv[3]  ->Fill((q1+q2).M());
      hTopPtv[3]       ->Fill(vtop.Pt());          hDijetPtv[3]   ->Fill(dijet.Pt());
      hTriPartonPtv[3] ->Fill((q1+q2+q3).Pt());    hDiPartonPtv[3]->Fill((q1+q2).Pt());

      float delpt[3]; delpt[0] = j1.Pt() - q1.Pt(); delpt[1] = j2.Pt() - q2.Pt(); delpt[2] = j3.Pt() - q3.Pt();
            
      delPt = maximum(fabs(delpt[0]),fabs(delpt[1]),fabs(delpt[2]));
      float signedDelPt=0;
      for(unsigned int i=0; i<3; ++i){
	if(delPt == fabs(delpt[i])){ signedDelPt = delpt[i];}
      }
      
      hDelPt_Mtopv[3]  ->Fill(vtop.M(),signedDelPt);
      
    }

    //if(dijet.Pt() < 160) continue;
    //if(dijet.M()<WMASSLOW || dijet.M()>WMASSHIGH) continue;
        
    if(vtop.M()<TOPMASSLOW && vtop.M()>TOPMASSHIGH) continue;

    bkgType = 4; // bkgType = 1: tthad combinatoric; bkgType = 2: Zjets; bkgType = 3: tt1l combinatoric; bkType = 4: tt2l combinatoric
    isSig   = wmatch && bmatch;
    b_mis   = wmatch && !bmatch;
    w_mis   = !wmatch && bmatch;
    wb_mis  = !wmatch && !bmatch;
    eventNum = evtNum;
    evtType = eventType;
    pttop   = vtop.Pt();
    mtop    = vtop.M();
    mdrop   = TMath::Max(j1.M(), j2.M())/(dijet.M())*(j1.DeltaR(j2));
    qgid1   = j1qgid;
    qgid2   = j2qgid;
    detaj1b = fabs(j1.Eta() - j3.Eta());
    detaj2b = fabs(j2.Eta() - j3.Eta());
    dphij1b = fabs(j1.DeltaPhi(j3));
    dphij2b = fabs(j2.DeltaPhi(j3));
    drj1b   = j1.DeltaR(j3);
    drj2b   = j2.DeltaR(j3);
    bjcsv    = j3csv;
    jet1csv  = j1csv;
    jet2csv  = j2csv;
    vjet1    = vjet1;
    vjet2    = vjet2;
    vjet3    = vjet3;
    q1vec    = q1vec;
    q2vec    = q2vec;
    q3vec    = q3vec;
    vpar1    = vpar1;
    vpar2    = vpar2;
    vpar3    = vpar3;

    thadpt   = thadpt;
    
    
    prob     = Prob;
    cost     = Cost;
    
    weight   = wgt;

    cosTS1    = costhetastar1;
    cosTS2    = costhetastar2;
      
    outTree->Fill();
    
    
  }
  
  delete infile;
  infile=0;
  intree=0; 

  //
  //Make plots
  //
  
  gStyle->SetTitleOffset(1.600,"Y");
  gStyle->SetPalette(1);
  TCanvas *c = MakeCanvas("c","c",800,600);
  c->SetTickx(1);
  c->SetTicky(1);
  c->SetGridx(1);
  c->SetGridy(1);
  
  vector<float> norm_TopMCombo, norm_dijet, norm_toppt, norm_dijetpt, norm_3parpt, norm_2parpt, norm_3parM, norm_2parM, norm_drjj;
  for(unsigned int k=0; k<4; ++k){
    norm_TopMCombo.push_back(1.0/hTopMassCombov[k]->Integral());
    hTopMassCombov[k]->Scale(norm_TopMCombo[k]);
    
    norm_dijet.push_back(1.0/hDijetv[k]->Integral());
    hDijetv[k]->Scale(norm_dijet[k]);

    norm_toppt.push_back(1.0/hTopPtv[k]->Integral());
    hTopPtv[k]->Scale(norm_toppt[k]);

    norm_dijetpt.push_back(1.0/hDijetPtv[k]->Integral());
    hDijetPtv[k]->Scale(norm_dijetpt[k]);

    norm_3parpt.push_back(1.0/hTriPartonPtv[k]->Integral());
    hTriPartonPtv[k]->Scale(norm_3parpt[k]);

    norm_2parpt.push_back(1.0/hDiPartonPtv[k]->Integral());
    hDiPartonPtv[k]->Scale(norm_2parpt[k]);

    norm_3parM.push_back(1.0/hTriPartonMv[k]->Integral());
    hTriPartonMv[k]->Scale(norm_3parM[k]);

    norm_2parM.push_back(1.0/hDiPartonMv[k]->Integral());
    hDiPartonMv[k]->Scale(norm_2parM[k]);

    norm_drjj.push_back(1.0/hdRj1j2v[k]->Integral());
    hdRj1j2v[k]->Scale(norm_drjj[k]);
  }

  char ylabel[100];

  sprintf(ylabel,"Fraction / %.2f",hTopPtv[0]->GetBinWidth(1));
  CPlot plotTopPt("TopPt_PartonPt_Combinations","","p_{T} [GeV]",ylabel);
  plotTopPt.AddHist1D(hTopPtv[0],     "B mis (jets)","hist",kGreen+2,1,0);
  plotTopPt.AddHist1D(hTriPartonPtv[0],"B mis (quarks)","hist",kGreen+2,7,0);
  plotTopPt.AddHist1D(hTopPtv[1],     "W mis (jets)","hist",kBlue,1,0);
  plotTopPt.AddHist1D(hTriPartonPtv[1],"W mis (quarks)","hist",kBlue,7,0);
  plotTopPt.AddHist1D(hTopPtv[2],     "B&W mis (jets)","hist",kOrange+7,1,0);
  plotTopPt.AddHist1D(hTriPartonPtv[2],"B&W mis (quarks)","hist",kOrange+7,7,0);
  plotTopPt.AddHist1D(hTopPtv[3],     "Signal (jets)","hist",kRed,1,0);
  plotTopPt.AddHist1D(hTriPartonPtv[3],"Signal (quarks)","hist",kRed,7,0);
  plotTopPt.TransLegend(-0.05,-0.02);
  plotTopPt.Draw(c,true,format.c_str());

  sprintf(ylabel,"Fraction / %.2f",hTopMassCombov[0]->GetBinWidth(1));
  CPlot plotTopMassCombo("TopMassCombinations","","Mass [GeV]",ylabel);
  plotTopMassCombo.AddHist1D(hTopMassCombov[0],"B mismatch","hist",kGreen+2,1,0);
  plotTopMassCombo.AddHist1D(hTopMassCombov[1],"W mismatch","hist",kBlue,1,0);
  plotTopMassCombo.AddHist1D(hTopMassCombov[2],"B & W mismatch","hist",kOrange+7,1,0);
  plotTopMassCombo.AddHist1D(hTopMassCombov[3],"Signal Top","hist",kRed,1,0);
  plotTopMassCombo.TransLegend(-0.05,-0.02);
  plotTopMassCombo.Draw(c,true,format.c_str());
  
  sprintf(ylabel,"Fraction / %.2f",hDijetPtv[0]->GetBinWidth(1));
  CPlot plotDijetPt("DijetPt","","Mass [GeV]",ylabel);
  plotDijetPt.AddHist1D(hDijetPtv[0],"B mismatch","hist",kGreen+2,1,0);
  plotDijetPt.AddHist1D(hDijetPtv[1],"W mismatch","hist",kBlue,1,0);
  plotDijetPt.AddHist1D(hDijetPtv[2],"B & W mismatch","hist",kOrange+7,1,0);
  plotDijetPt.AddHist1D(hDijetPtv[3],"Signal Top","hist",kRed,1,0);
  plotDijetPt.TransLegend(-0.05,-0.02);
  plotDijetPt.Draw(c,true,format.c_str());

  sprintf(ylabel,"Fraction / %.2f",hDijetv[0]->GetBinWidth(1));
  CPlot plotDijetMass("DijetMassCombinations","","Mass [GeV]",ylabel);
  plotDijetMass.AddHist1D(hDijetv[0],"B mismatch","hist",kGreen+2,1,0);
  plotDijetMass.AddHist1D(hDijetv[1],"W mismatch","hist",kBlue,1,0);
  plotDijetMass.AddHist1D(hDijetv[2],"B & W mismatch","hist",kOrange+7,1,0);
  plotDijetMass.AddHist1D(hDijetv[3],"Signal Top","hist",kRed,1,0);
  plotDijetMass.TransLegend(-0.05,-0.02);
  plotDijetMass.Draw(c,true,format.c_str());

  CPlot plotDipartonMass("DipartonMassCombinations","","Mass[GeV]",hDiPartonMv[3]->GetBinWidth(1));
  plotDipartonMass.AddHist1D(hDiPartonMv[3],"Signal","hist",kRed,1,0);
  plotDipartonMass.Draw(c,true,format.c_str());

  CPlot plotTripartonMass("TripartonMassCombinations","","Mass[GeV]",hTriPartonMv[3]->GetBinWidth(1));
  plotTripartonMass.AddHist1D(hTriPartonMv[3],"Signal","hist",kRed,1,0);
  plotTripartonMass.Draw(c,true,format.c_str());
  
  sprintf(ylabel,"Fraction / %.2f",hdRj1j2v[0]->GetBinWidth(1));
  CPlot plotdRj1j2("dRj1j2","","Mass [GeV]",ylabel);
  plotdRj1j2.AddHist1D(hdRj1j2v[0],"B mismatch","hist",kGreen+2,1,0);
  plotdRj1j2.AddHist1D(hdRj1j2v[1],"W mismatch","hist",kBlue,1,0);
  plotdRj1j2.AddHist1D(hdRj1j2v[2],"B & W mismatch","hist",kOrange+7,1,0);
  plotdRj1j2.AddHist1D(hdRj1j2v[3],"Signal Top","hist",kRed,1,0);
  plotdRj1j2.TransLegend(-0.05,-0.02);
  plotdRj1j2.Draw(c,true,format.c_str());

  //2D

  CPlot plotdelPtMtop("DelPt_Mtop","","Mass [GeV]","#delta p_{T}");
  plotdelPtMtop.AddHist2D(hDelPt_Mtopv[3],"COLZ");
  plotdelPtMtop.Draw(c,true,format.c_str());

  
  outFile->Write();
  outFile->Close();
}
