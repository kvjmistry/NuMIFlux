#ifndef NuMIFlux_h
#define NuMIFlux_h

#include <iostream>
#include <iomanip>
#include <string>

using namespace std;

#include "TChain.h"
#include "TH1.h"
#include "TH2D.h"
#include "TFile.h"
#include "TSystem.h"
#include "TVector3.h"

#include "dk2nu/tree/dk2nu.h"
#include "dk2nu/tree/dkmeta.h"
#include "dk2nu/tree/calcLocationWeights.h"

// Krishan: Updated this script to take in the dk2u file format instead of flugg 


class NuMIFlux {
public :

  int Nfiles = 0;

  static const int numu  =  14;
  static const int anumu = -14;
  static const int nue   =  12;
  static const int anue  = -12;

  int highest_evtno = 0;
  double NominalPOT = 6e20;
  bool debug = false;
  double fDefaultWeightCorrection = 1./(10000. * TMath::Pi());
  double Ntarget = 4.76e31/56.41e6* 256.35*233*1036.8; //TPC active!!!
  double AccumulatedPOT=0.;
  int treeNumber = -1;

  double histMin = 0;
  double histMax = 20;
  int histNbins = 4000;

  TChain *cflux;
  TChain *cflux_meta;

  TH1D* numuFluxHisto;
  TH1D* anumuFluxHisto;
  TH1D* nueFluxHisto;
  TH1D* anueFluxHisto;
  TH1D* numuCCHisto;
  TH1D* anumuCCHisto;
  TH1D* nueCCHisto;
  TH1D* anueCCHisto;

  // MIPP
  TH2D* pionplus_MIPP;
  TH2D* pionminus_MIPP;
  TH2D* Kplus_MIPP;
  TH2D* Kminus_MIPP;
  // NA49
  TH2D* pionplus_NA49;
  TH2D* pionminus_NA49;
  TH2D* Kplus_NA49;
  TH2D* Kminus_NA49;



  TGraph *genieXsecNumuCC;
  TGraph *genieXsecNumubarCC;
  TGraph *genieXsecNueCC;
  TGraph *genieXsecNuebarCC;
  TFile* f = new TFile("NuMIFlux.root", "RECREATE");


  NuMIFlux(string pattern="/pnfs/uboone/persistent/uboonebeam/numi_dk2nu_zero_threshold/FHC/g4numiv6_minervame_me000z200i_*.root");
  // NuMIFlux(string pattern="/pnfs/uboone/persistent/uboonebeam/numi_dk2nu_zero_threshold/FHC/g4numiv6_minervame_me000z200i_12*_0000.root");

  virtual ~NuMIFlux();

  void CalculateFlux();
  TVector3 RandomInTPC();
  TVector3 FromDetToBeam(const TVector3& det);
  double estimate_pots(int highest_potnum);
  int calcEnuWgt(bsim::Dk2Nu* decay, const TVector3& xyz, double& enu, double& wgt_xy);
  void GetConstraints_ThinTarg( bsim::Dk2Nu* fDk2Nu);
  void GetConstraints_ThickTarg( bsim::Dk2Nu* fDk2Nu);

};

#endif

#ifdef NuMIFlux_cxx

NuMIFlux::NuMIFlux(string pattern) {

  const char* path = "/uboone/app/users/kmistry/PPFX/numi-validation/NuMIFlux_dk2nu";
  if ( path ) {
    TString libs = gSystem->GetDynamicPath();
    libs += ":";
    libs += path;
    //libs += "/lib";
    gSystem->SetDynamicPath(libs.Data());       
    // gSystem->Load("FluxNtuple_C.so");
  }

  cflux = new TChain("dk2nuTree");
  cflux->Add(pattern.c_str());

  cflux_meta = new TChain("dkmetaTree");
  cflux_meta->Add(pattern.c_str());

  Nfiles = cflux->GetNtrees();
  cout << "Number of files: " << Nfiles << endl;

  //Inizialise histos
  TString titleBase1 = "Neutrino Flux;";
  TString titleBase2 = " Energy [GeV];";
  TString titleBase3 = " / cm^{2} / 6e20 POT";
  // numu
  numuFluxHisto = new TH1D("numuFluxHisto", (titleBase1 + "#nu_{#mu}" + titleBase2 +"#nu_{#mu}" + titleBase3),histNbins,histMin,histMax);
  numuCCHisto = new TH1D("numuCCHisto", "numu CC; #nu_{#mu} Energy [GeV]; #nu_{#mu} CC / 79 ton / 6e20 POT",histNbins,histMin,histMax);
  // anumu
  anumuFluxHisto = new TH1D("anumuFluxHisto", (titleBase1 + "#bar{#nu}_{#mu}" + titleBase2 +"#bar{#nu}_{#mu}" + titleBase3),histNbins,histMin,histMax);
  anumuCCHisto = new TH1D("anumuCCHisto", "numu bar CC; #bar{#nu}_{#mu} Energy [GeV]; #bar{#nu}_{#mu} CC / 79 ton / 6e20 POT",histNbins,histMin,histMax);
  // nue
  nueFluxHisto = new TH1D("nueFluxHisto", (titleBase1 + "#nu_{e}" + titleBase2 +"#nu_{e}" + titleBase3),histNbins,histMin,histMax);
  nueCCHisto = new TH1D("nueCCHisto", "nue CC; #nu_{e} Energy [GeV]; #nu_{e} CC / 79 ton / 6e20 POT",histNbins,histMin,histMax);
  // anue
  anueFluxHisto = new TH1D("anueFluxHisto", (titleBase1 + "#bar{#nu}_{e}" + titleBase2 + "#bar{#nu}_{e}" + titleBase3),histNbins,histMin,histMax);
  anueCCHisto = new TH1D("anueCCHisto", "nue bar CC; #bar{#nu}_{e} Energy [GeV]; #bar{#nu}_{#mu} CC / 79 ton / 6e20 POT",histNbins,histMin,histMax);
  


  // Constraint histograms
  // MIPP
  pionplus_MIPP  =  new TH2D("pionplus_MIPP", "pionplus_MIPP  ;p_{L} [GeV/c] ;p_{T} [GeV/c]", 100, 0, 100, 100, 0, 2 );
  pionminus_MIPP =  new TH2D("pionminus_MIPP","pionminus_MIPP ;p_{L} [GeV/c] ;p_{T} [GeV/c]", 100, 0, 100, 100, 0, 2 );
  Kplus_MIPP     =  new TH2D("kaonplus_MIPP", "Kplus_MIPP     ;p_{L} [GeV/c] ;p_{T} [GeV/c]", 100, 0, 100, 100, 0, 2 );
  Kminus_MIPP    =  new TH2D("kaonminus_MIPP","Kminus_MIPP    ;p_{L} [GeV/c] ;p_{T} [GeV/c]", 100, 0, 100, 100, 0, 2 );
  // NA49
  pionplus_NA49  =  new TH2D("pionplus_NA49",  "pionplus_NA49  ;xF  ;p_{T} [GeV/c]", 100, -1, 1, 100, 0, 2.0 );
  pionminus_NA49 =  new TH2D("pionminus_NA49", "pionminus_NA49 ;xF  ;p_{T} [GeV/c]", 100, -1, 1, 100, 0, 2.0 );
  Kplus_NA49     =  new TH2D("kaonplus_NA49",  "Kplus_NA49     ;xF  ;p_{T} [GeV/c]", 100, -1, 1, 100, 0, 2.0 );
  Kminus_NA49    =  new TH2D("kaonminus_NA49", "Kminus_NA49    ;xF  ;p_{T} [GeV/c]", 100, -1, 1, 100, 0, 2.0 );
}

NuMIFlux::~NuMIFlux() {

}

#endif // #ifdef NuMIFlux_cxx
