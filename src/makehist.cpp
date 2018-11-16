/**
 * Extract flux histograms from art ROOT files generateed with dk2nu,
 * with ppfx weights.
 *
 * Authors: J. Zennamo, A. Mastbaum
 */

#include <iostream>
#include <cstdlib>
#include <string>
#include <vector>
#include "TDirectory.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TH2D.h"
#include "TFile.h"
#include "canvas/Utilities/InputTag.h"
#include "gallery/Event.h"
#include "gallery/ValidHandle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "larsim/EventWeight/Base/MCEventWeight.h"
#include "geo/GeoVector.h"
#include "geo/GeoAABox.h"
#include "geo/GeoHalfLine.h"
#include "geo/GeoAlgo.h"

using namespace art;
using namespace std;

int main(int argc, char** argv) {
  
  std::pair<float, float>  _xRange;
  std::pair<float, float>  _yRange;
  std::pair<float, float>  _zRange;

  _xRange.first  = -176;
  _xRange.second = 177;
  _yRange.first  = -172;
  _yRange.second = 179;
  _zRange.first =  25;
  _zRange.second = 1150;

  geoalgo::GeoAlgo const _geo_algo_instance;

  geoalgo::AABox volAVTPC( _xRange.first, _yRange.first, _zRange.first, _xRange.second, _yRange.second, _zRange.second);

  InputTag mctruths_tag { "flux" };
  InputTag  evtwght_tag { "eventweight" };
 
  vector<string> filename;
  for (int i=1; i<argc; i++) { 
    std::cout << "FILE : " << argv[i] << std::endl; 
    filename.push_back(string(argv[i]));
  }

  // Histograms for each flavor
  std::vector<TH1D*> Enu_CV_Window;
  std::vector<TH1D*> Enu_CV_AV_TPC;
  std::vector<TH1D*> Enu_UW_Window;
  std::vector<TH1D*> Enu_UW_AV_TPC;

  // flavors - systematic - universe 
  std::vector<std::vector<std::vector<TH1D*> > > Enu_Syst_Window;
  std::vector<std::vector<std::vector<TH1D*> > > Enu_Syst_AV_TPC;

  std::vector<double> TotWeight;
  TotWeight.resize(100);

  // systematic - universe 
  std::vector< std::vector< double > > Weights;   

  std::vector<string> flav = { "numu", "nue", "numubar", "nuebar" };
  
  std::vector< std::vector<double> > bins; bins.resize(4);
  bins[0] = {
    0.0,  0.6,  0.8,  1.0,  1.2,  1.4,  1.6,  1.8,  2.0,  2.2,  2.4,  2.6,
    2.8,  3.0,  3.2,  3.4,  3.6,  3.8,  4.0,  4.2,  4.4,  4.6,  4.8,  5.0,
    6.0,  7.0,  8.0,  9.0,  10.0,  11.0,  12.0,  13.0,  14.0,  15.0,  16.0,
    17.0,  18.0,  19.0,  20.0
  };
  bins[1] = {  0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0};
  bins[2] = {  0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0};
  bins[3] = {  0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0};

  std::vector<string> labels = {"ms_PPFX","Total"};
  //labels = {"PPFXMIPPKaon","PPFXMIPPPion","PPFXOther","PPFXTargAtten","PPFXThinKaon","PPFXThinMeson","PPFXThinNeutron","PPFXThinNucA","PPFXThinNuc","PPFXThinPion","PPFXTotAbsorp","Total"};

  Weights.resize(labels.size());

  for (int i=0; i<labels.size(); i++) {
    Weights[i].resize(100);
  }

  Enu_CV_Window.resize(4);
  Enu_CV_AV_TPC.resize(4);
  Enu_UW_Window.resize(4);
  Enu_UW_AV_TPC.resize(4);
  Enu_Syst_Window.resize(4);
  Enu_Syst_AV_TPC.resize(4);
  std::vector<double> temp;

  // Flavors
  for(unsigned i=0; i<flav.size(); i++) {
    int const n = bins[i].size()-1;
    temp.clear();
    temp = bins[i];

    double* bin = &temp[0];

    Enu_CV_Window[i] = new TH1D(Form("%s_CV_Window",flav[i].c_str()),"",n, bin);
    Enu_CV_AV_TPC[i] = new TH1D(Form("%s_CV_AV_TPC",flav[i].c_str()),"",n, bin);

    Enu_UW_Window[i] = new TH1D(Form("%s_unweighted_Window",flav[i].c_str()),"",n, bin);
    Enu_UW_AV_TPC[i] = new TH1D(Form("%s_unweighted_AV_TPC",flav[i].c_str()),"",n, bin);

    Enu_Syst_Window[i].resize(labels.size());
    Enu_Syst_AV_TPC[i].resize(labels.size());

    // Labels
    for(unsigned j=0; j<labels.size(); j++) {
      Enu_Syst_Window[i][j].resize(100);
      Enu_Syst_AV_TPC[i][j].resize(100);

      // Universes
      for(int k=0; k<100; k++){
        Enu_Syst_Window[i][j][k] =  new TH1D(Form("%s_%s_Uni_%d_Window",flav[i].c_str(), labels[j].c_str(), k),"",n, bin);
        Enu_Syst_AV_TPC[i][j][k] =  new TH1D(Form("%s_%s_Uni_%d_AV_TPC",flav[i].c_str(), labels[j].c_str(), k),"",n, bin);
      }
    }
  }

  // Event loop
  int y = 0;
  int n = 0;

  for (gallery::Event ev(filename); !ev.atEnd(); ev.next()) {
    n++;

    auto const& mctruths = *ev.getValidHandle<vector<simb::MCTruth>>(mctruths_tag);   
    auto const& mcfluxs = *ev.getValidHandle<vector<simb::MCFlux>>(mctruths_tag);   
    bool EW = true;
    auto const& evtwghts = *ev.getValidHandle<vector<evwgh::MCEventWeight>>(evtwght_tag);  
    
    for (size_t i=0; i<mctruths.size(); i++) {
      auto const& mctruth = mctruths.at(i);
      auto const& mcflux = mcfluxs.at(i);
      evwgh::MCEventWeight evtwght;

      if (EW) {
        evtwght = evtwghts.at(i);
      }
      
      int pdg;
      if (mctruth.GetNeutrino().Nu().PdgCode() == 14) {pdg = 0;}
      else if(mctruth.GetNeutrino().Nu().PdgCode() == 12) {pdg = 1;}
      else if(mctruth.GetNeutrino().Nu().PdgCode() ==-14) {pdg = 2;}
      else if(mctruth.GetNeutrino().Nu().PdgCode() ==-12) {pdg = 3;}
      else {
        std::cout << "Unknown neutrino PDG: "
                  << mctruth.GetNeutrino().Nu().PdgCode()
                  << std::endl;
        continue;
      }

      geoalgo::HalfLine ray(mctruth.GetNeutrino().Nu().Vx(),
                            mctruth.GetNeutrino().Nu().Vy(),
                            mctruth.GetNeutrino().Nu().Vz(),
                            mctruth.GetNeutrino().Nu().Px(),
                            mctruth.GetNeutrino().Nu().Py(),
                            mctruth.GetNeutrino().Nu().Pz());
      
      auto vec = _geo_algo_instance.Intersection(volAVTPC, ray);
      bool intercept = false;

      if (vec.size() == 0) { intercept = false; }
      if (vec.size() == 2) { intercept = true; }
      if (vec.size() != 2 && vec.size() != 0) {
        std::cout << "Neutrino ray has " << vec.size()
                  << " intersection with the detector volume"
                  << std::endl;
      }

      double cv_weight = 1; 

      if(EW) {
        for (unsigned l=0; l<labels.size()-1; l++) {
          std::fill(Weights[l].begin(), Weights[l].end(), 1);
        }
        
        for (auto last : evtwght.fWeight) {
          if (last.first.find("PPFXCV") != std::string::npos) {
            cv_weight = last.second.at(0);
            if(last.second.at(0) > 90 || last.second.at(0) < 0){
              std::cout << "Bad CV weight: " << last.second.at(0) << std::endl;
            }
          }
        }

        cv_weight *= mcflux.fnimpwt * mcflux.fnwtfar;
        
        for (auto last : evtwght.fWeight) {
          for (unsigned l=0; l<labels.size()-1; l++) {
            if (last.first.find(labels[l].c_str()) != std::string::npos) {
              for (unsigned i=0; i<last.second.size(); i++) {
                if (last.second.at(i) > 0 && last.second.at(i) < 30)
                  Weights[l][i] *= last.second.at(i);               
                else
                  Weights[l][i] *= 1;
              }
            }
          }
        }
      }

      Enu_CV_Window[pdg]->Fill(mctruth.GetNeutrino().Nu().E(), cv_weight);
      Enu_UW_Window[pdg]->Fill(mctruth.GetNeutrino().Nu().E(), mcflux.fnimpwt * mcflux.fnwtfar);
      
      if (intercept) {
        Enu_CV_AV_TPC[pdg]->Fill(mctruth.GetNeutrino().Nu().E(), cv_weight);
        Enu_UW_AV_TPC[pdg]->Fill(mctruth.GetNeutrino().Nu().E(), mcflux.fnimpwt * mcflux.fnwtfar);
      }
     
      if (EW) {
        std::fill(TotWeight.begin(), TotWeight.end(), 1);
                
        for (unsigned l=0; l<labels.size()-1; l++) {
          for (unsigned i=0; i<Weights[l].size(); i++) {
            TotWeight[i] *= Weights[l][i];
            Enu_Syst_Window[pdg][l][i]->Fill(mctruth.GetNeutrino().Nu().E(), Weights[l][i]*cv_weight);
            
            if (intercept) {
              Enu_Syst_AV_TPC[pdg][l][i]->Fill(mctruth.GetNeutrino().Nu().E(), Weights[l][i]*cv_weight);
            }
          }
        }
        
        int full_label = labels.size() - 1;
        
        for (unsigned int i=0; i<Weights[0].size(); i++){
          Enu_Syst_Window[pdg][full_label][i]->Fill(mctruth.GetNeutrino().Nu().E(), TotWeight[i]*cv_weight);
          
          if (intercept) {
            Enu_Syst_AV_TPC[pdg][full_label][i]->Fill(mctruth.GetNeutrino().Nu().E(), TotWeight[i]*cv_weight);
          }        
        }
      }
    }
  }

  // Plotting 
  TFile* output = new TFile("output.root", "RECREATE");
  TDirectory* savdir = gDirectory;

  std::vector<std::vector<std::vector<TDirectory*> > > subdir(4); //flav //syst //cont 
  for (unsigned i=0; i<4; i++){
    subdir[i].resize(labels.size()+1);
    for(unsigned j=0; j<labels.size()+1; j++) {
      subdir[i][j].resize(3);
    }
  }

  std::vector<string> cont = { "Window", "Active_TPC_Volume", ""};

  for (int f=0; f<4; f++) {
    std::cout << flav[f] << std::endl;
    subdir[f][0][0] = savdir->mkdir(Form("%s",flav[f].c_str()));
    subdir[f][0][0]->cd();
    
    Enu_CV_Window[f]->Write();      
    Enu_CV_AV_TPC[f]->Write();
    Enu_UW_Window[f]->Write();      
    Enu_UW_AV_TPC[f]->Write();

    for (int s=1; s<labels.size()+1; s++) {
      std::cout << labels[s-1] << std::endl;
      subdir[f][s][0] = subdir[f][0][0]->mkdir(Form("%s",labels[s-1].c_str()));
      subdir[f][s][0]->cd();

      for (int c=1; c<3; c++) {
        std::cout << cont[c-1] << std::endl;
        subdir[f][s][c] = subdir[f][s][0]->mkdir(Form("%s",cont[c-1].c_str()));
        subdir[f][s][c]->cd();
        
        if (c == 1) {
          for(int i = 0; i < 100; i++){
            Enu_Syst_Window[f][s-1][i]->Write();
          }        
        }

        if (c == 2) {
          for(int i = 0; i < 100; i++){
            Enu_Syst_AV_TPC[f][s-1][i]->Write();
          }
        }
      }
    }

    savdir->cd();
  }

  output->Close();
  
  return 0;
}

