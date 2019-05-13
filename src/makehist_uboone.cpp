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
#include "dk2nu/tree/dk2nu.h"
#include "dk2nu/tree/dkmeta.h"
#include "dk2nu/tree/calcLocationWeights.h"

using namespace art;
using namespace std;

int main(int argc, char** argv) {

	std::pair<float, float>  _xRange;
	std::pair<float, float>  _yRange;
	std::pair<float, float>  _zRange;

	// MicroBooNE Fiducial Volume
	_xRange.first  =     -0;
	_xRange.second = 256.35;
	_yRange.first  = -116.5;
	_yRange.second =  116.5;
	_zRange.first  =      0;
	_zRange.second = 1036.8;

	geoalgo::GeoAlgo const _geo_algo_instance;

	geoalgo::AABox volAVTPC( _xRange.first, _yRange.first, _zRange.first, _xRange.second, _yRange.second, _zRange.second);

	InputTag mctruths_tag { "flux" };
	InputTag  evtwght_tag { "eventweight" };

	double totalPOT{0};
	double Kaontotal{0};


	vector<string> badfiles;
	vector<string> filename;
	for (int i=1; i<argc; i++) { 
		TFile *filein=new TFile(argv[i]);
		if (filein->IsZombie()) {
			std::cout << "ERROR: File is ZOMBIE:\t" << argv[i] << std::endl;
			badfiles.push_back(string(argv[i]));
			filein->Close();
		}
		else {
			std::cout << "FILE : " << argv[i] << std::endl; 
			filename.push_back(string(argv[i]));
			totalPOT+=100000; // 100 000 POT per dk2nu file 500 000 if running over nova files
			filein->Close();
		}
	}

	std::cout << "\nTotal POT read in:\t" << totalPOT << std::endl;

	std::cout << "\nUsing 1e5 POT per dk2nu file, the flux will be wrong if this is not the case!\n" << std::endl;

	// Histograms for each flavor
	std::vector<TH1D*> Enu_CV_Window; // Energy
	std::vector<TH1D*> Enu_CV_AV_TPC;
	std::vector<TH1D*> Enu_UW_Window;
	std::vector<TH1D*> Enu_UW_AV_TPC;
	
	std::vector<TH1D*> Th_CV_Window; // Theta
	std::vector<TH1D*> Th_CV_AV_TPC;
	std::vector<TH1D*> Th_UW_Window;
	std::vector<TH1D*> Th_UW_AV_TPC;

	// flavors - systematic - universe 
	std::vector<std::vector<std::vector<TH1D*> > > Enu_Syst_Window;
	std::vector<std::vector<std::vector<TH1D*> > > Enu_Syst_AV_TPC;
	std::vector<std::vector<std::vector<TH1D*> > > Th_Syst_Window;
	std::vector<std::vector<std::vector<TH1D*> > > Th_Syst_AV_TPC;

	// Tree for POT counting
	TTree* POTTree = new TTree("POT","Total POT");
	POTTree -> Branch("POT", &totalPOT);
	POTTree -> Fill();

	// systematic - universe 
	std::vector< std::vector< double > > Weights;   

	std::vector<string> flav = { "numu", "nue", "numubar", "nuebar" };

	// Bins for histograms
	std::vector< std::vector<double> > bins; bins.resize(5);
	bins[0] = { // numu
		0.00, 0.025, 0.03, 0.235 ,0.24, 0.50, 1.00, 1.25, 1.50, 1.75, 2.00, 2.25, 2.50, 3.00, 4.00, 5.00, 6.00, 7.00, 10.00 };

	bins[1] = {  // nue
		0.00 ,0.06, 0.125, 0.25, 0.5, 0.75, 1.00, 1.25, 1.50, 1.75, 2.00, 2.50, 3.00, 3.50, 4.00, 5.00 };

	bins[2] = {// numubar
		0.00, 0.025, 0.03, 0.235 ,0.24, 0.50, 1.00, 1.25, 1.50, 1.75, 2.00, 2.25, 2.50, 3.00, 4.00, 5.00, 6.00, 7.00, 10.00 };

	bins[3] = {  // nuebar
		0.00 ,0.06, 0.125,  0.25, 0.5, 0.75, 1.00, 1.25, 1.50, 1.75, 2.00, 2.50, 3.00, 3.50, 4.00, 5.00 };

	bins[4] = { 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180 }; // theta

	std::vector<string> labels;
	// labels = {"ms_PPFX","Total"};
	labels = {"PPFXMaster"};
	// labels = {"PPFXMIPPKaon","PPFXMIPPPion","PPFXOther","PPFXTargAtten",
	// 	"PPFXThinKaon","PPFXThinMeson","PPFXThinNeutron",
	// 	"PPFXThinNucA","PPFXThinNuc","PPFXThinPion","PPFXTotAbsorp",
	// 	"PPFXMaster", "ms_PPFX"};

	Weights.resize(labels.size());

	for (unsigned int i=0; i<labels.size(); i++) {
		Weights[i].resize(100);
	}

	Enu_CV_Window.resize(4); // Energy
	Enu_CV_AV_TPC.resize(4);
	Enu_UW_Window.resize(4);
	Enu_UW_AV_TPC.resize(4);
	Enu_Syst_Window.resize(4);
	Enu_Syst_AV_TPC.resize(4);
	
	Th_CV_Window.resize(4); // Theta 
	Th_CV_AV_TPC.resize(4);
	Th_UW_Window.resize(4);
	Th_UW_AV_TPC.resize(4);
	Th_Syst_Window.resize(4);
	Th_Syst_AV_TPC.resize(4);
	
	std::vector<double> temp, temp2;

	int const n_th = bins[4].size()-1; // theta bins
	temp2 = bins[4];

	// Flavors
	for(unsigned i=0; i<flav.size(); i++) {
		int const n = bins[i].size()-1;
		temp.clear();
		temp = bins[i];

		double* bin = &temp[0];
		double* bin_th = &temp2[0];

		// FLux histograms
		Enu_CV_Window[i] = new TH1D(Form("%s_CV_Window",flav[i].c_str()),"",n, bin); // Energy
		// Enu_CV_AV_TPC[i] = new TH1D(Form("%s_CV_AV_TPC",flav[i].c_str()),"",n, bin);

		Enu_UW_Window[i] = new TH1D(Form("%s_unweighted_Window",flav[i].c_str()),"",n, bin);
		// Enu_UW_AV_TPC[i] = new TH1D(Form("%s_unweighted_AV_TPC",flav[i].c_str()),"",n, bin);

		
		// Use these for a finer binning
		Enu_CV_AV_TPC[i] = new TH1D(Form("%s_CV_AV_TPC",flav[i].c_str()),"",4000, 0, 20);
		Enu_UW_AV_TPC[i] = new TH1D(Form("%s_unweighted_AV_TPC",flav[i].c_str()),"",4000, 0, 20);

		Th_CV_Window[i] = new TH1D(Form("Th_%s_CV_Window",flav[i].c_str()),"",n_th, bin_th); // Theta
		Th_CV_AV_TPC[i] = new TH1D(Form("Th_%s_CV_AV_TPC",flav[i].c_str()),"",n_th, bin_th);

		Th_UW_Window[i] = new TH1D(Form("Th_%s_unweighted_Window",flav[i].c_str()),"",n_th, bin_th);
		Th_UW_AV_TPC[i] = new TH1D(Form("Th_%s_unweighted_AV_TPC",flav[i].c_str()),"",n_th, bin_th);

		Enu_Syst_Window[i].resize(labels.size());
		Enu_Syst_AV_TPC[i].resize(labels.size());
		Th_Syst_Window[i].resize(labels.size());
		Th_Syst_AV_TPC[i].resize(labels.size());

		// Labels
		for(unsigned j=0; j<labels.size(); j++) {
			Enu_Syst_Window[i][j].resize(100);
			Enu_Syst_AV_TPC[i][j].resize(100);
			Th_Syst_Window[i][j].resize(100);
			Th_Syst_AV_TPC[i][j].resize(100);

			// Universes
			for(int k=0; k<100; k++){
				Enu_Syst_Window[i][j][k] =  new TH1D(Form("%s_%s_Uni_%d_Window",flav[i].c_str(), labels[j].c_str(), k),"",n, bin);
				Enu_Syst_AV_TPC[i][j][k] =  new TH1D(Form("%s_%s_Uni_%d_AV_TPC",flav[i].c_str(), labels[j].c_str(), k),"",n, bin);
				Th_Syst_Window[i][j][k] =   new TH1D(Form("Th_%s_%s_Uni_%d_Window",flav[i].c_str(), labels[j].c_str(), k),"",n_th, bin_th);
				Th_Syst_AV_TPC[i][j][k] =   new TH1D(Form("Th_%s_%s_Uni_%d_AV_TPC",flav[i].c_str(), labels[j].c_str(), k),"",n_th, bin_th);
			}
		}
	}

	// ++++++++++++++++++++++++++++++++
	// Event loop
	// ++++++++++++++++++++++++++++++++
	std::cout << "Starting Eventloop" << std::endl;
	bool disperrors{false}; // chose whether or not to display the unphysical nuetrino errors
	if (!disperrors){
		std::cout << "\n\033[1;31mUser has chosen to suppress the unknown pdg types!\033[0m" << std::endl;
		std::cout << "\033[1;33mChange \"disperrors\" -> true\033[0m\n" << std::endl;
	}

	int n = 0;

	// Loop over events
	for (gallery::Event ev(filename); !ev.atEnd(); ev.next()) {
		n++;

		// Alert the user
    	if (n % 1000000 == 0) std::cout << "On entry " << n/1000000.0 <<"M" << std::endl;

		auto const& mctruths = *ev.getValidHandle<vector<simb::MCTruth>>(mctruths_tag);   
		auto const& mcfluxs = *ev.getValidHandle<vector<simb::MCFlux>>(mctruths_tag);   
		bool EW = true;
		auto const& evtwghts = *ev.getValidHandle<vector<evwgh::MCEventWeight>>(evtwght_tag);  

		// Loop over MCTruths
		for (size_t i=0; i<mctruths.size(); i++) {
			auto const& mctruth = mctruths.at(i);
			auto const& mcflux = mcfluxs.at(i);
			evwgh::MCEventWeight evtwght;

			if (EW) {
				evtwght = evtwghts.at(i);
			}

			int pdg;
			if     (mctruth.GetNeutrino().Nu().PdgCode() == 14) pdg = 0; // numu
			else if(mctruth.GetNeutrino().Nu().PdgCode() == 12) pdg = 1; // nue
			else if(mctruth.GetNeutrino().Nu().PdgCode() ==-14) pdg = 2; // numubar
			else if(mctruth.GetNeutrino().Nu().PdgCode() ==-12) pdg = 3; // nuebar
			else {
				if (disperrors){
				std::cout << "Unknown neutrino PDG: "
					<< mctruth.GetNeutrino().Nu().PdgCode()
					<< std::endl;}
				continue;
			}

			geoalgo::HalfLine ray(mctruth.GetNeutrino().Nu().Vx(),
					mctruth.GetNeutrino().Nu().Vy(),
					mctruth.GetNeutrino().Nu().Vz(),
					mctruth.GetNeutrino().Nu().Px(),
					mctruth.GetNeutrino().Nu().Py(),
					mctruth.GetNeutrino().Nu().Pz());

			// Count nu intersections with tpc
			auto vec = _geo_algo_instance.Intersection(volAVTPC, ray); 
			bool intercept = false;

			if (vec.size() == 0) { intercept = false; } // no intersections
			if (vec.size() == 2) { intercept = true; }  // 2 intersections
			if (vec.size() != 2 && vec.size() != 0) {   // other intersection
				std::cout << "Neutrino ray has " << vec.size()
					<< " intersection with the detector volume"
					<< std::endl;
			}

			double cv_weight = 1; 
			
			// Now get the momentums to calculate theta
			double costheta = mctruth.GetNeutrino().Nu().Pz() / mctruth.GetNeutrino().Nu().P();
			double theta = std::acos(costheta) * 180 / 3.14159265;
			// Count the total Kaons
			if (mcflux.fptype == 321 || mcflux.fptype == -321) Kaontotal+=mcflux.fnimpwt;
			
			if(EW) {

				// Fill Weight vector with 1's to create size=labels
				for (unsigned l=0; l<labels.size(); l++) {
					std::fill(Weights[l].begin(), Weights[l].end(), 1);
				}

				// Get the cv_weight
				for (auto last : evtwght.fWeight) {

					if (last.first.find("PPFXCV") != std::string::npos) {

						if(last.second.at(0) > 30 || last.second.at(0) < 0){ // still fill even if bad weight, changed from >90 to >30
							std::cout << "Bad CV weight, setting to 1: " << last.second.at(0) << std::endl;
							cv_weight = 1;
							// cv_weight = last.second.at(0);
						}
						else {
							// std::cout << "CV weight:\t" << last.second.at(0) << std::endl;
							cv_weight = last.second.at(0);
						}  

					}

				} 
				 
				// Weight of neutrino parent (importance weight) * Neutrino weight for a decay forced at center of near detector 
				cv_weight *= mcflux.fnimpwt * mcflux.fnwtfar; // mcflux.fnwtfar == mcflux.fnwtnear

				// Error handling
				if (cv_weight < 0) cv_weight = 0; // get rid of them pesky negative weights
				if (std::isnan(cv_weight) == 1) { // catch NaN values
					std::cout << "got a nan:\t"<<cv_weight <<std::endl;
					cv_weight = 0;
				}

				// Loop over all event weight objs
				for (auto last : evtwght.fWeight) { 

					// Loop over all options
					for (unsigned l=0; l<labels.size(); l++) { 

						// Look for string name wishing to push back
						if (last.first.find(labels[l].c_str()) != std::string::npos) {

							// Loop over ms universes
							for (unsigned i=0; i<last.second.size(); i++) { 
								
								// Fill weights 0 < w < 30 otherwise fill 1's
								if (last.second.at(i) > 0 && last.second.at(i) < 30){ 
									Weights[l][i] *= last.second.at(i);          
								}     
								else {
									// std::cout << "Bad Univ weight, setting to 1: " << last.second.at(i) << std::endl;
									Weights[l][i] *= 1;
									// Weights[l][i] *= last.second.at(i);  
								}

							} // End loop over universes

						}

					}

				} // End loop over weights

			} 

			// ++++++++++++++++++++++++++++++++
			// Now got weights, fill histograms
			// ++++++++++++++++++++++++++++++++

			// Window
			Enu_CV_Window[pdg]->Fill(mctruth.GetNeutrino().Nu().E(), cv_weight);
			Enu_UW_Window[pdg]->Fill(mctruth.GetNeutrino().Nu().E(), mcflux.fnimpwt * mcflux.fnwtfar);
			Th_CV_Window[pdg]->Fill(theta, cv_weight);
			Th_UW_Window[pdg]->Fill(theta, mcflux.fnimpwt * mcflux.fnwtfar);

			// TPC AV
			if (intercept) {
				Enu_CV_AV_TPC[pdg]->Fill(mctruth.GetNeutrino().Nu().E(), cv_weight);
				Enu_UW_AV_TPC[pdg]->Fill(mctruth.GetNeutrino().Nu().E(), mcflux.fnimpwt * mcflux.fnwtfar);
				Th_CV_AV_TPC[pdg]->Fill(theta, cv_weight);
				Th_UW_AV_TPC[pdg]->Fill(theta, mcflux.fnimpwt * mcflux.fnwtfar);

				
			}

			// Now fill multisims
			if (EW) {

				// Options        
				for (unsigned l=0; l<labels.size(); l++) {

					// Universes
					for (unsigned i=0; i<Weights[l].size(); i++) {
						Enu_Syst_Window[pdg][l][i]->Fill(mctruth.GetNeutrino().Nu().E(), Weights[l][i]*cv_weight);
						Th_Syst_Window[pdg][l][i]->Fill(theta, Weights[l][i]*cv_weight);

						if (intercept) {
							Enu_Syst_AV_TPC[pdg][l][i]->Fill(mctruth.GetNeutrino().Nu().E(), Weights[l][i]*cv_weight);
							Th_Syst_AV_TPC[pdg][l][i]->Fill(theta, Weights[l][i]*cv_weight);

						}
					}
				}

			}

		} // End loop over mctruth

	} // End loop over events


	// ++++++++++++++++++++++++++++++++
	// Plotting 
	// ++++++++++++++++++++++++++++++++

	// Spit out the total number of Kaons in the file
	std::cout << "Kaon Total:\t" << Kaontotal << std::endl;

	TFile* output = new TFile("output.root", "RECREATE");
	TDirectory* savdir = gDirectory;

	std::cout << "flavour.size:\t" <<flav.size()<<std::endl;

	// Top Flav dir 
	std::vector<std::vector<std::vector<TDirectory*> > > subdir(flav.size()); 
	
	//Create label dirs
	for (unsigned i=0; i<flav.size(); i++){
		subdir[i].resize(labels.size()+1);
		
		// Create Win/AVTPC dirs
		for(unsigned j=0; j<labels.size()+1; j++) {
			subdir[i][j].resize(3);
		}
	
	}

	std::vector<string> cont = { "Window", "Active_TPC_Volume", ""};

	// Flavours
	for (unsigned int f=0; f<flav.size(); f++) {
	
		std::cout << "\n" <<flav[f] << std::endl;

		subdir[f][0][0] = savdir->mkdir(Form("%s",flav[f].c_str()));
		subdir[f][0][0]->cd();

		// Write CV fluxes
		Enu_CV_Window[f]->Write();      
		Enu_CV_AV_TPC[f]->Write();
		Enu_UW_Window[f]->Write();      
		Enu_UW_AV_TPC[f]->Write();
		Th_CV_Window[f]->Write();      
		Th_CV_AV_TPC[f]->Write();
		Th_UW_Window[f]->Write();      
		Th_UW_AV_TPC[f]->Write();

		// Labels
		for (unsigned int s=1; s<labels.size()+1; s++) {
			std::cout << labels[s-1] << std::endl;
			subdir[f][s][0] = subdir[f][0][0]->mkdir(Form("%s",labels[s-1].c_str()));
			subdir[f][s][0]->cd();

			// AV/TPC
			for (int c=1; c<3; c++) {
				std::cout << cont[c-1] << std::endl;
				subdir[f][s][c] = subdir[f][s][0]->mkdir(Form("%s",cont[c-1].c_str()));
				subdir[f][s][c]->cd();

				if (c == 1) {
					for(int i = 0; i < 100; i++){
						Enu_Syst_Window[f][s-1][i]->Write();
						Th_Syst_Window[f][s-1][i]->Write();
					}        
				}

				if (c == 2) {
					for(int i = 0; i < 100; i++){
						Enu_Syst_AV_TPC[f][s-1][i]->Write();
						Th_Syst_AV_TPC[f][s-1][i]->Write();
					}
				}
			}
		}

		savdir->cd();
	}

	POTTree->Write();

	output->Close();

	std::cout << "BADFiles:" << std::endl;
	for (unsigned int i=0;i < badfiles.size(); i++ ){
		std::cout << badfiles.at(i) << std::endl;
	}

	return 0;
}

