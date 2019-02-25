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
	_xRange.second =  177;
	_yRange.first  = -172;
	_yRange.second =  179;
	_zRange.first  =   25;
	_zRange.second = 1150;

	geoalgo::GeoAlgo const _geo_algo_instance;

	geoalgo::AABox volAVTPC( _xRange.first, _yRange.first, _zRange.first, _xRange.second, _yRange.second, _zRange.second);

	InputTag mctruths_tag { "flux" };
	InputTag  evtwght_tag { "eventweight" };

	double totalPOT{0};

	vector<string> filename;
	for (int i=1; i<argc; i++) { 
		std::cout << "FILE : " << argv[i] << std::endl; 
		filename.push_back(string(argv[i]));
		totalPOT+=5e5; // 5e5 POT per dk2nu file -- seems to be the standard
	}

	std::cout << "\nTotal POT read in:\t" << totalPOT << std::endl;

	std::cout << "\nUsing 500e3 POT per dk2nu file, the flux will be wrong if this is not the case!\n" << std::endl;

	// Histograms for each flavor
	std::vector<TH1D*> Enu_CV_Window;
	std::vector<TH1D*> Enu_CV_AV_TPC;
	std::vector<TH1D*> Enu_UW_Window;
	std::vector<TH1D*> Enu_UW_AV_TPC;

	// flavors - systematic - universe 
	std::vector<std::vector<std::vector<TH1D*> > > Enu_Syst_Window;
	std::vector<std::vector<std::vector<TH1D*> > > Enu_Syst_AV_TPC;

	// Tree for POT counting
	TTree* POTTree = new TTree("POT","Total POT");
	POTTree -> Branch("POT", &totalPOT);
	POTTree -> Fill();

	// Histogram for weight distributions
	// std::vector<TH2D*> Enu_Weight_CV; // Energy vs weight 2D Hist for CV
	// std::vector<std::vector<TH2D*>> Enu_Weight_MS; // Energy vs weight 2D Hist for Masterweight
	// std::vector<TH1D*> Weight_CV; // weight 1D Hist for CV
	// std::vector<std::vector<TH1D*>> Weight_MS; // weight 1D Hist for Masterweight

	// std::vector<double> TotWeight; // Product of weights in universe "i" for each label
	// TotWeight.resize(100);         // ** Might want to put this further down to make more universal **

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

	std::vector<string> labels;
	// labels = {"ms_PPFX","Total"};
	labels = {"PPFXMaster","ms_PPFX"};
	// labels = {"PPFXMIPPKaon","PPFXMIPPPion","PPFXOther","PPFXTargAtten",
	// 	"PPFXThinKaon","PPFXThinMeson","PPFXThinNeutron",
	// 	"PPFXThinNucA","PPFXThinNuc","PPFXThinPion","PPFXTotAbsorp",
	// 	"PPFXMaster", "ms_PPFX", "Total"};

	Weights.resize(labels.size());

	for (unsigned int i=0; i<labels.size(); i++) {
		Weights[i].resize(100);
	}

	Enu_CV_Window.resize(4);
	Enu_CV_AV_TPC.resize(4);
	Enu_UW_Window.resize(4);
	Enu_UW_AV_TPC.resize(4);
	Enu_Syst_Window.resize(4);
	Enu_Syst_AV_TPC.resize(4);
	// Enu_Weight_CV.resize(4);
	// Enu_Weight_MS.resize(4);
	// Weight_CV.resize(4);
	// Weight_MS.resize(4);
	
	std::vector<double> temp;


	// Flavors
	for(unsigned i=0; i<flav.size(); i++) {
		int const n = bins[i].size()-1;
		temp.clear();
		temp = bins[i];

		double* bin = &temp[0];

		// FLux histograms
		Enu_CV_Window[i] = new TH1D(Form("%s_CV_Window",flav[i].c_str()),"",n, bin);
		Enu_CV_AV_TPC[i] = new TH1D(Form("%s_CV_AV_TPC",flav[i].c_str()),"",n, bin);

		Enu_UW_Window[i] = new TH1D(Form("%s_unweighted_Window",flav[i].c_str()),"",n, bin);
		Enu_UW_AV_TPC[i] = new TH1D(Form("%s_unweighted_AV_TPC",flav[i].c_str()),"",n, bin);

		// Weight histograms
		// Enu_Weight_CV[i] = new TH2D(Form("%s_Enu_vs_CV_wght",flav[i].c_str()),";Enu; CV Weight",100, 0, 25, 200, -0.25, 5);
		// Weight_CV[i] = new TH1D(Form("%s_CV_wght",flav[i].c_str()),";CV Weight",200, -0.25, 5);
		

		Enu_Syst_Window[i].resize(labels.size());
		Enu_Syst_AV_TPC[i].resize(labels.size());
		// Enu_Weight_MS[i].resize(labels.size());
		// Weight_MS[i].resize(labels.size());

		// Labels
		for(unsigned j=0; j<labels.size(); j++) {
			Enu_Syst_Window[i][j].resize(100);
			Enu_Syst_AV_TPC[i][j].resize(100);

			// Enu_Weight_MS[i][j] = new TH2D(Form("%s_Enu_vs_MS_wght_%s",flav[i].c_str(), labels[j].c_str()), ";Enu; MS Weight",100, 0, 25, 200, -0.25, 5);
			// Weight_MS[i][j] = new TH1D(Form("%s_MS_wght_%s",flav[i].c_str(), labels[j].c_str()), ";MS Weight",200, -0.25, 5);

			// Universes
			for(int k=0; k<100; k++){
				Enu_Syst_Window[i][j][k] =  new TH1D(Form("%s_%s_Uni_%d_Window",flav[i].c_str(), labels[j].c_str(), k),"",n, bin);
				Enu_Syst_AV_TPC[i][j][k] =  new TH1D(Form("%s_%s_Uni_%d_AV_TPC",flav[i].c_str(), labels[j].c_str(), k),"",n, bin);
			}
		}
	}

	// ++++++++++++++++++++++++++++++++
	// Event loop
	// ++++++++++++++++++++++++++++++++
	std::cout << "Starting Eventloop" << std::endl;

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
			if (mctruth.GetNeutrino().Nu().PdgCode() == 14) {pdg = 0;}     // numu
			else if(mctruth.GetNeutrino().Nu().PdgCode() == 12) {pdg = 1;} // nue
			else if(mctruth.GetNeutrino().Nu().PdgCode() ==-14) {pdg = 2;} // numubar
			else if(mctruth.GetNeutrino().Nu().PdgCode() ==-12) {pdg = 3;} // nuebar
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

			
			if(EW) {

				// Fill Weight vector with 1's to create size=labels
				for (unsigned l=0; l<labels.size()-1; l++) {
					std::fill(Weights[l].begin(), Weights[l].end(), 1);
				}

				// Get the cv_weight
				for (auto last : evtwght.fWeight) {

					if (last.first.find("PPFXCV") != std::string::npos) {

						// Weights
						// Enu_Weight_CV[pdg]->Fill(mctruth.GetNeutrino().Nu().E(), last.second.at(0) );
						// Weight_CV[pdg]->Fill(last.second.at(0) );

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
				if (cv_weight < 0) cv_weight = 0; // get rid of them pesky negative weights
				if (std::isnan(cv_weight) == 1) { // catch NaN values
					std::cout << "got a nan:\t"<<cv_weight <<std::endl;
					cv_weight = 0;
				
				}

				// Loop over all event weight objs
				for (auto last : evtwght.fWeight) { 

					// Loop over all options
					for (unsigned l=0; l<labels.size()-1; l++) { 

						// Look for string name wishing to push back
						if (last.first.find(labels[l].c_str()) != std::string::npos) {

							// Loop over ms universes
							for (unsigned i=0; i<last.second.size(); i++) { 

								// Weight Hists
								// Enu_Weight_MS[pdg][l]->Fill(mctruth.GetNeutrino().Nu().E(),last.second.at(i));
								// Weight_MS[pdg][l]->Fill(last.second.at(i));
								

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

			// TPC AV
			if (intercept) {
				Enu_CV_AV_TPC[pdg]->Fill(mctruth.GetNeutrino().Nu().E(), cv_weight);
				Enu_UW_AV_TPC[pdg]->Fill(mctruth.GetNeutrino().Nu().E(), mcflux.fnimpwt * mcflux.fnwtfar);

				
			}

			// Now fill multisims
			if (EW) {
				// std::fill(TotWeight.begin(), TotWeight.end(), 1); 

				// Options        
				for (unsigned l=0; l<labels.size()-1; l++) {

					// Universes
					for (unsigned i=0; i<Weights[l].size(); i++) {
						// if (labels[l] != "PPFXMaster" || labels[l] != "ms_PPFX")  TotWeight[i] *= Weights[l][i]; // don't add masterweight to total 
						Enu_Syst_Window[pdg][l][i]->Fill(mctruth.GetNeutrino().Nu().E(), Weights[l][i]*cv_weight);

						if (intercept) {
							Enu_Syst_AV_TPC[pdg][l][i]->Fill(mctruth.GetNeutrino().Nu().E(), Weights[l][i]*cv_weight);

						}
					}
				}

				// int full_label = labels.size() - 1;

				// Now fill the total weights 
				// for (unsigned int i=0; i<Weights[0].size(); i++){
				// 	Enu_Syst_Window[pdg][full_label][i]->Fill(mctruth.GetNeutrino().Nu().E(), TotWeight[i]*cv_weight);

				// 	if (intercept) {
				// 		Enu_Syst_AV_TPC[pdg][full_label][i]->Fill(mctruth.GetNeutrino().Nu().E(), TotWeight[i]*cv_weight);

				// 	}        
				// }
			}

		} // End loop over mctruth

	} // End loop over events


	// ++++++++++++++++++++++++++++++++
	// Plotting 
	// ++++++++++++++++++++++++++++++++

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

		// Enu_Weight_CV[f]->Write();   
		// Weight_CV[f]->Write();  
		

		// Labels
		for (unsigned int s=1; s<labels.size()+1; s++) {
			std::cout << labels[s-1] << std::endl;
			subdir[f][s][0] = subdir[f][0][0]->mkdir(Form("%s",labels[s-1].c_str()));
			subdir[f][s][0]->cd();

			// Enu_Weight_MS[f][s-1]->Write();  
			// Weight_MS[f][s-1]->Write();  

			// AV/TPC
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

	POTTree->Write();

	output->Close();

	return 0;
}

