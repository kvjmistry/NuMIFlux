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
#include "TH2D.h"
#include "TGraph.h"
#include "TH2D.h"
#include "TFile.h"
#include "TRandom.h"
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
#include "functions_makehist.h"

using namespace art;
using namespace std;

// USAGE:
// makehist_uboone_2D <root_files>
// Make sure that the root files contain the correct flux reader module ran over the right geometry

//___________________________________________________________________________
int main(int argc, char** argv) {

	Detector Detector_;
	std::string detector_type = "uboone";
	Initialise(detector_type, Detector_);

	InputTag mctruths_tag { "flux" };
	InputTag  evtwght_tag { "eventweight" };

	double totalPOT{0};

	vector<string> badfiles;
	vector<string> filename;
	
	for (int i = 1; i < argc; i++) { 
		std::cout << "FILE : " << argv[i] << std::endl; 
		
		TFile *filein=new TFile(argv[i]);
		if (filein->IsZombie()) {
			std::cout << "ERROR: File is ZOMBIE:\t" << argv[i] << std::endl;
			badfiles.push_back(string(argv[i]));
			filein->Close();
		}
		else {
			
			filename.push_back(string(argv[i]));
			// totalPOT+=100000; // 100 000 POT per dk2nu file
			totalPOT+=500000;
			filein->Close();
		}
	}

	std::cout << "\nTotal POT read in:\t" << totalPOT << std::endl;

	std::cout << "\nUsing 5e5 POT per dk2nu file, the flux will be wrong if this is not the case!\n" << std::endl;

	// Histograms for each flavor
	std::vector<TH2D*> Enu_Th_CV_AV_TPC;
	std::vector<TH2D*> Enu_Th_UW_AV_TPC;

	// flavors - systematic - universe 
	std::vector<std::vector<std::vector<TH2D*> > > Enu_Th_Syst_AV_TPC;

	// Tree for POT counting
	TTree* POTTree = new TTree("POT","Total POT");
	POTTree -> Branch("POT", &totalPOT);
	POTTree -> Fill();

	// systematic - universe 
	std::vector< std::vector< double > > Weights;   

	std::vector<string> flav = { "numu", "nue", "numubar", "nuebar" };

	std::vector< std::vector<double> > bins; bins.resize(5);
	bins[0] = { // numu
		0.00, 0.025, 0.03, 0.235 ,0.24, 0.50, 1.00, 1.25, 1.50, 1.75, 2.00, 2.25, 2.50, 3.00, 4.00, 5.00, 6.00, 10.00 };

	bins[1] = {  // nue
		0.00 ,0.06, 0.125, 0.25, 0.5, 0.75, 1.00, 1.25, 1.50, 1.75, 2.00, 2.50, 3.00, 3.50, 5.00 };

	bins[2] = {// numubar
		0.00, 0.025, 0.03, 0.235 ,0.24, 0.50, 1.00, 1.25, 1.50, 1.75, 2.00, 2.25, 2.50, 3.00, 4.00, 5.00, 6.00,  10.00 };

	bins[3] = {  // nuebar
		0.00 ,0.06, 0.125,  0.25, 0.5, 0.75, 1.00, 1.25, 1.50, 1.75, 2.00, 2.50, 3.00, 3.50, 5.00 };

	bins[4] = {  0, 20, 110,  160 }; // theta -- removed edge theta bins where no events live and split into 3 bins for stats
	// bins[4] = {  20, 30, 40,50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 160 }; // theta -- removed edge theta bins where no events live OLD Bins

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

	Enu_Th_CV_AV_TPC.resize(4);
	Enu_Th_UW_AV_TPC.resize(4);
	Enu_Th_Syst_AV_TPC.resize(4);
	
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
		Enu_Th_CV_AV_TPC[i] = new TH2D(Form("%s_CV_AV_TPC",flav[i].c_str()),"",n, bin, n_th, bin_th);
		Enu_Th_UW_AV_TPC[i] = new TH2D(Form("%s_unweighted_AV_TPC",flav[i].c_str()),"",n, bin, n_th, bin_th);

		Enu_Th_Syst_AV_TPC[i].resize(labels.size());

		// Labels
		for(unsigned j=0; j<labels.size(); j++) {
			Enu_Th_Syst_AV_TPC[i][j].resize(100);

			// Universes
			for(int k=0; k<100; k++){
				Enu_Th_Syst_AV_TPC[i][j][k] =  new TH2D(Form("%s_%s_Uni_%d_AV_TPC",flav[i].c_str(), labels[j].c_str(), k),"",n, bin, n_th, bin_th);
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
		auto const& evtwghts = *ev.getValidHandle<vector<evwgh::MCEventWeight>>(evtwght_tag);  

		// Loop over MCTruths
		for (size_t i=0; i<mctruths.size(); i++) {
			auto const& mctruth = mctruths.at(i);
			auto const& mcflux = mcfluxs.at(i);
			evwgh::MCEventWeight evtwght;
			evtwght = evtwghts.at(i);
			

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

			double cv_weight = 1; 
			double detwgt; // New weight at a window value
			double tiltwght;

			double Enu = mctruth.GetNeutrino().Nu().E();

			// Now get the momentums to calculate theta
			TVector3 mom_det = {mctruth.GetNeutrino().Nu().Px(),mctruth.GetNeutrino().Nu().Py(),mctruth.GetNeutrino().Nu().Pz()};
			TVector3 mom_beam = FromDetToBeam(mom_det, true, Detector_);

			double costheta = mom_beam.Z() / std::sqrt( mom_beam.X()*mom_beam.X() + mom_beam.Y()*mom_beam.Y() + mom_beam.Z()*mom_beam.Z());
			double theta = std::acos(costheta) * 180 / 3.14159265;
			
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
			
			// Now re-calcuate weight at the window
			// Pick a random point in the TPC (in detector coordinates)
			TVector3 xyz_det = RandomInDet(Detector_);
			
			// From detector to beam coordinates
			TVector3 xyz_beam = FromDetToBeam(xyz_det, false, Detector_);

			// Get the new weight at the detector
			calcEnuWgt(mcflux, xyz_beam, Enu, detwgt);

			// Get the tiltweight -- unused here?
			// tiltwght = Get_tilt_wgt(xyz_beam, mcflux, Enu, Detector_);

			// Weight of neutrino parent (importance weight) * Neutrino weight for a decay forced at center of near detector 
			cv_weight *= mcflux.fnimpwt * detwgt / 3.1415926;

			check_weight(cv_weight);

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

			// ++++++++++++++++++++++++++++++++
			// Now got weights, fill histograms
			// ++++++++++++++++++++++++++++++++

			// TPC AV
			Enu_Th_CV_AV_TPC[pdg]->Fill(Enu, theta, cv_weight);
			Enu_Th_UW_AV_TPC[pdg]->Fill(Enu, theta, mcflux.fnimpwt * detwgt / 3.1415926);

			// Now fill multisims
			// Options        
			for (unsigned l=0; l<labels.size(); l++) {

				// Universes
				for (unsigned i=0; i<Weights[l].size(); i++) {
					Enu_Th_Syst_AV_TPC[pdg][l][i]->Fill(Enu,theta, Weights[l][i]*cv_weight);

				}
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
			subdir[i][j].resize(2);
		}
	
	}

	std::vector<string> cont = { "Active_TPC_Volume", ""};

	// Flavours
	for (unsigned int f=0; f<flav.size(); f++) {
	
		std::cout << "\n" <<flav[f] << std::endl;

		subdir[f][0][0] = savdir->mkdir(Form("%s",flav[f].c_str()));
		subdir[f][0][0]->cd();

		// Write CV fluxes   
		Enu_Th_CV_AV_TPC[f]->Write();
		Enu_Th_UW_AV_TPC[f]->Write();
		
		// Labels
		for (unsigned int s=1; s<labels.size()+1; s++) {
			std::cout << labels[s-1] << std::endl;
			subdir[f][s][0] = subdir[f][0][0]->mkdir(Form("%s",labels[s-1].c_str()));
			subdir[f][s][0]->cd();

			// AV/TPC
			for (int c=1; c<2; c++) {
				std::cout << cont[c-1] << std::endl;
				subdir[f][s][c] = subdir[f][s][0]->mkdir(Form("%s",cont[c-1].c_str()));
				subdir[f][s][c]->cd();

				if (c == 1) {
					for(int i = 0; i < 100; i++){
						Enu_Th_Syst_AV_TPC[f][s-1][i]->Write();
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

