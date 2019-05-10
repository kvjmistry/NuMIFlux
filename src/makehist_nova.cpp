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
#include "TRandom.h"
#include "TRandom3.h"
#include "functions_makehist.h"

using namespace art;
using namespace std;

//___________________________________________________________________________
int main(int argc, char** argv) {


	Detector Detector_;
	std::string detector_type = "nova";
	Initialise(detector_type, Detector_);

	geoalgo::GeoAlgo const _geo_algo_instance;

	geoalgo::AABox volAVTPC( Detector_.xRange.first, Detector_.yRange.first, Detector_.zRange.first, Detector_.xRange.second, Detector_.yRange.second, Detector_.zRange.second);

	InputTag mctruths_tag { "flux" };
	InputTag  evtwght_tag { "eventweight" };

	double totalPOT{0};

	vector<string> badfiles;
	vector<string> filename;
	for (int i=1; i<argc; i++) { 
		std::cout << "FILE : " << argv[i] << std::endl; 
		TFile *filein = TFile::Open(argv[i]);
		if (filein->IsZombie()) {
			std::cout << "ERROR: File is ZOMBIE:\t" << argv[i] << std::endl;
			badfiles.push_back(string(argv[i]));
			filein->Close();
		}
		else {
			// std::cout << "FILE : " << argv[i] << std::endl; 
			filename.push_back(string(argv[i]));
			totalPOT+=500000; // 50* 100 000 POT per dk2nu file
			filein->Close();
		}
	}

	std::cout << "\nTotal POT read in:\t" << totalPOT << std::endl;

	std::cout << "\nUsing 500e3 POT per dk2nu file, the flux will be wrong if this is not the case!\n" << std::endl;

	// Histograms for each flavor
	std::vector<TH1D*> Enu_CV_AV_TPC_intersection;
	std::vector<TH1D*> Enu_CV_AV_TPC_detweights;
	std::vector<TH1D*> Enu_CV_AV_TPC_detweights_notilt;


	// Tree for POT counting
	TTree* POTTree = new TTree("POT","Total POT");
	POTTree -> Branch("POT", &totalPOT);
	POTTree -> Fill();

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
	labels = {"PPFXMaster"};

	Enu_CV_AV_TPC_intersection.resize(4);
    Enu_CV_AV_TPC_detweights.resize(4);
    Enu_CV_AV_TPC_detweights_notilt.resize(4);
	
	std::vector<double> temp;

	// Flavors
	for(unsigned i=0; i<flav.size(); i++) {
		int const n = bins[i].size()-1;
		temp.clear();
		temp = bins[i];

		double* bin = &temp[0];

		// FLux histograms
		Enu_CV_AV_TPC_intersection[i] = new TH1D(Form("%s_CV_AV_TPC_intersection",flav[i].c_str()),"",n, bin);
        Enu_CV_AV_TPC_detweights[i] = new TH1D(Form("%s_CV_AV_TPC_detweights",flav[i].c_str()),"",n, bin);
        Enu_CV_AV_TPC_detweights_notilt[i] = new TH1D(Form("%s_CV_AV_TPC_detweights_notilt",flav[i].c_str()),"",n, bin);
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
            double cv_weight_notilt = 1; 
            double window_weight = 1;
            double detwgt; // New weight at a window value
            double tiltwght;
			double enu = mctruth.GetNeutrino().Nu().E();


			// Get the cv_weight
			for (auto last : evtwght.fWeight) {

				if (last.first.find("PPFXCV") != std::string::npos) {

					if(last.second.at(0) > 30 || last.second.at(0) < 0){ // still fill even if bad weight, changed from >90 to >30
						std::cout << "Bad CV weight, setting to 1: " << last.second.at(0) << std::endl;
						cv_weight = 1;
						window_weight = 1;
						cv_weight_notilt = 1; 
					}
					else {
						// std::cout << "CV weight:\t" << last.second.at(0) << std::endl;
						cv_weight        = last.second.at(0);
						window_weight    = last.second.at(0);
						cv_weight_notilt = last.second.at(0);
					}  

				}

			} 

			// Now re-calcuate weight at the window
			// Pick a random point in the TPC (in detector coordinates)
			TVector3 xyz_det = RandomInDet(Detector_);
			
			// From detector to beam coordinates
			TVector3 xyz_beam = FromDetToBeam(xyz_det, false, Detector_);

			// std::cout << "xyz_det:\t" << xyz_det.X() << ", " << xyz_det.Y() << ", " << xyz_det.Z() << "\n"
			// << "\txyz_beam:\t" << xyz_beam.X() << ", " << xyz_beam.Y() << ", " << xyz_beam.Z()<< std::endl;

			// Get the new weight at the detector
			calcEnuWgt(mcflux, xyz_beam, enu, detwgt);

			// Get the tiltweight
			tiltwght = Get_tilt_wgt(xyz_beam, mcflux, enu, Detector_);

			// std::cout << "tiltwgt:\t" << tiltwght << std::endl;

			// Weight of neutrino parent (importance weight) * Neutrino weight for a decay forced at center of near detector 
			cv_weight           *= mcflux.fnimpwt * detwgt / 3.1415926 * tiltwght ; // divide by area of circle equal to pi *r*r
			cv_weight_notilt    *= mcflux.fnimpwt * detwgt / 3.1415926; // for ppfx cases
			// window_weight       *= mcflux.fnimpwt * mcflux.fnwtfar; // mcflux.fnwtfar == mcflux.fnwtnear
			
			// window_weight                   *= mcflux.fnimpwt * Recalc_Intersection_wgt(_geo_algo_instance, volAVTPC, mcflux, mctruth, Detector_ ); // Recalculated for every event, already divide by Pi
			double window_weight_recalc;
			if (intercept) window_weight       *= mcflux.fnimpwt * mcflux.fnwtfar * tiltwght; // mcflux.fnwtfar == mcflux.fnwtnear
			else {
				window_weight_recalc           = Recalc_Intersection_wgt(_geo_algo_instance, volAVTPC, mcflux, mctruth, Detector_ );
				window_weight                 *= mcflux.fnimpwt * window_weight_recalc; // Recalculated for every event
			}
			
			intercept = true; // override above calculations
			
			// Error handling
			if (cv_weight < 0) cv_weight = 0; // get rid of them pesky negative weights
			if (std::isnan(cv_weight) == 1) { // catch NaN values
				std::cout << "got a nan:\t"<<cv_weight <<std::endl;
				cv_weight = 0;
			}

			if (cv_weight_notilt < 0) cv_weight_notilt = 0; // get rid of them pesky negative weights
			if (std::isnan(cv_weight_notilt) == 1) { // catch NaN values
				std::cout << "got a nan:\t"<<cv_weight_notilt <<std::endl;
				cv_weight_notilt = 0;
			}

			if (window_weight < 0) window_weight = 0; // get rid of them pesky negative weights
			if (std::isnan(window_weight) == 1) { // catch NaN values
				std::cout << "got a nan:\t"<<window_weight <<std::endl;
				window_weight = 0;
				
			}

			// Now fill the weights
			if (intercept) Enu_CV_AV_TPC_intersection[pdg]->Fill(enu, window_weight);
			// Enu_CV_AV_TPC_intersection[pdg]->Fill(Enu, window_weight);
			Enu_CV_AV_TPC_detweights[pdg]                 ->Fill(enu, cv_weight);
			Enu_CV_AV_TPC_detweights_notilt[pdg]          ->Fill(enu, cv_weight_notilt);


		} // End loop over mctruth

	} // End loop over events


	// ++++++++++++++++++++++++++++++++
	// Plotting 
	// ++++++++++++++++++++++++++++++++

	TFile* output = new TFile("output.root", "RECREATE");
	TDirectory* savdir = gDirectory;

	std::cout << "flavour.size:\t" <<flav.size()<<std::endl;

	// Top Flav dir 
	std::vector<TDirectory*> subdir(flav.size()); 
	
	// Flavours
	for (unsigned int f=0; f<flav.size(); f++) {
	
		std::cout << "\n" <<flav[f] << std::endl;

		subdir[f]= savdir->mkdir(Form("%s",flav[f].c_str()));
		subdir[f]->cd();

		// Write CV fluxes
		Enu_CV_AV_TPC_intersection[f]->Write();
		Enu_CV_AV_TPC_detweights[f]->Write();
		Enu_CV_AV_TPC_detweights_notilt[f]->Write();
	
		savdir->cd();
	}

	POTTree->Write();

	output->Close();

	return 0;
}

