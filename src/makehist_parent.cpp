/*
Galley script to take all the flux read, ppfx weighted art root files
and weights the flux based on the ppfx weights and importance weights etc.

This specific script overwrites the window flux method to compare the flux at
the detector (smeared). It also breaks the flux down by parent to give the
oppertunuty to investigate the flux by parent.

All the window method and detector smeared weights are preserved. 

* Authors: J. Zennamo, A. Mastbaum, K. Mistry
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
// makehist_parent <detector_type> <root_files>
// <detector_type>: uboone OR nova
// Make sure that the root files contain the correct flux reader module ran over the right geometry
//___________________________________________________________________________
int main(int argc, char** argv) {

	Detector Detector_;
	
	InputTag mctruths_tag { "flux" };
	InputTag  evtwght_tag { "eventweight" };

	double totalPOT{0};
	bool input_flag{false}; // flag to see if a detector has been specified

	vector<string> badfiles;
	vector<string> filename;
	for (int i = 1; i < argc; i++) {

		std::string input = string(argv[i]);

		// Initialise the detector type 
		if (input == "uboone" || input == "nova" ){
			std::string detector_type = string(argv[i]);
			Initialise(detector_type, Detector_);
			input_flag = true;
			continue;
		}

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

	std::cout << "\nUsing 5e5 POT per dk2nu file, the flux will be wrong if this is not the case!\n" << std::endl;

	// If no detector is specified then use uboone as default
	if (input_flag == false)
		Initialise("uboone", Detector_);

	// Set up intersection method stuff
	geoalgo::GeoAlgo const _geo_algo_instance;
	geoalgo::AABox volAVTPC( Detector_.xRange.first, Detector_.yRange.first, Detector_.zRange.first, Detector_.xRange.second, Detector_.yRange.second, Detector_.zRange.second);

	// Histograms for each flavor
	std::vector<TH1D*> Enu_CV_Window;
	std::vector<TH1D*> Enu_CV_AV_TPC;
	std::vector<TH1D*> Enu_UW_Window;
	std::vector<TH1D*> Enu_UW_AV_TPC;

	std::vector<TH1D*> Th_CV_AV_TPC;
	std::vector<TH1D*> Th_UW_AV_TPC;

	// 5Mev Bins
	std::vector<TH1D*> Enu_CV_Window_5MeV_bin;
	std::vector<TH1D*> Enu_CV_AV_TPC_5MeV_bin;
	std::vector<TH1D*> Enu_UW_Window_5MeV_bin;
	std::vector<TH1D*> Enu_UW_AV_TPC_5MeV_bin;

	// Detector intersection window method
	std::vector<TH1D*> Enu_CV_Window_5MeV_bin_intersect;

	// Flux by Parent
	std::vector<string> parent = {"PI_Plus", "PI_Minus", "Mu_Plus", "Mu_Minus", "Kaon_Plus", "Kaon_Minus" , "K0L"};
	std::vector<std::vector<TH1D*> > Enu_Parent_AV_TPC;		// Flux by Parent
	std::vector<std::vector<TH1D*> > Th_Parent_AV_TPC; 		// Flux by parent in theta
	std::vector<std::vector<TH1D*> > zpos_Parent_AV_TPC; 	// Flux by parent in z Pos at decay (production is unimplemented for now)
	std::vector<std::vector<TH1D*> > dk_Parent_mom; 		// Decay momentum by parent
	std::vector<std::vector<TH1D*> > impwght_Parent; 		// Importance weight by parent
	std::vector<std::vector<TH1D*> > Prod_energy_Parent; 	// Production energy by parent
	std::vector<std::vector<TH1D*> > Targ_mom_Parent; 	    // Momentum by parent as it leaves the target
	std::vector<std::vector<TH1D*> > DAR_Enu_Parent; 	    // Energy spectrum of decay at rest particles
	
	// Other hists
	TH1D* NuMu_PiDAR_zpos      = new TH1D("NuMu_PiDAR_zpos","", 400 , 0, 80000); // numu Pidar peak zpos
	TH1D* NuMu_KDAR_zpos       = new TH1D("NuMu_KDAR_zpos","", 400 , 0, 80000);  // numu kdar peak zpos
	TH1D* NuMu_peak_mom_muon   = new TH1D("NuMu_peak_mom_muon","", 100 , 0, 25);  // muon parent momentum at large enu peak
	TH1D* NuMu_peak_theta_muon = new TH1D("NuMu_peak_theta_muon","", 40 , 0, 180);  // muon parent thetaat large peak
	TH1D* NuMu_peak_zpos_muon  = new TH1D("NuMu_peak_zpos_muon","", 400 , 0, 80000);  // muon parent thetaat large peak

	// Tree for POT counting
	TTree* POTTree = new TTree("POT","Total POT");
	POTTree -> Branch("POT", &totalPOT);
	POTTree -> Fill();

	// systematic - universe 
	std::vector< std::vector< double > > Weights;   

	std::vector<string> flav = { "numu", "nue", "numubar", "nuebar" };

	std::vector< std::vector<double> > bins; bins.resize(4);
	
	if (Detector_.detector_name == "uboone"){
		bins[0] = { // numu
			0.00, 0.025, 0.03, 0.235 ,0.24, 0.50, 1.00, 1.25, 1.50, 1.75, 2.00, 2.25, 2.50, 3.00, 4.00, 5.00, 6.00, 7.00, 10.00 };

		bins[1] = {  // nue
			0.00 ,0.06, 0.125, 0.25, 0.5, 0.75, 1.00, 1.25, 1.50, 1.75, 2.00, 2.50, 3.00, 3.50, 4.00, 5.00 };

		bins[2] = {// numubar
			0.00, 0.025, 0.03, 0.235 ,0.24, 0.50, 1.00, 1.25, 1.50, 1.75, 2.00, 2.25, 2.50, 3.00, 4.00, 5.00, 6.00, 7.00, 10.00 };

		bins[3] = {  // nuebar
			0.00 ,0.06, 0.125,  0.25, 0.5, 0.75, 1.00, 1.25, 1.50, 1.75, 2.00, 2.50, 3.00, 3.50, 4.00, 5.00 };
	}
	else {
		bins[0] = {
		0.0,  0.6,  0.8,  1.0,  1.2,  1.4,  1.6,  1.8,  2.0,  2.2,  2.4,  2.6,
		2.8,  3.0,  3.2,  3.4,  3.6,  3.8,  4.0,  4.2,  4.4,  4.6,  4.8,  5.0,
		6.0,  7.0,  8.0,  9.0,  10.0,  11.0,  12.0,  13.0,  14.0,  15.0,  16.0,
		17.0,  18.0,  19.0,  20.0
		};
		bins[1] = {  0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0};
		bins[2] = {  0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0};
		bins[3] = {  0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0};
	}

	std::vector<string> labels;
	// labels = {"ms_PPFX","Total"};
	labels = {"PPFXMaster"};

	Weights.resize(labels.size());

	for (unsigned int i=0; i<labels.size(); i++) {
		Weights[i].resize(100);
	}

	Enu_CV_Window.resize(4);
	Enu_CV_AV_TPC.resize(4);
	Enu_UW_Window.resize(4);
	Enu_UW_AV_TPC.resize(4);

	Th_CV_AV_TPC.resize(4);
	Th_UW_AV_TPC.resize(4);
	
	Enu_CV_Window_5MeV_bin.resize(4);
	Enu_CV_AV_TPC_5MeV_bin.resize(4);
	Enu_UW_Window_5MeV_bin.resize(4);
	Enu_UW_AV_TPC_5MeV_bin.resize(4);

	Enu_CV_Window_5MeV_bin_intersect.resize(4);

	Enu_Parent_AV_TPC.resize(4);
	Th_Parent_AV_TPC.resize(4);
	zpos_Parent_AV_TPC.resize(4);
	impwght_Parent.resize(4);
	Targ_mom_Parent.resize(4);
	DAR_Enu_Parent.resize(4);
	
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
		Enu_UW_Window[i] = new TH1D(Form("%s_UW_Window",flav[i].c_str()),"",n, bin);
		Enu_UW_AV_TPC[i] = new TH1D(Form("%s_UW_AV_TPC",flav[i].c_str()),"",n, bin);
		Th_CV_AV_TPC [i] = new TH1D(Form("Th_%s_CV_TPC", flav[i].c_str()), "", 18, 0, 180);
		Th_UW_AV_TPC [i] = new TH1D(Form("Th_%s_UW_TPC", flav[i].c_str()), "", 18, 0, 180);

		// new binning schmeme to be the same as marcos
		Enu_CV_Window_5MeV_bin[i] = new TH1D(Form("%s_CV_Window_5MeV_bin",flav[i].c_str()),"",4000, 0, 20);
		Enu_CV_AV_TPC_5MeV_bin[i] = new TH1D(Form("%s_CV_AV_TPC_5MeV_bin",flav[i].c_str()),"",4000, 0, 20);
		Enu_UW_Window_5MeV_bin[i] = new TH1D(Form("%s_UW_Window_5MeV_bin",flav[i].c_str()),"",4000, 0, 20);
		Enu_UW_AV_TPC_5MeV_bin[i] = new TH1D(Form("%s_UW_AV_TPC_5MeV_bin",flav[i].c_str()),"",4000, 0, 20);

		Enu_CV_Window_5MeV_bin_intersect[i] = new TH1D(Form("%s_CV_Window_5MeV_bin_intersect",flav[i].c_str()),"",4000, 0, 20);

		Enu_Parent_AV_TPC[i].resize(parent.size());
		Th_Parent_AV_TPC[i].resize(parent.size());
		zpos_Parent_AV_TPC[i].resize(parent.size());
		impwght_Parent[i].resize(parent.size());
		Targ_mom_Parent[i].resize(parent.size());
		DAR_Enu_Parent[i].resize(parent.size());
		
		// Parent
		for(unsigned k = 0; k < parent.size(); k++){
			Enu_Parent_AV_TPC[i][k]  = new TH1D(Form("Enu_%s_%s_AV_TPC",     flav[i].c_str(), parent[k].c_str()),"", 4000, 0, 20);
			Th_Parent_AV_TPC[i][k]   = new TH1D(Form("Th_%s_%s_AV_TPC",      flav[i].c_str(), parent[k].c_str()),"", 18 , 0, 180);
			zpos_Parent_AV_TPC[i][k] = new TH1D(Form("zpos_%s_%s_AV_TPC",    flav[i].c_str(), parent[k].c_str()),"", 400 , 0, 80000);
			impwght_Parent[i][k]     = new TH1D(Form("impwght_Parent_%s_%s", flav[i].c_str(), parent[k].c_str()),"", 1 , 0, -1);
			Targ_mom_Parent[i][k]    = new TH1D(Form("Targ_mom_Parent_%s_%s",flav[i].c_str(), parent[k].c_str()),"", 4000, 0, 20);
			DAR_Enu_Parent[i][k]     = new TH1D(Form("DAR_Enu_%s_%s_AV_TPC", flav[i].c_str(), parent[k].c_str()),"", 4000, 0, 20);
		}
		
	}

	// ++++++++++++++++++++++++++++++++
	// Event loop
	// ++++++++++++++++++++++++++++++++
	std::cout << "Starting Eventloop" << std::endl;
	bool disperrors{false}; // chose whether or not to display the unphysical nuetrino errors

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
						<< std::endl;
				}
				continue;
			}

			// Do the neutrino ray calculation for flux at the window
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


			double cv_weight = 1;    // PPFX CV
			double dk2nu_weight = 1; // UW 
			double detwgt; // New weight at a window value
			double tiltwght;
			double window_weight = 1;
			double dk2nu_window_weight = 1;
			
			// Now get the momentums to calculate theta
			TVector3 mom_det = {mctruth.GetNeutrino().Nu().Px(),mctruth.GetNeutrino().Nu().Py(),mctruth.GetNeutrino().Nu().Pz()};
			TVector3 mom_beam = FromDetToBeam(mom_det, true, Detector_);

			double costheta = mom_beam.Z() / std::sqrt( mom_beam.X()*mom_beam.X() + mom_beam.Y()*mom_beam.Y() + mom_beam.Z()*mom_beam.Z());
			double theta = std::acos(costheta) * 180 / 3.14159265;
			// double costheta = mctruth.GetNeutrino().Nu().Pz() / mctruth.GetNeutrino().Nu().P(); // I believe these are not correct
			// double theta = std::acos(costheta) * 180 / 3.14159265;

			double Enu = mctruth.GetNeutrino().Nu().E();
			double Pmom_dk = std::sqrt( mcflux.fpdpz*mcflux.fpdpz + mcflux.fpdpy*mcflux.fpdpy + mcflux.fpdpx*mcflux.fpdpx ); // Parent momentum at decay
			double Pmom_tg = std::sqrt(mcflux.ftpx*mcflux.ftpx + mcflux.ftpy*mcflux.ftpy + mcflux.ftpz*mcflux.ftpz);         // Parent moment

			// Get the ppfx cv_weight
			for (auto last : evtwght.fWeight) {

				if (last.first.find("PPFXCV") != std::string::npos) {

					if(last.second.at(0) > 30 || last.second.at(0) < 0){ // still fill even if bad weight, changed from >90 to >30
					// if(last.second.at(0) < 0){ // change this to only throw out negative weights
						std::cout << "Bad CV weight, setting to 1: " << last.second.at(0) << std::endl;
						cv_weight = 1;
						window_weight = 1;
						// cv_weight = last.second.at(0);
					}
					else {
						// std::cout << "CV weight:\t" << last.second.at(0) << std::endl;
						cv_weight     = last.second.at(0);
						window_weight = last.second.at(0);
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

			// Get the tiltweight
			tiltwght = Get_tilt_wgt(xyz_beam, mcflux, Enu, Detector_);

			// Weight of neutrino parent (importance weight) * Neutrino weight for a decay forced at center of near detector 
			cv_weight        *= mcflux.fnimpwt * detwgt / 3.1415926; // for ppfx cases
			dk2nu_weight     *= mcflux.fnimpwt * detwgt / 3.1415926; // for UW cases 
			
			window_weight       *= mcflux.fnimpwt * mcflux.fnwtfar * tiltwght ; // mcflux.fnwtfar == mcflux.fnwtnear
			dk2nu_window_weight *= mcflux.fnimpwt * mcflux.fnwtfar * tiltwght ; // mcflux.fnwtfar == mcflux.fnwtnear
			
			// Error handling, sets to zero if bad
			check_weight(cv_weight);
			check_weight(window_weight);
			check_weight(dk2nu_weight);
			check_weight(dk2nu_window_weight);
			
			// ++++++++++++++++++++++++++++++++
			// Now got weights, fill histograms
			// ++++++++++++++++++++++++++++++++

			// Window
			Enu_CV_Window[pdg]      ->Fill(Enu, window_weight);
			Enu_UW_Window[pdg]      ->Fill(Enu, dk2nu_window_weight);
			Enu_CV_Window_5MeV_bin[pdg]->Fill(Enu, window_weight);
			Enu_UW_Window_5MeV_bin[pdg]->Fill(Enu, dk2nu_window_weight);

			if (intercept) Enu_CV_Window_5MeV_bin_intersect[pdg] ->Fill(Enu, window_weight); // Flux that intersects the detector

			// TPC AV
			Enu_CV_AV_TPC[pdg]      ->Fill(Enu, cv_weight);
			Enu_UW_AV_TPC[pdg]      ->Fill(Enu, dk2nu_weight);
			Enu_CV_AV_TPC_5MeV_bin[pdg]->Fill(Enu, cv_weight);
			Enu_UW_AV_TPC_5MeV_bin[pdg]->Fill(Enu, dk2nu_weight);
			Th_CV_AV_TPC[pdg]->Fill(theta, cv_weight);
			Th_UW_AV_TPC[pdg]->Fill(theta, dk2nu_weight);

			// INDEXING: 0: PI_Plus 1: PI_Minus 2: Mu Plus 3: Mu_Minus 4: Kaon_Plus 5: Kaon_Minus 6: K0L 
			if (mcflux.fptype == 211){ // pi plus
				Enu_Parent_AV_TPC[pdg][0]  ->Fill(Enu, cv_weight);
				Th_Parent_AV_TPC[pdg][0]   ->Fill(theta, cv_weight);
				zpos_Parent_AV_TPC[pdg][0] ->Fill(mcflux.fvz, cv_weight);
				impwght_Parent[pdg][0]     ->Fill(mcflux.fnimpwt);
				Targ_mom_Parent[pdg][0]    ->Fill(Pmom_tg, cv_weight);

				// Fill DAR Energy spectrum
				if (Pmom_dk == 0 ) DAR_Enu_Parent[pdg][0]->Fill(Enu, cv_weight);
				
			}
			else if (mcflux.fptype == -211){ // pi minus
				Enu_Parent_AV_TPC[pdg][1]  ->Fill(Enu, cv_weight);
				Th_Parent_AV_TPC[pdg][1]   ->Fill(theta, cv_weight);
				zpos_Parent_AV_TPC[pdg][1] ->Fill(mcflux.fvz, cv_weight);
				impwght_Parent[pdg][1]     ->Fill(mcflux.fnimpwt);
				Targ_mom_Parent[pdg][1]    ->Fill(Pmom_tg, cv_weight);

				// Fill DAR Energy spectrum
				if (Pmom_dk == 0 ) DAR_Enu_Parent[pdg][1]->Fill(Enu, cv_weight);
				
			}
			else if (mcflux.fptype == -13){ // mu plus
				Enu_Parent_AV_TPC[pdg][2]  ->Fill(Enu, cv_weight);
				Th_Parent_AV_TPC[pdg][2]   ->Fill(theta, cv_weight);
				zpos_Parent_AV_TPC[pdg][2] ->Fill(mcflux.fvz, cv_weight);
				impwght_Parent[pdg][2]     ->Fill(mcflux.fnimpwt);
				Targ_mom_Parent[pdg][2]    ->Fill(Pmom_tg, cv_weight);
	
				// Fill DAR Energy spectrum
				if (Pmom_dk == 0 ) DAR_Enu_Parent[pdg][2]->Fill(Enu, cv_weight);

			} 
			else if (mcflux.fptype == 13){ // mu minus
				Enu_Parent_AV_TPC[pdg][3]  ->Fill(Enu, cv_weight);
				Th_Parent_AV_TPC[pdg][3]   ->Fill(theta, cv_weight);
				zpos_Parent_AV_TPC[pdg][3] ->Fill(mcflux.fvz, cv_weight);
				impwght_Parent[pdg][3]     ->Fill(mcflux.fnimpwt);
				Targ_mom_Parent[pdg][3]    ->Fill(Pmom_tg, cv_weight);

				// Fill DAR Energy spectrum
				if (Pmom_dk == 0 ) DAR_Enu_Parent[pdg][3]->Fill(Enu, cv_weight);

			} 
			else if (mcflux.fptype == 321){ // K+
				Enu_Parent_AV_TPC[pdg][4]  ->Fill(Enu, cv_weight);
				Th_Parent_AV_TPC[pdg][4]   ->Fill(theta, cv_weight);
				zpos_Parent_AV_TPC[pdg][4] ->Fill(mcflux.fvz, cv_weight);
				impwght_Parent[pdg][4]     ->Fill(mcflux.fnimpwt);
				Targ_mom_Parent[pdg][4]    ->Fill(Pmom_tg, cv_weight);

				// Fill DAR Energy spectrum
				if (Pmom_dk == 0 ) DAR_Enu_Parent[pdg][4]->Fill(Enu, cv_weight);
				
			}
			else if ( mcflux.fptype == -321){ // K-
				Enu_Parent_AV_TPC[pdg][5]  ->Fill(Enu, cv_weight);
				Th_Parent_AV_TPC[pdg][5]   ->Fill(theta, cv_weight);
				zpos_Parent_AV_TPC[pdg][5] ->Fill(mcflux.fvz, cv_weight);
				impwght_Parent[pdg][5]     ->Fill(mcflux.fnimpwt);
				Targ_mom_Parent[pdg][5]    ->Fill(Pmom_tg, cv_weight);

				// Fill DAR Energy spectrum
				if (Pmom_dk == 0 ) DAR_Enu_Parent[pdg][5]->Fill(Enu, cv_weight);
				
			}
			else if (mcflux.fptype == 130 ){ // K0L
				Enu_Parent_AV_TPC[pdg][6]  ->Fill(Enu, cv_weight);
				Th_Parent_AV_TPC[pdg][6]   ->Fill(theta, cv_weight);
				zpos_Parent_AV_TPC[pdg][6] ->Fill(mcflux.fvz, cv_weight);
				impwght_Parent[pdg][6]     ->Fill(mcflux.fnimpwt);
				Targ_mom_Parent[pdg][6]    ->Fill(Pmom_tg, cv_weight);

				// Fill DAR Energy spectrum
				if (Pmom_dk == 0 ) DAR_Enu_Parent[pdg][6]->Fill(Enu, cv_weight);

			}
			
			if (Enu > 0.025 && Enu < 0.03 && mcflux.fptype == 211 ){ // Pidar peak
				NuMu_PiDAR_zpos->Fill(mcflux.fvz, cv_weight);
			}
			if (Enu > 0.235 && Enu < 0.24 && (mcflux.fptype == 321 || mcflux.fptype == -321)){ // kdar peak
				NuMu_KDAR_zpos->Fill(mcflux.fvz, cv_weight);
			}

			// Look in the energy peak for muons
			if (Enu > 0.074 && Enu < 0.082 && mcflux.fptype == 13 ){
				NuMu_peak_mom_muon   ->Fill(Pmom_dk,cv_weight );
				NuMu_peak_theta_muon ->Fill(theta, cv_weight );
				NuMu_peak_zpos_muon  ->Fill(mcflux.fvz,cv_weight );
				// std::cout << "E:\t" << Enu << "Parent:\t" << mcflux.fptype <<"  theta:\t" <<theta<<"   ntype:\t"<< mcflux.fndecay<<   std::endl;
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
	std::vector<std::vector<TDirectory*>> subdir(flav.size()); 
	
	//Create label dirs
	for (unsigned i=0; i<flav.size(); i++){
		subdir[i].resize(parent.size()+2);
	
	}

	// Flavours
	for (unsigned int f=0; f<flav.size(); f++) {
	
		std::cout << "\n" <<flav[f] << std::endl;

		subdir[f][0] = savdir->mkdir(Form("%s",flav[f].c_str()));
		subdir[f][0]->cd();

		// Write CV fluxes
		Enu_CV_Window[f]->Write();      
		Enu_CV_AV_TPC[f]->Write();
		Enu_UW_Window[f]->Write();      
		Enu_UW_AV_TPC[f]->Write();
		Th_CV_AV_TPC[f]->Write();
		Th_UW_AV_TPC[f]->Write();

		Enu_CV_Window_5MeV_bin[f]->Write();      
		Enu_CV_AV_TPC_5MeV_bin[f]->Write();
		Enu_UW_Window_5MeV_bin[f]->Write();      
		Enu_UW_AV_TPC_5MeV_bin[f]->Write();
		
		Enu_CV_Window_5MeV_bin_intersect[f]->Write();

		// Parent
		// INDEXING: 1: PI_Plus 2: PI_Minus 3: Mu Plus 4: Mu_Minus 5: Kaon_Plus 6: Kaon_Minus 7: K0L 
		for(unsigned int k = 1; k < parent.size()+1; k++){
			std::cout << parent[k-1] << std::endl;
			subdir[f][k] = subdir[f][0]->mkdir(Form("%s",parent[k-1].c_str()));
			subdir[f][k]->cd();

			Enu_Parent_AV_TPC[f][k-1]->Write();
			Th_Parent_AV_TPC[f][k-1]->Write();
			zpos_Parent_AV_TPC[f][k-1]->Write();
			impwght_Parent[f][k-1]->Write();
			Targ_mom_Parent[f][k-1]->Write();
			DAR_Enu_Parent[f][k-1]->Write();
	
		}

		// Make other plots folder for miscalanious variables
		std::cout << "OtherPlots" << std::endl;
		subdir.at(f).at(parent.size()+1) = subdir[f][0]->mkdir("OtherPlots");
		subdir.at(f).at(parent.size()+1)->cd();
		
		NuMu_PiDAR_zpos->Write();
		NuMu_KDAR_zpos->Write();
		NuMu_peak_mom_muon->Write();
		NuMu_peak_theta_muon->Write();
		NuMu_peak_zpos_muon->Write();

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
