/*
Galley script to take all the flux read, ppfx weighted art root files
and weights the flux based on the ppfx weights and importance weights etc.

This specific script overwrites the window flux method to compare the flux at
the detector (smeared). It also breaks the flux down by parent to give the
oppertunuty to investigate the flux by parent.

All the window method and detector smeared weights are preserved. 

This script still hardcodes the POT. It is locked away in the subrun
sumdata::POTsummary dataproduct which I have no idea how to access in 
gallery...

* Authors: J. Zennamo, A. Mastbaum, K. Mistry
*/

// Helper functions and some histogram defintions are found in functions_makehist.h
#include "functions_makehist.h"

using namespace art;
using namespace std;
using namespace std::chrono;

// USAGE:
// makehist_parent <detector_type> <root_files>
// <detector_type>: uboone OR nova DEFAULT IS uboone
// Make sure that the root files contain the correct flux reader module ran over the right geometry
//___________________________________________________________________________
int main(int argc, char** argv) {

	auto start = high_resolution_clock::now(); // Start time of script

	Detector Detector_;
	InputTag mctruths_tag { "flux" };
	InputTag  evtwght_tag { "eventweight" };

	// labels = {"ms_PPFX","Total"};
	//labels = {"PPFXMaster"};
	labels = {"PPFXMIPPKaon","PPFXMIPPPion","PPFXOther","PPFXTargAtten","PPFXThinKaon","PPFXThinMeson","PPFXThinNeutron","thinna_PPFXThinNucA","thinn_PPFXThinNuc","PPFXThinPion","PPFXTotAbsorp","PPFXMaster"};

	// Loop over input arguments
	for (int i = 1; i < argc; i++) {

		std::string input = string(argv[i]);

		// Initialise the detector type 
		if (input == "uboone" || input == "nova" || input == "minerva" ){
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

	// Resizing of histograms
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

	// 2D
	Enu_Th_CV_AV_TPC.resize(4);
	Enu_Th_UW_AV_TPC.resize(4);

	// Other
	parent_mom.resize(4);
	parent_angle.resize(4);
	parent_zpos.resize(4); 
	PiDAR_zpos.resize(4);
	KDAR_zpos.resize(4); 
	MuDAR_zpos.resize(4);
	parent_zpos_angle.resize(4);
	parent_zpos_angle_energy.resize(4);
	flux_targ.resize(4);
	flux_pipe.resize(4);
	flux_dump.resize(4);
	flux_targ_theta.resize(4);
	flux_pipe_theta.resize(4);
	flux_dump_theta.resize(4);

	// Weighted
	Enu_Syst_AV_TPC.resize(4);    // 1D
	Enu_Th_Syst_AV_TPC.resize(4); // 2D
	
	// Tree for POT counting
	TTree* POTTree = new TTree("POT","Total POT");
	POTTree -> Branch("POT", &totalPOT);
	POTTree -> Fill();

	
	// Resize the weight labels
	Weights.resize(labels.size());
	for (unsigned int i=0; i<labels.size(); i++) Weights[i].resize(100);

	int const n_th = Detector_.bins.at(4).size()-1; // theta bins
	temp2 = Detector_.bins.at(4); 

	// Flavors
	for(unsigned i=0; i<flav.size(); i++) {
		
		int const n = Detector_.bins.at(i).size()-1;
		temp.clear();
		temp = Detector_.bins.at(i);

		double* bin = &temp[0];
		double* bin_th = &temp2[0];

		// FLux histograms
		Enu_CV_Window[i] = new TH1D(Form("%s_CV_Window",flav[i].c_str()),"",n, bin);
		Enu_CV_AV_TPC[i] = new TH1D(Form("%s_CV_AV_TPC",flav[i].c_str()),"",n, bin);
		Enu_UW_Window[i] = new TH1D(Form("%s_UW_Window",flav[i].c_str()),"",n, bin);
		Enu_UW_AV_TPC[i] = new TH1D(Form("%s_UW_AV_TPC",flav[i].c_str()),"",n, bin);
		Th_CV_AV_TPC [i] = new TH1D(Form("Th_%s_CV_TPC", flav[i].c_str()), "", 18, 0, 180);
		Th_UW_AV_TPC [i] = new TH1D(Form("Th_%s_UW_TPC", flav[i].c_str()), "", 18, 0, 180);

		// Other Histograms
		parent_mom[i]   = new TH1D(Form("%s_parent_mom",   flav[i].c_str()), "", 100, 0, 25);    // momentum distribution of nu parent
		parent_angle[i] = new TH1D(Form("%s_parent_angle", flav[i].c_str()), "", 180, 0, 180);   // angle distribution of nu parent
		parent_zpos[i]  = new TH1D(Form("%s_parent_zpos",  flav[i].c_str()), "", 400, 0, 80000); // zpos distribution of nu parent
		PiDAR_zpos[i]   = new TH1D(Form("%s_PiDAR_zpos",   flav[i].c_str()),"",  400, 0, 80000); // Pidar peak zpos
		KDAR_zpos[i]    = new TH1D(Form("%s_KDAR_zpos",    flav[i].c_str()), "", 400, 0, 80000); // Kdar peak zpos 
		MuDAR_zpos[i]   = new TH1D(Form("%s_MuDAR_zpos",   flav[i].c_str()),"",  400, 0, 80000); // Mudar peak zpos
		parent_zpos_angle[i]        = new TH2D(Form("%s_parent_zpos_angle", flav[i].c_str()), "", 200, 0, 80000, 180 , 0, 180);
		parent_zpos_angle_energy[i] = new TH2D(Form("%s_parent_zpos_angle_energy", flav[i].c_str()), "", 200, 0, 80000, 180 , 0, 180);
		flux_targ[i] = new TH1D( Form("%s_flux_targ", flav[i].c_str()),"", 4000, 0, 20 );
		flux_pipe[i] = new TH1D( Form("%s_flux_pipe", flav[i].c_str()),"", 4000, 0, 20 );
		flux_dump[i] = new TH1D( Form("%s_flux_dump", flav[i].c_str()),"", 4000, 0, 20 );
		flux_targ_theta[i] = new TH1D( Form("%s_flux_targ_theta", flav[i].c_str()),"", 4000, 0, 20 );
		flux_pipe_theta[i] = new TH1D( Form("%s_flux_pipe_theta", flav[i].c_str()),"", 4000, 0, 20 );
		flux_dump_theta[i] = new TH1D( Form("%s_flux_dump_theta", flav[i].c_str()),"", 4000, 0, 20 );

		// 2D
		Enu_Th_CV_AV_TPC[i] = new TH2D(Form("%s_CV_AV_TPC_2D",flav[i].c_str()),"",n, bin, n_th, bin_th);
		Enu_Th_UW_AV_TPC[i] = new TH2D(Form("%s_unweighted_AV_TPC_2D",flav[i].c_str()),"",n, bin, n_th, bin_th);

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
			Th_Parent_AV_TPC[i][k]   = new TH1D(Form("Th_%s_%s_AV_TPC",      flav[i].c_str(), parent[k].c_str()),"", 18,   0, 180);
			zpos_Parent_AV_TPC[i][k] = new TH1D(Form("zpos_%s_%s_AV_TPC",    flav[i].c_str(), parent[k].c_str()),"", 400,  0, 80000);
			impwght_Parent[i][k]     = new TH1D(Form("impwght_Parent_%s_%s", flav[i].c_str(), parent[k].c_str()),"", 1000, 0, 1000);
			Targ_mom_Parent[i][k]    = new TH1D(Form("Targ_mom_Parent_%s_%s",flav[i].c_str(), parent[k].c_str()),"", 4000, 0, 20);
			DAR_Enu_Parent[i][k]     = new TH1D(Form("DAR_Enu_%s_%s_AV_TPC", flav[i].c_str(), parent[k].c_str()),"", 4000, 0, 20);
		}
		
		// Weighted Stuff
		Enu_Syst_AV_TPC[i].resize(labels.size());
		Enu_Th_Syst_AV_TPC[i].resize(labels.size());
		
		// Labels
		for(unsigned j=0; j<labels.size(); j++) {
			Enu_Syst_AV_TPC[i][j].resize(100);
			Enu_Th_Syst_AV_TPC[i][j].resize(100);

			// Universes
			for(int k=0; k<100; k++){
				Enu_Syst_AV_TPC[i][j][k]    =  new TH1D(Form("%s_%s_Uni_%d_AV_TPC",flav[i].c_str(), labels[j].c_str(), k),"",n, bin);
				Enu_Th_Syst_AV_TPC[i][j][k] =  new TH2D(Form("%s_%s_Uni_%d_AV_TPC_2D",flav[i].c_str(), labels[j].c_str(), k),"",n, bin, n_th, bin_th);
			}
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
		auto const& mcfluxs  = *ev.getValidHandle<vector<simb::MCFlux>>(mctruths_tag);   
		auto const& evtwghts = *ev.getValidHandle<vector<evwgh::MCEventWeight>>(evtwght_tag);

		// Loop over MCTruths
		for (size_t i=0; i<mctruths.size(); i++) {
			auto const& mctruth = mctruths.at(i);
			auto const& mcflux  = mcfluxs.at(i);
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

			// Check if the ray intercepted with the detector
			bool intercept = check_intercept(_geo_algo_instance, volAVTPC, mctruth);

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
			double KRDET       =  100.0; // Radius of circle in cm
			double KRDET_Area  = 3.1415926*KRDET*KRDET/10000.0; // Area of circle in m2
			calcEnuWgt(mcflux, xyz_beam, Enu, detwgt, KRDET);

			// Get the tiltweight
			tiltwght = Get_tilt_wgt(xyz_beam, mcflux, Enu, Detector_);

			// if (tiltwght > 0.95) std::cout << "tiltweight:\t" << std::setprecision(5) << tiltwght << "\tzpos:\t" << mcflux.fvz <<"\ttheta:\t" << theta << std::endl;

			// Window weight recalculations
			double window_weight_recalc;
			if (intercept){
				window_weight       *= mcflux.fnimpwt * mcflux.fnwtfar * tiltwght / KRDET_Area; // mcflux.fnwtfar == mcflux.fnwtnear
				dk2nu_window_weight *= mcflux.fnimpwt * mcflux.fnwtfar * tiltwght / KRDET_Area; // mcflux.fnwtfar == mcflux.fnwtnear
			} 
			else {
				// Only do this for nova since it doesnt really work that well for uboone with the tiltweight and stuff
				if (Detector_.detector_name == "nova"){
					window_weight_recalc           = Recalc_Intersection_wgt(_geo_algo_instance, volAVTPC, mcflux, mctruth, Detector_, KRDET, Enu );
					window_weight                 *= mcflux.fnimpwt * window_weight_recalc / KRDET_Area; // Recalculated for every event
					dk2nu_window_weight           *= mcflux.fnimpwt * window_weight_recalc / KRDET_Area; // Recalculated for every event
					intercept = true; // override above calculations
				}
			}
			// Uncomment this if wanting to re-write all the window weights
			// window_weight_recalc           = Recalc_Intersection_wgt(_geo_algo_instance, volAVTPC, mcflux, mctruth, Detector_, KRDET, Enu );
			// window_weight                 *= mcflux.fnimpwt * window_weight_recalc / KRDET_Area; // Recalculated for every event
			// dk2nu_window_weight           *= mcflux.fnimpwt * window_weight_recalc / KRDET_Area; // Recalculated for every event
			// intercept = true; // override above calculations


			// Weight of neutrino parent (importance weight) * Neutrino weight for a decay forced at center of near detector 
			cv_weight        *= mcflux.fnimpwt * detwgt / KRDET_Area; // for ppfx cases
			dk2nu_weight     *= mcflux.fnimpwt * detwgt / KRDET_Area; // for UW cases 
			
			// Error handling, sets to zero if bad
			check_weight(cv_weight);
			check_weight(window_weight);
			check_weight(dk2nu_weight);
			check_weight(dk2nu_window_weight);

			// Fill Weight vector with 1's to create size=labels
			for (unsigned l=0; l<labels.size(); l++) {
				std::fill(Weights[l].begin(), Weights[l].end(), 1);
			}

			// Get the PPFX weights
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

			// Window
			Enu_CV_Window[pdg]              ->Fill(Enu, window_weight);
			Enu_UW_Window[pdg]              ->Fill(Enu, dk2nu_window_weight);
			Enu_CV_Window_5MeV_bin[pdg]     ->Fill(Enu, window_weight);
			Enu_UW_Window_5MeV_bin[pdg]     ->Fill(Enu, dk2nu_window_weight);

			if (intercept) Enu_CV_Window_5MeV_bin_intersect[pdg] ->Fill(Enu, window_weight); // Flux that intersects the detector

			// TPC AV
			Enu_CV_AV_TPC[pdg]              ->Fill(Enu, cv_weight);
			Enu_UW_AV_TPC[pdg]              ->Fill(Enu, dk2nu_weight);
			Enu_CV_AV_TPC_5MeV_bin[pdg]     ->Fill(Enu, cv_weight);
			Enu_UW_AV_TPC_5MeV_bin[pdg]     ->Fill(Enu, dk2nu_weight);
			Th_CV_AV_TPC[pdg]               ->Fill(theta, cv_weight);
			Th_UW_AV_TPC[pdg]               ->Fill(theta, dk2nu_weight);

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
			
			if (Pmom_dk == 0 && (mcflux.fptype == 211 || mcflux.fptype == -211) ){ // Pidar peak
				PiDAR_zpos[pdg]->Fill(mcflux.fvz, cv_weight);
			}
			if (Pmom_dk == 0 && (mcflux.fptype == 321 || mcflux.fptype == -321)){  // Kdar peak
				KDAR_zpos[pdg]->Fill(mcflux.fvz, cv_weight);
			}
			if (Pmom_dk == 0 && (mcflux.fptype == 13 || mcflux.fptype == -13)){    // Mudar peak
				MuDAR_zpos[pdg]->Fill(mcflux.fvz, cv_weight);
			}

			// Other plots
			parent_mom[pdg]  ->Fill(Pmom_tg,    cv_weight);
			parent_angle[pdg]->Fill(theta,      cv_weight);
			parent_zpos[pdg] ->Fill(mcflux.fvz, cv_weight);
			parent_zpos_angle[pdg]        ->Fill(mcflux.fvz, theta, cv_weight);
			parent_zpos_angle_energy[pdg]->Fill(mcflux.fvz, theta, Enu);
		
			if ( mcflux.fvz < 10000) flux_targ[pdg] ->Fill(Enu, cv_weight);
			if ( mcflux.fvz > 10000 &&  mcflux.fvz < 72000) flux_pipe[pdg] ->Fill(Enu, cv_weight);
			if ( mcflux.fvz > 72000) flux_dump[pdg] ->Fill(Enu, cv_weight);

			if ( theta < 20) flux_targ_theta[pdg] ->Fill(Enu, cv_weight);
			if ( theta > 20 &&  theta < 110) flux_pipe_theta[pdg] ->Fill(Enu, cv_weight);
			if ( theta > 110) flux_dump_theta[pdg] ->Fill(Enu, cv_weight);

			// 2D Stuff
			Enu_Th_CV_AV_TPC[pdg]->Fill(Enu, theta, cv_weight);
			Enu_Th_UW_AV_TPC[pdg]->Fill(Enu, theta, dk2nu_weight);


			// Now fill multisims
			for (unsigned l=0; l<labels.size(); l++) {

				// Universes
				for (unsigned i=0; i<Weights[l].size(); i++) {
					
					Enu_Syst_AV_TPC[pdg][l][i]->Fill(Enu, Weights[l][i]*cv_weight);
					Enu_Th_Syst_AV_TPC[pdg][l][i]->Fill(Enu, theta, Weights[l][i]*cv_weight);
				}
			}



		} // End loop over mctruth

	} // End loop over events


	// ++++++++++++++++++++++++++++++++
	// Plotting 
	// ++++++++++++++++++++++++++++++++

	TFile* output = new TFile("output.root", "RECREATE");

	POTTree->SetDirectory(output);

	TDirectory* savdir = gDirectory;

	std::cout << "flavour.size:\t" <<flav.size()<<std::endl;

	// Top Flav dir 
	std::vector<std::vector<TDirectory*>> subdir(flav.size());

	//Create label dirs
	for (unsigned i=0; i<flav.size(); i++){
		subdir[i].resize(parent.size()+5);
	
	}

	// Flavours
	for (unsigned int f=0; f<flav.size(); f++) {
	
		std::cout << "\n" <<flav[f] << std::endl;

		subdir[f][0] = savdir->mkdir(Form("%s",flav[f].c_str()));
		subdir[f][0]->cd();

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
		
		PiDAR_zpos[f]->Write();
		KDAR_zpos[f]->Write();
		MuDAR_zpos[f]->Write();
		parent_mom[f]->Write();
		parent_angle[f]->Write();
		parent_zpos[f]->Write();
		parent_zpos_angle[f]->Write();
		parent_zpos_angle_energy[f]->Write();
		flux_targ[f]->Write();
		flux_pipe[f]->Write();
		flux_dump[f]->Write();
		flux_targ_theta[f]->Write();
		flux_pipe_theta[f]->Write();
		flux_dump_theta[f]->Write();

		// Bolted this stuff in as an afterthought, could improve but it works ;)
		std::cout << "Window" << std::endl;
		subdir.at(f).at(parent.size()+2) = subdir[f][0]->mkdir("Window");
		subdir.at(f).at(parent.size()+2)->cd();
		Enu_CV_Window[f]->Write(); 
		Enu_UW_Window[f]->Write();     
		Enu_CV_Window_5MeV_bin[f]->Write();    
		Enu_UW_Window_5MeV_bin[f]->Write();    
		Enu_CV_Window_5MeV_bin_intersect[f]->Write();

		std::cout << "Detsmear" << std::endl;
		subdir.at(f).at(parent.size()+3) = subdir[f][0]->mkdir("Detsmear");
		subdir.at(f).at(parent.size()+3)->cd();
		Enu_CV_AV_TPC[f]->Write();
		Enu_UW_AV_TPC[f]->Write();
		Th_CV_AV_TPC[f]->Write();
		Th_UW_AV_TPC[f]->Write();
		Enu_CV_AV_TPC_5MeV_bin[f]->Write();
		Enu_UW_AV_TPC_5MeV_bin[f]->Write();
		Enu_Th_CV_AV_TPC[f]->Write();
		Enu_Th_UW_AV_TPC[f]->Write();

		std::cout << "Multisims" << std::endl;
		subdir.at(f).at(parent.size()+4) = subdir[f][0]->mkdir("Multisims");
		subdir.at(f).at(parent.size()+4)->cd();
		
		for (unsigned p=0; p < labels.size(); p++){
			
			for(int i = 0; i < 100; i++){
				Enu_Syst_AV_TPC.at(f).at(p).at(i)->Write();
				Enu_Th_Syst_AV_TPC.at(f).at(p).at(i)->Write();
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

	auto stop = high_resolution_clock::now();  // end time of script
	auto duration_sec = duration_cast<seconds>(stop - start); // time taken to run script
	auto duration_min = duration_cast<minutes>(stop - start); // time taken to run script
	std::cout << "Time taken by function: " << duration_sec.count() << " seconds" << std::endl; 
	std::cout << "Time taken by function: " << duration_min.count() << " minutes" << std::endl; 

	return 0;
}
