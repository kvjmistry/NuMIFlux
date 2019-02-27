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
			totalPOT+=100000*50; // 50* 100 000 POT per dk2nu file
			filein->Close();
		}
	}

	std::cout << "\nTotal POT read in:\t" << totalPOT << std::endl;

	std::cout << "\nUsing 5e6 POT per dk2nu file, the flux will be wrong if this is not the case!\n" << std::endl;

	// Histograms for each flavor
	std::vector<TH1D*> Enu_CV_Window;
	std::vector<TH1D*> Enu_CV_AV_TPC;
	std::vector<TH1D*> Enu_UW_Window;
	std::vector<TH1D*> Enu_UW_AV_TPC;

	// Flux by Parent
	std::vector<string> parent = {"PI_Plus", "PI_Minus", "Mu_Plus", "Mu_Minus", "KaonPM", "K0L"};
	std::vector<std::vector<TH1D*> > Enu_Parent_AV_TPC;		// Flux by Parent
	std::vector<std::vector<TH1D*> > Th_Parent_AV_TPC; 		// Flux by parent in theta
	std::vector<std::vector<TH1D*> > zpos_Parent_AV_TPC; 	// Flux by parent in z Pos at decay (production is unimplemented for now)
	std::vector<std::vector<TH1D*> > dk_Parent_mom; 		// Decay momentum by parent
	std::vector<std::vector<TH1D*> > impwght_Parent; 		// Importance weight by parent
	std::vector<std::vector<TH1D*> > Prod_energy_Parent; 	// Production energy by parent
	std::vector<std::vector<TH1D*> > Targ_mom_Parent; 	    // Momentum by parent as it leaves the target
	
	// Other hists
	TH1D* NuMu_PiDAR_zpos = new TH1D("NuMu_PiDAR_zpos","", 400 , 0, 80000); // numu Pidar peak zpos
	TH1D* NuMu_KDAR_zpos =  new  TH1D("NuMu_KDAR_zpos","", 400 , 0, 80000);  // numu kdar peak zpos
	TH1D* NuMu_peak_mom_muon =   new  TH1D("NuMu_peak_mom_muon","", 100 , 0, 25);  // muon parent momentum at 550m peak

	// Tree for POT counting
	TTree* POTTree = new TTree("POT","Total POT");
	POTTree -> Branch("POT", &totalPOT);
	POTTree -> Fill();

	// systematic - universe 
	std::vector< std::vector< double > > Weights;   

	std::vector<string> flav = { "numu", "nue", "numubar", "nuebar" };

	std::vector< std::vector<double> > bins; bins.resize(4);
	bins[0] = {
		0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50,
		0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15, 
		1.20, 1.25, 1.30, 1.35, 1.40, 1.45, 1.50, 1.75, 1.90, 2.10, 2.50,
		3.00, 3.50, 4.00, 4.50, 5.00, 5.50, 6.00, 6.50, 7.00, 7.50, 8.00,
		9.00, 10.0, 11.0, 12.0, 13.0, 14.0};

	bins[1] = {  
		0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50,
		0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15, 
		1.20, 1.25, 1.30, 1.35, 1.40, 1.45, 1.50, 1.75, 1.90, 2.10, 2.50,
		3.00, 3.50, 4.00, 4.50, 5.00, 6.00, 7.00, 8.00, 9.00, 10.0};

	bins[2] = {
		0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50,
		0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15, 
		1.20, 1.25, 1.30, 1.35, 1.40, 1.45, 1.50, 1.75, 1.90, 2.10, 2.50,
		3.00, 3.50, 4.00, 4.50, 5.00, 5.50, 6.00, 6.50, 7.00, 7.50, 8.00,
		9.00, 10.0, 11.0, 12.0, 13.0, 14.0};

	bins[3] = {  
		0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50,
		0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15, 
		1.20, 1.25, 1.30, 1.35, 1.40, 1.45, 1.50, 1.75, 1.90, 2.10, 2.50,
		3.00, 3.50, 4.00, 4.50, 5.00, 6.00, 7.00, 8.00, 9.00, 10.0};

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

	Enu_CV_Window.resize(4);
	Enu_CV_AV_TPC.resize(4);
	Enu_UW_Window.resize(4);
	Enu_UW_AV_TPC.resize(4);
	Enu_Parent_AV_TPC.resize(4);
	Th_Parent_AV_TPC.resize(4);
	zpos_Parent_AV_TPC.resize(4);
	impwght_Parent.resize(4);
	Targ_mom_Parent.resize(4);

	
	
	std::vector<double> temp;


	// Flavors
	for(unsigned i=0; i<flav.size(); i++) {
		int const n = bins[i].size()-1;
		temp.clear();
		temp = bins[i];

		double* bin = &temp[0];

		// FLux histograms
		Enu_CV_Window[i] = new TH1D(Form("%s_CV_Window",flav[i].c_str()),"",n, bin);
		// Enu_CV_AV_TPC[i] = new TH1D(Form("%s_CV_AV_TPC",flav[i].c_str()),"",n, bin);

		// new binning schmeme to be the same as marcos
		// Enu_CV_Window[i] = new TH1D(Form("%s_CV_Window",flav[i].c_str()),"",4000, 0, 20);
		Enu_CV_AV_TPC[i] = new TH1D(Form("%s_CV_AV_TPC",flav[i].c_str()),"",4000, 0, 20);

		Enu_UW_Window[i] = new TH1D(Form("%s_unweighted_Window",flav[i].c_str()),"",n, bin);
		Enu_UW_AV_TPC[i] = new TH1D(Form("%s_unweighted_AV_TPC",flav[i].c_str()),"",n, bin);


		Enu_Parent_AV_TPC[i].resize(parent.size());
		Th_Parent_AV_TPC[i].resize(parent.size());
		zpos_Parent_AV_TPC[i].resize(parent.size());
		impwght_Parent[i].resize(parent.size());
		Targ_mom_Parent[i].resize(parent.size());
		

		// Parent
		for(unsigned k = 0; k < parent.size(); k++){
			Enu_Parent_AV_TPC[i][k]  = new TH1D(Form("Enu_%s_%s_AV_TPC",flav[i].c_str(), parent[k].c_str()),"", 4000, 0, 20);
			Th_Parent_AV_TPC[i][k]   = new TH1D(Form("Th_%s_%s_AV_TPC",flav[i].c_str(), parent[k].c_str()),"", 40 , 0, 180);
			zpos_Parent_AV_TPC[i][k] = new TH1D(Form("zpos_%s_%s_AV_TPC",flav[i].c_str(), parent[k].c_str()),"", 400 , 0, 80000);
			impwght_Parent[i][k]     = new TH1D(Form("impwght_Parent_%s_%s",flav[i].c_str(), parent[k].c_str()),"", 1 , 0, -1);
			Targ_mom_Parent[i][k]    = new TH1D(Form("Targ_mom_Parent_%s_%s",flav[i].c_str(), parent[k].c_str()),";Parent P at Target [GeV/c]", 4000, 0, 20);
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
		
		// Loop over MCTruths
		for (size_t i=0; i<mctruths.size(); i++) {
			auto const& mctruth = mctruths.at(i);
			auto const& mcflux = mcfluxs.at(i);

			int pdg;
			if     (mctruth.GetNeutrino().Nu().PdgCode() == 14) pdg = 0; // numu
			else if(mctruth.GetNeutrino().Nu().PdgCode() == 12) pdg = 1; // nue
			else if(mctruth.GetNeutrino().Nu().PdgCode() ==-14) pdg = 2; // numubar
			else if(mctruth.GetNeutrino().Nu().PdgCode() ==-12) pdg = 3; // nuebar
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
			
			// Now get the momentums to calculate theta
			double costheta = mctruth.GetNeutrino().Nu().Pz() / mctruth.GetNeutrino().Nu().P();
			double theta = std::acos(costheta) * 180 / 3.14159265;
			double Enu = mctruth.GetNeutrino().Nu().E();
			double Pmom_dk = std::sqrt( mcflux.fpdpz*mcflux.fpdpz + mcflux.fpdpy*mcflux.fpdpy + mcflux.fpdpx*mcflux.fpdpx ); // Parent momentum at decay
			double Pmom_tg = std::sqrt(mcflux.ftpx*mcflux.ftpx + mcflux.ftpy*mcflux.ftpy + mcflux.ftpz*mcflux.ftpz); // Parent moment

			// Weight of neutrino parent (importance weight) * Neutrino weight for a decay forced at center of near detector 
			cv_weight *= mcflux.fnimpwt * mcflux.fnwtfar; // mcflux.fnwtfar == mcflux.fnwtnear
			
			if (cv_weight < 0) cv_weight = 0; // get rid of them pesky negative weights
			if (std::isnan(cv_weight) == 1) { // catch NaN values
				std::cout << "got a nan:\t"<<cv_weight <<std::endl;
				cv_weight = 0;
				
			}
			
			if (Enu > 0.074 && Enu < 0.082 ){
				// std::cout << "E:\t" << Enu << "Parent:\t" << mcflux.fptype <<"  theta:\t" <<theta<<"   ntype:\t"<< mcflux.fndecay<<   std::endl;
			}

			// ++++++++++++++++++++++++++++++++
			// Now got weights, fill histograms
			// ++++++++++++++++++++++++++++++++

			// Window
			Enu_CV_Window[pdg]->Fill(Enu, cv_weight);
			Enu_UW_Window[pdg]->Fill(Enu, mcflux.fnimpwt * mcflux.fnwtfar);

			// TPC AV
			if (intercept) {
				Enu_CV_AV_TPC[pdg]->Fill(Enu, cv_weight);
				Enu_UW_AV_TPC[pdg]->Fill(Enu, mcflux.fnimpwt * mcflux.fnwtfar);

				// INDEXING: 0: PI_Plus 1: PI_Minus 2: Mu Plus 3: Mu_Minus 4: Kaon+/- 5: K0L 
				if (mcflux.fptype == 211){ // pi plus
					Enu_Parent_AV_TPC[pdg][0]->Fill(Enu, cv_weight);
					Th_Parent_AV_TPC[pdg][0]->Fill(theta, cv_weight);
					zpos_Parent_AV_TPC[pdg][0]->Fill(mcflux.fvz, cv_weight);
					impwght_Parent[pdg][0]->Fill(mcflux.fnimpwt);
					Targ_mom_Parent[pdg][0]->Fill(Pmom_tg, cv_weight);

				}
				else if (mcflux.fptype == -211){ // pi minus
					Enu_Parent_AV_TPC[pdg][1]->Fill(Enu, cv_weight);
					Th_Parent_AV_TPC[pdg][1]->Fill(theta, cv_weight);
					zpos_Parent_AV_TPC[pdg][1]->Fill(mcflux.fvz, cv_weight);
					impwght_Parent[pdg][1]->Fill(mcflux.fnimpwt);
					Targ_mom_Parent[pdg][1]->Fill(Pmom_tg, cv_weight);


				}
				else if (mcflux.fptype == -13){ // mu plus
					Enu_Parent_AV_TPC[pdg][2]->Fill(Enu, cv_weight);
					Th_Parent_AV_TPC[pdg][2]->Fill(theta, cv_weight);
					zpos_Parent_AV_TPC[pdg][2]->Fill(mcflux.fvz, cv_weight);
					impwght_Parent[pdg][2]->Fill(mcflux.fnimpwt);
					Targ_mom_Parent[pdg][2]->Fill(Pmom_tg, cv_weight);


				} 
				else if (mcflux.fptype == 13){ // mu minus
					Enu_Parent_AV_TPC[pdg][3]->Fill(Enu, cv_weight);
					Th_Parent_AV_TPC[pdg][3]->Fill(theta, cv_weight);
					zpos_Parent_AV_TPC[pdg][3]->Fill(mcflux.fvz, cv_weight);
					impwght_Parent[pdg][3]->Fill(mcflux.fnimpwt);
					Targ_mom_Parent[pdg][3]->Fill(Pmom_tg, cv_weight);


				} 
				else if (mcflux.fptype == 321 || mcflux.fptype == -321){ // K+/-
					Enu_Parent_AV_TPC[pdg][4]->Fill(Enu, cv_weight);
					Th_Parent_AV_TPC[pdg][4]->Fill(theta, cv_weight);
					zpos_Parent_AV_TPC[pdg][4]->Fill(mcflux.fvz, cv_weight);
					impwght_Parent[pdg][4]->Fill(mcflux.fnimpwt);
					Targ_mom_Parent[pdg][4]->Fill(Pmom_tg, cv_weight);

				}
				else if (mcflux.fptype == 130 ){ // K0L
					Enu_Parent_AV_TPC[pdg][5]->Fill(Enu, cv_weight);
					Th_Parent_AV_TPC[pdg][5]->Fill(theta, cv_weight);
					zpos_Parent_AV_TPC[pdg][5]->Fill(mcflux.fvz, cv_weight);
					impwght_Parent[pdg][5]->Fill(mcflux.fnimpwt);
					Targ_mom_Parent[pdg][5]->Fill(Pmom_tg, cv_weight);

				}
				
				if (Enu > 0.025 && Enu < 0.03 && mcflux.fptype == 211 ){ // Pidar peak
					NuMu_PiDAR_zpos->Fill(mcflux.fvz, cv_weight);
				}
				if (Enu > 0.235 && Enu < 0.24 && (mcflux.fptype == 321 || mcflux.fptype == -321)){ // kdar peak
					NuMu_KDAR_zpos->Fill(mcflux.fvz, cv_weight);
				}

				// Beam Pipe Peak at 550 m
				if (mcflux.fvz> 54000 && mcflux.fvz < 56000 && mcflux.fptype == 13 ){ 
					NuMu_peak_mom_muon->Fill(Pmom_dk,cv_weight );
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
	std::vector<std::vector<TDirectory*>> subdir(flav.size()); 
	
	//Create label dirs
	for (unsigned i=0; i<flav.size(); i++){
		subdir[i].resize(parent.size()+1);
		
		// Create Win/AVTPC dirs
		// for(unsigned j=0; j<labels.size()+1; j++) {
		// 	subdir[i][j].resize(3);
		// }
	
	}

	// std::vector<string> cont = { "Window", "Active_TPC_Volume", ""};

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

		// Parent
		for(unsigned int k = 1; k < parent.size()+1; k++){
			std::cout << parent[k-1] << std::endl;
			subdir[f][k] = subdir[f][0]->mkdir(Form("%s",parent[k-1].c_str()));
			subdir[f][k]->cd();

			Enu_Parent_AV_TPC[f][k-1]->Write();
			Th_Parent_AV_TPC[f][k-1]->Write();
			zpos_Parent_AV_TPC[f][k-1]->Write();
			impwght_Parent[f][k-1]->Write();
			Targ_mom_Parent[f][k-1]->Write();
		}

		// Make other plots folder for miscalanious variables
		std::cout << "OtherPlots" << std::endl;
		subdir[f][parent.size()+2] = subdir[f][0]->mkdir("OtherPlots");
		subdir[f][parent.size()+2]->cd();
		NuMu_PiDAR_zpos->Write();
		NuMu_KDAR_zpos->Write();
		NuMu_peak_mom_muon->Write();

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