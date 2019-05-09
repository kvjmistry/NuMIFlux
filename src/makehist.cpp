/**
 * Script will take the output files of that have been ppfx weighted and fill the flux histograms
 * This script contains a unified approach for code for each of the methods to calculate the flux  
 *
 * Authors: J. Zennamo, A. Mastbaum, K Mistry
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
#include "functions_makehist.h"

using namespace art;
using namespace std;
//___________________________________________________________________________
// Contains hard coded routines because it is a hacky test function
TVector3 RandomInDet() {

	// Randomly choose point in microboone
	double x = gRandom->Uniform(0.0    , 256.35);    //cm
	double y = gRandom->Uniform(-116.5 , 116.5);     //cm
	double z = gRandom->Uniform(0.0    , 1036.8);    //cm

	return TVector3(x, y, z);
}
//___________________________________________________________________________
TVector3 FromDetToBeam( const TVector3 det, bool rotate_only ) {

    TVector3 beam;
    TRotation R;
	bool debug{false};

    // Rotation matrix using the 0,0,0 position for MicroBooNE (det to beam input)
    TVector3 newX(0.92103853804025682, 0.0000462540012621546684, -0.38947144863934974);
    TVector3 newY(0.0227135048039241207, 0.99829162468141475, 0.0538324139386641073);
    TVector3 newZ(0.38880857519374290, -0.0584279894529063024, 0.91946400794392302);

    R.RotateAxes(newX,newY,newZ); // Also inverts to beam to det
    if (debug) {
        cout << "R_{beam to det} = " << endl;
        cout << " [ " << R.XX() << " " << R.XY() << " " << R.XZ() << " ] " << endl;
        cout << " [ " << R.YX() << " " << R.YY() << " " << R.YZ() << " ] " << endl;
        cout << " [ " << R.ZX() << " " << R.ZY() << " " << R.ZZ() << " ] " << endl;
        cout << endl;
    }
    R.Invert(); // R is now the inverse
    if (debug) {
        cout << "R_{det to beam} = " << endl;
        cout << " [ " << R.XX() << " " << R.XY() << " " << R.XZ() << " ] " << endl;
        cout << " [ " << R.YX() << " " << R.YY() << " " << R.YZ() << " ] " << endl;
        cout << " [ " << R.ZX() << " " << R.ZY() << " " << R.ZZ() << " ] " << endl;
        cout << endl;
    }
    
	// Now R allows to go from detector to beam coordinates.
    // NuMIDet is vector from NuMI target to uB detector (in beam coordinates)
    TVector3 NuMIDet (55.02, 72.59,  672.70); //m
    NuMIDet *= 100.; // To have NuMIDet in cm

	if (rotate_only) beam = R * det; // Only rotate the vector
    else beam = R * det + NuMIDet;

    return beam;
}
//___________________________________________________________________________
// Get the window normal for the tiltweight
double Get_tilt_wgt( const TVector3& detxyz, auto const& mcflux, double enu){

	TVector3 xyzDk(mcflux.fvx,mcflux.fvy,mcflux.fvz);  // origin of decay
	
	TVector3 p3beam = enu * ( detxyz - xyzDk ).Unit();

	// Hardcoded for testing, but in priciple would want
	TVector3 windownorm = { 0.528218812, 0.8308577046, 0.175101003}; 

	double tiltweight =  p3beam.Unit().Dot( windownorm );

	return tiltweight;

}
//___________________________________________________________________________
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

	InputTag mctruths_tag { "flux" };
	InputTag  evtwght_tag { "eventweight" };

	double totalPOT{0};
	double Kaontotal{0};


	vector<string> badfiles;
	vector<string> filename;
	
	// Load in the input arguments
	for (int i=1; i<argc; i++) { 
		
		
		
		
		std::cout << "FILE : " << argv[i] << std::endl; 
		TFile *filein=new TFile(argv[i]);
		
		if (filein->IsZombie()) {
			std::cout << "ERROR: File is ZOMBIE:\t" << argv[i] << std::endl;
			badfiles.push_back(string(argv[i]));
			filein->Close();
		}
		else {
			filename.push_back(string(argv[i]));
			totalPOT+=500000; // 100 000 POT per dk2nu file 500 000 if running over nova files
			filein->Close();
		}
	}

	std::cout << "\nTotal POT read in:\t" << totalPOT << std::endl;

	std::cout << "\nUsing 5e5 POT per dk2nu file, the flux will be wrong if this is not the case!\n" << std::endl;

	// Histograms for each flavor
	std::vector<TH1D*> Enu_CV_AV_TPC;
	std::vector<TH1D*> Enu_UW_AV_TPC;
	
	std::vector<TH1D*> Th_CV_AV_TPC;
	std::vector<TH1D*> Th_UW_AV_TPC;

	// flavors - systematic - universe 
	std::vector<std::vector<std::vector<TH1D*> > > Enu_Syst_AV_TPC;
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

	Weights.resize(labels.size());

	for (unsigned int i=0; i<labels.size(); i++) {
		Weights[i].resize(100);
	}

	Enu_CV_AV_TPC.resize(4);
	Enu_UW_AV_TPC.resize(4);
	Enu_Syst_AV_TPC.resize(4);
	
	Th_CV_AV_TPC.resize(4);
	Th_UW_AV_TPC.resize(4);
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
		// Enu_CV_AV_TPC[i] = new TH1D(Form("%s_CV_AV_TPC",flav[i].c_str()),"",n, bin);
		// Enu_UW_AV_TPC[i] = new TH1D(Form("%s_unweighted_AV_TPC",flav[i].c_str()),"",n, bin);

		
		// Use these for a finer binning
		Enu_CV_AV_TPC[i] = new TH1D(Form("%s_CV_AV_TPC",flav[i].c_str()),"",4000, 0, 20);
		Enu_UW_AV_TPC[i] = new TH1D(Form("%s_unweighted_AV_TPC",flav[i].c_str()),"",4000, 0, 20);

		Th_CV_AV_TPC[i] = new TH1D(Form("Th_%s_CV_AV_TPC",flav[i].c_str()),"",n_th, bin_th);
		Th_UW_AV_TPC[i] = new TH1D(Form("Th_%s_unweighted_AV_TPC",flav[i].c_str()),"",n_th, bin_th);

		Enu_Syst_AV_TPC[i].resize(labels.size());
		Th_Syst_AV_TPC[i].resize(labels.size());

		// Labels
		for(unsigned j=0; j<labels.size(); j++) {
			Enu_Syst_AV_TPC[i][j].resize(100);
			Th_Syst_AV_TPC[i][j].resize(100);

			// Universes
			for(int k=0; k<100; k++){			
				Enu_Syst_AV_TPC[i][j][k] =  new TH1D(Form("%s_%s_Uni_%d_AV_TPC",flav[i].c_str(), labels[j].c_str(), k),"",n, bin);
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


			double cv_weight = 1; 
			double detwgt; // New weight at a window value
			double tiltwght;
			double dk2nu_cv = 1;

			double enu = mctruth.GetNeutrino().Nu().E();
			
			// Now get the momentum to calculate theta
			TVector3 mom_det = {mctruth.GetNeutrino().Nu().Px(),mctruth.GetNeutrino().Nu().Py(),mctruth.GetNeutrino().Nu().Pz()};
			TVector3 mom_beam = FromDetToBeam(mom_det, true);

			double costheta_beam = mom_beam.Z() / std::sqrt(mom_beam.X()*mom_beam.X() + mom_beam.Y()*mom_beam.Y() + mom_beam.Z()*mom_beam.Z());
			double theta_beam = std::acos(costheta_beam) * 180 / 3.14159265;

			double costheta = mctruth.GetNeutrino().Nu().Pz() / mctruth.GetNeutrino().Nu().P();
			double theta = std::acos(costheta) * 180 / 3.14159265;

			std::cout << "theta_current:\t" << theta << "\ttheta_new:\t" << theta_beam << "\tdecay pos:\t"<< mcflux.fvz/100.0 << std::endl;
			std::cout << "mctruth px:\t" << mctruth.GetNeutrino().Nu().Px() << "\tnew px:\t" << mom_beam.X() << std::endl;
			
			// Count the total Kaons
			if (mcflux.fptype == 321 || mcflux.fptype == -321) Kaontotal+=mcflux.fnimpwt;
			
			
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
			TVector3 xyz_det = RandomInDet();
			
			// From detector to beam coordinates
			TVector3 xyz_beam = FromDetToBeam(xyz_det,false);

			// Get the new weight at the detector
			calcEnuWgt(mcflux, xyz_beam, enu,  detwgt);

			// Get the tiltweight
			tiltwght = Get_tilt_wgt(xyz_beam, mcflux, enu);

			// Weight of neutrino parent (importance weight) * Neutrino weight for a decay forced at center of near detector 
			cv_weight *= mcflux.fnimpwt * detwgt * tiltwght; 
			dk2nu_cv   = mcflux.fnimpwt * detwgt * tiltwght; 

			// Error handling
			if (cv_weight < 0) cv_weight = 0; // get rid of them pesky negative weights
			if (std::isnan(cv_weight) == 1) { // catch NaN values
				std::cout << "got a nan:\t"<<cv_weight <<std::endl;
				cv_weight = 0;
			}

			if (dk2nu_cv < 0) dk2nu_cv = 0; // get rid of them pesky negative weights
			if (std::isnan(dk2nu_cv) == 1) { // catch NaN values
				std::cout << "got a nan:\t"<<dk2nu_cv <<std::endl;
				dk2nu_cv = 0;
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

			

			// ++++++++++++++++++++++++++++++++
			// Now got weights, fill histograms
			// ++++++++++++++++++++++++++++++++

			// TPC AV
			Enu_CV_AV_TPC[pdg]->Fill(enu, cv_weight);
			Enu_UW_AV_TPC[pdg]->Fill(enu, dk2nu_cv);
			Th_CV_AV_TPC[pdg]->Fill(theta, cv_weight);
			Th_UW_AV_TPC[pdg]->Fill(theta, dk2nu_cv);

			// Now fill multisims
			
			// Options
			for (unsigned l=0; l<labels.size(); l++) {

				// Universes
				for (unsigned i=0; i<Weights[l].size(); i++) {
				
					Enu_Syst_AV_TPC[pdg][l][i]->Fill(enu, Weights[l][i]*cv_weight);
					Th_Syst_AV_TPC[pdg][l][i]->Fill(theta, Weights[l][i]*cv_weight);
					
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

	std::vector<string> cont = {"Active_TPC_Volume", ""};

	// Flavours
	for (unsigned int f=0; f<flav.size(); f++) {
	
		std::cout << "\n" <<flav[f] << std::endl;

		subdir[f][0][0] = savdir->mkdir(Form("%s",flav[f].c_str()));
		subdir[f][0][0]->cd();

		// Write CV fluxes
		Enu_CV_AV_TPC[f]->Write();
		Enu_UW_AV_TPC[f]->Write();
		Th_CV_AV_TPC[f]->Write();    
		Th_UW_AV_TPC[f]->Write();

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
