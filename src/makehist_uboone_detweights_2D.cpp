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

using namespace art;
using namespace std;

//___________________________________________________________________________
int calcEnuWgt(auto const& mcflux, const TVector3& xyz, double& wgt_xy){
	// Neutrino Energy and Weight at arbitrary point
	// Arguments:
	//    dk2nu    :: contains current decay information
	//    xyz      :: 3-vector of position to evaluate
	//                in *beam* frame coordinates  (cm units)
	//    enu      :: resulting energy
	//    wgt_xy   :: resulting weight
	// Return:
	//    (int)    :: error code
	// Assumptions:
	//    Energies given in GeV
	//    Particle codes have been translated from GEANT into PDG codes
	const double kPIMASS = 0.13957;
	const double kKMASS  = 0.49368;
	const double kK0MASS = 0.49767;
	const double kMUMASS = 0.105658389;
	const double kOMEGAMASS = 1.67245;

	const int kpdg_nue       =   12;  // extended Geant 53
	const int kpdg_nuebar    =  -12;  // extended Geant 52
	const int kpdg_numu      =   14;  // extended Geant 56
	const int kpdg_numubar   =  -14;  // extended Geant 55

	const int kpdg_muplus     =   -13;  // Geant  5
	const int kpdg_muminus    =    13;  // Geant  6
	const int kpdg_pionplus   =   211;  // Geant  8
	const int kpdg_pionminus  =  -211;  // Geant  9
	const int kpdg_k0long     =   130;  // Geant 10  ( K0=311, K0S=310 )
	const int kpdg_k0short    =   310;  // Geant 16
	const int kpdg_k0mix      =   311;  
	const int kpdg_kaonplus   =   321;  // Geant 11
	const int kpdg_kaonminus  =  -321;  // Geant 12
	const int kpdg_omegaminus =  3334;  // Geant 24
	const int kpdg_omegaplus  = -3334;  // Geant 32

	const double kRDET = 100.0;   // set to flux per 100 cm radius

	double xpos = xyz.X();
	double ypos = xyz.Y();
	double zpos = xyz.Z();

	double enu    = 0.0;  // don't know what the final value is
	wgt_xy = 0.0;  // but set these in case we return early due to error

	// in principle we should get these from the particle DB
	// but for consistency testing use the hardcoded values
	double parent_mass = kPIMASS;
	switch ( mcflux.fptype ) {
		case kpdg_pionplus:
		case kpdg_pionminus:
			parent_mass = kPIMASS;
			break;
		case kpdg_kaonplus:
		case kpdg_kaonminus:
			parent_mass = kKMASS;
			break;
		case kpdg_k0long:
		case kpdg_k0short:
		case kpdg_k0mix:
			parent_mass = kK0MASS;
			break;
		case kpdg_muplus:
		case kpdg_muminus:
			parent_mass = kMUMASS;
			break;
		case kpdg_omegaminus:
		case kpdg_omegaplus:
			parent_mass = kOMEGAMASS;
			break;
		default:
		std::cerr << "bsim::calcEnuWgt unknown particle type " << mcflux.fptype << std::endl << std::flush;
		assert(0);
		return 1;
	}

	double parentp2 = (mcflux.fpdpx*mcflux.fpdpx +
			  mcflux.fpdpy*mcflux.fpdpy +
			  mcflux.fpdpz*mcflux.fpdpz );
	double parent_energy = TMath::Sqrt( parentp2 +
					parent_mass*parent_mass);
	double parentp = TMath::Sqrt( parentp2 );

	double gamma     = parent_energy / parent_mass;
	double gamma_sqr = gamma * gamma;
	double beta_mag  = TMath::Sqrt( ( gamma_sqr - 1.0 )/gamma_sqr );

	// Get the neutrino energy in the parent decay CM
	double enuzr = mcflux.fnecm;
	// Get angle from parent line of flight to chosen point in beam frame
	double rad = TMath::Sqrt( (xpos-mcflux.fvx)*(xpos-mcflux.fvx) +
				(ypos-mcflux.fvy)*(ypos-mcflux.fvy) +
				(zpos-mcflux.fvz)*(zpos-mcflux.fvz) );

	double emrat = 1.0;
	double costh_pardet = -999.;
	// double theta_pardet = -999.;

	// boost correction, but only if parent hasn't stopped
	if ( parentp > 0. ) {
		
		costh_pardet = ( mcflux.fpdpx*(xpos-mcflux.fvx) +
				mcflux.fpdpy*(ypos-mcflux.fvy) +
				mcflux.fpdpz*(zpos-mcflux.fvz) ) 
				/ ( parentp * rad);
		
		if ( costh_pardet >  1.0 ) costh_pardet =  1.0;
		if ( costh_pardet < -1.0 ) costh_pardet = -1.0;
		// theta_pardet = TMath::ACos(costh_pardet);

		// Weighted neutrino energy in beam, approx, good for small theta
		emrat = 1.0 / ( gamma * ( 1.0 - beta_mag * costh_pardet ));
	}

	enu = emrat * enuzr;  // the energy ... normally

	// Get solid angle/4pi for detector element
	// small angle approximation, fixed by Alex Radovic
	//SAA//  double sangdet = ( kRDET*kRDET / 
	//SAA//                   ( (zpos-fDk2Nu->decay.vz)*(zpos-fDk2Nu->decay.vz) ) ) / 4.0;
	double sanddetcomp = TMath::Sqrt( ( (xpos-mcflux.fvx)*(xpos-mcflux.fvx) ) +
					( (ypos-mcflux.fvy)*(ypos-mcflux.fvy) ) +
					( (zpos-mcflux.fvz)*(zpos-mcflux.fvz) )   );
	double sangdet = (1.0-TMath::Cos(TMath::ATan( kRDET / sanddetcomp )))/2.0;

	// Weight for solid angle and lorentz boost
	wgt_xy = sangdet * ( emrat * emrat );  // ! the weight ... normally

	
	// Done for all except polarized muon decay
	// in which case need to modify weight 
	// (must be done in double precision)
	if ( mcflux.fptype == kpdg_muplus || mcflux.fptype == kpdg_muminus) {
		double beta[3], p_dcm_nu[4], p_nu[3], p_pcm_mp[3], partial;

		// Boost neu neutrino to mu decay CM
		beta[0] = mcflux.fpdpx / parent_energy;
		beta[1] = mcflux.fpdpy / parent_energy;
		beta[2] = mcflux.fpdpz / parent_energy;
		p_nu[0] = (xpos-mcflux.fvx)*enu/rad;
		p_nu[1] = (ypos-mcflux.fvy)*enu/rad;
		p_nu[2] = (zpos-mcflux.fvz)*enu/rad;
		
		partial = gamma * 
		(beta[0]*p_nu[0] + beta[1]*p_nu[1] + beta[2]*p_nu[2] );
		
		partial = enu - partial/(gamma+1.0);
		// the following calculation is numerically imprecise
		// especially p_dcm_nu[2] leads to taking the difference of numbers 
		//  of order ~10's and getting results of order ~0.02's
		// for g3numi we're starting with floats (ie. good to ~1 part in 10^7)
		p_dcm_nu[0] = p_nu[0] - beta[0]*gamma*partial;
		p_dcm_nu[1] = p_nu[1] - beta[1]*gamma*partial;
		p_dcm_nu[2] = p_nu[2] - beta[2]*gamma*partial;
		p_dcm_nu[3] = TMath::Sqrt( p_dcm_nu[0]*p_dcm_nu[0] +
					p_dcm_nu[1]*p_dcm_nu[1] +
					p_dcm_nu[2]*p_dcm_nu[2] );

		// Boost parent of mu to mu production CM
		double particle_energy = mcflux.fppenergy;
	
		gamma = particle_energy/parent_mass;
	
		beta[0] = mcflux.fppdxdz * mcflux.fpppz / particle_energy;
		beta[1] = mcflux.fppdydz * mcflux.fpppz / particle_energy;
		beta[2] =                    mcflux.fpppz / particle_energy;

		partial = gamma * ( beta[0]*mcflux.fmuparpx + 
				beta[1]*mcflux.fmuparpy + 
				beta[2]*mcflux.fmuparpz );
		partial = mcflux.fmupare - partial/(gamma+1.0);
	  
		p_pcm_mp[0] = mcflux.fmuparpx - beta[0]*gamma*partial;
		p_pcm_mp[1] = mcflux.fmuparpy - beta[1]*gamma*partial;
		p_pcm_mp[2] = mcflux.fmuparpz - beta[2]*gamma*partial;
		double p_pcm = TMath::Sqrt ( p_pcm_mp[0]*p_pcm_mp[0] +
					p_pcm_mp[1]*p_pcm_mp[1] +
					p_pcm_mp[2]*p_pcm_mp[2] );

		const double eps = 1.0e-30;  // ? what value to use
		
		if ( p_pcm < eps || p_dcm_nu[3] < eps ) {
			return 3; // mu missing parent info?
	  	
		}
		
		// Calc new decay angle w.r.t. (anti)spin direction
		double costh = ( p_dcm_nu[0]*p_pcm_mp[0] +
				p_dcm_nu[1]*p_pcm_mp[1] +
				p_dcm_nu[2]*p_pcm_mp[2] ) /
				( p_dcm_nu[3]*p_pcm );
	
		if ( costh >  1.0 ) costh =  1.0;
		if ( costh < -1.0 ) costh = -1.0;
		// Calc relative weight due to angle difference
		double wgt_ratio = 0.0;
		double xnu = 0.0;
		switch ( mcflux.fntype ) {
			case kpdg_nue:
			case kpdg_nuebar:
				wgt_ratio = 1.0 - costh;
				break;
			case kpdg_numu:
			case kpdg_numubar:
			{
				xnu = 2.0 * enuzr / kMUMASS;
				wgt_ratio = ( (3.0-2.0*xnu )  - (1.0-2.0*xnu)*costh ) / (3.0-2.0*xnu);
				if (xnu > 1) wgt_ratio = 0; // cut out unphysical numu decays
				break;
			}
			default:
			return 2; // bad neutrino type
		}

		wgt_xy = wgt_xy * wgt_ratio;
	
	} // ptype is muon

	return 0;
}
//___________________________________________________________________________
// Contains hard coded routines because it is a hacky test function
TVector3 RandomInDet() {

	// Randomly choose point in microboone
	double x = gRandom->Uniform(0.0    , 256.35);    //cm
	double y = gRandom->Uniform(-116.5 , 116.5); //cm
	double z = gRandom->Uniform(0.0    , 1036.8);    //cm

	return TVector3(x, y, z);
}
//___________________________________________________________________________
TVector3 FromDetToBeam( const TVector3& det ) {

    TVector3 beam;
    TRotation R;
	bool debug{false};

    //corrected rotation matrix using the 0,0,0 position for MicroBooNE
    //Previous matrix is calculated relative to MiniBooNE, which is not in the centre of the BNB!

    TVector3 newX(0.92103853804025682, 0.0000462540012621546684, -0.38947144863934974);
    TVector3 newY(0.0227135048039241207, 0.99829162468141475, 0.0538324139386641073);
    TVector3 newZ(0.38880857519374290, -0.0584279894529063024, 0.91946400794392302);
    //old matrix
    /*
    TVector3 newX(0.921228671,   0.00136256111, -0.389019125);
    TVector3 newY(0.0226872648,  0.998103714,    0.0572211871);
    TVector3 newZ(0.388359401,  -0.061539578,    0.919450845);
    */

    R.RotateAxes(newX,newY,newZ);
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
    // Updated position - leaving old positions here (July 2018)
    //TVector3 NuMIDet (54.499, 74.461,  677.611); // m
    TVector3 NuMIDet (55.02, 72.59,  672.70); //m
    NuMIDet *= 100.; // To have NuMIDet in cm

    beam = R * det + NuMIDet;

    return beam;
}
//___________________________________________________________________________
// Get the window normal for the tiltweight
double Get_tilt_wgt( const TVector3& detxyz, auto const& mcflux, auto const& mctruth){

	TVector3 xyzDk(mcflux.fvx,mcflux.fvy,mcflux.fvz);  // origin of decay
	
	TVector3 p3beam = mctruth.GetNeutrino().Nu().E() * ( detxyz - xyzDk ).Unit();

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

	bins[4] = {  20, 40, 120,  160 }; // theta -- removed edge theta bins where no events live and split into 3 bins for stats
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

			// Now get the momentums to calculate theta
			double costheta = mctruth.GetNeutrino().Nu().Pz() / mctruth.GetNeutrino().Nu().P();
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
			TVector3 xyz_det = RandomInDet();
			
			// From detector to beam coordinates
			TVector3 xyz_beam = FromDetToBeam(xyz_det);

			// Get the new weight at the detector
			calcEnuWgt(mcflux, xyz_beam, detwgt);

			// Get the tiltweight
			tiltwght = Get_tilt_wgt(xyz_beam, mcflux, mctruth);

			// Weight of neutrino parent (importance weight) * Neutrino weight for a decay forced at center of near detector 
			cv_weight *= mcflux.fnimpwt * detwgt * tiltwght;

			// if (cv_weight < 0) std::cout << "Still got a negative weight:\t" << cv_weight << " pdg" << pdg << " Enu:\t" << mctruth.GetNeutrino().Nu().E()<< std::endl;
			if (cv_weight < 0) cv_weight = 0; // get rid of them pesky negative weights
			if ( std::isnan(cv_weight) == 1) { // catch NaN values
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

			

			// ++++++++++++++++++++++++++++++++
			// Now got weights, fill histograms
			// ++++++++++++++++++++++++++++++++

			// TPC AV
			Enu_Th_CV_AV_TPC[pdg]->Fill(mctruth.GetNeutrino().Nu().E(),theta, cv_weight);
			Enu_Th_UW_AV_TPC[pdg]->Fill(mctruth.GetNeutrino().Nu().E(),theta, mcflux.fnimpwt * detwgt * tiltwght);

			// Now fill multisims
			// Options        
			for (unsigned l=0; l<labels.size(); l++) {

				// Universes
				for (unsigned i=0; i<Weights[l].size(); i++) {
					Enu_Th_Syst_AV_TPC[pdg][l][i]->Fill(mctruth.GetNeutrino().Nu().E(),theta, Weights[l][i]*cv_weight);

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

