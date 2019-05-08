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
	double y = gRandom->Uniform(-116.5 , 116.5);     //cm
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

	geoalgo::GeoAlgo const _geo_algo_instance;

	geoalgo::AABox volAVTPC( _xRange.first, _yRange.first, _zRange.first, _xRange.second, _yRange.second, _zRange.second);

	InputTag mctruths_tag { "flux" };
	InputTag  evtwght_tag { "eventweight" };

	double totalPOT{0};
	double Kaontotal{0};


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

	std::cout << "\nUsing 5e5 POT per dk2nu file, the flux will be wrong if this is not the case!\n" << std::endl;

	// Histograms for each flavor
	std::vector<TH1D*> Enu_CV_Window;
	std::vector<TH1D*> Enu_CV_AV_TPC;
	std::vector<TH1D*> Enu_UW_Window;
	std::vector<TH1D*> Enu_UW_AV_TPC;

	std::vector<TH1D*> Th_CV_AV_TPC;
	std::vector<TH1D*> Th_UW_AV_TPC;

	// 5Mev Bins
	std::vector<TH1D*> Enu_CV_Window_rebin;
	std::vector<TH1D*> Enu_CV_AV_TPC_rebin;
	std::vector<TH1D*> Enu_UW_Window_rebin;
	std::vector<TH1D*> Enu_UW_AV_TPC_rebin;

	// Detector intersection window method
	std::vector<TH1D*> Enu_CV_Window_rebin_intersect;

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
	bins[0] = { // numu
		0.00, 0.025, 0.03, 0.235 ,0.24, 0.50, 1.00, 1.25, 1.50, 1.75, 2.00, 2.25, 2.50, 3.00, 4.00, 5.00, 6.00, 7.00, 10.00 };

	bins[1] = {  // nue
		0.00 ,0.06, 0.125, 0.25, 0.5, 0.75, 1.00, 1.25, 1.50, 1.75, 2.00, 2.50, 3.00, 3.50, 4.00, 5.00 };

	bins[2] = {// numubar
		0.00, 0.025, 0.03, 0.235 ,0.24, 0.50, 1.00, 1.25, 1.50, 1.75, 2.00, 2.25, 2.50, 3.00, 4.00, 5.00, 6.00, 7.00, 10.00 };

	bins[3] = {  // nuebar
		0.00 ,0.06, 0.125,  0.25, 0.5, 0.75, 1.00, 1.25, 1.50, 1.75, 2.00, 2.50, 3.00, 3.50, 4.00, 5.00 };

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
	
	Enu_CV_Window_rebin.resize(4);
	Enu_CV_AV_TPC_rebin.resize(4);
	Enu_UW_Window_rebin.resize(4);
	Enu_UW_AV_TPC_rebin.resize(4);

	Enu_CV_Window_rebin_intersect.resize(4);


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
		Enu_CV_Window_rebin[i] = new TH1D(Form("%s_CV_Window_rebin",flav[i].c_str()),"",4000, 0, 20);
		Enu_CV_AV_TPC_rebin[i] = new TH1D(Form("%s_CV_AV_TPC_rebin",flav[i].c_str()),"",4000, 0, 20);
		Enu_UW_Window_rebin[i] = new TH1D(Form("%s_UW_Window_rebin",flav[i].c_str()),"",4000, 0, 20);
		Enu_UW_AV_TPC_rebin[i] = new TH1D(Form("%s_UW_AV_TPC_rebin",flav[i].c_str()),"",4000, 0, 20);

		Enu_CV_Window_rebin_intersect[i] = new TH1D(Form("%s_CV_Window_rebin_intersect",flav[i].c_str()),"",4000, 0, 20);

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
			double costheta = mctruth.GetNeutrino().Nu().Pz() / mctruth.GetNeutrino().Nu().P();
			double theta = std::acos(costheta) * 180 / 3.14159265;
			double Enu = mctruth.GetNeutrino().Nu().E();
			double Pmom_dk = std::sqrt( mcflux.fpdpz*mcflux.fpdpz + mcflux.fpdpy*mcflux.fpdpy + mcflux.fpdpx*mcflux.fpdpx ); // Parent momentum at decay
			double Pmom_tg = std::sqrt(mcflux.ftpx*mcflux.ftpx + mcflux.ftpy*mcflux.ftpy + mcflux.ftpz*mcflux.ftpz); // Parent moment

			// Count the total Kaons
			if (mcflux.fptype == 321 || mcflux.fptype == -321) Kaontotal+=mcflux.fnimpwt;

			// Get the ppfx cv_weight
			for (auto last : evtwght.fWeight) {

				if (last.first.find("PPFXCV") != std::string::npos) {

					// if(last.second.at(0) > 30 || last.second.at(0) < 0){ // still fill even if bad weight, changed from >90 to >30
					if(last.second.at(0) < 0){ // change this to only throw out negative weights
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
			TVector3 xyz_det = RandomInDet();
			
			// From detector to beam coordinates
			TVector3 xyz_beam = FromDetToBeam(xyz_det);

			// Get the new weight at the detector
			calcEnuWgt(mcflux, xyz_beam, detwgt);

			// Get the tiltweight
			tiltwght = Get_tilt_wgt(xyz_beam, mcflux, mctruth);

			// Weight of neutrino parent (importance weight) * Neutrino weight for a decay forced at center of near detector 
			cv_weight        *= mcflux.fnimpwt * detwgt * tiltwght; // for ppfx cases
			dk2nu_weight     *= mcflux.fnimpwt * detwgt * tiltwght; // for UW cases 
			
			window_weight       *= mcflux.fnimpwt * mcflux.fnwtfar; // mcflux.fnwtfar == mcflux.fnwtnear
			dk2nu_window_weight *= mcflux.fnimpwt * mcflux.fnwtfar; // mcflux.fnwtfar == mcflux.fnwtnear
			
			// Error handling
			if (cv_weight < 0) cv_weight = 0; // get rid of them pesky negative weights
			if (std::isnan(cv_weight) == 1) { // catch NaN values
				std::cout << "got a nan:\t"<<cv_weight <<std::endl;
				cv_weight = 0;
			}

			if (window_weight < 0) window_weight = 0; // get rid of them pesky negative weights
			if (std::isnan(window_weight) == 1) { // catch NaN values
				std::cout << "got a nan:\t"<<window_weight <<std::endl;
				window_weight = 0;
				
			}

			if (dk2nu_weight < 0) dk2nu_weight = 0; // get rid of them pesky negative weights
			if (std::isnan(dk2nu_weight) == 1) { // catch NaN values
				std::cout << "got a nan:\t"<<dk2nu_weight <<std::endl;
				dk2nu_weight = 0;
				
			}

			if (dk2nu_window_weight < 0) dk2nu_window_weight = 0; // get rid of them pesky negative weights
			if (std::isnan(dk2nu_window_weight) == 1) { // catch NaN values
				std::cout << "got a nan:\t"<<dk2nu_window_weight <<std::endl;
				dk2nu_window_weight = 0;
				
			}
			
			// ++++++++++++++++++++++++++++++++
			// Now got weights, fill histograms
			// ++++++++++++++++++++++++++++++++

			// Window
			Enu_CV_Window[pdg]      ->Fill(Enu, window_weight);
			Enu_UW_Window[pdg]      ->Fill(Enu, dk2nu_window_weight);
			Enu_CV_Window_rebin[pdg]->Fill(Enu, window_weight);
			Enu_UW_Window_rebin[pdg]->Fill(Enu, dk2nu_window_weight);

			if (intercept) Enu_CV_Window_rebin_intersect[pdg] ->Fill(Enu, window_weight); // Flux that intersects the detector

			// TPC AV
			Enu_CV_AV_TPC[pdg]      ->Fill(Enu, cv_weight);
			Enu_UW_AV_TPC[pdg]      ->Fill(Enu, dk2nu_weight);
			Enu_CV_AV_TPC_rebin[pdg]->Fill(Enu, cv_weight);
			Enu_UW_AV_TPC_rebin[pdg]->Fill(Enu, dk2nu_weight);
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

	// Spit out the total Kaon
	std::cout << "Total Kaons:\t" << Kaontotal << std::endl;

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

		Enu_CV_Window_rebin[f]->Write();      
		Enu_CV_AV_TPC_rebin[f]->Write();
		Enu_UW_Window_rebin[f]->Write();      
		Enu_UW_AV_TPC_rebin[f]->Write();
		
		Enu_CV_Window_rebin_intersect[f]->Write();

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
