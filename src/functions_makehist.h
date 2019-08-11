#ifndef FUNCTIONS_MAKEHIST_H
#define FUNCTIONS_MAKEHIST_H

// Functions file for make hist script
// Contains the predetermined detector geometries too

// Std Includes
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <stdlib.h>
#include <string>
#include <vector>
#include <chrono> 

// ROOT Includes
#include "TDirectory.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TH2D.h"
#include "TFile.h"
#include "TRandom.h"
#include "TVector3.h"
#include "TRotation.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "canvas/Utilities/InputTag.h"

// LArsoft Includes
#include "geo/GeoVector.h"
#include "geo/GeoAABox.h"
#include "geo/GeoHalfLine.h"
#include "geo/GeoAlgo.h"
#include "gallery/Event.h"
#include "gallery/ValidHandle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "larsim/EventWeight/Base/MCEventWeight.h"
#include "dk2nu/tree/dk2nu.h"
#include "dk2nu/tree/dkmeta.h"
#include "dk2nu/tree/calcLocationWeights.h"

double totalPOT{0};
bool input_flag{false}; // flag to see if a detector has been specified

std::vector<std::string> badfiles;
std::vector<std::string> filename;

// systematic - universe 
std::vector< std::vector< double > > Weights;   

std::vector<std::string> flav = { "numu", "nue", "numubar", "nuebar" };

std::vector<std::string> labels;

std::vector<double> temp, temp2; // For the bins

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
std::vector<std::string> parent = {"PI_Plus", "PI_Minus", "Mu_Plus", "Mu_Minus", "Kaon_Plus", "Kaon_Minus" , "K0L"};
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

// 2D histograms
std::vector<TH2D*> Enu_Th_CV_AV_TPC;
std::vector<TH2D*> Enu_Th_UW_AV_TPC;


// Weighted Histograms
std::vector<std::vector<std::vector<TH1D*>>> Enu_Syst_AV_TPC;     //1D
std::vector<std::vector<std::vector<TH2D*>>> Enu_Th_Syst_AV_TPC;  //2D

//___________________________________________________________________________
// Container class for a detector
class Detector {
	public:
		Detector()=default;
		
		// Constructor
		Detector(std::string detector_name_,
			std::pair<float, float>  xRange_,
			std::pair<float, float>  yRange_,
			std::pair<float, float>  zRange_,
			TVector3 Trans_Targ2Det_beam_,
			TVector3 Trans_Targ2Det_det_,
			TVector3 Rot_row_x_,
			TVector3 Rot_row_y_,
			TVector3 Rot_row_z_, 
			TVector3 Win_Base_,
			TVector3 Win_pt1_,
			TVector3 Win_pt2_,
			std::vector<std::vector<double>> bins_  ) {
				detector_name = detector_name_;
				xRange = xRange_;
				yRange = yRange_;
				zRange = zRange_;
				Trans_Targ2Det_beam = Trans_Targ2Det_beam_;
				Trans_Targ2Det_det  = Trans_Targ2Det_det_;
				Rot_row_x = Rot_row_x_;
				Rot_row_y = Rot_row_y_;
				Rot_row_z = Rot_row_z_;
				Win_Base = Win_Base_;
				Win_pt1 = Win_pt1_;
				Win_pt2 = Win_pt2_;
				bins.resize(5);
				bins = bins_;
			};
		
		std::string detector_name;

		// Fiducial Volume Definition
		std::pair<float, float>  xRange; 
		std::pair<float, float>  yRange;
		std::pair<float, float>  zRange;

		// Translation vector from target to detector in beam coords [cm]
		TVector3 Trans_Targ2Det_beam;

		// Translation vector from target to detector in det coords [cm]
		TVector3 Trans_Targ2Det_det;

		// bool to decide whether to use the translation in beam or detector coords
		// bool useBeam = true; 

		// Rotation matrix from beam to detector
		TVector3 Rot_row_x; // Row x of rotation matrix
		TVector3 Rot_row_y; // Row y of rotation matrix
		TVector3 Rot_row_z; // Row z of rotation matrix

		// Window defintion in detector coordinates [cm]
		TVector3 Win_Base;
		TVector3 Win_pt1;
		TVector3 Win_pt2;

		// Histogram Bins | 1 for each flavour + theta
		std::vector<std::vector<double>> bins;

};
//___________________________________________________________________________
// Create detector definition
void Initialise(std::string detector_type, Detector &Detector_){
	std::cout << "Initialising detector for:\t" << detector_type << std::endl;

	// Fiducial Volume Definition
	std::pair<float, float>  xRange; 
	std::pair<float, float>  yRange;
	std::pair<float, float>  zRange;

	// Translation vector from target to detector in beam coords [cm]
	TVector3 Trans_Targ2Det_beam;

	// Translation vector from target to detector in det coords [cm]
	TVector3 Trans_Targ2Det_det;

	// Rotation matrix from beam to detector
	TVector3 Rot_row_x; // Row x of rotation matrix
	TVector3 Rot_row_y; // Row y of rotation matrix
	TVector3 Rot_row_z; // Row z of rotation matrix

	// Window defintion in detector coordinates [cm]
	TVector3 Win_Base;
	TVector3 Win_pt1;
	TVector3 Win_pt2;

	// Histogram Bins | 1 for each flavour + theta
	std::vector<std::vector<double>> bins(5);

	// MicroBooNE is the Input
	if (detector_type == "uboone"){

		xRange.first  =     -0; // cm
		xRange.second = 256.35;
		yRange.first  = -116.5;
		yRange.second =  116.5;
		zRange.first  =      0;
		zRange.second = 1036.8;

		Trans_Targ2Det_det = { -31387.58422, -3316.402543, -60100.2414}; //cm -- detector coords
		Trans_Targ2Det_beam = { 5502, 7259, 67270}; //cm -- in beam coords

		// Rotation matrix using the 0,0,0 position for MicroBooNE (beam to det input)
		Rot_row_x = { 0.92103853804025681562, 0.022713504803924120662, 0.38880857519374290021  };
		Rot_row_y = { 4.6254001262154668408e-05, 0.99829162468141474651, -0.058427989452906302359 };
		Rot_row_z = { -0.38947144863934973769, 0.053832413938664107345, 0.91946400794392302291   };

		Win_Base = { 500, -500, -3500 };
		Win_pt1  = {-500,  200, -3500 };
		Win_pt2  = { 500, -500,  2000 };

		bins[0] = { // numu
			0.00, 0.025, 0.03, 0.235 ,0.24, 0.50, 1.00, 1.25, 1.50, 1.75, 2.00, 2.25, 2.50, 3.00, 4.00, 5.00, 6.00, 7.00, 10.00 };

		bins[1] = {  // nue
			0.00 ,0.06, 0.125, 0.25, 0.5, 0.75, 1.00, 1.25, 1.50, 1.75, 2.00, 2.50, 3.00, 3.50, 4.00, 5.00 };

		bins[2] = {// numubar
			0.00, 0.025, 0.03, 0.235 ,0.24, 0.50, 1.00, 1.25, 1.50, 1.75, 2.00, 2.25, 2.50, 3.00, 4.00, 5.00, 6.00, 7.00, 10.00 };

		bins[3] = {  // nuebar
			0.00 ,0.06, 0.125,  0.25, 0.5, 0.75, 1.00, 1.25, 1.50, 1.75, 2.00, 2.50, 3.00, 3.50, 4.00, 5.00 };
		
		// Theta
		bins[4] = {  0, 20, 110,  160 }; // theta -- removed edge theta bins where no events live and split into 3 bins for stats


	}
	else if (detector_type == "nova"){

		xRange.first  = -176; // cm
		xRange.second =  177;
		yRange.first  = -172;
		yRange.second =  179;
		zRange.first  =   25;
		zRange.second = 1150;

		Trans_Targ2Det_det = {226.9447, 6100.1882, -99113.1313}; //cm -- detector coords
		// Trans_Targ2Det_beam = {1171.74545 ,       -331.51325 ,      99293.47347}; // beam coords
		Trans_Targ2Det_beam = {1150.170113 ,      -280.0752339 ,    100099.1001}; // beam coords -- with updated attempt
		// Trans_Targ2Det_beam = { 1150.170113 ,     -280.0752339 ,      100099.1001}; // new test with beam coords from genie page


		// Rotation matrix using the 0,0,0 position for NOvA (beam to det input)
		Rot_row_x = {9.9990e-01, -8.2300e-04, -1.4088e-02 };
		Rot_row_y = {3.0533e-06, 9.9831e-01,  -5.8103e-02 };
		Rot_row_z = {1.4112e-02, 5.8097e-02,  9.9821e-01  };

		Win_Base  = { 500, -250,  -500 };
		Win_pt1   = { 500,  500,  -500 };
		Win_pt2   = { -500, -250, -500 };

		// Window normal for nova (for debug comparisons)
		// Window Normal:	 X 0.01411202435 Y 0.05809720554 Z 0.9982111828
		// Window area is 74

		// numu
		bins[0] = { 0.0, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0};
		
		// nue
		bins[1] = {  0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0};
		
		// numubar
		bins[2] = {  0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0};
		
		// nuebar
		bins[3] = {  0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0};
		
		// Theta
		bins[4] = {0.0, 180.0};

	}
	else{
		std::cout << "Unknown detector type given, EXITING....." << std::endl;
		exit(1);
	}
	std::cout << "\nFiducial Volume:\n"
		"X:\t(" << xRange.first << ", "<< xRange.second << ") cm" << "\n" <<
		"Y:\t(" << yRange.first << ", "<< yRange.second << ") cm" << "\n" <<
		"Z:\t(" << zRange.first << ", "<< zRange.second << ") cm\n" << "\n"<<
		"R_{beam to det} = \n" << 
		" [ " << Rot_row_x.X() << " " << Rot_row_x.Y() << " " << Rot_row_x.Z() << " ] " << "\n" <<
		" [ " << Rot_row_y.X() << " " << Rot_row_y.Y() << " " << Rot_row_y.Z() << " ] " << "\n" <<
		" [ " << Rot_row_z.X() << " " << Rot_row_z.Y() << " " << Rot_row_z.Z() << " ] " << "\n\n" <<
		"Translation from target to detector in beam coords [cm] = \n" <<
		" [ " << Trans_Targ2Det_beam.X() << ", " << Trans_Targ2Det_beam.Y() << ", " << Trans_Targ2Det_beam.Z() << " ] " << "\n\n" <<
		"Translation from target to detector in detector coords [cm] = \n" <<
		" [ " << Trans_Targ2Det_det.X() << ", " << Trans_Targ2Det_det.Y() << ", " << Trans_Targ2Det_det.Z() << " ] " << "\n\n" <<
		"Window (in det coords) [cm] = \n" << 
		" [ " << Win_Base.X() << " " << Win_Base.Y() << " " << Win_Base.Z() << " ] " << "\n" <<
		" [ " << Win_pt1.X()  << " " << Win_pt1.Y()  << " " << Win_pt1.Z()  << " ] " << "\n" <<
		" [ " << Win_pt2.X()  << " " << Win_pt2.Y()  << " " << Win_pt2.Z()  << " ] " << "\n" <<
		std::endl; 

	Detector_ = Detector(detector_type, xRange, yRange, zRange, Trans_Targ2Det_beam, Trans_Targ2Det_det, Rot_row_x, Rot_row_y, Rot_row_z, Win_Base, Win_pt1, Win_pt2, bins );

}
//___________________________________________________________________________
int calcEnuWgt(auto const& mcflux, const TVector3& xyz, double& enu, double& wgt_xy, double kRDET) {
	// Neutrino Energy and Weight at arbitrary point
	// Based on:
	// NuMI-NOTE-BEAM-0109 (MINOS DocDB # 109)
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
	// for now ... these masses _should_ come from TDatabasePDG
	// but use these hard-coded values to "exactly" reproduce old code
	//
	// old mass values are v01_07_* and before
	// new mass values (v01_08_* and after) are Geant4 v4.10.3 values
	//
#ifdef HISTORIC_MASS
	const double kPIMASS      = 0.13957;
	const double kKMASS       = 0.49368;
	const double kK0MASS      = 0.49767;
	const double kMUMASS      = 0.105658389;
	const double kOMEGAMASS   = 1.67245;
#else
	const double kPIMASS      = 0.1395701;     // 0.13957;
	const double kKMASS       = 0.493677;      // 0.49368;
	const double kK0MASS      = 0.497614;      // 0.49767;
	const double kMUMASS      = 0.1056583715;  // 0.105658389;
	const double kOMEGAMASS   = 1.67245;       // 1.67245;
#endif
	
	// from CLHEP/Units/PhysicalConstants.h
	// used by Geant as CLHEP::neutron_mass_c2
	const double kNEUTRONMASS = 0.93956536;
	const int kpdg_nue       =   12;  // extended Geant 53
	const int kpdg_nuebar    =  -12;  // extended Geant 52
	const int kpdg_numu      =   14;  // extended Geant 56
	const int kpdg_numubar   =  -14;  // extended Geant 55
	const int kpdg_muplus      =   -13;  // Geant  5
	const int kpdg_muminus     =    13;  // Geant  6
	const int kpdg_pionplus    =   211;  // Geant  8
	const int kpdg_pionminus   =  -211;  // Geant  9
	const int kpdg_k0long      =   130;  // Geant 10  ( K0=311, K0S=310 )
	const int kpdg_k0short     =   310;  // Geant 16
	const int kpdg_k0mix       =   311;
	const int kpdg_kaonplus    =   321;  // Geant 11
	const int kpdg_kaonminus   =  -321;  // Geant 12
	const int kpdg_omegaminus  =  3334;  // Geant 24
	const int kpdg_omegaplus   = -3334;  // Geant 32
	const int kpdg_neutron     =  2112;
	const int kpdg_antineutron = -2112;
	
	// const double kRDET = 100.0;   // set to flux per 100 cm radius
	
	double xpos = xyz.X();
	double ypos = xyz.Y();
	double zpos = xyz.Z();
	
	enu    = 0.0;  // don't know what the final value is
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
		case kpdg_neutron:
		case kpdg_antineutron:
			parent_mass = kNEUTRONMASS;
			break;
		default:
			std::cerr << "bsim::calcEnuWgt unknown particle type " << mcflux.fptype << std::endl << std::flush;
			enu    = 0.0;
			wgt_xy = 0.0;
			return 1;
	}

	double parentp2 = ( mcflux.fpdpx*mcflux.fpdpx +
						mcflux.fpdpy*mcflux.fpdpy +
						mcflux.fpdpz*mcflux.fpdpz );
	double parent_energy = TMath::Sqrt( parentp2 + parent_mass * parent_mass);
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
	//SAA//                   ( (zpos-mcflux.fvz)*(zpos-mcflux.fvz) ) ) / 4.0;
	double sanddetcomp = TMath::Sqrt( ( (xpos-mcflux.fvx)*(xpos-mcflux.fvx) ) +
									  ( (ypos-mcflux.fvy)*(ypos-mcflux.fvy) ) +
									  ( (zpos-mcflux.fvz)*(zpos-mcflux.fvz) ) );
	double sangdet = (1.0-TMath::Cos(TMath::ATan( kRDET / sanddetcomp )))/2.0;
	
	// Weight for solid angle and lorentz boost
	wgt_xy = sangdet * ( emrat * emrat );  // ! the weight ... normally
	
	// Done for all except polarized muon decay in which case need to modify weight
	// (must be done in double precision)
	
	// BUT do this only for case of muon decay, not muon capture
	// until beamline simulation code gets updated these generally show up as
	// mcflux.fndecay == 0, but certainly not dkp_mup_nusep or dkp_mum_nusep
	
	// so was:
	// if ( mcflux.fptype  == kpdg_muplus      ||
	//      mcflux.fptype  == kpdg_muminus        ) {
	
	// now:
	if ( mcflux.fndecay == bsim::dkp_mup_nusep || mcflux.fndecay == bsim::dkp_mum_nusep ) {
		double beta[3], p_dcm_nu[4], p_nu[3], p_pcm_mp[3], partial;
		
		// Boost neu neutrino to mu decay CM
		beta[0] = mcflux.fpdpx / parent_energy;
		beta[1] = mcflux.fpdpy / parent_energy;
		beta[2] = mcflux.fpdpz / parent_energy;
		p_nu[0] = (xpos-mcflux.fvx)*enu/rad;
		p_nu[1] = (ypos-mcflux.fvy)*enu/rad;
		p_nu[2] = (zpos-mcflux.fvz)*enu/rad;
		
		partial = gamma * (beta[0]*p_nu[0] + beta[1]*p_nu[1] + beta[2]*p_nu[2] );
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
		
		// protect against small excursions
		if ( costh >  1.0 ) costh =  1.0;
		if ( costh < -1.0 ) costh = -1.0;
		
		// Calc relative weight due to angle difference
		double wgt_ratio = 0.0;
		switch ( mcflux.fntype ) {
			// Nue/nuebar
			case kpdg_nue:
			case kpdg_nuebar:
				wgt_ratio = 1.0 - costh;
				break;
			
			// Numu/numubar
			case kpdg_numu:
			case kpdg_numubar:
			{
				double xnu = 2.0 * enuzr / kMUMASS;
				wgt_ratio = ( (3.0-2.0*xnu )  - (1.0-2.0*xnu)*costh ) / (3.0-2.0*xnu);
				
				if ( wgt_ratio < 0.0 ) {
					std::cerr << "bsim::calcEnuWgt encountered serious problem: "
							<< " wgt_ratio " << wgt_ratio
							<< " enu " << enu << " costh " << costh << " xnu " << xnu
							<< " enuzr=mcflux.fnecm " << enuzr << " kMUMASS " << kMUMASS
							<< " norig " << mcflux.fnorig
							<< " ndecay " << mcflux.fndecay
							<< " ntype " << mcflux.fntype
							<< " ptype " << mcflux.fptype
							<< std::endl;
					enu    = 0;
					wgt_xy = 0;
					return 4; // bad, bad, bad calculation
				}
				break;
			}
			
			default:
				enu    = 0.0;
				wgt_xy = 0.0;
				return 2; // bad neutrino type for muon decay
		}
		
		wgt_xy = wgt_xy * wgt_ratio;
	
	} // ptype is muon
	
	return 0;
}
//___________________________________________________________________________
// Picks a random point in the detector
TVector3 RandomInDet(Detector Detector_) {

	// Randomly choose point in microboone
	double x = gRandom->Uniform(Detector_.xRange.first, Detector_.xRange.second); //cm
	double y = gRandom->Uniform(Detector_.yRange.first, Detector_.yRange.second); //cm
	double z = gRandom->Uniform(Detector_.zRange.first, Detector_.zRange.second); //cm

	return TVector3(x, y, z);
}
//___________________________________________________________________________
TVector3 FromDetToBeam( const TVector3 det, bool rotate_only, Detector Detector_ ) {

	TVector3 beam;
	TRotation R;
	bool debug{false};
	
	R.RotateAxes(Detector_.Rot_row_x, Detector_.Rot_row_y, Detector_.Rot_row_z); // Also inverts so now to det to beam
	R.Invert(); // go back to beam to det
	if (debug) {
		std::cout << "R_{beam to det} = " << std::endl;
		std::cout << " [ " << R.XX() << " " << R.XY() << " " << R.XZ() << " ] " << std::endl;
		std::cout << " [ " << R.YX() << " " << R.YY() << " " << R.YZ() << " ] " << std::endl;
		std::cout << " [ " << R.ZX() << " " << R.ZY() << " " << R.ZZ() << " ] " << std::endl;
		std::cout << std::endl;
	}
	R.Invert(); // R is now the inverse
	if (debug) {
		std::cout << "R_{det to beam} = " << std::endl;
		std::cout << " [ " << R.XX() << " " << R.XY() << " " << R.XZ() << " ] " << std::endl;
		std::cout << " [ " << R.YX() << " " << R.YY() << " " << R.YZ() << " ] " << std::endl;
		std::cout << " [ " << R.ZX() << " " << R.ZY() << " " << R.ZZ() << " ] " << std::endl;
		std::cout << std::endl;
	}

	if (rotate_only) beam = R * det;                         // Only rotate the vector
	// else beam = R * (det - Detector_.Trans_Targ2Det_det); // for when transaltion is given in det coords
	else beam = R * det + Detector_.Trans_Targ2Det_beam;     // for when translation is given in beam coords

	return beam;
}
//___________________________________________________________________________
// Get the window normal for the tiltweight --check
double Get_tilt_wgt( const TVector3& detxyz, auto const& mcflux, double enu, Detector Detector_){

	TVector3 xyzDk(mcflux.fvx, mcflux.fvy, mcflux.fvz);  // origin of decay in beam coords
	
	TVector3 p3beam = enu  * ( detxyz - xyzDk ).Unit(); // Momentum in beam coords

	// Convert from user to beam coord and from 3 points to base + 2 directions
	TVector3 fWin_Base_beam = FromDetToBeam( Detector_.Win_Base, false,  Detector_ );
	TVector3 fWin_pt1_beam  = FromDetToBeam( Detector_.Win_pt1,  false,  Detector_ );
	TVector3 fWin_pt2_beam  = FromDetToBeam( Detector_.Win_pt2,  false,  Detector_ );

	// Define direction vectors
	TVector3 fFluxWindowDir1 = fWin_pt1_beam - fWin_Base_beam;
	TVector3 fFluxWindowDir2 = fWin_pt2_beam - fWin_Base_beam;

	// Get the window normal
	TVector3 fWindowNormal = fFluxWindowDir1.Cross(fFluxWindowDir2).Unit();

	// Window normal debug
	// std::cout << "window normal:" <<std::endl;
	// std::cout << " [ " << fWindowNormal.X() << " " << fWindowNormal.Y() << " " << fWindowNormal.Z() << " ] " << std::endl;

	double tiltweight =  p3beam.Unit().Dot( fWindowNormal );

	return tiltweight;

}
//___________________________________________________________________________
// This function will recaclulate the missed nu rays for the intersection method
// The intension is that this will fix the normalistion problems
// Returns the weight/pi
double Recalc_Intersection_wgt(geoalgo::GeoAlgo const _geo_algo_instance, geoalgo::AABox volAVTPC, auto const& mcflux, auto const& mctruth, Detector Detector_, double KRDET, double &enu ){

	TRotation R;
	TVector3 x3beam, x3beam_det;
	int retries{0};
	double weight = 0;
	// bool debug = true;
	bool debug = true;
	
	R.RotateAxes(Detector_.Rot_row_x, Detector_.Rot_row_y, Detector_.Rot_row_z); // R is now det to beam
	TRotation R_Beam_2_Det = R.Invert();

	// Convert from user to beam coord and from 3 points to base + 2 directions
	TVector3 fWin_Base_beam = FromDetToBeam( Detector_.Win_Base, false,  Detector_ );
	TVector3 fWin_pt1_beam  = FromDetToBeam( Detector_.Win_pt1,  false,  Detector_ );
	TVector3 fWin_pt2_beam  = FromDetToBeam( Detector_.Win_pt2,  false,  Detector_ );

	// Define direction vectors
	TVector3 fFluxWindowDir1 = fWin_pt1_beam - fWin_Base_beam;
	TVector3 fFluxWindowDir2 = fWin_pt2_beam - fWin_Base_beam;

	// Now see if this neutrino vector is going to intersect with the detector
	// If it doesn't then recalculate
	while(true){
		retries++;
		
		// Extra precautions to make sure we are definately overwriting -- can remove
		x3beam.SetXYZ(0,0,0);
		x3beam_det.SetXYZ(0,0,0);
		weight=0;

		double random = ((double) rand() / (RAND_MAX)); // Get a random number 0 to 1

		if (debug)std::cout <<"r uniform:\t" << random << std::endl;

		// Pick a new point on the window in beam coords
		x3beam = fWin_Base_beam + (random * fFluxWindowDir1) + (random * fFluxWindowDir2);

		// Get the weight
		calcEnuWgt(mcflux, x3beam, enu, weight, KRDET);

		// Get the nu ray direction in beam coords
		TVector3 xyzDk(mcflux.fvx,mcflux.fvy,mcflux.fvz);  // Origin of decay in beam coords
		TVector3 NuRay_Dir =  enu * (x3beam - xyzDk).Unit();
		
		// Convert to detector coordinates
		NuRay_Dir = R_Beam_2_Det * NuRay_Dir;

		// Now convert xbeam to detector coordinates too
		x3beam_det = R_Beam_2_Det * (x3beam - Detector_.Trans_Targ2Det_beam);

		// debug, mom *1k to improve readability
		if (debug) std::cout << "retry no:\t" << retries << std::endl;
		if (debug) std::cout << "px:\t" <<mctruth.GetNeutrino().Nu().Px()*1000 << "  " << "  " << NuRay_Dir.X()*1000 << std::endl; 
		if (debug) std::cout << "py:\t" <<mctruth.GetNeutrino().Nu().Py()*1000 << "  " << "  " << NuRay_Dir.Y()*1000 << std::endl;
		if (debug) std::cout << "pz:\t" <<mctruth.GetNeutrino().Nu().Pz()*1000 << "  " << "  " << NuRay_Dir.Z()*1000 <<  std::endl;
		if (debug) std::cout << "vx:\t" <<mctruth.GetNeutrino().Nu().Vx() << "  " << "  " << x3beam_det.X() << std::endl;
		if (debug) std::cout << "vy:\t" <<mctruth.GetNeutrino().Nu().Vy() << "  " << "  " << x3beam_det.Y() << std::endl;
		if (debug) std::cout << "vz:\t" <<mctruth.GetNeutrino().Nu().Vz() << "  " << "  " << x3beam_det.Z() << std::endl;
		if (debug) std::cout << "Enu:\t" << enu << std::endl;
		if (debug) std::cout << "Weight:\t" << weight << std::endl;

		// Make the neutrino ray
		geoalgo::HalfLine ray(x3beam_det.X(), // point on window in detector coordinates units have to be m (maybe divide by 100??)
					x3beam_det.Y(),
					x3beam_det.Z(),
					NuRay_Dir.X(),  // px in detector coordinates
					NuRay_Dir.Y(),  // py
					NuRay_Dir.Z()); // pz

		// Count nu intersections with tpc
		auto vec = _geo_algo_instance.Intersection(volAVTPC, ray); 
		bool intercept = false;

		if (debug) std::cout << "vec size:\t" << vec.size()<< std::endl;

		if (vec.size() == 0) { intercept = false; } // no intersections
		if (vec.size() == 2) { intercept = true; }  // 2 intersections
		if (vec.size() != 2 && vec.size() != 0) {   // other intersection
			// std::cout << "Neutrino ray has " << vec.size()
			// 	<< " intersection with the detector volume"
			// 	<< std::endl;
		}
		
		// Got an interception, so break!
		if (intercept) {
			if (debug) std::cout << "Passed with " << retries << " retries"<< "\n" << std::endl;
			break; 
		}
		if (retries > 200) {
			// If there are more than 1000 attempts and still no intersection then give up!
			// uboone has cases where they never intersect -- could be due to window not big enough to catch them? For now supress these
			if (Detector_.detector_name == "nova") std::cout << "Recalculation failed due to > 200 tries and no interception" << std::endl;
			return 0.0; // throw the event away
		} 

	}
	
	double tiltweight = Get_tilt_wgt( x3beam, mcflux, enu, Detector_);

	return weight * tiltweight ;
}
//___________________________________________________________________________
// Function to check the status of weights, if they are bad, set to zero
void check_weight(double &weight){
	
	if (weight < 0) weight = 0; // get rid of them pesky negative weights
	
	if (std::isnan(weight) == 1) { // catch NaN values
		std::cout << "got a nan:\t" << weight << std::endl;
		weight = 0;
	}
}
//___________________________________________________________________________
// Function to check interceptions with the detector
bool check_intercept(geoalgo::GeoAlgo const _geo_algo_instance, geoalgo::AABox volAVTPC, auto const& mctruth){

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
		// std::cout << "Neutrino ray has " << vec.size()
		// 	<< " intersection with the detector volume"
		// 	<< std::endl;
	}
	return intercept;
}
//___________________________________________________________________________
#endif