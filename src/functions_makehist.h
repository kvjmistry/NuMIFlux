// Functions file for make hist script
// Intention is to unify the makehist scripts into one

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
#include "TVector3.h"
#include "TRotation.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "geo/GeoVector.h"
#include "geo/GeoAABox.h"
#include "geo/GeoHalfLine.h"
#include "geo/GeoAlgo.h"
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
			TVector3 Trans_Det2Beam_,
			TVector3 Rot_row_x_,
			TVector3 Rot_row_y_,
			TVector3 Rot_row_z_, 
			TVector3 Win_Base_,
			TVector3 Win_pt1_,
			TVector3 Win_pt2_  ) {
				detector_name = detector_name_;
				xRange = xRange_;
				yRange = yRange_;
				zRange = zRange_;
				Trans_Det2Beam = Trans_Det2Beam_;
				Rot_row_x = Rot_row_x_;
				Rot_row_y = Rot_row_y_;
				Rot_row_z = Rot_row_z_;
				Win_Base = Win_Base_;
				Win_pt1 = Win_pt1_;
				Win_pt2 = Win_pt2_;
			};
		
		std::string detector_name;

		// Fiducial Volume Definition
		std::pair<float, float>  xRange; 
		std::pair<float, float>  yRange;
		std::pair<float, float>  zRange;

		// Translation vector from beam origin to detector origin [cm]
		TVector3 Trans_Det2Beam;

		// Rotation matrix from beam to detector
		TVector3 Rot_row_x; // Row x of rotation matrix
		TVector3 Rot_row_y; // Row y of rotation matrix
		TVector3 Rot_row_z; // Row z of rotation matrix

		// Window defintion in detector coordinates [cm]
		TVector3 Win_Base;
		TVector3 Win_pt1;
		TVector3 Win_pt2;

};
//___________________________________________________________________________
// Create detector definition
void Initialise(std::string detector_type, Detector &Detector_){
	std::cout << "Initialising detector for:\t" << detector_type << std::endl;

	// Fiducial Volume Definition
	std::pair<float, float>  xRange; 
	std::pair<float, float>  yRange;
	std::pair<float, float>  zRange;

	// Translation vector from det origin to beam origin in det coords [cm]
	TVector3 Trans_Det2Beam;

	// Rotation matrix from beam to detector
	TVector3 Rot_row_x; // Row x of rotation matrix
	TVector3 Rot_row_y; // Row y of rotation matrix
	TVector3 Rot_row_z; // Row z of rotation matrix

	// Window defintion in detector coordinates [cm]
	TVector3 Win_Base;
	TVector3 Win_pt1;
	TVector3 Win_pt2;

	// MicroBooNE is the Input
	if (detector_type == "uboone"){

		xRange.first  =     -0; // cm
		xRange.second = 256.35;
		yRange.first  = -116.5;
		yRange.second =  116.5;
		zRange.first  =      0;
		zRange.second = 1036.8;

		// Trans_Det2Beam = { -31387.58422, -3316.402543, -60100.2414}; //cm
		Trans_Det2Beam = { 5502, 7259, 67270}; //cm in beam coords

		// Rotation matrix using the 0,0,0 position for MicroBooNE (beam to det input)
		Rot_row_x = { 0.92103853804025681562, 0.022713504803924120662, 0.38880857519374290021  };
		Rot_row_y = { 4.6254001262154668408e-05, 0.99829162468141474651, -0.058427989452906302359 };
		Rot_row_z = { -0.38947144863934973769, 0.053832413938664107345, 0.91946400794392302291   };

		Win_Base = { 500, -500, -3500 };
		Win_pt1  = {-500,  200, -3500 };
		Win_pt2  = { 500, -500,  2000 };

	}
	else if (detector_type == "nova"){

		xRange.first  = -176; // cm
		xRange.second =  177;
		yRange.first  = -172;
		yRange.second =  179;
		zRange.first  =   25;
		zRange.second = 1150;

		// Trans_Det2Beam = {1226.9447, 6100.1882, -99113.1313}; //cm
		 Trans_Det2Beam = {1171.74545 ,       -331.51325 ,      99293.47347}; // new test with beam coords

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

	}
	else{
		std::cout << "Unknown detector type given, EXITING....." << std::endl;
		exit(1);
	}
	std::cout << "\nFiducial Volume:\n"
		"X:\t(" << xRange.first << ", "<< xRange.second << ")" << "\n" <<
		"Y:\t(" << yRange.first << ", "<< yRange.second << ")" << "\n" <<
		"Z:\t(" << zRange.first << ", "<< zRange.second << ")\n" << "\n"<<
		"R_{beam to det} = \n" << 
		" [ " << Rot_row_x.X() << " " << Rot_row_x.Y() << " " << Rot_row_x.Z() << " ] " << "\n" <<
		" [ " << Rot_row_y.X() << " " << Rot_row_y.Y() << " " << Rot_row_y.Z() << " ] " << "\n" <<
		" [ " << Rot_row_z.X() << " " << Rot_row_z.Y() << " " << Rot_row_z.Z() << " ] " << "\n\n" <<
		"Translation from target to detector in detector coords = \n" <<
		" [ " << Trans_Det2Beam.X() << ", " << Trans_Det2Beam.Y() << ", " << Trans_Det2Beam.Z() << " ] " << "\n\n" <<
		"Window (in det coords) = \n" << 
		" [ " << Win_Base.X() << " " << Win_Base.Y() << " " << Win_Base.Z() << " ] " << "\n" <<
		" [ " << Win_pt1.X()  << " " << Win_pt1.Y()  << " " << Win_pt1.Z()  << " ] " << "\n" <<
		" [ " << Win_pt2.X()  << " " << Win_pt2.Y()  << " " << Win_pt2.Z()  << " ] " << "\n" <<
		std::endl; 

	Detector_ = Detector(detector_type, xRange, yRange, zRange, Trans_Det2Beam, Rot_row_x, Rot_row_y, Rot_row_z, Win_Base, Win_pt1, Win_pt2 );

}
//___________________________________________________________________________
int calcEnuWgt(auto const& mcflux, const TVector3& xyz, double &enu, double& wgt_xy){
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

	if (rotate_only) beam = R * det;              // Only rotate the vector
	// else beam = R * (det - Detector_.Trans_Det2Beam);
	else beam = R * det + Detector_.Trans_Det2Beam; // new test

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
double Recalc_Intersection_wgt(geoalgo::GeoAlgo const _geo_algo_instance, geoalgo::AABox volAVTPC, auto const& mcflux, auto const& mctruth, Detector Detector_ ){

	TRotation R;
	int retries{0};
	double enu = mctruth.GetNeutrino().Nu().E();
	TVector3 x3beam;
	TRandom3 fRnd;
	
	// Now get the weight
	double weight;

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
		// Pick a new point on the window in beam coords
		x3beam = fWin_Base_beam + (fRnd.Uniform() * fFluxWindowDir1) + (fRnd.Uniform() * fFluxWindowDir2);

		calcEnuWgt(mcflux, x3beam, enu, weight);

		// Get the nu ray direction in beam coords
		TVector3 xyzDk(mcflux.fvx,mcflux.fvy,mcflux.fvz);  // Origin of decay in beam coords
		TVector3 NuRay_Dir =  enu * (x3beam - xyzDk).Unit();
		
		// Convert to detector coordinates
		NuRay_Dir = R_Beam_2_Det * NuRay_Dir;

		// Now convert xbeam to detector coordinates too
		TVector3 x3beam_det = R_Beam_2_Det * (x3beam - Detector_.Trans_Det2Beam);
	
		// std::cout << "vx:\t" <<mctruth.GetNeutrino().Nu().Vx() << "  "<<  mcflux.fvx<<  "   " <<x3beam_det.X()/100.0 << std::endl;
		// std::cout << "vy:\t" <<mctruth.GetNeutrino().Nu().Vy() << "  "<<  mcflux.fvy<<  "   " <<x3beam_det.Y()/100.0 << std::endl;
		// std::cout << "vz:\t" <<mctruth.GetNeutrino().Nu().Vz() << "  "<<  mcflux.fvz<<  "   " <<x3beam_det.Z()/100.0 << std::endl;

		// Make the neutrino ray
		geoalgo::HalfLine ray(x3beam_det.X()/100.0, // point on window in detector coordinates units have to be m
					x3beam_det.Y()/100.0,
					x3beam_det.Z()/100.0,
					NuRay_Dir.X(),  // px in detector coordinates
					NuRay_Dir.Y(),  // py
					NuRay_Dir.Z()); // pz

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
		
		// Got an interception, so break!
		if (intercept) {
			// std::cout << "Passed with " << retries << " retries"<<  std::endl;
			break;
		}
		if (retries > 1000) {
			// If there are more than 1000 attempts and still no intersection then give up!
			std::cout << "Recalculation failed due to > 1000 tries and no interception" << std::endl;
			return 0.0; // throw the event away
		} 
	}
	
	return weight/3.1415926;
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