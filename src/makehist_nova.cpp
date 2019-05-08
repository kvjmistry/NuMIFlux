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
	// const double kRDET = 10.0;       // set to flux per 10 cm radius

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
	double x = gRandom->Uniform(-176    , 177);    //cm
	double y = gRandom->Uniform(-172    ,179);     //cm
	double z = gRandom->Uniform(25      ,1150);    //cm

	return TVector3(x, y, z);
}
//___________________________________________________________________________
TVector3 FromDetToBeam( const TVector3& det ) {

    TVector3 beam;
    TRotation R;
	bool debug{false};

	// rotation from beam to det
    TVector3 newX(9.9990e-01, -8.2300e-04, -1.4088e-02);
    TVector3 newY(3.0533e-06, 9.9831e-01, -5.8103e-02);
    TVector3 newZ(1.4112e-02, 5.8097e-02, 9.9821e-01);

    R.RotateAxes(newX,newY,newZ);
	R.Invert();
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
    // NuMIDet is vector from NuMI target to detector (in beam coordinates)
    TVector3 NuMIDet ( 1171.74545 ,       -331.51325 ,      99293.47347); // cm

    beam = R * det + NuMIDet;

    return beam;
}
//___________________________________________________________________________
// Get the window normal for the tiltweight
double Get_tilt_wgt( const TVector3& detxyz, auto const& mcflux, auto const& mctruth){

	TVector3 xyzDk(mcflux.fvx,mcflux.fvy,mcflux.fvz);  // origin of decay
	
	TVector3 p3beam = mctruth.GetNeutrino().Nu().E() * ( detxyz - xyzDk ).Unit();

	// Hardcoded for testing, but in priciple would want
	TVector3 windownorm = { 0.01411202435, 0.05809720554, 0.9982111828}; 

	double tiltweight =  p3beam.Unit().Dot( windownorm );

	return tiltweight;

}
//___________________________________________________________________________
// This function will recaclulate the missed nu rays for the intersection method
// The intension is that this will fix the normalistion problems
double Recalc_Intersection_wgt(geoalgo::GeoAlgo const _geo_algo_instance, geoalgo::AABox volAVTPC, auto const& mcflux, auto const& mctruth ){
    TRotation R;
	int retries{0};
	double Enu = mctruth.GetNeutrino().Nu().E();
	TVector3 x3beam;
	TRandom3 fRnd;

	// rotation from beam to det
    TVector3 newX(9.9990e-01, -8.2300e-04, -1.4088e-02);
    TVector3 newY(3.0533e-06, 9.9831e-01, -5.8103e-02);
    TVector3 newZ(1.4112e-02, 5.8097e-02, 9.9821e-01);
    
	// R is now  *** det to beam ***
	R.RotateAxes(newX,newY,newZ); 
	TRotation R_Beam_2_Det = R.Inverse();

	// Define the window in det coords
	TVector3 fWBase  = { 500, -250, -500 };
    TVector3 fWpt1  = { 500,  500,  -500 };
    TVector3 fWpt2  = { -500, -250, -500 };
	TVector3 NuMIDet { 226.9447, 6100.1882, -99113.1313}; // cm in detector coords

	// Convert to beam coordinates and define directions by subtracting the base
	fWBase = R * (fWBase - NuMIDet);
	TVector3 fWdir1 = R * (fWpt1 - NuMIDet) - fWBase;
	TVector3 fWdir2 = R * (fWpt2 - NuMIDet) - fWBase;

	// Now see if this neutrino vector is going to intersect with the detector
	// If it doesn't then recalculate
	
	while(true){

		retries++;
		// Pick a new point on the window in beam coords
		x3beam = fWBase + (fRnd.Uniform() * fWdir1) + (fRnd.Uniform() * fWdir2);

		// Get the nu ray direction in beam coords
		TVector3 NuRay_Dir = { (x3beam.X() - mcflux.fvx), (x3beam.Y() - mcflux.fvy), (x3beam.Z() - mcflux.fvz)  };
		
		// Convert to a momentum and detector coordinates
		NuRay_Dir = R_Beam_2_Det * (Enu * NuRay_Dir.Unit());

		// Now convert xbeam to detector coordinates too
		TVector3 x3beam_det = R_Beam_2_Det * x3beam + NuMIDet;

		// std::cout << "vx:\t" <<mctruth.GetNeutrino().Nu().Vx() << "  "<<  mcflux.fvx<<  "   " <<x3beam_det.X() << std::endl;
		// std::cout << "vy:\t" <<mctruth.GetNeutrino().Nu().Vy() << "  "<<  mcflux.fvy<<  "   " <<x3beam_det.Y() << std::endl;
		// std::cout << "vz:\t" <<mctruth.GetNeutrino().Nu().Vz() << "  "<<  mcflux.fvz<<  "   " <<x3beam_det.Z() << std::endl;

		// Make the neutrino ray
		geoalgo::HalfLine ray(x3beam_det.X(), // point on window in detector coordinates
					x3beam_det.Y(),
					x3beam_det.Z(),
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
			std::cout << "Passed with " << retries << " retries"<<  std::endl;
			break;
			
		}
		if (retries > 1000) {
			// If there are more than 1000 attempts and still no intersection then give up!
			std::cout << "Recalculation failed due to > 1000 tries and no interception" << std::endl;
			return 0.0; // throw the event away
		} 
	}
	
	// Now get the weight
	double weight;
	calcEnuWgt(mcflux, x3beam, weight);
	return weight/3.1415926;
}

//___________________________________________________________________________
int main(int argc, char** argv) {

	std::pair<float, float>  _xRange;
	std::pair<float, float>  _yRange;
	std::pair<float, float>  _zRange;

	_xRange.first  = -176;
	_xRange.second =  177;
	_yRange.first  = -172;
	_yRange.second =  179;
	_zRange.first  =   25;
	_zRange.second = 1150;

	geoalgo::GeoAlgo const _geo_algo_instance;

	geoalgo::AABox volAVTPC( _xRange.first, _yRange.first, _zRange.first, _xRange.second, _yRange.second, _zRange.second);

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
			double Enu = mctruth.GetNeutrino().Nu().E();


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
			TVector3 xyz_det = RandomInDet();
			
			// From detector to beam coordinates
			TVector3 xyz_beam = FromDetToBeam(xyz_det);

			// Get the new weight at the detector
			calcEnuWgt(mcflux, xyz_beam, detwgt);

			// Get the tiltweight
			tiltwght = Get_tilt_wgt(xyz_beam, mcflux, mctruth);

			// Weight of neutrino parent (importance weight) * Neutrino weight for a decay forced at center of near detector 
			cv_weight           *= mcflux.fnimpwt * detwgt / 3.1415926 * tiltwght ; // divide by area of circle equal to pi *r*r
			cv_weight_notilt    *= mcflux.fnimpwt * detwgt / 3.1415926; // for ppfx cases
			// window_weight       *= mcflux.fnimpwt * mcflux.fnwtfar; // mcflux.fnwtfar == mcflux.fnwtnear
			
			double window_weight_recalc         = Recalc_Intersection_wgt(_geo_algo_instance, volAVTPC, mcflux, mctruth );
			// window_weight                   *= mcflux.fnimpwt * window_weight_recalc; // Recalculated for every event, already divide by Pi
			if (intercept) window_weight       *= mcflux.fnimpwt * mcflux.fnwtfar; // mcflux.fnwtfar == mcflux.fnwtnear
			else window_weight                 *= mcflux.fnimpwt * window_weight_recalc; // Recalculated for every event
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
			if (intercept) Enu_CV_AV_TPC_intersection[pdg]->Fill(Enu, window_weight);
			// Enu_CV_AV_TPC_intersection[pdg]->Fill(Enu, window_weight);
			Enu_CV_AV_TPC_detweights[pdg]                 ->Fill(Enu, cv_weight);
			Enu_CV_AV_TPC_detweights_notilt[pdg]          ->Fill(Enu, cv_weight_notilt);


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

