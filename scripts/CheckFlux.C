// Script to make sure that the fluxes are being calculated properly
#include "/uboone/app/users/kmistry/PPFX/numi-validation/scripts/plot_comp_functions.h"

void CheckFlux(){

	TH1D* h_g_simp_nuebar, *hdk2nu_nuebar, *hdk2nu_nue, *hppfx_nue, *hppfx_nuebar, *hppfx_mod_nue, *hppfx_mod_nuebar;
	TFile* f_ppfx2d;
	TH2D* hppfx2d, * hppfx2d_nuebar;
    bool boolfile, boolhist;

    TH1D *h_g_simp, *hnovafileflux, *hppfx; 
	TFile *f_gsimple, *f_novafiles, *f_ppfx, *f_dk2nu, *f_ppfx_mod;
	
	// ----------------- GSimple ----------------------
	boolfile  = GetFile(f_gsimple , "/uboone/data/users/kmistry/work/PPFX/uboone/NuMIFlux_update_morebins.root"); if (boolfile == false) gSystem->Exit(0);
    // double POT_gsimp = GetPOT(f_gsimple);
    boolhist  = GetHist(f_gsimple, h_g_simp, "nueFluxHisto"); if (boolhist == false) gSystem->Exit(0);
    boolhist  = GetHist(f_gsimple, h_g_simp_nuebar, "anueFluxHisto"); if (boolhist == false) gSystem->Exit(0);
	// Normalise(h_g_simp_nuebar);
	double gsimp_flux_nuebar = IntegrateHist1D(h_g_simp_nuebar);
	double gsimp_flux_nue    = IntegrateHist1D(h_g_simp);	

	// ----------------- Dk2nu ----------------------
    boolfile  = GetFile(f_dk2nu , "/uboone/data/users/kmistry/work/PPFX/gsimple/gsimple_dk2nu_bugfix/output.root"); if (boolfile == false) gSystem->Exit(0);
	double POT_dk2nu = GetPOT(f_dk2nu);
    boolhist = GetHist(f_dk2nu, hdk2nu_nue, "nue/nue_CV_AV_TPC"); if (boolhist == false) gSystem->Exit(0);
	boolhist = GetHist(f_dk2nu, hdk2nu_nuebar, "nuebar/nuebar_CV_AV_TPC"); if (boolhist == false) gSystem->Exit(0);
	// Normalise(hdk2nu_nue);
    // Normalise(hdk2nu_nuebar);
	double flux_dk2nu_nue    = IntegrateHist1D(hdk2nu_nue); 
	double flux_dk2nu_nuebar = IntegrateHist1D(hdk2nu_nuebar);

	// ----------------- PPFX ----------------------
    boolfile  = GetFile(f_ppfx , "/uboone/data/users/kmistry/work/PPFX/uboone/bugfix_release_notilt/1D_finebinning/output.root"); if (boolfile == false) gSystem->Exit(0);
	double POT_ppfx = GetPOT(f_ppfx);
    boolhist  = GetHist(f_ppfx, hppfx_nue, "nue/nue_CV_AV_TPC"); if (boolhist == false) gSystem->Exit(0);
	boolhist  = GetHist(f_ppfx, hppfx_nuebar, "nuebar/nuebar_CV_AV_TPC"); if (boolhist == false) gSystem->Exit(0);
	// Normalise(hppfx_nue);
    // Normalise( hppfx_nuebar);
	double flux_ppfx_nue    = IntegrateHist1D(hppfx_nue); 
	double flux_ppfx_nuebar = IntegrateHist1D(hppfx_nuebar);

    // ----------------- PPFX_modified weights +tilt----------------------
    boolfile  = GetFile(f_ppfx_mod , "/uboone/data/users/kmistry/work/PPFX/uboone/DetectorWeights_withtilt/output.root"); if (boolfile == false) gSystem->Exit(0);
	double POT_ppfx_mod = GetPOT(f_ppfx_mod);
    boolhist  = GetHist(f_ppfx_mod, hppfx_mod_nue, "nue/nue_CV_AV_TPC"); if (boolhist == false) gSystem->Exit(0);
	boolhist  = GetHist(f_ppfx_mod, hppfx_mod_nuebar, "nuebar/nuebar_CV_AV_TPC"); if (boolhist == false) gSystem->Exit(0);
	// Normalise(hppfx_nue);
    // Normalise( hppfx_nuebar);
	double flux_ppfx_mod_nue    = IntegrateHist1D(hppfx_mod_nue); 
	double flux_ppfx_mod_nuebar = IntegrateHist1D(hppfx_mod_nuebar);
	
	// ----------------- PPFX 2D with tilt ----------------------
	boolfile  = GetFile(f_ppfx2d , "/uboone/data/users/kmistry/work/PPFX/uboone/DetectorWeights_withtilt/2D/output.root"); if (boolfile == false) gSystem->Exit(0); // with tilt
	// boolfile  = GetFile(f_ppfx2d , "/uboone/data/users/kmistry/work/PPFX/uboone/beamline_modified/run15/output.root"); if (boolfile == false) gSystem->Exit(0); // run15 beamline nominal
	double POT_ppfx2d = GetPOT(f_ppfx2d);
    boolhist  = GetHist(f_ppfx2d, hppfx2d, "nue/nue_CV_AV_TPC"); if (boolhist == false) gSystem->Exit(0);
	boolhist  = GetHist(f_ppfx2d, hppfx2d_nuebar, "nuebar/nuebar_CV_AV_TPC"); if (boolhist == false) gSystem->Exit(0);
	// Normalise(hppfx2d);
	// Normalise(hppfx2d_nuebar);
	TH2D* hppfx2d_clone = (TH2D*)  hppfx2d->Clone("hCV2d"); // Clone for alternative way of normalising the flux
	hppfx2d_clone->Add(hppfx2d_nuebar); // Combine the fluxes

	double flux_ppfx_nue_2D    = IntegrateHist2D(hppfx2d); // divide out by the gsimple POT to scale to coltons POT
	double flux_ppfx_nuebar_2D = IntegrateHist2D(hppfx2d_nuebar);
	double flux_ppfx_clone_2D  = IntegrateHist2D(hppfx2d_clone);

	// Print
	std::cout << "gsimple flux:\t\t"    << (gsimp_flux_nue    + gsimp_flux_nuebar)    * (2.334e+20 / 6e20)                     << std::endl; // 6e20 POT becasue thats how much pot is in the gsimple, we then scale this to data POT
	std::cout << "dk2nu flux:\t\t"      << (flux_dk2nu_nue    + flux_dk2nu_nuebar)    * ((2.334e+20 ) / (POT_dk2nu*1.0e4) )    << std::endl; 
	std::cout << "ppfx flux:\t\t"       << (flux_ppfx_nue     + flux_ppfx_nuebar)     * ((2.334e+20 ) / (POT_ppfx*1.0e4) )     << std::endl; 
	std::cout << "ppfx flux 2d:\t\t"    << (flux_ppfx_nue_2D  + flux_ppfx_nuebar_2D)  * (2.334e+20 / (POT_ppfx2d*1e4) )        << std::endl; 
	std::cout << "ppfx flux 2d swich:\t"  << (flux_ppfx_clone_2D)                     * (2.334e+20 / (POT_ppfx2d*1e4) )        << std::endl; 
    std::cout << "ppfx flux modifed:\t" << (flux_ppfx_mod_nue + flux_ppfx_mod_nuebar) * ((2.334e+20 ) / (POT_ppfx_mod*1.0e4) ) << std::endl; 

    // Plot the flux as a function of energy
    for (double E_th = 0; E_th < 5; E_th+=0.05){
		gsimp_flux_nuebar = IntegrateHist1D(h_g_simp_nuebar, E_th);
		gsimp_flux_nue    = IntegrateHist1D(h_g_simp, E_th);
		double gsimp_flux = (gsimp_flux_nue + gsimp_flux_nuebar) * (2.334e+20 / 6e20);

		flux_ppfx_mod_nue    = IntegrateHist1D(hppfx_mod_nue, E_th); 
	    flux_ppfx_mod_nuebar = IntegrateHist1D(hppfx_mod_nuebar, E_th);

		double ppfx_mod_flux = ((flux_ppfx_mod_nue + flux_ppfx_mod_nuebar) * ((2.334e+20 ) / (POT_ppfx_mod*1.0e4)) );

		double flux_rat = gsimp_flux / ppfx_mod_flux ;

		// std::cout << "E_th[GeV]:\t"<< E_th << "\tGsimp/PPFX:\t"<< flux_rat << "   "<< ppfx_mod_flux << "   " << gsimp_flux << std::endl;	

	}


    gSystem->Exit(0);
}

