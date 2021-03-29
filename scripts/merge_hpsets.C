
#include "functions.h"
// Root script to get the different sets of HP Universes and merge them into 1 file for easier use.
// Will also scale the histograms to per pot which is easier to work with


std::vector<std::string> cv_names= {
    "CV_AV_TPC",
    "UW_AV_TPC",
    "CV_AV_TPC_5MeV_bin",
    "UW_AV_TPC_5MeV_bin",
    "CV_AV_TPC_2D",
    "unweighted_AV_TPC_2D"
};

std::vector<std::string> flav = {
    "numu",
    "nue",
    "numubar",
    "nuebar"
};

std::vector<std::string> inputmode = {
        "ppfx_mippk_PPFXMIPPKaon",
        "ppfx_mipppi_PPFXMIPPPion",
        "ppfx_other_PPFXOther",
        "ppfx_targatt_PPFXTargAtten",
        "ppfx_think_PPFXThinKaon",
        "ppfx_thinmes_PPFXThinMeson",
        "ppfx_thinnpi_PPFXThinNeutronPion",
        "ppfx_thinna_PPFXThinNucA",
        "ppfx_thinn_PPFXThinNuc",
        "ppfx_thinpi_PPFXThinPion",
        "ppfx_totabs_PPFXTotAbsorp",
        "ppfx_ms_UBPPFX"};


void merge_hpsets(){

    int n_uni = 600;

    // Open the TFiles
    TFile *f_set1 = TFile::Open("/uboone/data/users/kmistry/work/PPFX/uboone/beamline_zero_threshold_v46/FHC/sept/output_uboone_fhc_run0_set1.root");
    TFile *f_set2 = TFile::Open("/uboone/data/users/kmistry/work/PPFX/uboone/beamline_zero_threshold_v46/FHC/sept/output_uboone_fhc_run0_set2.root");
    TFile *f_set3 = TFile::Open("/uboone/data/users/kmistry/work/PPFX/uboone/beamline_zero_threshold_v46/FHC/sept/output_uboone_fhc_run0_set3.root");

    // Get the POT for each file
    double pot_set1 =  GetPOT(f_set1);
    double pot_set2 =  GetPOT(f_set2);
    double pot_set3 =  GetPOT(f_set3);

    // Now we want to define the histograms
    std::vector<std::vector<TH1D*>> h_cv_vec;
    h_cv_vec.resize(flav.size());
    for (unsigned int f = 0; f < h_cv_vec.size(); f++){
        h_cv_vec.at(f).resize(cv_names.size());
    }

    
    // flavour -- inputmode -- universe
    std::vector<std::vector<std::vector<TH1D*>>> h_univ_1D;
    std::vector<std::vector<std::vector<TH2D*>>> h_univ_2D;

    h_univ_1D.resize(flav.size());
    h_univ_2D.resize(flav.size());
    
    for (unsigned int f = 0; f < h_univ_1D.size(); f++){
        h_univ_1D.at(f).resize(inputmode.size());
        h_univ_2D.at(f).resize(inputmode.size());
    }

    
    // Loop over flavours
    for (unsigned int f = 0; f < h_univ_1D.size(); f++){

        // Loop over input modes
        for (unsigned int t = 0; t < h_univ_1D.at(f).size(); t++){
            h_univ_1D.at(f).at(t).resize(n_uni);
            h_univ_2D.at(f).at(t).resize(n_uni);
        }
    }
    

    // Now lets get the histgrams from the file
    bool boolhist;
    
    // CV histograms -- we only need 1 hist, so get it from set 1

    // Loop over flavours
    for (unsigned int f = 0; f < h_cv_vec.size(); f++){
        
        // loop over cv types
        for (unsigned int t = 0; t < h_cv_vec.at(f).size(); t++){
            boolhist = GetHist(f_set1, h_cv_vec.at(f).at(t), Form("%s/Detsmear/%s_%s", flav.at(f).c_str(), flav.at(f).c_str(), cv_names.at(t).c_str())); 
    
            if (boolhist == false) gSystem->Exit(0);
            
            // Scale its POT
            h_cv_vec.at(f).at(t)->Scale(1.0 / pot_set1);

        }
    }

    // Now lets get the multisims

    // Loop over the neutrino flavours
    for (unsigned int f = 0; f < h_univ_1D.size(); f++){
        
        // loop over input modes
        for (unsigned int t = 0; t < h_univ_1D.at(f).size(); t++){
            
            // Loop over the universes
            for (unsigned int uni = 0; uni < h_univ_1D.at(f).at(t).size(); uni++){
                
                if (uni < 200) {
                    
                    boolhist = GetHist(f_set1, h_univ_1D.at(f).at(t).at(uni), Form("%s/Multisims/%s_%s_Uni_%i_AV_TPC",    flav.at(f).c_str(), flav.at(f).c_str(), inputmode.at(t).c_str(), uni)); 
                    boolhist = GetHist(f_set1, h_univ_2D.at(f).at(t).at(uni), Form("%s/Multisims/%s_%s_Uni_%i_AV_TPC_2D", flav.at(f).c_str(), flav.at(f).c_str(), inputmode.at(t).c_str(), uni));

                    if (boolhist == false) gSystem->Exit(0);
                
                    // Scale its POT
                    h_univ_1D.at(f).at(t).at(uni)->Scale(1.0 / pot_set1);
                    h_univ_2D.at(f).at(t).at(uni)->Scale(1.0 / pot_set1); 
                }
                // Second set of universes
                else if (uni >= 200 && uni < 400) {
                    
                    boolhist = GetHist(f_set2, h_univ_1D.at(f).at(t).at(uni), Form("%s/Multisims/%s_%s_Uni_%i_AV_TPC",    flav.at(f).c_str(), flav.at(f).c_str(), inputmode.at(t).c_str(), uni)); 
                    boolhist = GetHist(f_set2, h_univ_2D.at(f).at(t).at(uni), Form("%s/Multisims/%s_%s_Uni_%i_AV_TPC_2D", flav.at(f).c_str(), flav.at(f).c_str(), inputmode.at(t).c_str(), uni)); 
                    
                    if (boolhist == false) gSystem->Exit(0);
                
                    // Scale its POT
                    h_univ_1D.at(f).at(t).at(uni)->Scale(1.0 / pot_set2);
                    h_univ_2D.at(f).at(t).at(uni)->Scale(1.0 / pot_set2); 
                }
                else {
                    
                    boolhist = GetHist(f_set3, h_univ_1D.at(f).at(t).at(uni), Form("%s/Multisims/%s_%s_Uni_%i_AV_TPC",    flav.at(f).c_str(), flav.at(f).c_str(), inputmode.at(t).c_str(), uni)); 
                    boolhist = GetHist(f_set3, h_univ_2D.at(f).at(t).at(uni), Form("%s/Multisims/%s_%s_Uni_%i_AV_TPC_2D", flav.at(f).c_str(), flav.at(f).c_str(), inputmode.at(t).c_str(), uni)); 
                    
                    if (boolhist == false) gSystem->Exit(0);
                
                    // Scale its POT
                    h_univ_1D.at(f).at(t).at(uni)->Scale(1.0 / pot_set3);
                    h_univ_2D.at(f).at(t).at(uni)->Scale(1.0 / pot_set3); 

                }
                
            } // end loop over universes

        } // end loop over input modes
    
    } // end loop over flavours

    // Now we have the histograms we need to create a TFile, the directories and write the histograms
    std::cout << "Got the histograms so now writing to file" << std::endl;
    TFile *f_out = new TFile("f_out.root", "UPDATE");
    TTree* POTTree = new TTree("POT","Total POT");
    double totalPOT = 1.0;
    POTTree -> Branch("POT", &totalPOT);
    POTTree -> Fill();
    POTTree->Write("",TObject::kOverwrite);

    std::vector<std::string> folders = {"Detsmear", "Multisims"};

    // Create subdirectory for each flavour
    TDirectory *dir_flav[flav.size()];

    // Create subdirectory for either the CV or multisims
    TDirectory *dir_labels_var[folders.size()];

    // Loop over the flavours
    for (unsigned int f = 0; f < flav.size(); f++) {
        
        // See if the directory already exists
        bool bool_dir = GetDirectory(f_out, dir_flav[f], flav.at(f).c_str());

        // If it doesnt exist then create it
        if (!bool_dir) dir_flav[f] = f_out->mkdir(flav.at(f).c_str());

        // Go into the directory
        dir_flav[f]->cd();

        // Loop over the directories
        for (unsigned int fold = 0; fold < folders.size(); fold ++){

            // See if the directory already exists
            bool bool_dir = GetDirectory(f_out, dir_labels_var[fold], Form("%s/%s", flav.at(f).c_str(), folders.at(fold).c_str()));

            // If it doesnt exist then create it
            if (!bool_dir) dir_labels_var[fold] = dir_flav[f]->mkdir(folders.at(fold).c_str());

            // Go into the directory
            dir_labels_var[fold]->cd();
        
            // Now write the histograms, 
            if (folders.at(fold) == "Detsmear"){
                for (unsigned int p = 0; p < h_cv_vec.at(f).size(); p++){
                    h_cv_vec.at(f).at(p)->Write("",TObject::kOverwrite);
                }
            }
            else {

                // loop over input modes
                for (unsigned int t = 0; t < h_univ_1D.at(f).size(); t++){
                    
                    // Loop over the universes
                    for (unsigned int uni = 0; uni < h_univ_1D.at(f).at(t).size(); uni++){
                        h_univ_1D.at(f).at(t).at(uni)->Write("",TObject::kOverwrite);
                        h_univ_2D.at(f).at(t).at(uni)->Write("",TObject::kOverwrite);

                    } // End loop over universes
                
                } // End loop over input modes
                
            }
        
        } // End loop over the variables

        f_out->cd();    // change current directory to top

    }
    f_out->cd();    // change current directory to top

    f_out->Close();

}
