
// Macro will loop over files in a directory and see if a specific tree exists. If it doesn't then dont want to run over these files.
// Creates a new file with the goodruns
#include <fstream>

// ------------------------------------------------------------------------------------------------------------
bool GetFile(TFile* &f , TString string){
	f = TFile::Open(string);
	
	if (f == NULL) {
		// std::cout << "failed to get:\t" << string << "\tThis file might not exist in the file" << std::endl;
		return false;
	}
	else {
		return true;
	}
}
// ------------------------------------------------------------------------------------------------------------
void checkfiles(){
    bool boolfile{false};
    std::string run{"12"}; // 09 or 10 etc.
    int numfiles{1000}; // change this to 5000, 1000 etc.
    std::ofstream outfile (Form("novafilelist_Run00%s_goodruns.txt", run.c_str()));
    std::ofstream badfile (Form("badfiles_run00%s.txt", run.c_str()));

    for (int i = 0; i < numfiles; i++){ 

        std::string inname =   Form("/pnfs/numix/flux/g4numi_syst_3rd_ana/test_g4numi_v6/me000z200i/minervame/run00%s/g4numiv6_minervame_me000z200i_%i_00%s.root", run.c_str(), i, run.c_str());
        
        std::cout << "FILE:\t" << inname<< std::endl;
        
        // Get the file
        TFile *fIn;
        boolfile  = GetFile(fIn , inname); 
        if (boolfile == false) continue;

        // File was Zombie
        if (fIn->IsZombie()) {
			std::cout << "\033[1;31mERROR: File is ZOMBIE\033[0m\t" << std::endl;
            badfile << inname << std::endl;// file was bad so add to bad runs
			fIn->Close();
		}
        // File was ran through a recovery procedure so skip as its unreadable
		else if (fIn->TestBit(TFile::kRecovered)){
			std::cout << "\033[1;31mERROR: File Needed Recovery\033[0m\t" << std::endl;
            badfile << inname << std::endl;// file was bad so add to bad runs
			fIn->Close();
		}	
        else {
            TTree* TMeta = (TTree*) fIn->Get("dkmetaTree"); // Try getting the tree
            if (TMeta == NULL) {
                std::cout << "\033[1;31mError cant get Dk2Nu Meta Info\033[0m" << std::endl;
                badfile << inname << std::endl;// file was bad so add to bad runs
            }
            else outfile << inname << std::endl;// file was good so add to good runs

            fIn->Close();
        }
    
    }
}
// ------------------------------------------------------------------------------------------------------------