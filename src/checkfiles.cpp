/**
 * Check files by opening them and seeing if they are zombie
 * Authors: J. Zennamo, A. Mastbaum
 */

#include <iostream>
#include <cstdlib>
#include <string>
#include <vector>
#include <fstream>
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
		else if (filein->TestBit(TFile::kRecovered)){
			// File was ran through a recovery procedure so skip as its unreadable
			std::cout << "ERROR: File Needed Recovery:\t" << argv[i] << std::endl;
			badfiles.push_back(string(argv[i]));
			filein->Close();
		}		
		else {
			std::cout << "FILE : " << argv[i] << std::endl; 
			filename.push_back(string(argv[i]));
			totalPOT+=100000; // 100 000 POT per dk2nu file 500 000 if running over nova files
			filein->Close();
		}
	}

	std::cout << "\nTotal POT read in:\t" << totalPOT << std::endl;

	std::cout << "\nUsing 1e5 POT per dk2nu file, the flux will be wrong if this is not the case!\n" << std::endl;


	std::cout << "BADFiles:" << std::endl;
	for (unsigned int i=0;i < badfiles.size(); i++ ){
		std::cout << badfiles.at(i) << std::endl;
	}

	std::ofstream outfile ("goodruns.txt");
	std::cout << "Now creating a file with the good runs..."<<std::endl;
	for (unsigned int i=0;i < filename.size(); i++ ){
		outfile << filename[i] << std::endl;
	}

	outfile.close();

	return 0;
}

