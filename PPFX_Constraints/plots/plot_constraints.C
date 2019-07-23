// A script that makes all the hadron constraint plots in one
// This script adapts the current pyroot scripts 

#include "functions.h"

void plot_constraints(){
    gStyle->SetOptStat(0); // say no to stats box

    bool boolfile, boolhist;
    TFile *f;
    
    // Get the file
    boolfile  = GetFile(f , "../NuMIFlux.root"); 
    if (boolfile == false) gSystem->Exit(0);

   

    std::vector<std::string> pars =  {
        "pionplus_NA49",
        "pionminus_NA49",
        //"kaonplus_NA49",
        //"kaonminus_NA49",
        "pionplus_MIPP",
        "pionminus_MIPP"
        //"kaonplus_MIPP",
        //"kaonminus_MIPP"
    };

    for (int k=0; k < pars.size();k++){
       
        TCanvas *c = new TCanvas();
        TH2D *hist;

        // Get the histogram
        boolhist  = GetHist(f, hist, pars.at(k)); 
        if (boolhist == false) gSystem->Exit(0);

        hist->RebinX(4);
        hist->RebinY(4);
        
        DrawSpecifiers(hist, pars.at(k));
        IncreaseLabelSize(hist);
        gStyle->SetTitleH(0.07);
        
        
        std::vector<TLine*> lineVector = GetLineVector(pars.at(k));

        for ( int i=0; i<lineVector.size(); i++){
            lineVector.at(i)->SetLineColor(kBlack);
            lineVector.at(i)->SetLineWidth(3);
            lineVector.at(i)->Draw();
        }

        c->Print(Form("%s.pdf", pars.at(k).c_str() ));

    }

}
