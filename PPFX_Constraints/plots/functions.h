// -----------------------------------------------------------------------------
// Function to increase to axes labels
void IncreaseLabelSize(TH1D* h){

    // h->GetXaxis()->SetRangeUser(0,3.5);
    h->GetXaxis()->SetLabelSize(0.05);
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetYaxis()->SetLabelSize(0.05);
    h->GetYaxis()->SetTitleSize(0.05);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.12);
}
// -----------------------------------------------------------------------------
void IncreaseLabelSize(TH2D* h){

    // h->GetXaxis()->SetRangeUser(0,3.5);
    h->GetXaxis()->SetLabelSize(0.05);
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetYaxis()->SetLabelSize(0.05);
    h->GetYaxis()->SetTitleSize(0.05);
    h->GetZaxis()->SetLabelSize(0.05);
    h->GetZaxis()->SetTitleSize(0.05);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.2);
    gPad->SetBottomMargin(0.14);
    h->SetMarkerSize(1.8);
    h->GetYaxis()->SetTitleOffset(1);
    // gPad->SetGridx(); 
}
// -----------------------------------------------------------------------------
bool GetHist(TFile* f, TH1D* &h, TString string){
    h = (TH1D*) f->Get(string);
    if (h == NULL) {
        std::cout << "\nfailed to get:\t" << string <<
         "\tThis histogram might not exist in the file\n" << std::endl;
        return false;
    }
    else {
        return true;
    }
}
// -----------------------------------------------------------------------------
bool GetHist(TFile* f, TH2D* &h, TString string){
    h = (TH2D*) f->Get(string);
    if (h == NULL) {
        std::cout << "\nfailed to get:\t" << string <<
         "\tThis histogram might not exist in the file\n" << std::endl;
        return false;
    }
    else {
        return true;
    }
}
// -----------------------------------------------------------------------------
bool GetTree(TFile* f, TTree* &T, TString string){
    T = (TTree*) f->Get(string);
    if (T == NULL) {
        std::cout << "\nfailed to get:\t" << string <<
         "\tThis tree might not exist in the file\n" << std::endl;
        return false;
    }
    else {
        return true;
    }
}
// -----------------------------------------------------------------------------
bool GetFile(TFile* &f , TString string){
    f = TFile::Open(string);
    
    if (f == NULL) {
        std::cout << "failed to get:\t" << string << 
        "\tThis file might not exist in the file" << std::endl;
        return false;
    }
    else {
        return true;
    }
}
// -----------------------------------------------------------------------------
void DrawSpecifiers(TH2D* hist, std::string histname){

    
    TLegend *leg;
    TLatex *text;
    
   
   

    if (histname == "pionplus_NA49"){
        hist->SetTitle("#pi^{+}; x_{F}; p_{T} [GeV/c]");
        hist->GetXaxis()->SetRangeUser(-0.15,1);
        hist->GetYaxis()->SetRangeUser(0,3.5);
        leg = new TLegend(.64, .82 , .68, .86);
        text = new TLatex(.41, .825, "NA49 Coverage");
    }
    else if (histname == "pionminus_NA49"){
        hist->SetTitle("#pi^{-}; x_{F}; p_{T} [GeV/c]");
        hist->GetXaxis()->SetRangeUser(-0.15,1);
        hist->GetYaxis()->SetRangeUser(0,3.5);
        leg = new TLegend(.64, .82 , .68, .86);
        text = new TLatex(.41, .825, "NA49 Coverage");
    }
    else if (histname == "kaonplus_NA49"){
        hist->SetTitle("K^{+}; x_{F}; p_{T} [GeV/c]");
        hist->GetXaxis()->SetRangeUser(-0.04,1);
        hist->GetYaxis()->SetRangeUser(0,3.5);
        leg = new TLegend(.40, .75 , .44 , .79);
        text = new TLatex(.17, .755, "NA49 Coverage");
    }
    else if (histname == "kaonminus_NA49"){
        hist->SetTitle("K^{-}; x_{F}; p_{T} [GeV/c]");
        hist->GetXaxis()->SetRangeUser(-0.04,1);
        hist->GetYaxis()->SetRangeUser(0,3.5);
        leg = new TLegend(.40, .75 , .44 , .79);
        text = new TLatex(.17, .755, "NA49 Coverage");
    }
    else if (histname == "pionplus_MIPP"){
        hist->SetTitle("#pi^{+}; p_{Z}; p_{T} [GeV/c]");
         hist->GetXaxis()->SetRangeUser(0,130);
        hist->GetYaxis()->SetRangeUser(0,5);
        leg = new TLegend(.44, .55, .48, .59);
        text = new TLatex(.5, .55, "MIPP Coverage");
    }
    else if (histname == "pionminus_MIPP"){
        hist->SetTitle("#pi^{-}; p_{Z}; p_{T} [GeV/c]");
        hist->GetXaxis()->SetRangeUser(0,130);
        hist->GetYaxis()->SetRangeUser(0,5);
        leg = new TLegend(.44, .55, .48, .59);
        text = new TLatex(.5, .55, "MIPP Coverage");
    }
    else if (histname == "kaonplus_MIPP"){
        hist->SetTitle("K^{+}; p_{Z}; p_{T} [GeV/c]");
        //hist->GetXaxis()->SetRangeUser(-0.15,1);
        //hist->GetYaxis()->SetRangeUser(0,3.5);
        leg = new TLegend(.44, .55, .48, .59);
        text = new TLatex(.5, .55, "MIPP Coverage");
    }
    else if (histname == "kaonminus_MIPP"){
        hist->SetTitle("K^{-}; p_{Z}; p_{T} [GeV/c]");
        //hist->GetXaxis()->SetRangeUser(-0.15,1);
        //hist->GetYaxis()->SetRangeUser(0,3.5);
        leg = new TLegend(.44, .55, .48, .59);
        text = new TLatex(.5, .55, "MIPP Coverage");
    }
    else return;

    hist->Draw("colz");

    leg->SetFillStyle(0);
    leg->SetLineWidth(4);
    
    text->SetTextColor(kBlack);
    text->SetNDC();
    text->SetTextSize(1.4/30.);
    text->SetTextAlign(11);
    text->Draw();
    leg->Draw();


}
// -----------------------------------------------------------------------------
std::vector<TLine*> GetLineVector(std::string histname){

     std::vector<TLine*> lineVector;

    if (histname == "pionplus_NA49" || histname == "pionminus_NA49"){
        TLine *line1  = new TLine( -0.1125 ,   0.35 , -0.0875  ,  0.35);
        TLine *line2  = new TLine( -0.0875 ,   0.35   ,  -0.0875  ,  0.25);
        TLine *line3  = new TLine( -0.0875 ,   0.25 ,   -0.0625 ,   0.25);
        TLine *line4  = new TLine( -0.0625 ,   0.25,    -0.0625  ,  0.75);
        TLine *line5  = new TLine( -0.0625  ,  0.75,    -0.055  ,  0.75);
        TLine *line6  = new TLine( -0.055  ,  0.75,  -0.055   ,  0.35);
        TLine *line7  = new TLine( -0.055   ,  0.35 ,    0.17  ,  0.35);
        TLine *line8  = new TLine( 0.17  ,  0.35,   0.17  ,  0.65);
        TLine *line9  = new TLine( 0.17  ,  0.65 ,   0.175  ,  0.65);
        TLine *line10 = new TLine( 0.175  ,  0.65   ,0.175  ,  0.35);
        TLine *line11 = new TLine( 0.175  ,  0.35, 0.225  ,  0.35);
        TLine *line12 = new TLine( 0.225  ,  0.35,  0.225  ,  0.325);
        TLine *line13 = new TLine( 0.225  ,  0.325, 0.175 ,   0.325);
        TLine *line14 = new TLine( 0.175 ,   0.325  ,   0.175  ,  0.025);
        TLine *line15 = new TLine( 0.175  ,  0.025, 0.225  ,  0.025);
        TLine *line16 = new TLine(0.225  ,  0.025  ,   0.225  ,  0.045);
        TLine *line17 = new TLine(0.225  ,  0.045,    0.325  ,  0.045);
        TLine *line18 = new TLine(0.325  ,  0.045, 0.325  ,  0.65);
        TLine *line19 = new TLine(0.325  ,  0.65,   0.225  ,  0.65);
        TLine *line20 = new TLine(0.225  ,  0.65  ,   0.225   , 0.7);
        TLine *line21 = new TLine(0.225   , 0.7  ,   0.325  ,  0.7);
        TLine *line22 = new TLine(0.325  ,  0.7   ,   0.325  ,  1.3);
        TLine *line23 = new TLine(0.325  ,  1.3,  0.35  ,  1.3);
        TLine *line24 = new TLine(0.35  ,  1.3    ,0.35  ,  0.7);
        TLine *line25 = new TLine(0.35  ,  0.7 ,  0.55 ,   0.7);
        TLine *line26 = new TLine( 0.55 ,   0.7    ,   0.55  ,  1.5);
        TLine *line27 = new TLine( 0.55  ,  1.5    ,   0.45 ,   1.5);
        TLine *line28 = new TLine( 0.45 ,   1.5,  0.45 ,   1.7);
        TLine *line29 = new TLine( 0.45 ,   1.7    ,0.25  ,  1.7);
        TLine *line30 = new TLine( 0.25  ,  1.7    ,   0.25 ,   1.9);
        TLine *line31 = new TLine( 0.25 ,   1.9    ,   0.15  ,  1.9);
        TLine *line32 = new TLine( 0.15  ,  1.9    ,0.15   , 1.3);
        TLine *line33 = new TLine( 0.15   , 1.3    ,0.125  ,  1.3);
        TLine *line34 = new TLine( 0.125  ,  1.3   ,   0.125  ,  1.9);
        TLine *line35 = new TLine( 0.125  ,  1.9   ,   -0.125  ,  1.9);
        TLine *line36 = new TLine( -0.125   , 1.9  ,   -0.125  ,  1.3);
        TLine *line37 = new TLine( -0.125  ,  1.3  ,   -0.1125 ,   1.3);
        TLine *line38 = new TLine( -0.1125 ,   1.3 ,   -0.1125    ,1.1);
        TLine *line39 = new TLine( -0.1125    ,1.1    ,0.1125   , 1.1);
        TLine *line40 = new TLine( 0.1125    ,1.1 ,   0.1125   , 1.3);
        TLine *line41 = new TLine( 0.1125   , 1.3 ,   0.125  ,  1.3);
        TLine *line42 = new TLine(0.125  ,  1.3   ,   0.125  ,  1.1);
        TLine *line43 = new TLine( 0.125  ,  1.1   ,   0.225   , 1.1);
        TLine *line44 = new TLine( 0.225   , 1.1   ,   0.225  ,  1.05);
        TLine *line45 = new TLine( 0.225  ,  1.05   ,   0.125  ,  1.05);
        TLine *line46 = new TLine( 0.125  ,  1.05  ,   0.125  ,  0.65);
        TLine *line47 = new TLine( 0.125  ,  0.65  ,   0.1125  ,  0.65);
        TLine *line48 = new TLine( 0.1125  ,  0.65 ,   0.1125  ,  1.05);
        TLine *line49 = new TLine( 0.1125  ,  1.05 ,   -0.1125 ,   1.05);
        TLine *line50 = new TLine(-0.1125 ,   1.05,    -0.1125, 0.35);

        // box 2
        TLine *line51 = new TLine(-0.055 ,   0.325  ,  -0.055, 0.175);
        TLine *line52 = new TLine(-0.055, 0.175,    -0.045, 0.175);
        TLine *line53 = new TLine(-0.045, 0.175,  -0.045  ,  0.12);
        TLine *line54 = new TLine(-0.045  ,  0.12, -0.035  ,  0.125);
        TLine *line55 = new TLine( -0.035  ,  0.125,    -0.035 ,   0.075);
        TLine *line56 = new TLine( -0.035  ,  0.075,   -0.025  ,  0.075);
        TLine *line57 = new TLine( -0.025  ,  0.075,    -0.025  ,  0.025);
        TLine *line58 = new TLine( -0.025 ,   0.025,   0.1625  ,  0.025);
        TLine *line59 = new TLine( 0.1625  ,  0.025,   0.1625  ,  0.325);
        TLine *line60 = new TLine(  0.1625  ,  0.325, -0.055 ,   0.325 );

        // box 3

        TLine *line61 = new TLine( 0.35, 0.65, 0.35, 0.045 );
        TLine *line62 = new TLine( 0.35, 0.045, 0.55, 0.045 );
        TLine *line63 = new TLine( 0.55, 0.045, 0.55, 0.65 );
        TLine *line64 = new TLine( 0.55, 0.65, 0.35, 0.65 );
        
        lineVector = { line1, line2, line3, line4, line5, line6, line7, line8, line9,
                    line10, line11, line12, line13, line14, line15, line16, line17,
                    line18, line19, line20, line21, line22, line23, line24, line25,
                    line26, line27, line28, line29, line30, line31, line32, line33,
                    line34, line35, line36, line37, line38, line39, line40, line41,
                    line42, line43, line44, line45, line46, line47, line48, line49,
                    line50, line51, line52, line53, line54, line55, line56, line57,
                    line58, line59,line60, line61, line62, line63, line64 };
    }
    else if (histname == "kaonplus_NA49" || histname == "kaonminus_NA49"){
        //box 1
        TLine *line1 = new TLine( -0.0125, 0.05, 0.225, 0.05 );
        TLine *line2 = new TLine( 0.225, 0.05, 0.225, 0.75 );
        TLine *line3 = new TLine( 0.225, 0.75 , -0.0125, 0.75 );
        TLine *line4 = new TLine( -0.0125, 0.75 , -0.0125 ,0.05 );

        //box 2
        TLine *line5 = new TLine( -0.0125, 0.8, 0.225, 0.8 );
        TLine *line6 = new TLine( 0.225, 0.8, 0.225, 1 );
        TLine *line7 = new TLine( 0.225, 1 , -0.0125, 1 );
        TLine *line8 = new TLine( -0.0125, 1 , -0.0125 ,0.8 );

        lineVector= { line1, line2, line3, line4, line5, line6, line7, line8 };
    }
    else if (histname == "pionplus_MIPP" || histname == "pionminus_MIPP" ||
             histname == "kaonplus_MIPP" || histname == "kaonminus_MIPP"){

            TLine *line1  = new TLine( 20, 0, 20, 1 );
            TLine *line2  = new TLine( 20, 1, 24, 1 );
            TLine *line3  = new TLine( 24, 1, 24, 1.2 );
            TLine *line4  = new TLine( 24, 1.2, 31, 1.2 );
            TLine *line5  = new TLine( 31, 1.2, 31, 1.55 );
            TLine *line6  = new TLine( 31, 1.55, 42, 1.55 );
            TLine *line7  = new TLine( 42, 1.55, 42, 2 );
            TLine *line8  = new TLine( 42, 2, 60, 2 );
            TLine *line9  = new TLine( 60, 2, 90, 2 );
            TLine *line10 = new TLine( 90, 2, 90, 0 );
            TLine *line11 = new TLine( 90, 0, 20, 0 );

            lineVector = { line1, line2, line3, line4, line5, line6, line7,
                           line8, line9, line10, line11 };
    }
    else return lineVector;


    return lineVector;
}
