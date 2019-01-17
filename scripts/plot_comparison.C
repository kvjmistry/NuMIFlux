/**
 * NuMI Validation Plotting
 *
 * A. Mastbaum <mastbaum@uchicago.edu> 2018/11
 */

void plot_comparison(TString  inputfile) {
  gStyle->SetOptStat(0);

  // Choose one mode
  // std::string inputmode = "ms_PPFX"; 
  std::string inputmode = "PPFXMaster"; 
  // std::string inputmode = "Total"; 

  std::cout << "\nUsing "<< inputmode << " Mode!" << std::endl;
  std::cout << "If the code breaks after the second plot, might want to chnage this...\n" << std::endl;

  enum type {knumu, knue};

  type mode = knumu; // select whether plotting numu or nue
 
  if (mode == knumu) std::cout << "\nUsing NuMu Mode!\n" << std::endl;
  if (mode == knue) std::cout << "\nUsing Nue Mode!\n" << std::endl;

  // Root is dumb and so need to pre-decalre here
  TDirectory* d;
  TH1D* h; 
  TH1D* h2;

  // Overlay of output plot with official NOvA FHC numu flux
  TFile* f1 = TFile::Open(inputfile);
  TCanvas* c1 = new TCanvas();

  if (mode == knumu) d = (TDirectory*)f1->Get("numu");
  if (mode == knue) d = (TDirectory*)f1->Get("nue");

  d->cd();
  if (mode == knumu) h = (TH1D*) (gDirectory->Get("numu_CV_AV_TPC")->Clone("fx"));
  if (mode == knue) h = (TH1D*) (gDirectory->Get("nue_CV_AV_TPC")->Clone("fx"));
  h->SetDirectory(0);
  for (int i=0;i<h->GetNbinsX()+1;i++) {
    h->SetBinContent(i, h->GetBinContent(i)/h->GetBinWidth(i));
  }
  h->Sumw2();
  h->SetLineColor(kRed);
  h->SetLineWidth(2);
  h->SetTitle(";E_{#nu} (GeV);Fraction/GeV");
  h->Draw("");
  TH1D* horig = (TH1D*) h->Clone("horig");

  TFile* f2 = TFile::Open("/uboone/app/users/kmistry/PPFX/numi-validation/nova_flux/FHC_Flux_NOvA_ND_2017.root");
  f2->ls();
  if (mode == knumu) h2 = (TH1D*) f2->Get("flux_numu");
  if (mode == knue) h2 = (TH1D*) f2->Get("flux_nue");
  h2->SetLineColor(kBlue);
  h2->SetLineWidth(2);

  h2->Draw("same");
  h->GetYaxis()->SetTitle(h2->GetYaxis()->GetTitle());
  h->Scale(h2->Integral(2,-1)/h->Integral(2,-1));
  // h->Scale(1.0/(0.356e4)); // Try calculating norm


  TH1D* hratio = (TH1D*) h->Clone("hratio");
  hratio->SetDirectory(0);

  TLegend* l = new TLegend(0.5, 0.6, 0.7, 0.8);
  l->AddEntry(h2, "NOvA","l");
  l->AddEntry(h, "This work","l");
  l->Draw();

  gPad->SetLogy();
  gPad->Update();

  // Ratio of ours to NOvA
  TCanvas* c2 = new TCanvas();
  hratio->Divide(h2);
  hratio->SetLineColor(kBlack);
  hratio->SetMarkerStyle(7);
  hratio->Draw("e1");
  hratio->SetTitle(";E_{#nu} (GeV);Ratio to NOvA");
  hratio->GetYaxis()->SetRangeUser(0.85, 1.15);
  TLine* flat = new TLine(0, 1, 20, 1);
  flat->SetLineStyle(7);
  flat->Draw();

  // Correlations & uncertainties
  f1->cd();
  const int nuni = 100;
  const int nbins = h->GetNbinsX();
  double* edges = new double[nbins+1];

  for (int i=1; i<nbins+1; i++) {
    edges[i-1] = h->GetBinLowEdge(i);
  }
  edges[nbins] = h->GetBinLowEdge(nbins-1) + h->GetBinWidth(nbins-1);

  // Covariance matrix
  TH2D* cov = new TH2D("cov", ";E_{#nu} (GeV);E_{#nu} (GeV)", nbins, edges, nbins, edges);
  for (int k=0; k<nuni; k++) {
    char name[500];
    if (mode == knumu) snprintf(name, 500, "numu/%s/Active_TPC_Volume/numu_%s_Uni_%i_AV_TPC",inputmode.c_str(),inputmode.c_str(), k); 
    if (mode == knue) snprintf(name, 500, "nue/%s/Active_TPC_Volume/nue_%s_Uni_%i_AV_TPC",inputmode.c_str(),inputmode.c_str(), k);    
    
    TH1D* hu = (TH1D*) f1->Get(name);
    for (int m=0; m<nbins; m++) {
      hu->SetBinContent(m, hu->GetBinContent(m) / hu->GetBinWidth(m));
    }

    for (int i=1; i<nbins+1; i++) {
      double cvi = horig->GetBinContent(i);

      double uvi = hu->GetBinContent(i);

      for (int j=1; j<nbins+1; j++) {
        double cvj = horig->GetBinContent(j);
        double uvj = hu->GetBinContent(j);

        double c = (uvi - cvi) * (uvj - cvj);

        // std::cout << c << "\t" << uvi << "\t" << cvi << "\t" << uvj << "\t" << cvj << "\t"<< std::endl;
        
        cov->SetBinContent(i, j, cov->GetBinContent(i, j) + c / nuni);
      }
    }
  }

  double cii{0}, cjj{0}, n{1}, temp{0}, horig_cont{0};

  // Correlation matrix & fractional errors
  TH2D* cor = (TH2D*) cov->Clone("cor");
  TH1D* herr = new TH1D("herr", ";E_{#nu};Fractional Uncertainty", nbins, edges);
  
  // loop over rows
  for (int i=1; i<nbins+1; i++) {
    cii = cov->GetBinContent(i, i);
    
    // Catch zeros
    if (horig->GetBinContent(i) <= 0) horig_cont = 0.5;
    else horig_cont = horig->GetBinContent(i);

    herr->SetBinContent(i, sqrt(cii) / horig_cont);
    
    // Loop over columns
    for (int j=1; j<nbins+1; j++) {
      cjj = cov->GetBinContent(j, j);
      n = sqrt(cii * cjj);
      
      // Catch Zeros
      if (n == 0) temp=0.5;
      else temp = cov->GetBinContent(i, j) / n;
      
      cor->SetBinContent(i, j, temp );
    }
  }

  // Plot correlations
  TCanvas* c3 = new TCanvas();
  cor->Draw("colz");
  gStyle->SetPalette(55); // kRainbow

  // Plot fractional errors overlaid with official NOvA plot
  TCanvas* c4 = new TCanvas();
  herr->SetLineColor(kRed);
  herr->SetLineWidth(2);
  herr->Draw("hist");

  f2->cd();
  TH1D* herr2 = (TH1D*) f2->Get("fractional_uncertainty_numu");
  herr2->SetLineColor(kBlue);
  herr2->SetLineWidth(2);
  herr2->Draw("same");

  TLegend* lfrac = new TLegend(0.2, 0.65, 0.4, 0.8);
  lfrac->AddEntry(herr2, "NOvA","l");
  lfrac->AddEntry(herr, "This work", "l");
  lfrac->Draw();

  // create plots folder if it does not exist
  gSystem->Exec("if [ ! -d \"plots\" ]; then echo \"\nPlots folder does not exist... creating\"; mkdir plots; fi"); 
  
  // Save the plots as pdfs in the plots folder
  c1->Print("plots/CV_Flux_Prediction.pdf");
  c2->Print("plots/Ratio_FLux_Prediction.pdf");
  c3->Print("plots/Correlation_Matrix.pdf");
  c4->Print("plots/Fractional_Uncertainties.pdf");
  std::cout << "\n"<< std::endl;
}
