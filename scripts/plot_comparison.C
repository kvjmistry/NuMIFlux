/**
 * NuMI Validation Plotting
 *
 * A. Mastbaum <mastbaum@uchicago.edu> 2018/11
 */

void draw(TString inputfile="output.root") {
  gStyle->SetOptStat(0);

  // Overlay of output plot with official NOvA FHC numu flux
  TFile* f1 = TFile::Open(inputfile);
  TDirectory* d = (TDirectory*)f1->Get("numu");
  d->cd();
  TH1D* h = (TH1D*) (gDirectory->Get("numu_CV_AV_TPC")->Clone("fx"));
  h->SetDirectory(0);
  for (int i=0;i<h->GetNbinsX()+1;i++) {
    h->SetBinContent(i, h->GetBinContent(i)/h->GetBinWidth(i));
  }
  h->Sumw2();
  h->SetLineColor(kRed);
  h->SetTitle(";E_{#nu} (GeV);Fraction/GeV");
  h->Draw("");
  TH1D* horig = (TH1D*) h->Clone("horig");

  TFile* f2 = TFile::Open("FHC_Flux_NOvA_ND_2017.root");
  f2->ls();
  TH1D* h2 = (TH1D*) f2->Get("flux_numu");
  h2->SetLineColor(kBlue);

  h2->Draw("same");
  h->GetYaxis()->SetTitle(h2->GetYaxis()->GetTitle());
  h->Scale(h2->Integral(2,-1)/h->Integral(2,-1));
  TH1D* hratio = (TH1D*) h->Clone("hratio");
  hratio->SetDirectory(0);

  TLegend* l = new TLegend(0.5, 0.5, 0.88, 0.7);
  l->AddEntry(h2, "NOvA");
  l->AddEntry(h, "This work");
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
    snprintf(name, 500, "numu/ms_PPFX/Active_TPC_Volume/numu_ms_PPFX_Uni_%i_AV_TPC", k);
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
        cov->SetBinContent(i, j, cov->GetBinContent(i, j) + c / nuni);
      }
    }
  }

  // Correlation matrix & fractional errors
  TH2D* cor = (TH2D*) cov->Clone("cor");
  TH1D* herr = new TH1D("herr", ";E_{#nu};Fractional Uncertainty", nbins, edges);
  for (int i=1; i<nbins+1; i++) {
    double cii = cov->GetBinContent(i, i);
    herr->SetBinContent(i, sqrt(cii) / horig->GetBinContent(i));
    for (int j=1; j<nbins+1; j++) {
      double cjj = cov->GetBinContent(j, j);
      double n = sqrt(cii * cjj);
      cor->SetBinContent(i, j, cov->GetBinContent(i, j) / n);
    }
  }

  // Plot correlations
  TCanvas* c3 = new TCanvas();
  cor->Draw("colz");

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
}

