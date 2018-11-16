/**
 * Draw a comparison between our flux and the standard NOvA one.
 *
 * A. Mastbaum <mastbaum@uchicago.edu>, 2018/11/12
 */

void plot_comparison(TString filename="output.root") {
  TFile* f1 = TFile::Open(filename);
  TDirectory* d = (TDirectory*)f1->Get("numu");
  d->cd();
  TH1D* h = (TH1D*) (gDirectory->Get("numu_CV_AV_TPC")->Clone("fx"));
  h->SetDirectory(0);

  // Normalize by bin width
  for (int i=0;i<h->GetNbinsX()+1;i++) {
    h->SetBinContent(i, h->GetBinContent(i)/h->GetBinWidth(i));
  }

  h->SetLineColor(kRed);
  h->SetTitle(";E_{#nu} (GeV);Fraction/GeV");
  h->Draw("");

  TFile* f2 = TFile::Open("nova_flux/FHC_Flux_NOvA_ND_2017.root");
  f2->ls();
  TH1D* h2 = (TH1D*) f2->Get("flux_numu");
  h2->SetLineColor(kBlue);

  h2->Draw("same");
  h->GetYaxis()->SetTitle(h2->GetYaxis()->GetTitle());

  // Shape-only comparison: scale to NOvA
  h->Scale(h2->Integral(2,-1)/h->Integral(2,-1));

  TLegend* l = new TLegend(0.5, 0.5, 0.88, 0.7);
  l->AddEntry(h2, "NOvA");
  l->AddEntry(h, "This");
  l->Draw();

  gPad->SetLogy();
  gPad->Update();
}

