const double liveTimePerEvent = 1024. / 3.2e9;
//const double pCPerPE = 0.56;
const double pCPerPE = 1.0; // durp

TGraphErrors * integralHistogram(TTree * tree) {

  const int n = 100;
  const double xmax = 30.;

  Double_t thresholds[n];
  Double_t thresholdErrors[n];
  Double_t rates[n];
  Double_t rateErrors[n];

  double totalLiveTime = liveTimePerEvent * tree->GetEntries();

  cout << "For " << tree->GetEntries() << " triggers, a live time of " 
                                       << totalLiveTime << " seconds" << endl;

  for(int i = 0; i < n; i++) {

    thresholds[i] = i * xmax / n / pCPerPE;
    thresholdErrors[i] = 0.;

    TString cut = "-1.*Min > ";
    cut += Form("%f", thresholds[i]);
    cut += "&& -1.*Min < 1200";

    rates[i] = tree->Draw("-1.*Min", cut, "goff");
    rateErrors[i] = TMath::Sqrt(rates[i]) / totalLiveTime;
    rates[i] /= totalLiveTime;

  }

  TGraphErrors * gr = new TGraphErrors(n, thresholds, rates, thresholdErrors, rateErrors);
  return gr;
}

void darkRate() {

  //gROOT->SetBatch(true);
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0000);
  gStyle->SetOptTitle(0);

  TFile * f1300 = new TFile("Spectra/08_12_16_1700V_LightsOff_Short.root", "READ");
  TFile * f1500 = new TFile("Spectra/08_12_16_1800V_LightsOff_Short.root", "READ");
  TFile * f1700 = new TFile("Spectra/08_12_16_1900V_LightsOff_Short.root", "READ");

  TTree * t1300 = (TTree*)f1300->Get("Channel_1");
  TTree * t1500 = (TTree*)f1500->Get("Channel_1");
  TTree * t1700 = (TTree*)f1700->Get("Channel_1");

  TGraphErrors * h1300 = integralHistogram(t1300);
  TGraphErrors * h1500 = integralHistogram(t1500);
  TGraphErrors * h1700 = integralHistogram(t1700);

  h1300->GetYaxis()->SetTitle("Rate [Hz]");
  h1300->GetYaxis()->SetRangeUser(0.5, 10e6);

  h1300->GetXaxis()->SetTitle("Threshold [mV]");
  h1300->GetXaxis()->SetRangeUser(0.5, 1000);

  h1300->SetLineColor(kBlack);
  h1500->SetLineColor(kRed);
  h1700->SetLineColor(kBlue);

  TLegend * leg = new TLegend(0.60, 0.70, 0.88, 0.88);
  leg->SetHeader("PMT Voltage");
  leg->AddEntry(h1300, "1700V", "l");
  leg->AddEntry(h1500, "1800V", "l");
  leg->AddEntry(h1700, "1900V", "l");

  h1300->Draw("ALP");
  h1500->Draw("LP");
  h1700->Draw("LP");
  leg->Draw();


}
