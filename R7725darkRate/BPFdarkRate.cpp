const double liveTimePerEvent = 1024. / 3.2e9;
//const double pCPerPE = 0.56;
const double pCPerPE = 1.0; // durp

TGraphErrors * integralHistogram(TTree * tree) {

  const int n = 125;
  const double xmax = 25.;

  Double_t thresholds[n];
  Double_t thresholdErrors[n];
  Double_t rates[n];
  Double_t rateErrors[n];

  double totalLiveTime = liveTimePerEvent * tree->GetEntries();

  cout << "For " << tree->GetEntries() << " triggers, a live time of " << totalLiveTime << " seconds" << endl;

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

  TFile * fLights = new TFile("/data/users/bfrancis/ExtTrigger_SG_1375.root", "READ");
  TFile * fNoLights = new TFile("/data/users/bfrancis/ExtTrigger_SG_1375_lightsoff.root", "READ");

  TFile * f1250 = new TFile("/data/users/MilliQan/caen/DarkRate/ExtTrigger_1250.root", "READ");
  TFile * f1325 = new TFile("/data/users/MilliQan/caen/DarkRate/ExtTrigger_1325.root", "READ");
  TFile * f1400 = new TFile("/data/users/MilliQan/caen/DarkRate/ExtTrigger_1400.root", "READ");

  TTree * tLights = (TTree*)fLights->Get("Channel_6");
  TTree * tNoLights = (TTree*)fNoLights->Get("Channel_6");

  TTree * t1250 = (TTree*)f1250->Get("Channel_13");
  TTree * t1325 = (TTree*)f1325->Get("Channel_13");
  TTree * t1400 = (TTree*)f1400->Get("Channel_13");

  TGraphErrors * hLights = integralHistogram(tLights);
  TGraphErrors * hNoLights = integralHistogram(tNoLights);

  TGraphErrors * h1250 = integralHistogram(t1250);
  TGraphErrors * h1325 = integralHistogram(t1325);
  TGraphErrors * h1400 = integralHistogram(t1400);

  hLights->GetYaxis()->SetTitle("Rate [Hz]");
  hLights->GetYaxis()->SetRangeUser(93-1, 6e5);

  hLights->GetXaxis()->SetTitle("Threshold [mV]");
  hLights->GetXaxis()->SetRangeUser(0, 1000);

  hLights->SetLineColor(kRed);
  hLights->SetLineColor(kBlack);
  hNoLights->SetLineColor(kBlack);

  hLights->Draw("ALP");
  //hNoLights->Draw("LP");

  h1250->GetYaxis()->SetTitle("Rate [Hz]");
  h1250->GetYaxis()->SetRangeUser(9e-2, 6e6);

  h1250->GetXaxis()->SetTitle("Threshold [mV]");
  h1250->GetXaxis()->SetRangeUser(0.5, 50);

  h1250->SetLineColor(kBlack);
  h1325->SetLineColor(kRed);
  h1400->SetLineColor(kBlue);

  //h1250->Draw("ALP");
  //h1325->Draw("LP");
  //h1400->Draw("LP");


}
