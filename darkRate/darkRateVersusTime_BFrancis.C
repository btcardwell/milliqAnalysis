const double liveTimePerEvent = 1024. / 0.4e9;

TGraphErrors * integralHistogram(TTree * tree) {

  const int n = 200;
  const double xmax = 20.;

  Double_t thresholds[n];
  Double_t thresholdErrors[n];
  Double_t rates[n];
  Double_t rateErrors[n];

  double totalLiveTime = liveTimePerEvent * tree->GetEntries();

  cout << endl << "Triggers: " << tree->GetEntries() << endl;
  cout << "Live time: " << totalLiveTime << " seconds" << endl;

  // For bare R7725 specifically 5 mV
  double nAt5_bare = tree->Draw("-1. * min_2", "-1. * min_2 > 5.0", "goff");
  double rateAt5_bare = nAt5_bare / totalLiveTime;
  double rateErrorAt5_bare = TMath::Sqrt(nAt5_bare) / totalLiveTime;
  cout << endl << endl << "Overall bare R7725 rate above 5 mV = " << rateAt5_bare << " +/- " << rateErrorAt5_bare << endl;

  // For Saint Gobain module specifically 5 mV
  double nAt5_SG = tree->Draw("-1. * min_4", "-1. * min_4 > 5.0", "goff");
  double rateAt5_SG = nAt5_SG / totalLiveTime;
  double rateErrorAt5_SG = TMath::Sqrt(nAt5_SG) / totalLiveTime;
  cout << endl << endl << "Overall Saint Gobain rate above 5 mV = " << rateAt5_SG << " +/- " << rateErrorAt5_SG << endl << endl;

  for(int i = 0; i < n; i++) {
    thresholds[i] = i * xmax / n;
    thresholdErrors[i] = 0.;

    TString cut = "-1. * min_2 > ";
    cut += Form("%f", thresholds[i]);

    rates[i] = tree->Draw("-1. * min_2", cut, "goff");
    rateErrors[i] = TMath::Sqrt(rates[i]) / totalLiveTime;
    rates[i] /= totalLiveTime;

  }

  TGraphErrors * gr = new TGraphErrors(n, thresholds, rates, thresholdErrors, rateErrors);
  return gr;
}

void drawTimeSlice(TTree * tree, int sliceNumber,
                   unsigned int eventStart, unsigned int nEvents) {

   TString cut = "Entry$ >= ";
   cut += Form("%u", eventStart);

   cut += " && Entry$ < ";
   cut += Form("%u", eventStart + nEvents);

   tree->SetLineColor(sliceNumber + 1);

   if(sliceNumber == 0) tree->Draw("-1. * min_2 >> hist(200,0,400)", cut, "");
   else tree->Draw("-1. * min_2", cut, "same");

}

void calculateDarkRate(TTree * tree,
                       float threshold,
                       unsigned int eventStart, unsigned int nEvents,
                       double &rate, double &rateError, double &meanTime, double &meanTimeError,
                       int &nTDCRollovers) {

  double totalLiveTime = liveTimePerEvent * nEvents;

  float amplitude;
  ULong64_t tdc;
  tree->SetBranchAddress("min_4", &amplitude);
  tree->SetBranchAddress("TDC_4", &tdc);

  tree->GetEntry(0);
  ULong64_t initialTDC = tdc;

  if(eventStart == 0) {
    cout << endl << "initial TDC counter = " << initialTDC << " (" << initialTDC / 200.e6 / 60.0 << ") minutes" << endl << endl;
  }

  tree->GetEntry(eventStart);
  ULong64_t startTDC = tdc + nTDCRollovers*TMath::Power(2.0, 40.0);

  tree->GetEntry(eventStart + nEvents - 1);
  ULong64_t endTDC = tdc + nTDCRollovers*TMath::Power(2.0, 40.0);

  // The TDC clock is a 40-bit counter running at 200 MHz
  // So it resets every ~90 minutes
  if(endTDC <= startTDC) {
    endTDC += TMath::Power(2.0, 40.0);
    nTDCRollovers++;
  }

  float startTime = (startTDC - initialTDC) / 200.e6 / 60.0;
  float endTime = (endTDC - initialTDC) / 200.e6 / 60.0;

  meanTime = (endTime + startTime) / 2.0;
  meanTimeError = (endTime - startTime) / 2.0;

  rate = 0.;

  for(unsigned int i = eventStart; i < eventStart + nEvents; i++) {
    tree->GetEntry(i);
    if(-1. * amplitude > threshold) rate++;
  }

  cout << "rate = " << rate << " / " << totalLiveTime << endl;

  rateError = TMath::Sqrt(rate) / totalLiveTime;
  rate /= totalLiveTime;

  cout << "At time " << meanTime << " +/- " << meanTimeError << ", rate = "<< rate << " +/- " << rateError << " Hz" << endl;

  tree->ResetBranchAddresses();

  return;
}



void darkRateVersusTime() {

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  TTree * tree = (TTree*)_file0->Get("data");

  ///////////////////////////
  // Overall dark rate
  ///////////////////////////

  double rateAt5;
  double rateErrorAt5;

  TCanvas * canOverallRate = new TCanvas("overallRate", "Overall dark rate", 10, 10, 800, 600);
  canOverallRate->SetLogy(true);

  TGraphErrors * gr_overall = integralHistogram(tree);
  gr_overall->GetYaxis()->SetRangeUser(0.1, 8.e6);
  gr_overall->GetYaxis()->SetTitle("Rate above threshold [Hz]");
  gr_overall->GetXaxis()->SetTitle("Threshold [mV]");
  gr_overall->Draw("ALP");

  ///////////////////////////
  // Dark rate versus time
  ///////////////////////////

  const int nDivisions = 20;

  double rates[nDivisions];
  double rateErrors[nDivisions];
  double tdcs[nDivisions];
  double tdcErrors[nDivisions];

  unsigned int nEvents = tree->GetEntries() / nDivisions;

  int nTDCRollovers = 0;

  TCanvas * canTimeSlices = new TCanvas("timeSlices", "Amplitude Distributions in sub-runs", 10, 10, 800, 600);
  canOverallRate->SetLogy(true);

  for(int i = 0; i < nDivisions; i++) {
    calculateDarkRate(tree, 10.0, i*nEvents, nEvents, rates[i], rateErrors[i], tdcs[i], tdcErrors[i], nTDCRollovers);
    if(i < 5) drawTimeSlice(tree, i, i*nEvents, nEvents);
  }

  TCanvas * canVersusTime = new TCanvas("versusTime", "Dark rate versus time", 10, 10, 800, 600);

  TGraphErrors * gr = new TGraphErrors(nDivisions, tdcs, rates, tdcErrors, rateErrors);
  gr->GetXaxis()->SetTitle("Time in run [minutes]");
  gr->GetYaxis()->SetTitle("Event rate above 5 mV [Hz]");
  gr->Draw("ALP");

  TH1D * hist_ratesVersusTime = new TH1D("ratesVersusTime", "ratesVersusTime;Dark rate in sub-run [Hz];Runs", 500, 0, 200000);
  for(int i = 0; i < nDivisions; i++) {
    hist_ratesVersusTime->Fill(rates[i]);
  }

  TCanvas * canHist = new TCanvas("canHist", "Rates in sub-runs", 10, 10, 800, 600);
  hist_ratesVersusTime->Draw();

}
