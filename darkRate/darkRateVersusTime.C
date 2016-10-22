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
                       int &nTDCRollovers, float &meanTemp) {

  double totalLiveTime = liveTimePerEvent * nEvents;

  float amplitude;
  float temp = 0;
  ULong64_t tdc;
  tree->SetBranchAddress("min_2", &amplitude);
  tree->SetBranchAddress("TDC_2", &tdc);
  tree->SetBranchAddress("temperature", &temp);

  tree->GetEntry(0);
  ULong64_t initialTDC = tdc;

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
    //if(-1. * amplitude > threshold) rate++;
    meanTemp += temp;
  }

  rateError = TMath::Sqrt(rate) / totalLiveTime;
  rate /= totalLiveTime;

  meanTemp /= nEvents;

  tree->ResetBranchAddresses();

  return;
}

void darkRateVersusTime() {

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  vector<TTree*> trees;
  
  // Read all files
  const char * dirName = "/data/users/MilliQan/caen/PeriodicDarkRate/";
  const char * ext     = ".root";

  TSystemDirectory dir(dirName, dirName);
  TList * files = dir.GetListOfFiles();

  if (files) {
    cout << endl << "Reading the following files:" << endl << endl;
    TSystemFile * file;
    TString fName;
    TIter next(files);
    while ((file = (TSystemFile*)next())) {
      fName = file->GetName();
      if (!file->IsDirectory() && fName.EndsWith(ext)) {
        cout << fName.Data() << endl;
        TFile * tempFile = new TFile(dirName + fName);
        trees.push_back((TTree*)tempFile->Get("data"));
      }
    }
    cout << endl;
  }

  ///////////////////////////
  // Dark rate versus time
  ///////////////////////////

  const int nDivisions = 150;

  double rates[nDivisions][trees.size()];
  double rateErrors[nDivisions][trees.size()];
  double tdcs[nDivisions][trees.size()];
  double tdcErrors[nDivisions][trees.size()];
  int nTDCRollovers = 0;
  unsigned int nEvents[trees.size()];
  float temps;
  float meanTemps[nDivisions][trees.size()];

  for(int j = 0; j < trees.size(); j++) {
    nEvents[j] = trees[j]->GetEntries() / nDivisions;
    cout << trees[j]->GetEntries() << " " << nEvents[j] << endl;
   
    for(int i = 0; i < nDivisions; i++) {
      calculateDarkRate(trees[j], 5.0, i*nEvents[j], nEvents[j], rates[i][j], rateErrors[i][j],
                        tdcs[i][j], tdcErrors[i][j], nTDCRollovers, meanTemps[i][j]);
    }
  }

  TH1D * h_rates = new TH1D("ratesVersusTime", "rates;Dark rate in 6s sub-run [Hz];Sub-Runs", 50, 10000, 50000);
  TH2D * h_ratesVsEventNum = new TH2D("ratesVersusEventNum", "ratesVersusEventNum", nDivisions/2, 0, trees[0]->GetEntries(), 25, 10000, 50000);
  TH1D * h_temps = new TH1D("ratesVersusTemp", "ratesVersusTemp", 30, 0, 29); 
  TH2D * h_ratesVsTemp = new TH2D("ratesVersusTemp", "ratesVersusTemp", 2,
                                       21, 23, 25, 10000, 50000);
  TH2D * h_tempVsEventNum = new TH2D("TempVersusEventNum", "TempVersusEventNum", nDivisions/2,
                                       0, trees[0]->GetEntries(), 30, 0, 29);

  for (int j = 0; j < trees.size(); j++) {
    trees[j]->SetBranchAddress("temperature", &temps);
    for (int i = 0; i < nDivisions; i++) {
      if (nEvents[j] < 1360000) {
        //h_rates->Fill(rates[i][j]);
        h_temps->Fill(meanTemps[i][j]);
        //h_ratesVsEventNum->Fill(i*nEvents[j], rates[i][j]);
        //h_ratesVsTemp->Fill(meanTemps[i][j], rates[i][j]);
        h_tempVsEventNum->Fill(i*nEvents[j], meanTemps[i][j]);
      }
    }
  }

  TCanvas * c_rates = new TCanvas("c_rates", "Rates in sub-runs", 10, 10, 800, 600);
  h_rates->Draw();

  TCanvas * c_ratesVsEventNum = new TCanvas("c_ratesVsEventNum", "Rates Vs Event Number", 10, 10, 800, 600);
  h_ratesVsEventNum->Draw("COLZ");

  TCanvas * c_temps = new TCanvas("c_temps", "Temps in all events", 10, 10, 800, 600);
  h_temps->Draw();

  TCanvas * c_ratesVsTemp = new TCanvas("c_ratesVsTemp", "Rates Vs Temp", 10, 10, 800, 600);
  h_ratesVsTemp->Draw("COLZ");

  TCanvas * c_tempVsEventNum = new TCanvas("c_TempVsEventNum", "Temp Vs Event Number", 10, 10, 800, 600);
  h_tempVsEventNum->Draw("COLZ");
}
