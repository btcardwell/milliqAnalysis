const double liveTimePerEvent = 1024. / 0.4e9;

TGraphErrors * integralHistogram(TTree * tree, TString branchName, double &rateAt10, double &rateErrorAt10) {
  const int n = 200;
  const double xmax = 20.;

  Double_t thresholds[n];
  Double_t thresholdErrors[n];
  Double_t rates[n];
  Double_t rateErrors[n];

  double totalLiveTime = liveTimePerEvent * tree->GetEntries();

  cout << endl << "Triggers: " << tree->GetEntries() << endl;
  cout << "Live time: " << totalLiveTime << " seconds" << endl;

  rateAt10 = tree->Draw("-1. * " + branchName, "-1. * " + branchName + " > 5.0", "goff");
  rateErrorAt10 = TMath::Sqrt(rateAt10) / totalLiveTime;
  rateAt10 /= totalLiveTime;

  for(int i = 0; i < n; i++) {
    thresholds[i] = i * xmax / n;
    thresholdErrors[i] = 0.;

    TString cut = "-1. * " + branchName + " > ";
    cut += Form("%f", thresholds[i]);

    rates[i] = tree->Draw("-1. * " + branchName, cut, "goff");
    rateErrors[i] = TMath::Sqrt(rates[i]) / totalLiveTime;
    rates[i] /= totalLiveTime;

  }

  TGraphErrors * gr = new TGraphErrors(n, thresholds, rates, thresholdErrors, rateErrors);
  return gr;
}

void darkRate() {
  TTree * tree = (TTree*)_file0->Get("data");

  double rateBare;
  double rateBare;

  double rateSG;
  double rateSG;

  TGraphErrors * gr = integralHistogram(tree, "min_2", rateBare, rateBare);
  cout << endl << endl << "Bare R7725 Rate above 5 mV = " << rateBare 
       << " +/- " << rateBare << endl << endl;

  TGraphErrors * gr = integralHistogram(tree, "min_4", rateSG, rateSG);
  cout << endl << endl << "Saint Gobain Rate above 5 mV = " << rateSG 
       << " +/- " << rateSG << endl << endl;
  
  gr->GetYaxis()->SetRangeUser(1, 8e6);
  gr->GetYaxis()->SetTitle("Rate [Hz]");
  gr->GetXaxis()->SetTitle("Threshold [mV]");
  gr->Draw("ALP");

  TCanvas * can = new TCanvas("can", "can", 10, 10, 800, 600);
  can->SetLogy(true);
  tree->Draw("-1.*min_2>>hist(200,0,20)");
}

void darkRate_fake() {

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TString dir = "C:/Users/milliqan/Documents/Visual Studio 2015/Projects/ConsoleApplication1/Release/sep13/";

  const unsigned int n = 7;
  TString fileNames[n] = {"DarkRate1100v.root", "DarkRate1200v.root",
                          "DarkRate1375v.root", "DarkRate1475v_v2.root",
                          "DarkRate1575v.root", "DarkRate1675v.root",
                          "DarkRate1775v.root"};

  int colors[10] = {kBlack, kRed, kBlue, kViolet, 8, kMagenta, kOrange, kCyan, kPink, kYellow+3};

  Double_t voltages[n] = {1100, 1200, 1375, 1475, 1575, 1675, 1775};
  Double_t voltageErrors[n] = {0, 0, 0, 0, 0, 0, 0};

  TFile * files[n];
  TTree * trees[n];
  TGraphErrors * graphs[n];

  Double_t rateAt10[n];
  Double_t rateErrorsAt10[n];

  for(unsigned int i = 0; i < n; i++) {
    files[i] = new TFile(dir + fileNames[i]);
    trees[i] = (TTree*)files[i]->Get("data");
    graphs[i] = integralHistogram(trees[i], "min_2", rateAt10[i], rateErrorsAt10[i]);
  }

  TGraphErrors * rateVersusVoltage = new TGraphErrors(n, voltages, rateAt10, voltageErrors, rateErrorsAt10);

  for(unsigned int i = 0; i < n; i++) {
    if(i == 0) {
      graphs[i]->GetYaxis()->SetRangeUser(1, 8e6);
      graphs[i]->GetYaxis()->SetTitle("Rate [Hz]");
      graphs[i]->GetXaxis()->SetTitle("Threshold [mV]");
    }

    graphs[i]->SetLineColor(colors[i]);

    if(i == 0) graphs[i]->Draw("ALP");
    else graphs[i]->Draw("LP");

  }

  TCanvas * can = new TCanvas("can", "can", 10, 10, 800, 600);
  can->SetLogy(true);

  for(unsigned int i = 0; i < n; i++) {
    trees[i]->SetLineColor(i+1);
    if(i == 0) trees[i]->Draw("-1.*min_2>>hist(625,0,1250)");
    else trees[i]->Draw("-1.*min_2", "", "same");
  }

  TCanvas * can2 = new TCanvas("can2", "can2", 10, 10, 800, 600);
  can->SetLogy(true);
  can->SetLogx(true);
  rateVersusVoltage->Draw("ALP");

}

  void LightvsDark() {

    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    TString dir = "C:/Users/milliqan/Documents/Visual Studio 2015/Projects/ConsoleApplication1/Release/work/sep21/";

    TFile * fLight = new TFile(dir + "v3_dark_lightsOn_1750v.root");
    TFile * fDark = new TFile(dir + "v3_dark_lightsOff_1750v.root");

    TTree * tLight = (TTree *)fLight->Get("data");
    TTree * tDark = (TTree *)fDark ->Get("data");

    Double_t rateLight;
    Double_t rateErrorLight;
    Double_t rateDark;
    Double_t rateErrorDark;

    TGraphErrors * graphLight = integralHistogram(tLight, "min_2", rateLight, rateErrorLight);
    TGraphErrors * graphDark = integralHistogram(tDark, "min_2", rateDark, rateErrorDark);

    graphLight->SetLineColor(kRed);
    graphDark->SetLineColor(kBlack);

    graphLight->GetYaxis()->SetRangeUser(1., 1.e8);

    TCanvas * can = new TCanvas("can", "can", 10, 10, 800, 600);
    can->SetLogy(true);

    graphLight->Draw("ALP");
    graphDark->Draw("LP");

    TCanvas * can2 = new TCanvas("can2", "can2", 10, 10, 800, 600);
    tLight->SetLineColor(kRed);

    tLight->Draw("-1.*min_2>>hist(200,0,500)");
    tDark->Draw("-1.*min_2", "", "same");


    cout << "For light: Rate Above 10 mV is " << rateLight << " +/- " << rateErrorLight << endl;
    cout << "For dark: Rate Above 10 mV is " << rateDark << " +/- " << rateErrorDark << endl;


  }
