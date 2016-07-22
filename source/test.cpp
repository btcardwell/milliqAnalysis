



const vector<string> fileNames = {"spectra/07_12_16_10cm.root",
                                  "spectra/07_12_16_20cm.root",
                                  "spectra/07_12_16_30cm.root",
                                  "spectra/07_12_16_40cm.root",
                                  "spectra/07_12_16_50cm.root",
                                  "spectra/07_12_16_60cm.root",
                                /*"spectra/07_12_16_70cm.root",*/
                                  "spectra/07_12_16_80cm.root"};

const vector<int> x = {10, 20, 30, 40, 50, 60, /*70,*/ 80};


class Spectrum {
  public:
    int x;         // Source distance from PMT in cm
    double peak;   // Peak of Cesium Spectrum in pC
    TH1F * h;      // Histogram

    Spectrum(string,int);
};

Spectrum::Spectrum(string f, int d) {
    
    // Set distance
    x = d;

    // Get hist
    TFile * file = TFile::Open((TString)f);
    TTree * tree = (TTree*)file->Get("Channel_7");
    tree->Draw("Charge*(-1)>>h_temp(200, -10, 300)");
    h = (TH1F*)gDirectory->Get("h_temp");

    // Find peak
    TSpectrum * s = new TSpectrum(2);
    s->Search(h, .9);
    peak = s->GetPositionX()[1];
}

Double_t func(Double_t * x, Double_t * par) {
  Double_t xx = x[0];
  Double_t C = par[0];
  Double_t lambda = par[1];

  return C * TMath::Exp(-1. * xx / lambda);
}

void plotPeaks(vector<Spectrum> s) {
    
    TGraphErrors * plotAL = new TGraphErrors();
    plotAL->SetMarkerColor(4);
    plotAL->SetMarkerStyle(21);
 
    for (int i = 0; i < s.size(); i++) {
        plotAL->SetPoint(i, s[i].x, s[i].peak);
        plotAL->SetPointError(i, 1., s[i].peak*.1);
    }
 
    plotAL->Draw("AP");

    TF1 * myfunc = new TF1("myfunc", func, -10, 300, 2);
    plotAL->Fit("myfunc");
}


void test() {
    vector<Spectrum> spectra;

    for (int i = 0; i < fileNames.size(); i++) { 
        spectra.push_back(Spectrum(fileNames[i], x[i]));
    }

    plotPeaks(spectra);
}
