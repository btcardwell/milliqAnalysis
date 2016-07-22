

const string bg_fileName = "spectra/07_12_16_OutBox_NoSource.root";

const vector<string> fileNames = {"spectra/07_12_16_10cm.root",  
                                  "spectra/07_12_16_20cm.root",
                                  "spectra/07_12_16_30cm.root",
                                  "spectra/07_12_16_40cm.root",
                                  "spectra/07_12_16_50cm.root",
                                  "spectra/07_12_16_60cm.root",
                                  "spectra/07_12_16_70cm.root",
                                  "spectra/07_12_16_80cm.root"};

const vector<int> x = {10, 20, 30, 40, 50, 60, 70, 80};


class Spectrum {
  public:
    int x;         // Source distance from PMT in cm
    double peak;   // Peak of Cesium Spectrum in pC
    TH1F * h;      // Histogram
    Spectrum * bg; // Associated background

    Spectrum(string);
    Spectrum(string,int,Spectrum&);
    void Fit();
};

Spectrum::Spectrum(string f) {
    
    // Get hist
    TFile * file = TFile::Open((TString)f);
    TTree * tree = (TTree*)file->Get("Channel_7");
    tree->Draw("Charge*(-1)>>h_temp(200, -10, 300)");
    h = (TH1F*)gDirectory->Get("h_temp");
}

Spectrum::Spectrum(string f, int d, Spectrum &b) {
    
    // Set distance
    x = d;

    // Set Background
    bg = &b;

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

void Spectrum::Fit() {

    RooRealVar q("q", "q", -10, 300); 
    RooRealVar w("w", "weight", 0.5, 0., 1.);

    // Make Voigtian peak
    RooRealVar  mean("mean", "mean", 150, -10, 300);
    RooRealVar  width("width", "width", 1, 0., 100.);
    RooRealVar  sigma("sigma", "sigma", 1, 0., 10.);
    RooVoigtian peak("peak", "peak", q, mean, width, sigma);

    // Make background pdf from histogram
    RooDataHist bg_h("bg", "bg", q, bg->h);
    RooHistPdf  bg_pdf("bg", "bg", q, bg_h, 2);

    // Fit to signal
    RooDataHist sig_h("sig", "sig", q, h);
    RooAddPdf   model("model", "model", RooArgList(peak, bg_pdf), w);
    model.fitTo(sig_h);

    RooPlot * frame = q.frame();
    sig_h.plotOn(frame);
    model.plotOn(frame);
    
    frame->Draw();

}



void plotPeaks(vector<Spectrum> s) {
    
    // Plot Peak Charge vs Source Distance
    TGraphErrors* plotAL = new TGraphErrors();
    plotAL->SetMarkerColor(4);
    plotAL->SetMarkerStyle(21);
 
    for (int i = 0; i < s.size(); i++) {
        plotAL->SetPoint(i, s[i].x, s[i].peak);
        plotAL->SetPointError(i, 1., s[i].peak*.1);
    }
 
    plotAL->Draw("AP");
    plotAL->SetTitle("Location of Cesium Peak as a Function of Source Placement");
    plotAL->GetXaxis()->SetTitle("Horizontal Distance Between Source and PMT / (cm)");
    plotAL->GetYaxis()->SetTitle("Cesium Peak / (pC)");

    // Fit with exponential
    TF1 * exp = new TF1("exp", "[0]*exp(-1*x/[1])", 0, 80);
    exp->SetParameters(160, 210);
    exp->SetParName(0, "Normalization / (pC)");
    exp->SetParName(1, "Attenuation Length / (cm)");
    plotAL->Fit("exp");

    TLegend * leg = new TLegend(0.60, 0.70, 0.88, 0.88);
    leg->AddEntry(plotAL, "Data", "P");
    leg->AddEntry(exp, "Exponential Fit", "L");
    leg->Draw();
}


void plotSpectra(vector<Spectrum> s) {
    
    TCanvas * c = new TCanvas("c", "c", 600, 400); 

    for(int i = 0; i < s.size(); i+=2) {
        s[i].h->SetLineColor(i+2);
        s[i].h->Draw("SAME"); 
    }

    s[0].h->SetTitle("Cesium Spectra with Marked Peaks while Varying Source Distance from PMT");
    s[0].h->GetXaxis()->SetTitle("Integrated Charge / (pC)");
    s[0].h->GetYaxis()->SetTitle("Events");
    gStyle->SetOptStat(0);
    c->Update();

    TLegend * leg = new TLegend(0.60, 0.70, 0.88, 0.88);
    leg->AddEntry(s[0].h, "10cm");
    leg->AddEntry(s[2].h, "30cm");
    leg->AddEntry(s[4].h, "50cm");
    leg->AddEntry(s[6].h, "70cm");
    leg->Draw();

}


void attenuationLength() {

    Spectrum bg(bg_fileName);
    vector<Spectrum> spectra;

    for (int i = 0; i < fileNames.size(); i++) { 
        if (i == 0) { 
            spectra.push_back(Spectrum(fileNames[i], x[i], bg));
            spectra[i].Fit();
        }
    }

    //plotPeaks(spectra);
    //plotSpectra(spectra);

}
