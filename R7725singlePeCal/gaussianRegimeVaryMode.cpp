#include "TH1.h"
#include "TTree.h"
#include "TFile.h"
#include "TGraph.h"

#include <vector>
#include <stdlib.h>

const vector<string> fileNames = {"spectra/06_28_16_NormalMode_9.0ns.root",
                                  "spectra/06_28_16_NormalMode_9.5ns.root",
                                  "spectra/06_28_16_NormalMode_10.0ns.root",
                                  "spectra/06_28_16_NormalMode_10.5ns.root",
                                  "spectra/06_28_16_ChargeMode_9.0ns.root",
                                  "spectra/06_28_16_ChargeMode_9.5ns.root",
                                  "spectra/06_28_16_ChargeMode_10.0ns.root",
                                  "spectra/06_28_16_ChargeMode_10.5ns.root"};
 
const vector<double> pulseWidths = {9.0, 9.5, 10.0, 10.5, 9.0, 9.5, 10.0, 10.5}; 
const vector<bool> chargeMode = {false, false, false, false, true, true, true, true};

class inputValues {
  public:
    double pw, // LED pulse width
    	   v,  // Distribution variance
    	   e;  // Distribution mean
    bool   cm;
    TH1D * h;  // Distribution histogram
};


inputValues getInputValues(TString fileName, double pw, bool cm) {

    TFile * file = new TFile(fileName, "READ");
    TTree * tree = (TTree*)file->Get("Channel_8");
    tree->Draw("Charge>>hist");
    TH1D * hist = (TH1D*)gDirectory->Get("hist");
    
    double h_mean = hist->GetMean();
    double h_sigma = hist->GetStdDev();
    hist->Fit("gaus","","",h_mean-5*h_sigma,h_mean+5*h_sigma);
    TF1 *gauss = hist->GetFunction("gaus");

    // Draw histogram to show fit
    hist->GetXaxis()->SetTitle("Integrated Charge / pC");
    hist->GetYaxis()->SetTitle("Events");
    hist->SetTitle("PMT Spectrum in Gaussian Regime");
    hist->Draw();
   
    inputValues inVal;
   
    // Get mean and variance of fit
    inVal.pw = pw;
    inVal.cm = cm;
    inVal.v = gauss->CentralMoment(2, -120, 100);
    inVal.e = -gauss->Moment(1, -120, 100); // Flip sign
    inVal.h = hist;

    return inVal;
}


void gaussianRegimeVaryMode() {

    vector<inputValues> t;

    for (int i = 0; i < fileNames.size(); i++) {
        t.push_back(getInputValues(fileNames[i], pulseWidths[i], chargeMode[i]));
    }

    // Plot charge per photoelectron vs pulse width
    TMultiGraph * qPerPEPlot = new TMultiGraph();
    vector<TGraph*> qPerPEPlots;

    for (int i = 0; i < 2; i++) {
        qPerPEPlots.push_back(new TGraph());
        qPerPEPlots[i]->SetMarkerStyle(21);
        qPerPEPlots[i]->SetMarkerColor(1+i);
    
        for (int j = 0; j < t.size() / 2; j++) {
            int t_ix = j + i*(t.size()/2);
            qPerPEPlots[i]->SetPoint(j, t[t_ix].pw, t[t_ix].v / t[t_ix].e);
        }
        qPerPEPlot->Add(qPerPEPlots[i]);
    }
    
    TCanvas * qPerPECanvas = new TCanvas("qPerPE", "qPerPE Canvas", 200, 10, 700, 500);
    qPerPEPlot->Draw("AP");
    qPerPEPlot->SetTitle("Average Charge per Photoelectron in Gaussian Regime");
    qPerPEPlot->GetXaxis()->SetTitle("LED Pulse Width / (ns)");
    qPerPEPlot->GetYaxis()->SetTitle("pC");
    qPerPEPlot->GetYaxis()->SetRangeUser(0.4, 1.0);
    legend = new TLegend(0.55, 0.70, 0.87, 0.87);
    legend->SetHeader("   Integration Window: 78.125ns");
    legend->AddEntry(qPerPEPlots[0], "Normal Mode", "P");
    legend->AddEntry(qPerPEPlots[1], "Charge Mode", "P");
    legend->Draw();
    gPad->Modified();

    // Debugging
    //for (int i = 0; i < t.size(); i++) {
    //    cout << t[i].pw << "  " << t[i].e << "  " << t[i].v << "  " << t[i].v / t[i].e << endl;
    //}
}
