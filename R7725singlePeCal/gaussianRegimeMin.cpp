#include "TH1.h"
#include "TTree.h"
#include "TFile.h"
#include "TGraph.h"

#include <vector>
#include <stdlib.h>

const vector<string> fileNames = {"spectra/1700V_09.5ns.root",
                                  "spectra/1700V_10.0ns.root",
                                  "spectra/1700V_10.5ns.root",
                                  "spectra/1700V_11.0ns.root"};
 
const vector<double> pulseWidths = {9.5, 10.0, 10.5, 11.0}; 


class inputValues {
  public:
    double pw, // LED pulse width
    	   v,  // Distribution variance
    	   e;  // Distribution mean
    TH1D * h;  // Distribution histogram
};


inputValues getInputValues(TString fileName, double pw) {

    TFile * file = new TFile(fileName, "READ");
    TTree * tree = (TTree*)file->Get("Channel_1");
    tree->Draw("Min>>hist");
    TH1D * hist = (TH1D*)gDirectory->Get("hist");
    
    double h_mean = hist->GetMean();
    double h_sigma = hist->GetStdDev();
    hist->Fit("gaus","","",h_mean-5*h_sigma,h_mean+5*h_sigma);
    TF1 *gauss = hist->GetFunction("gaus");

    // Draw histogram to show fit
    hist->GetXaxis()->SetTitle("Amplitude / mV");
    hist->GetYaxis()->SetTitle("Events");
    hist->SetTitle("PMT Spectrum in Gaussian Regime");
    hist->Draw();
   
    inputValues inVal;
   
    // Get mean and variance of fit
    inVal.pw = pw;
    inVal.v = gauss->CentralMoment(2, -120, 100);
    inVal.e = -gauss->Moment(1, -120, 100); // Flip sign
    inVal.h = hist;

    return inVal;
}


void gaussianRegimeMin() {

    vector<inputValues> t;

    for (int i = 0; i < fileNames.size(); i++) {
        t.push_back(getInputValues(fileNames[i], pulseWidths[i]));
    }

    // Plot charge per photoelectron vs pulse width
    TGraph * qPerPEPlot = new TGraph();
    qPerPEPlot->SetTitle("Average Amplitude per Photoelectron in Gaussian Regime");
    qPerPEPlot->SetMarkerColor(4);
    qPerPEPlot->SetMarkerStyle(21);
    
    for (int i = 0; i < t.size(); i++) {
        qPerPEPlot->SetPoint(i, t[i].pw, t[i].v / t[i].e);
    }
    
    TCanvas * qPerPECanvas = new TCanvas("mVperPE", "mVperPE Canvas", 200, 10, 700, 500);
    qPerPEPlot->GetXaxis()->SetTitle("LED Pulse Width / (ns)");
    qPerPEPlot->GetYaxis()->SetTitle("mV");
    //qPerPEPlot->GetYaxis()->SetRangeUser(0.2, 1.0);
    qPerPEPlot->Draw("AP");
    TF1 * fit = new TF1("fit", "[0]", 0, 10);   
    qPerPEPlot->Fit(fit);
}
