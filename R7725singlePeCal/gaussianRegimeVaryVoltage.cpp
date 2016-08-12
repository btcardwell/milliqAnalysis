#include "TH1.h"
#include "TTree.h"
#include "TFile.h"
#include "TGraph.h"

#include <vector>
#include <stdlib.h>

const vector<string> fileNames = {"spectra/06_14_16_VaryVoltage_1250V_11.0ns.root",
                                  "spectra/06_14_16_VaryVoltage_1300V_11.0ns.root",
                                  "spectra/06_14_16_VaryVoltage_1350V_11.0ns.root",
                                  "spectra/06_14_16_VaryVoltage_1400V_11.0ns.root"};
 
const vector<double> voltages = {1250, 1300, 1350, 1400};
const vector<double> pulseWidths(fileNames.size(), 11.0); 


class inputValues {
  public:
    double pw, // LED pulse width
           hv, // PMT voltage
    	   v,  // Distribution variance
    	   e;  // Distribution mean
    TH1D * h;  // Distribution histogram
};


inputValues getInputValues(TString fileName, double pw, double hv) {

    TFile * file = new TFile(fileName, "READ");
    TTree * tree = (TTree*)file->Get("Channel_8");
    tree->Draw("Charge>>hist");
    // tree->Draw("Charge>>hist", "Charge < 0"); 
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
    inVal.hv = hv;
    inVal.v = gauss->CentralMoment(2, -120, 100);
    inVal.e = -gauss->Moment(1, -120, 100); // Flip sign
    inVal.h = hist;

    return inVal;
}


void gaussianRegimeVaryVoltage() {

    vector<inputValues> t;

    for (int i = 0; i < fileNames.size(); i++) {
        t.push_back(getInputValues(fileNames[i], pulseWidths[i], voltages[i]));
    }

    // Plot charge per photoelectron vs PMT voltage
    TGraph * qPerPEPlot = new TGraph();
    qPerPEPlot->SetTitle("Single-photoelectron Response Mean in Few-PE Regime");
    qPerPEPlot->SetMarkerColor(4);
    qPerPEPlot->SetMarkerStyle(21);
    
    for (int i = 0; i < t.size(); i++) {
        qPerPEPlot->SetPoint(i, t[i].hv, t[i].v / t[i].e);
    }
    
    TCanvas * qPerPECanvas = new TCanvas("qPerPE", "qPerPE Canvas", 200, 10, 700, 500);
    qPerPEPlot->GetXaxis()->SetTitle("PMT Voltage / (V)");
    qPerPEPlot->GetYaxis()->SetTitle("Average Charge per PE / (pC)");
    qPerPEPlot->GetYaxis()->SetRangeUser(0.2, 1.0);
    qPerPEPlot->Draw("AP");

    // Debugging
    for (int i = 0; i < t.size(); i++) {
        cout << t[i].pw << "  " << t[i].e << "  " << t[i].v << "  " << t[i].v / t[i].e << endl;
    }
}
