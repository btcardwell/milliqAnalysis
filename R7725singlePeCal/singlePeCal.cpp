#include "TH1.h"
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TGraph.h"

#include <vector>
#include <string>
#include <stdlib.h>

const vector<string> fileNames_b = {"spectra/1300V_Blank.root"};

const vector<string> fileNames_t = {
                                    "spectra/1300V_06.8ns.root",
                                    "spectra/1300V_07.0ns.root",
                                    "spectra/1300V_07.2ns.root",
                                    "spectra/1300V_07.4ns.root",
                                    "spectra/1300V_07.6ns.root",
                                    "spectra/1300V_07.8ns.root",
                                    "spectra/1300V_08.0ns.root",
                                    "spectra/1300V_08.2ns.root"
                                   };


const vector<double> pulseWidths_b(fileNames_b.size(), 7.5); // ns
const vector<double> pulseWidths_t = {6.8, 7.0, 7.2, 7.4, 7.6, 7.8, 8.0, 8.2}; // ns
const vector<double> cut_b = {0.48}; // pC

class inputValues {
  public:
    double pw, // LED pulse width
           n,  // Number of entries
	       v,  // Distribution variance
	       e,  // Distribution mean
	       a;  // Number of entries below cut

    TH1D * h;  // Distribution histogram
};
 
class outputValues {
  public:
    double n0,     // Number of zero-pe triggers 
           eL,     // Average number of pe per trigger
	       ePsi,   // Mean of single-pe response
	       vPsi;   // Variance of single-pe response

    double v_eL,   // Variance in eL 
           v_ePsi; // Variance in ePsi
};

inputValues getInputValues(TString fileName, double pw, string cutString) {

    TFile * file = new TFile(fileName, "READ");
    TTree * tree = (TTree*)file->Get("Channel_1");
    tree->Draw("(-1)*Charge>>hist(300, -2, 100)"); 
    TH1D * hist = (TH1D*)gDirectory->Get("hist");
    
    inputValues inVal;
    
    inVal.pw = pw;
    inVal.n = hist->GetEntries();
    inVal.v = pow(hist->GetStdDev(), 2);
    inVal.e = hist->GetMean(); 
    inVal.a = tree->Draw("Charge", (TString)cutString);
    inVal.h = hist;
    
    return inVal;
}

outputValues getOutputValues(inputValues b, inputValues t) {

    outputValues outVal;

    // Calculate output values
    outVal.n0   = t.a*b.n/b.a;
    outVal.eL   = -log(outVal.n0/t.n);
    outVal.ePsi = (t.e - b.e)/outVal.eL;
    outVal.vPsi = (t.v - b.v)/outVal.eL - pow(outVal.ePsi, 2);

    // Calculate uncertainties in output values
    outVal.v_eL   = ( exp(outVal.eL) + 1 - 2*b.a/b.n ) / b.a;
    outVal.v_ePsi = (outVal.eL * (pow(outVal.ePsi, 2) + outVal.vPsi) + 2*b.v) / (b.n * pow(outVal.eL, 2)) + pow(outVal.ePsi, 2) * (outVal.v_eL) / pow(outVal.eL, 2);

    return outVal;
}

void makePlots(vector<inputValues> b, vector<inputValues> t, vector<outputValues> o, double cut) {

    // Plot output values vs pulse width
    TGraphErrors * eLPlot   = new TGraphErrors();
    TGraphErrors  * ePsiPlot = new TGraphErrors ();
    eLPlot->SetTitle("Average Number of LED-induced Photoelectrons per Trigger");
    ePsiPlot->SetTitle("Single-photoelectron Response Mean");

    eLPlot->SetMarkerColor(4);
    ePsiPlot->SetMarkerColor(4);
    eLPlot->SetMarkerStyle(21);
    ePsiPlot->SetMarkerStyle(21);

    for (int i = 0; i < b.size(); i++) {
        for (int j = 0; j < t.size(); j++) {
            eLPlot->SetPoint((i+1)*j, t[j].pw, o[(i+1)*j].eL);
            eLPlot->SetPointError((i+1)*j, 0.02, sqrt(o[(i+1)*j].v_eL));
            ePsiPlot->SetPoint((i+1)*j, t[j].pw, o[(i+1)*j].ePsi);
            ePsiPlot->SetPointError((i+1)*j, 0.02, sqrt(o[(i+1)*j].v_ePsi));
	    }
    }

    TCanvas * eL_c = new TCanvas("eL_c", "eL Canvas", 200, 10, 700, 500);
    eLPlot->GetXaxis()->SetTitle("LED Pulse Width / (ns)");
    eLPlot->GetYaxis()->SetTitle("Photoelectrons / Trigger");
    eLPlot->Draw("AP");

    TCanvas * ePsi_c = new TCanvas("ePsi_c", "ePsi Canvas", 200, 10, 700, 500);
    ePsiPlot->GetXaxis()->SetTitle("LED Pulse Width / (ns)");
    ePsiPlot->GetYaxis()->SetTitle("Mean / (pC)");
    ePsiPlot->GetYaxis()->SetRangeUser(0.3, 0.8);
    ePsiPlot->Draw("AP");
    TF1 * fit = new TF1("fit", "[0]", 0, 10);
    ePsiPlot->Fit(fit);

    // Plot spectra
    TCanvas * spectra_c = new TCanvas("spectra_c", "Spectra Canvas", 200, 10, 700, 500);
    spectra_c->SetLogy();
    for (int i = t.size() - 1; i >= 0; i--) {
        t[i].h->SetLineColor(4);
        t[i].h->SetLineWidth(2);
        t[i].h->Draw("SAME");
        t[i].h->SetStats(0);
        t[i].h->SetTitle("PMT Spectra of Pulsed LED w/ Varying Pulse Width");
        t[i].h->GetXaxis()->SetTitle("Integrated Charge / (pC)");
        t[i].h->GetYaxis()->SetTitle("Events");
        t[i].h->GetYaxis()->SetRangeUser(0.5, 20000);
    }
    b[0].h->SetLineColor(1);
    b[0].h->SetLineWidth(2);
    b[0].h->Draw("SAME");
    b[0].h->SetStats(0);
    TLine * cutLine_t = new TLine(-1*cut, 0.0, -1*cut, 20000.0);
    cutLine_t->SetLineWidth(2);
    cutLine_t->SetLineColor(1);
    cutLine_t->SetLineStyle(2);
    cutLine_t->Draw();
    TLine * meanLine_t = new TLine(0.66, 0.0, 0.66, 20000.0);
    meanLine_t->SetLineWidth(2);
    meanLine_t->SetLineColor(2);
    meanLine_t->SetLineStyle(2);
    meanLine_t->Draw();

    // Add legend
    //TLegend * leg = new TLegend(0.1, 0.7, 0.48, 0.9);
    //leg->SetHeader("Legend");
    //leg->AddEntry(t[0].h, "

    // Plot background
    //TCanvas * b_c = new TCanvas("b_c", "Background Canvas", 200, 10, 700, 500);
    //b_c->SetLogy();
    //b[0].h->Sumw2();
    //b[0].h->Draw();
    //b[0].h->GetYaxis()->SetRangeUser(1.0, 100000.0);
    //b[0].h->SetTitle("Background PMT Spectrum");
    //b[0].h->GetXaxis()->SetTitle("Integrated Charge / (pC)");
    //b[0].h->GetYaxis()->SetTitle("Events");
    //TLine * cutLine_b = new TLine(cut, 0.0, cut, 1000000.0);
    //cutLine_b->SetLineWidth(2);
    //cutLine_b->Draw();
}

void singlePeCal() {

    vector<inputValues> b; // Blank runs
    vector<inputValues> t; // Total runs
    vector<outputValues> o;
    vector<string> cutString(cut_b.size(), "Charge > ");
    
    // Create strings to perform cuts - I would be able to use to_string if C++11 worked 
    for (int i = 0; i < cut_b.size(); i++) {
        ostringstream ss;
        ss << cut_b[i];
        cutString[i].append(ss.str());
    }
    
    // Get input and output values from .root spectra
    for (int i = 0; i < fileNames_b.size(); i++) {
        b.push_back(getInputValues(fileNames_b[i], pulseWidths_b[i], cutString[i]));
        
        for (int j = 0; j < fileNames_t.size(); j++) {
                t.push_back(getInputValues(fileNames_t[j], pulseWidths_t[j], cutString[i]));
                o.push_back(getOutputValues(b[i], t[j]));
            
        }
    }
    makePlots(b, t, o, cut_b[0]);

}





