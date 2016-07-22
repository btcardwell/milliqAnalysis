#include "TH1.h"
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TGraph.h"

#include <vector>
#include <string>
#include <stdlib.h>

const vector<string> fileNames_b = {"spectra/06_16_16_ChargeLength100_Blank.root",
                                    "spectra/06_16_16_ChargeLength150_Blank.root",
                                    "spectra/06_16_16_ChargeLength200_Blank.root",
                                    "spectra/06_16_16_ChargeLength250_Blank.root"};

const vector<string> fileNames_t = {"spectra/06_16_16_ChargeLength100_7.2ns.root",
                                    "spectra/06_16_16_ChargeLength100_7.5ns.root",
                                    "spectra/06_16_16_ChargeLength100_8.0ns.root",
                                    "spectra/06_16_16_ChargeLength150_7.2ns.root",
                                    "spectra/06_16_16_ChargeLength150_7.5ns.root",
                                    "spectra/06_16_16_ChargeLength150_8.0ns.root",
                                    "spectra/06_16_16_ChargeLength200_7.2ns.root",
                                    "spectra/06_16_16_ChargeLength200_7.5ns.root",
                                    "spectra/06_16_16_ChargeLength200_8.0ns.root",
                                    "spectra/06_16_16_ChargeLength250_7.2ns.root",
                                    "spectra/06_16_16_ChargeLength250_7.5ns.root",
                                    "spectra/06_16_16_ChargeLength250_8.0ns.root"};

const vector<double> pulseWidths_b(fileNames_b.size(), 7.5); // ns
const vector<double> pulseWidths_t = {7.2, 7.5, 8.0, 7.2, 7.5, 8.0, 7.2, 7.5, 8.0, 7.2, 7.5, 8.0}; // ns
const vector<double> cl_b = {31.250, 46.875, 62.500, 78.125}; //ns
const vector<double> cl_t = {31.250, 31.250, 31.250, 46.875, 46.875, 46.875, 62.500, 62.500, 62.500, 78.125, 78.125, 78.125}; //ns
const vector<double> cut_b = {-0.008, -0.016, -0.025, -0.033}; //33% blank below cut

class inputValues {
    public:
    double pw, // LED pulse width
           n,  // Number of entries
    	   v,  // Distribution variance
    	   e,  // Distribution mean
    	   a,  // Number of entries below cut
           cl; // Charge length in ns
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

inputValues getInputValues(TString fileName, double pw, string cutString, int cl) {

    TFile * file = new TFile(fileName, "READ");
    TTree * tree = (TTree*)file->Get("Channel_8");
    tree->Draw("Charge>>hist"); 
    TH1D * hist = (TH1D*)gDirectory->Get("hist");
    
    inputValues inVal;

    inVal.pw = pw;
    inVal.cl = cl;
    inVal.n = hist->GetEntries();
    inVal.v = pow(hist->GetStdDev(), 2);
    inVal.e = -hist->GetMean(); // Flip sign
    inVal.a = tree->Draw("Charge", (TString)cutString);
    inVal.h = hist;

    // Debugging
    //cout << inVal.pw << " " << inVal.cl << " " << inVal.a << endl;

    return inVal;
}

outputValues getOutputValues(inputValues b, inputValues t) {

    outputValues outVal;

    // Calculate output values
    outVal.n0   = t.a*b.n/b.a;
    outVal.eL   = -log(outVal.n0/t.n);
    outVal.ePsi = (t.e - b.e)/outVal.eL;
    outVal.vPsi = (t.v - b.v)/outVal.eL - pow(outVal.ePsi, 2);

    // Debugging
    //cout << t.cl << " " << t.pw << " " << outVal.n0 << " " << outVal.eL << " " << outVal.ePsi<< endl;

    // Calculate uncertainties in output values
    outVal.v_eL   = ( exp(outVal.eL) + 1 - 2*b.a/b.n ) / b.a;
    outVal.v_ePsi = (outVal.eL * (pow(outVal.ePsi, 2) + outVal.vPsi) + 2*b.v) / (b.n * pow(outVal.eL, 2)) + pow(outVal.ePsi, 2) * (outVal.v_ePsi) / pow(outVal.eL, 2);

    // Debugging
    // cout << sqrt(outVal.v_eL) << " " << sqrt(outVal.v_ePsi) << endl;

    return outVal;
}

void makePlots(vector<inputValues> b, vector<inputValues> t, vector<outputValues> o) {

    // Plot output values vs pulse width
    TMultiGraph * eLPlot   = new TMultiGraph();
    TMultiGraph * ePsiPlot = new TMultiGraph();
    TMultiGraph * vPsiPlot = new TMultiGraph();

    vector<TGraph*> eLPlots;
    vector<TGraph*> ePsiPlots;
    vector<TGraph*> vPsiPlots;

    int pointNum = 0;
    for (int i = 0; i < b.size(); i++) {
        eLPlots.push_back(new TGraph());
        ePsiPlots.push_back(new TGraph());
        vPsiPlots.push_back(new TGraph());
        eLPlots[i]->SetMarkerStyle(21);
        ePsiPlots[i]->SetMarkerStyle(21);
        vPsiPlots[i]->SetMarkerStyle(21);
        eLPlots[i]->SetMarkerColor(i+1);
        ePsiPlots[i]->SetMarkerColor(i+1);
        vPsiPlots[i]->SetMarkerColor(i+1);
        //eLPlots[i]->SetLineColor(i+1);
        //ePsiPlots[i]->SetLineColor(i+1);
        //vPsiPlots[i]->SetLineColor(i+1);

        for (int j = 0; j < t.size(); j++) {
            if (t[j].cl == b[i].cl) {
                eLPlots[i]->SetPoint(pointNum, t[j].pw, o[j].eL);
                ePsiPlots[i]->SetPoint(pointNum, t[j].pw, o[j].ePsi);
                vPsiPlots[i]->SetPoint(pointNum, t[j].pw, o[j].vPsi);
                
                pointNum++;

                // Debugging
                // cout << pointNum << " " << " " << t[j].pw << " " << o[j].eL << " " << o[j].ePsi << " " << o[j].vPsi << endl;
            }
	    }
        eLPlot->Add(eLPlots[i]);
        ePsiPlot->Add(ePsiPlots[i]);
        vPsiPlot->Add(vPsiPlots[i]);
    }

    // Make pe / trig plot
    TCanvas * eL_c = new TCanvas("eL_c", "eL Canvas", 200, 10, 700, 500);
    eLPlot->Draw("AP");
    eLPlot->SetTitle("Average Number of LED-induced Photoelectrons per Trigger");
    eLPlot->GetXaxis()->SetTitle("LED Pulse Width / (ns)");
    eLPlot->GetYaxis()->SetTitle("Photoelectrons / Trigger");
    eLPlot->GetXaxis()->SetRangeUser(7.0, 8.2);
    eLPlot->GetYaxis()->SetRangeUser(0.0, 6.0);
    eL_leg = new TLegend(0.15, 0.65, 0.5, 0.85);
    eL_leg->AddEntry(eLPlots[0], "Integration window = 31.250 ns", "P");
    eL_leg->AddEntry(eLPlots[1], "Integration window = 46.875 ns", "P");
    eL_leg->AddEntry(eLPlots[2], "Integration window = 62.500 ns", "P");
    eL_leg->AddEntry(eLPlots[3], "Integration window = 78.125 ns", "P");
    eL_leg->Draw();
    gPad->Modified();

    // Make mean plot
    TCanvas * ePsi_c = new TCanvas("ePsi_c", "ePsi Canvas", 200, 10, 700, 500);
    ePsiPlot->Draw("AP");
    ePsiPlot->SetTitle("Single-photoelectron Response Mean in Few-PE Regime");
    ePsiPlot->GetXaxis()->SetTitle("LED Pulse Width / (ns)");
    ePsiPlot->GetYaxis()->SetTitle("pC");
    ePsiPlot->GetXaxis()->SetRangeUser(7.0, 8.2);
    ePsiPlot->GetYaxis()->SetRangeUser(0.0, 1.0);
    ePsi_leg = new TLegend(0.50, 0.65, 0.85, 0.85);
    ePsi_leg->AddEntry(eLPlots[0], "Integration window = 31.250 ns", "P");
    ePsi_leg->AddEntry(eLPlots[1], "Integration window = 46.875 ns", "P");
    ePsi_leg->AddEntry(eLPlots[2], "Integration window = 62.500 ns", "P");
    ePsi_leg->AddEntry(eLPlots[3], "Integration window = 78.125 ns", "P");
    ePsi_leg->Draw();
    gPad->Modified();

    // Make variance plot
    //TCanvas * vPsi_c = new TCanvas("vPsi_c", "vPsi Canvas", 200, 10, 700, 500);
    //vPsiPlot->Draw("AP");
    //vPsiPlot->SetTitle("Single-photoelectron Response Variance");
    //vPsiPlot->GetXaxis()->SetTitle("LED Pulse Width / (ns)");
    //vPsiPlot->GetYaxis()->SetTitle("pC^2");
    //vPsiPlot->GetXaxis()->SetRangeUser(7.2, 8.4);
    //vPsi_leg = new TLegend(0.51, 0.65, 0.85, 0.85);
    //vPsi_leg->AddEntry(eLPlots[0], "Integration window = 31.250 ns", "P");
    //vPsi_leg->AddEntry(eLPlots[1], "Integration window = 46.875 ns", "P");
    //vPsi_leg->AddEntry(eLPlots[2], "Integration window = 62.500 ns", "P");
    //vPsi_leg->AddEntry(eLPlots[3], "Integration window = 78.125 ns", "P");
    //vPsi_leg->Draw();
    //gPad->Modified();

}

void singlePeVaryChargeLength() {

    vector<inputValues> b; // Blank runs
    vector<inputValues> t; // Total runs
    vector<outputValues> o;
    vector<string> cutString(cut_b.size(), "Charge > ");
    
    // Create strings to perform cuts
    for (int i = 0; i < cut_b.size(); i++) {
        ostringstream ss;
        ss << cut_b[i];
        cutString[i].append(ss.str());
    }
    
    // Get input and output values from .root spectra
    for (int i = 0; i < fileNames_b.size(); i++) {
        b.push_back(getInputValues(fileNames_b[i], pulseWidths_b[i], cutString[i], cl_b[i]));
        
        for (int j = 0; j < fileNames_t.size(); j++) {
            if (cl_t[j] == cl_b[i]) {
                t.push_back(getInputValues(fileNames_t[j], pulseWidths_t[j], cutString[i], cl_t[j]));
                o.push_back(getOutputValues(b[i], t[j]));
            }
        }
    }
    makePlots(b, t, o);
}





