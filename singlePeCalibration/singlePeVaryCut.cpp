#include "TH1.h"
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TGraph.h"

#include <vector>
#include <string>
#include <stdlib.h>

const vector<string> fileNames_b = {"spectra/06_28_16_NormalMode_ChargeLength500_Blank.root"};

const vector<string> fileNames_t = {"spectra/06_28_16_NormalMode_ChargeLength500_6.8ns.root",
                                    "spectra/06_28_16_NormalMode_ChargeLength500_7.0ns.root",
                                    //"spectra/06_28_16_NormalMode_ChargeLength500_7.2ns.root",
                                    "spectra/06_28_16_NormalMode_ChargeLength500_7.4ns.root",
                                    //"spectra/06_28_16_NormalMode_ChargeLength500_6ns.root",
                                    "spectra/06_28_16_NormalMode_ChargeLength500_7.8ns.root",
                                    //"spectra/06_28_16_NormalMode_ChargeLength500_8.0ns.root",
                                    "spectra/06_28_16_NormalMode_ChargeLength500_8.2ns.root"};

const vector<double> pulseWidths_b(fileNames_b.size(), 7.5); // ns
const vector<double> pulseWidths_t = {6.8, 7.0, 7.4, 7.8, 8.2}; // ns
const vector<double> cuts = {0.2, 0.1, -0.01, -0.2}; // pC

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
    TTree * tree = (TTree*)file->Get("Channel_8");
    tree->Draw("Charge>>hist"); 
    TH1D * hist = (TH1D*)gDirectory->Get("hist");
    
    inputValues inVal;
    
    inVal.pw = pw;
    inVal.n  = hist->GetEntries();
    inVal.v  = pow(hist->GetStdDev(), 2);
    inVal.e  = -hist->GetMean(); // Flip sign
    inVal.a  = tree->Draw("Charge", (TString)cutString);
    inVal.h  = hist;
    
    // Debugging
    if (inVal.pw == 7.5) { cout << inVal.a/inVal.n << endl; }

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
    //cout << setw(8) << outVal.n0 << " " << setw(8) << outVal.eL << " " << outVal.ePsi << " " << outVal.vPsi << endl;

    // Calculate uncertainties in output values
    outVal.v_eL   = ( exp(outVal.eL) + 1 - 2*b.a/b.n ) / b.a;
    outVal.v_ePsi = (outVal.eL * (pow(outVal.ePsi, 2) + outVal.vPsi) + 2*b.v) / (b.n * pow(outVal.eL, 2)) + pow(outVal.ePsi, 2) * (outVal.v_ePsi) / pow(outVal.eL, 2);

    // Debugging
    //cout << sqrt(outVal.v_eL) << " " << sqrt(outVal.v_ePsi) << endl;

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

    for (int i = 0; i < cuts.size(); i++) {
        eLPlots.push_back(new TGraph());
        ePsiPlots.push_back(new TGraph());
        vPsiPlots.push_back(new TGraph());
        eLPlots[i]->SetMarkerStyle(21);
        ePsiPlots[i]->SetMarkerStyle(21);
        vPsiPlots[i]->SetMarkerStyle(21);
        eLPlots[i]->SetMarkerColor(i+1);
        ePsiPlots[i]->SetMarkerColor(i+1);
        vPsiPlots[i]->SetMarkerColor(i+1);

        for (int j = 0; j < fileNames_t.size(); j++) {
            int o_ix = j + (i*(fileNames_t.size()));
            eLPlots[i]->SetPoint(j, t[j].pw, o[o_ix].eL);
            ePsiPlots[i]->SetPoint(j, t[j].pw, o[o_ix].ePsi);
            vPsiPlots[i]->SetPoint(j, t[j].pw, o[o_ix].vPsi);
	    }

        eLPlot->Add(eLPlots[i]);
        ePsiPlot->Add(ePsiPlots[i]);
        vPsiPlot->Add(vPsiPlots[i]);
    }

    eLPlot->SetTitle("Average Number of LED-induced Photoelectrons per Trigger in Normal Mode");
    ePsiPlot->SetTitle("Single-photoelectron Response Mean in Few-PE Regime in Normal Mode");
    vPsiPlot->SetTitle("Single-photoelectron Response Variance");


    TCanvas * eL_c = new TCanvas("eL_c", "eL Canvas", 200, 10, 700, 500);
    eLPlot->Draw("AP");
    eLPlot->GetXaxis()->SetTitle("LED Pulse Width / (ns)");
    eLPlot->GetYaxis()->SetTitle("Photoelectrons / Trigger");
    eL_leg = new TLegend(0.15, 0.6, 0.55, 0.85);
    eL_leg->SetHeader("   Integration Window: 156.25ns");
    eL_leg->AddEntry(eLPlots[0], "Cut at  0.2pC", "P");
    eL_leg->AddEntry(eLPlots[1], "Cut at  0.1pC", "P");
    eL_leg->AddEntry(eLPlots[2], "Cut at -0.1pC", "P");
    eL_leg->AddEntry(eLPlots[3], "Cut at -0.2pC", "P");
    eL_leg->Draw();
    gPad->Modified();


    TCanvas * ePsi_c = new TCanvas("ePsi_c", "ePsi Canvas", 200, 10, 700, 500);
    ePsiPlot->Draw("AP");
    ePsiPlot->GetXaxis()->SetTitle("LED Pulse Width / (ns)");
    ePsiPlot->GetYaxis()->SetTitle("pC");
    ePsiPlot->GetYaxis()->SetRangeUser(0.2, 1.4);
    ePsi_leg = new TLegend(0.45, 0.18, 0.85, 0.43);
    ePsi_leg->SetHeader("   Integration Window: 156.25ns");
    ePsi_leg->AddEntry(eLPlots[0], "Cut at  0.2pC", "P");
    ePsi_leg->AddEntry(eLPlots[1], "Cut at  0.1pC", "P");
    ePsi_leg->AddEntry(eLPlots[2], "Cut at -0.1pC", "P");
    ePsi_leg->AddEntry(eLPlots[3], "Cut at -0.2pC", "P");
    ePsi_leg->Draw();
    gPad->Modified();

    //TCanvas * vPsi_c = new TCanvas("vPsi_c", "vPsi Canvas", 200, 10, 700, 500);
    //vPsiPlot->Draw("AP");
    //vPsiPlot->GetXaxis()->SetTitle("LED Pulse Width / (ns)");
    //vPsiPlot->GetYaxis()->SetTitle("pC^2");
    //vPsi_leg = new TLegend(0.15, 0.6, 0.4, 0.85);
    //vPsi_leg->AddEntry(eLPlots[0], "Cut at 0.00pC", "P");
    //vPsi_leg->AddEntry(eLPlots[1], "Cut at 0.01pC", "P");
    //vPsi_leg->AddEntry(eLPlots[2], "Cut at 0.02pC", "P");
    //vPsi_leg->AddEntry(eLPlots[3], "Cut at 0.03pC", "P");
    //vPsi_leg->AddEntry(eLPlots[4], "Cut at 0.04pC", "P");
    //vPsi_leg->Draw();
    //gPad->Modified();
}

void singlePeVaryCut() {

    vector<inputValues> b; // Blank runs
    vector<inputValues> t; // Total runs
    vector<outputValues> o;
    vector<string> cutString(cuts.size(), "Charge > ");
    
    for (int i = 0; i < cuts.size(); i++) {

        // Create strings to perform cuts - I would be able to use to_string if C++11 worked 
        ostringstream ss;
        ss << cuts[i];
        cutString[i].append(ss.str());

        // Get input and output values from .root spectra
        b.push_back(getInputValues(fileNames_b[0], pulseWidths_b[0], cutString[i]));
        
        for (int j = 0; j < fileNames_t.size(); j++) {
            int t_ix = j + i*fileNames_t.size();
            t.push_back(getInputValues(fileNames_t[j], pulseWidths_t[j], cutString[i]));
            o.push_back(getOutputValues(b[i], t[t_ix]));
        }
    }
    makePlots(b, t, o);
}

