#include "TH1.h"
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TGraph.h"

#include <vector>
#include <string>
#include <stdlib.h>

const vector<string> fileNames_b = {"spectra/06_14_16_VaryVoltage_1250V_Blank_v2.root",
                                    "spectra/06_14_16_VaryVoltage_1300V_Blank_v2.root",
                                    "spectra/06_14_16_VaryVoltage_1350V_Blank_v2.root",
                                    "spectra/06_14_16_VaryVoltage_1400V_Blank_v2.root"};

const vector<string> fileNames_t = {"spectra/06_14_16_VaryVoltage_1250V_7.5ns.root",
                                    "spectra/06_14_16_VaryVoltage_1250V_8.0ns.root",
                                    "spectra/06_14_16_VaryVoltage_1250V_8.2ns.root",
                                    "spectra/06_14_16_VaryVoltage_1300V_7.5ns.root",
                                    "spectra/06_14_16_VaryVoltage_1300V_8.0ns.root",
                                    "spectra/06_14_16_VaryVoltage_1300V_8.2ns.root",
                                    "spectra/06_14_16_VaryVoltage_1350V_7.5ns.root",
                                    "spectra/06_14_16_VaryVoltage_1350V_8.0ns.root",
                                    "spectra/06_14_16_VaryVoltage_1350V_8.2ns.root",
                                    "spectra/06_14_16_VaryVoltage_1400V_7.5ns.root",
                                    "spectra/06_14_16_VaryVoltage_1400V_8.0ns.root",
                                    "spectra/06_14_16_VaryVoltage_1400V_8.2ns.root"};

const vector<double> pulseWidths_b(fileNames_b.size(), 7.5); // ns
const vector<double> pulseWidths_t = {7.5, 8.0, 8.2, 7.5, 8.0, 8.2, 7.5, 8.0, 8.2, 7.5, 8.0, 8.2}; // ns
const vector<int> hv_b = {1250, 1300, 1350, 1400}; // Volts
const vector<int> hv_t = {1250,1250, 1250, 1300, 1300, 1300, 1350, 1350, 1350, 1400, 1400, 1400}; // Volts
const vector<double> cut_b = {-0.023, -0.022, -0.022, -0.022};

class inputValues {
    public:
    double pw, // LED pulse width
           n,  // Number of entries
    	   v,  // Distribution variance
    	   e,  // Distribution mean
    	   a,  // Number of entries below cut
           hv; // PMT voltage setting

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

inputValues getInputValues(TString fileName, double pw, string cutString, int hv) {

    TFile * file = new TFile(fileName, "READ");
    TTree * tree = (TTree*)file->Get("Channel_8");
    tree->Draw("Charge>>hist"); 
    TH1D * hist = (TH1D*)gDirectory->Get("hist");
    
    inputValues inVal;

    inVal.pw = pw;
    inVal.hv = hv;
    inVal.n = hist->GetEntries();
    inVal.v = pow(hist->GetStdDev(), 2);
    inVal.e = -hist->GetMean(); // Flip sign
    inVal.a = tree->Draw("Charge", (TString)cutString);
    inVal.h = hist;

    // Debugging
    //cout << inVal.pw << " " << inVal.hv << " " << inVal.a << endl;

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
    //cout << t.hv << " " << t.pw << " " << outVal.n0 << " " << outVal.eL << " " << outVal.ePsi<< endl;

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
            if (t[j].hv == b[i].hv) {
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

    eLPlot->SetTitle("Average Number of LED-induced Photoelectrons per Trigger");
    ePsiPlot->SetTitle("Single-photoelectron Response Mean in Few-PE Regime");
    vPsiPlot->SetTitle("Single-photoelectron Response Variance");


    TCanvas * eL_c = new TCanvas("eL_c", "eL Canvas", 200, 10, 700, 500);
    eLPlot->Draw("AP");
    eLPlot->GetXaxis()->SetTitle("LED Pulse Width / (ns)");
    eLPlot->GetYaxis()->SetTitle("Photoelectrons / Trigger");
    eLPlot->GetXaxis()->SetRangeUser(7.2, 8.4);
    eL_leg = new TLegend(0.15, 0.65, 0.5, 0.85);
    eL_leg->AddEntry(eLPlots[0], "PMT Voltage = 1250V", "P");
    eL_leg->AddEntry(eLPlots[1], "PMT Voltage = 1300V", "P");
    eL_leg->AddEntry(eLPlots[2], "PMT Voltage = 1350V", "P");
    eL_leg->AddEntry(eLPlots[3], "PMT Voltage = 1400V", "P");
    eL_leg->Draw();
    gPad->Modified();


    TCanvas * ePsi_c = new TCanvas("ePsi_c", "ePsi Canvas", 200, 10, 700, 500);
    ePsiPlot->Draw("AP");
    ePsiPlot->GetXaxis()->SetTitle("LED Pulse Width / (ns)");
    ePsiPlot->GetYaxis()->SetTitle("pC");
    ePsiPlot->GetXaxis()->SetRangeUser(7.2, 8.4);
    ePsiPlot->GetYaxis()->SetRangeUser(0.1, 0.6);
    ePsi_leg = new TLegend(0.15, 0.65, 0.5, 0.85);
    ePsi_leg->AddEntry(eLPlots[0], "PMT Voltage = 1250V", "P");
    ePsi_leg->AddEntry(eLPlots[1], "PMT Voltage = 1300V", "P");
    ePsi_leg->AddEntry(eLPlots[2], "PMT Voltage = 1350V", "P");
    ePsi_leg->AddEntry(eLPlots[3], "PMT Voltage = 1400V", "P");
    ePsi_leg->Draw();
    gPad->Modified();

    TCanvas * vPsi_c = new TCanvas("vPsi_c", "vPsi Canvas", 200, 10, 700, 500);
    vPsiPlot->Draw("AP");
    vPsiPlot->GetXaxis()->SetTitle("LED Pulse Width / (ns)");
    vPsiPlot->GetYaxis()->SetTitle("pC^2");
    vPsiPlot->GetXaxis()->SetRangeUser(7.2, 8.4);
    vPsi_leg = new TLegend(0.15, 0.65, 0.5, 0.85);
    vPsi_leg->AddEntry(eLPlots[0], "PMT Voltage = 1250V", "P");
    vPsi_leg->AddEntry(eLPlots[1], "PMT Voltage = 1300V", "P");
    vPsi_leg->AddEntry(eLPlots[2], "PMT Voltage = 1350V", "P");
    vPsi_leg->AddEntry(eLPlots[3], "PMT Voltage = 1400V", "P");
    vPsi_leg->Draw();
    gPad->Modified();

    // Plot spectra
    //TCanvas * spectra_c = new TCanvas("spectra_c", "Spectra Canvas", 200, 10, 700, 500);
    //spectra_c->SetLogy();
    //for (int i = t.size() - 1; i >= 0; i--) {
    //    //t[i].h->Sumw2();
    //    t[i].h->SetLineColor(i+1);
    //    t[i].h->Draw("SAME");
    //    t[i].h->SetStats(0);
    //    t[i].h->GetYaxis()->SetRangeUser(1.0, 500000.0);
    //    t[i].h->SetTitle("PMT Spectra of Pulsed LED w/ Varying Pulse Width");
    //    t[i].h->GetXaxis()->SetTitle("Integrated Charge / (pC)");
    //    t[i].h->GetYaxis()->SetTitle("Events");
    //}
    //b[0].h->Sumw2();
    //b[0].h->Draw("SAME");
    //b[0].h->SetStats(0);
    //TLine * cutLine_t = new TLine(cut, 0.0, cut, 500000.0);
    //cutLine_t->SetLineWidth(2);
    //cutLine_t->Draw();

    // Add legend
    //TLegend * leg = new TLegend(0.1, 0.7, 0.48, 0.9);
    //leg->SetHeader("Legend");
    //leg->AddEntry(t[0].h, "

    // Plot background
    //TCanvas * b_c = new TCanvas("b_c", "Background Canvas", 200, 10, 700, 500);
    //b_c->SetLogy();
    ////b[0].h->Sumw2();
    //b[0].h->Draw();
    //b[0].h->GetYaxis()->SetRangeUser(1.0, 100000.0);
    //b[0].h->SetTitle("Background PMT Spectrum");
    //b[0].h->GetXaxis()->SetTitle("Integrated Charge / (pC)");
    //b[0].h->GetYaxis()->SetTitle("Events");
    //TLine * cutLine_b = new TLine(cut, 0.0, cut, 1000000.0);
    //cutLine_b->SetLineWidth(2);
    //cutLine_b->Draw();
}

void singlePeVaryVoltage() {

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
        b.push_back(getInputValues(fileNames_b[i], pulseWidths_b[i], cutString[i], hv_b[i]));
        
        for (int j = 0; j < fileNames_t.size(); j++) {
            if (hv_t[j] == hv_b[i]) {
                t.push_back(getInputValues(fileNames_t[j], pulseWidths_t[j], cutString[i], hv_t[j]));
                o.push_back(getOutputValues(b[i], t[j]));
            }
        }
    }
    makePlots(b, t, o);
}





