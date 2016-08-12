


const vector<string> fileNames = {"spectra/06_28_16_NormalMode_ChargeLength500_Blank.root",
                                  //"spectra/06_28_16_NormalMode_ChargeLength500_6.6ns.root",
                                  "spectra/06_28_16_NormalMode_ChargeLength500_6.8ns.root",
                                  //"spectra/06_28_16_NormalMode_ChargeLength500_7.0ns.root",
                                  //"spectra/06_28_16_NormalMode_ChargeLength500_7.2ns.root",
                                  "spectra/06_28_16_NormalMode_ChargeLength500_7.4ns.root",
                                  //"spectra/06_28_16_NormalMode_ChargeLength500_7.6ns.root",
                                  "spectra/06_28_16_NormalMode_ChargeLength500_7.8ns.root",
                                  //"spectra/06_28_16_NormalMode_ChargeLength500_8.0ns.root",
                                  //"spectra/06_28_16_NormalMode_ChargeLength500_8.2ns.root",
                                  "spectra/06_28_16_NormalMode_ChargeLength500_9.0ns.root"};

const vector<string> leg_entries = {"Blank Run", "6.8ns", "7.4ns", "7.8ns", "9.0ns"};

TH1F* getHist(string fileName) {

  TFile * file = TFile::Open((TString)fileName);
  TTree * tree = (TTree*)file->Get("Channel_8");

  tree->Draw("Charge*(-1)>>h(500, -2, 40)");
  
  return (TH1F*)gDirectory->Get("h");
}

void plotSpectra() {

    vector<TH1F*> hists;

    for(int i = 0; i < fileNames.size(); i++) {
        hists.push_back(getHist(fileNames[i]));
    }

    TLegend * leg = new TLegend(0.60, 0.70, 0.88, 0.88);
    leg->SetHeader("LED Pulse Width");
    
    for(int i =0; i < hists.size(); i++) {
        hists[i]->SetLineColor(i+1);
        hists[i]->SetLineWidth(2);
        //hists[i]->Sumw2();
        hists[i]->Draw("SAME"); 
        leg->AddEntry(hists[i], (TString)leg_entries[i]);
        
    }

    hists[hists.size()-1]->SetLineColor(hists.size()+1);
    hists[hists.size()-1]->SetTitle("PMT Spectra with Varying LED Pulse Width");
    hists[hists.size()-1]->GetXaxis()->SetTitle("Integrated Charge / (pC)");
    hists[hists.size()-1]->GetYaxis()->SetTitle("Events");
    hists[hists.size()-1]->GetYaxis()->SetRangeUser(0.5, 20000);
    gStyle->SetOptStat(0);
    leg->Draw();

}

