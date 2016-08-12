


const vector<string> fileNames = {
                                  "spectra/1700V_Blank.root",
                                  "spectra/1700V_07.4ns.root",
                                  "spectra/1700V_08.0ns.root"
                                 }; 

const vector<string> leg_entries = { "1300V", "1500V", "1700V" };

TH1F* getHist(string fileName) {

  TFile * file = TFile::Open((TString)fileName);
  TTree * tree = (TTree*)file->Get("Channel_1");

  tree->Draw("Min*(-1)>>h(200, 0, 200)");
  
  return (TH1F*)gDirectory->Get("h");
}

void plotSpectraMin() {

    vector<TH1F*> hists;

    for(int i = 0; i < fileNames.size(); i++) {
        hists.push_back(getHist(fileNames[i]));
    }

    TLegend * leg = new TLegend(0.60, 0.70, 0.88, 0.88);
    leg->SetHeader("PMT Voltage");
    
    for(int i =0; i < hists.size(); i++) {
        hists[i]->SetLineColor(i+1);
        hists[i]->SetLineWidth(2);
        hists[i]->Draw("SAME"); 
        leg->AddEntry(hists[i], (TString)leg_entries[i]);
    }

    hists[hists.size()-1]->SetLineColor(hists.size()+1);
    hists[hists.size()-1]->SetTitle(
            "PMT Spectra from 7.2ns LED Pulse with Varying PMT Supply Voltage");
    hists[hists.size()-1]->GetXaxis()->SetTitle("Maximum Amplitude / (mV)");
    hists[hists.size()-1]->GetYaxis()->SetTitle("Events");
    hists[hists.size()-1]->GetYaxis()->SetRangeUser(0.5, 100000);
    gStyle->SetOptStat(0);
    leg->Draw();

}

