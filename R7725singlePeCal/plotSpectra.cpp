


const vector<string> fileNames = {
                                  "spectra/1700V_07.2ns.root",
                                  };

//const vector<string> leg_entries = {"Blank Run", "6.8ns", "7.4ns", "7.8ns", "9.0ns"};

TH1F* getHist(string fileName) {

  TFile * file = TFile::Open((TString)fileName);
  TTree * tree = (TTree*)file->Get("Channel_1");

  tree->Draw("Charge*(-1)>>h(300, -10, 200)");
  
  return (TH1F*)gDirectory->Get("h");
}

void plotSpectra() {

    vector<TH1F*> hists;

    for(int i = 0; i < fileNames.size(); i++) {
        hists.push_back(getHist(fileNames[i]));
    }

    //TLegend * leg = new TLegend(0.60, 0.70, 0.88, 0.88);
    //leg->SetHeader("LED Pulse Width");
    
    for(int i =0; i < hists.size(); i++) {
        hists[i]->SetLineColor(i);
        hists[i]->SetLineWidth(2);
        //hists[i]->Sumw2();
        hists[i]->Draw("SAME"); 
        //leg->AddEntry(hists[i], (TString)leg_entries[i]);
        
    }

    hists[hists.size()-1]->SetLineColor(hists.size()+1);
    hists[hists.size()-1]->SetTitle("R7725 PMT Spectrum from a 7.2ns LED Pulse");
    hists[hists.size()-1]->GetXaxis()->SetTitle("Integrated Charge / (pC)");
    hists[hists.size()-1]->GetYaxis()->SetTitle("Events");
    hists[hists.size()-1]->GetYaxis()->SetRangeUser(0.5, 30000);
    gStyle->SetOptStat(0);
    //leg->Draw();

}

