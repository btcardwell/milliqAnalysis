


const vector<string> fileNames = {"spectra/07_12_16_10cm.root",
                                  "spectra/07_12_16_50cm.root",
                                  "spectra/07_12_16_80cm.root"};

TH1F* getHist(string fileName) {

  TFile * file = TFile::Open((TString)fileName);
  TTree * tree = (TTree*)file->Get("Channel_7");

  tree->Draw("Charge*(-1)>>h(200, -10, 300)");
  
  return (TH1F*)gDirectory->Get("h");
}

void plotSpectra() {

    vector<TH1F*> hists;

    for(int i = 0; i < fileNames.size(); i++) {
        hists.push_back(getHist(fileNames[i]));
    }
    

    for(int i =0; i < hists.size(); i++) {
        hists[i]->SetLineColor(i+2);
        hists[i]->Draw("SAME"); 
        
    }

    hists[2]->SetTitle("Cesium Spectra while Varying Source Distance from PMT");
    hists[2]->GetXaxis()->SetTitle("Integrated Charge / (pC)");
    hists[2]->GetYaxis()->SetTitle("Events");
    gStyle->SetOptStat(0);

    TLegend * leg = new TLegend(0.60, 0.70, 0.88, 0.88);
    leg->AddEntry(hists[0], "10cm");
    leg->AddEntry(hists[1], "50cm");
    leg->AddEntry(hists[2], "80cm");
    leg->Draw();
}

