


const vector<string> fileNames = {"spectra/07_15_16_Source_Sig3.root",
                                  "spectra/07_15_16_Source_Sig5.root",
                                  "spectra/07_15_16_Source_Sig7.root"}; 

TH1F* getHist(string fileName) {

  TFile * file = TFile::Open((TString)fileName);
  TTree * tree = (TTree*)file->Get("Channel_7");

  tree->Draw("PeakDeltaT>>h(180, 0, 60)");
  
  return (TH1F*)gDirectory->Get("h");
}

void plotSpectra() {

    vector<TH1F*> hists;

    for(int i = 0; i < fileNames.size(); i++) {
        hists.push_back(getHist(fileNames[i]));
    }
    

    for(int i =0; i < hists.size(); i++) {
        hists[i]->SetLineColor(i+2);
        hists[i]->SetLineWidth(2);
        hists[i]->Draw("SAME"); 
        
    }

    hists[hists.size() - 1]->SetTitle("Time Between First Two Detected Minima");
    hists[hists.size() - 1]->GetXaxis()->SetTitle("Time / (ns)");
    hists[hists.size() - 1]->GetXaxis()->SetLabelSize(0.025);
    hists[hists.size() - 1]->GetYaxis()->SetTitle("Events");
    hists[hists.size() - 1]->GetYaxis()->SetLabelSize(0.025);
    gStyle->SetOptStat(0);

    TLegend * leg = new TLegend(0.60, 0.70, 0.88, 0.88);
    leg->SetHeader("Ascending/Descending Samples Required");
    leg->AddEntry(hists[0], "3");
    leg->AddEntry(hists[1], "5");
    leg->AddEntry(hists[2], "7");
    leg->Draw();
}

