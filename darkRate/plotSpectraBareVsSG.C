const char * dirName = "/data/users/MilliQan/caen/PeriodicDarkRate/";
const char * ext     = ".root";

void plotSpectraBareVsSG() {

    vector<TTree*> trees; 

    TSystemDirectory dir(dirName, dirName);
    TList * files = dir.GetListOfFiles();

    // Read all files that contain Saint-Gobain data
    if (files) {
        cout << endl << "Files w/ Saint-Gobain Data:" << endl << endl;
        TSystemFile * file;
        TString fName;
        TIter next(files);
        while ((file = (TSystemFile*)next())) {
            fName = file->GetName();
            if (!file->IsDirectory() && fName.EndsWith(ext) && fName != "bfrancis_dark_10_17_16.root") {
                TFile * tempFile = new TFile(dirName + fName);
                TTree * tempTree = (TTree*)tempFile->Get("data");
                float temp;
                if (tempTree->SetBranchAddress("min_4", &temp) == 0) {
                    cout << fName.Data() << endl;
                    trees.push_back(tempTree);
                }
            }
        }
        cout << endl;
    }

    TH1F * h_BT = new TH1F("BareTubeChargeSpectrum", "BareTubeChargeSpectrum", 500, -10, 100);
    TH1F * h_SG = new TH1F("SaintGobChargeSpectrum", "SaintGobChargeSpectrum", 500, -10, 100);

    float btMin;
    float sgMin;
    
    for (int j = 0; j < trees.size(); j++) {
        trees[j]->SetBranchAddress("min_2", &btMin);
        trees[j]->SetBranchAddress("min_4", &sgMin);
        
        for (int i = 0; i < trees[j]->GetEntries(); i++) {
            trees[j]->GetEntry(i);
            h_BT->SetLineColor(kRed);
            h_BT->Fill(-1*btMin);
            h_SG->SetLineColor(kBlue);
            h_SG->Fill(-1*sgMin);
        }
    }

    TCanvas * c = new TCanvas("c", "Charge Spectra", 10, 10, 800, 600);
    gStyle->SetOptStat(0);

    h_BT->SetTitle("Dark Run Amplitude Spectra");
    h_BT->GetXaxis()->SetTitle("Pulse Amplitude [mV]");
    h_BT->GetYaxis()->SetTitle("Events");
    h_BT->SetLineWidth(2);
    h_SG->SetLineWidth(2);
    h_BT->Draw();
    h_SG->Draw("SAME");

    TLegend * leg = new TLegend(0.6, 0.7, 0.88, 0.88);
    leg->AddEntry(h_BT, "Bare R7725", "l");
    leg->AddEntry(h_SG, "Saint-Gobain Assembly", "l");
    leg->Draw();
}
