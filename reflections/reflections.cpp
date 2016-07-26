

const vector<string> fileNames = {
                                    //"spectra/07_15_16_Source_Sig3.root",
                                    //"spectra/07_15_16_Source_Sig4.root",
                                    "spectra/07_15_16_Source_Sig5.root",
                                    "spectra/07_15_16_Source_Sig6.root",
                                    "spectra/07_15_16_Source_Sig7.root",
                                    "spectra/07_15_16_Source_Sig8.root"
                                 };

const vector<int> significances = {/*3, 4,*/ 5, 6, 7, 8};


class Spectrum {
  public:
    int sig;
    vector<vector<float>> peakLocations;
    vector<vector<float>> adjPeakLocations;
    vector<float> h_adjPeakLocations;
    TH1F * h;
    TH1F * adj_h;
    

    Spectrum(string, int);
};


Spectrum::Spectrum(string f, int s) {

    sig = s;
    h = new TH1F("h", "h", 300, 0, 300);
    adj_h = new TH1F("adj_h", "adj_h", 300, 0, 300);

    TFile * file = TFile::Open((TString)f);
    TTree * tree = (TTree*)file->Get("Channel_7");

    vector<float> *  pL = 0;
    tree->SetBranchAddress("PeakLocations", &pL);
    for (int i = 0; i < tree->GetEntries(); i++) {
        tree->GetEntry(i);
        peakLocations.push_back(*pL);
        adjPeakLocations.push_back(*pL); // Copy for now, adjust below
    }

    // Fill histograms
    for (int i = 0; i < peakLocations.size(); i++) {
        for (int j =0; j < peakLocations[i].size(); j++) {
            
            adjPeakLocations[i][j] -= peakLocations[i][0];

            h->Fill(peakLocations[i][j]);
            adj_h->Fill(adjPeakLocations[i][j]);
        }
    }

    // Find peaks in adjusted histogram
    TSpectrum * spec = new TSpectrum();
    spec->Search(adj_h, 2);
    h_adjPeakLocations.push_back(0);
    for (int i = 0; i < spec->GetNPeaks(); i++) {
        h_adjPeakLocations.push_back(spec->GetPositionX()[i]);
    }
    sort(h_adjPeakLocations.begin(), h_adjPeakLocations.end());

}


void drawHists(vector<Spectrum> spectra) {
    
    const vector<string> legLabels = {"5", "6", "7", "8"};

    gPad->GetCanvas()->Clear();

    spectra[0].h->SetTitle("All Found Minima in PMT Traces");
    spectra[0].h->GetXaxis()->SetRangeUser(0, 300);
    spectra[0].h->GetYaxis()->SetRangeUser(0, 13000);
    spectra[0].h->GetXaxis()->SetTitle("Sample Cell");
    spectra[0].h->GetYaxis()->SetTitle("Events");

    TLegend * leg = new TLegend(0.1, 0.7, 0.5, 0.9);
    leg->SetHeader("Number of consecutive ascending/descending samples to define a minimum");
 
    for (int i = 0; i < spectra.size(); i++) { 
        spectra[i].h->SetLineWidth(2);
        spectra[i].h->SetLineColor(i+6);
        spectra[i].h->Draw("SAME");
        leg->AddEntry(spectra[i].h, (TString)legLabels[i]);
    }

    leg->Draw();
}


void drawAdjHists(vector<Spectrum> spectra) {
    
    const vector<string> legLabels = {"5", "6", "7", "8"};

    gPad->GetCanvas()->Clear();

    spectra[0].adj_h->SetTitle("All Found Minima in PMT Traces Relative to First Found Minima");
    spectra[0].adj_h->GetXaxis()->SetRangeUser(0, 300);
    spectra[0].adj_h->GetYaxis()->SetRangeUser(0, 13000);
    spectra[0].adj_h->GetXaxis()->SetTitle("Sample Cell");
    spectra[0].adj_h->GetYaxis()->SetTitle("Events");

    TLegend * leg = new TLegend(0.1, 0.7, 0.5, 0.9);
    leg->SetHeader("Number of consecutive ascending/descending samples to define a minimum");

    for (int i = 0; i < spectra.size(); i++) { 
        spectra[i].adj_h->SetLineWidth(2);
        spectra[i].adj_h->SetLineColor(i+6);
        spectra[i].adj_h->Draw("SAME");
        leg->AddEntry(spectra[i].adj_h, (TString)legLabels[i]);
    }
    
    leg->Draw();

}


TH1F * drawDiffHist(vector<Spectrum> spectra) {

    gPad->GetCanvas()->Clear();

    TH1F * h_diff = new TH1F("h_diff", "h_diff", 30, 0, 30);

    h_diff->SetTitle("Distances Between Consecutive Peaks in Relative Minima");
    h_diff->GetXaxis()->SetTitle("Difference in Peak Sample Cell");
    h_diff->GetYaxis()->SetTitle("Events");
    h_diff->GetYaxis()->SetLabelSize(0.025);
    h_diff->SetLineWidth(2);

    for (int i = 0; i < spectra.size(); i++) {
        for(int j = 1; j < spectra[i].h_adjPeakLocations.size(); j++) {
            h_diff->Fill(spectra[i].h_adjPeakLocations[j] - spectra[i].h_adjPeakLocations[j-1]);
        }
    }
    h_diff->Draw();

    return h_diff;
}


THStack * drawStackedBySigHist(vector<Spectrum> spectra) {

    const vector<string> legLabels = {"5", "6", "7", "8"};

    gPad->GetCanvas()->Clear();

    vector<TH1F*> h_diffs;
    THStack  * h_diff_stack = new THStack("h_diff_stack", "");
    TLegend * leg = new TLegend(0.1, 0.7, 0.5, 0.9);
    leg->SetHeader("Number of consecutive ascending/descending samples to define a minimum");

    for (int i = 0; i < spectra.size(); i++) {

        h_diffs.push_back(new TH1F("h_diff", "h_diff", 30, 0, 30));
        h_diffs[i]->SetLineColor(1);
        h_diffs[i]->SetFillColor(i+6);
        h_diffs[i]->SetLineWidth(2);

        for(int j = 1; j < spectra[i].h_adjPeakLocations.size(); j++) {
            h_diffs[i]->Fill(spectra[i].h_adjPeakLocations[j] - spectra[i].h_adjPeakLocations[j-1]);
        }

        h_diff_stack->Add(h_diffs[i]);
        leg->AddEntry(h_diffs[i], (TString)legLabels[i], "f");
    }

    h_diff_stack->Draw();
    leg->Draw();

    h_diff_stack->SetTitle("Distances Between Consecutive Peaks in Relative Minima");
    h_diff_stack->GetXaxis()->SetTitle("Difference in Peak Sample Cell");
    h_diff_stack->GetYaxis()->SetTitle("Events");
    h_diff_stack->GetYaxis()->SetLabelSize(0.025);

    gPad->GetCanvas()->Modified();

    return h_diff_stack;
}


THStack * drawStackedByNpeaksHist(vector<Spectrum> spectra) {

    THStack  * h_diff_stack = new THStack("h_diff_stack", "");
    vector<TH1F*> h_diffs;
    TLegend * leg = new TLegend(0.1, 0.7, 0.5, 0.9);
    leg->SetHeader("Peaks used to define difference");

    for (int i = 0; i < spectra.size(); i++) {
        for(int j = 1; j < spectra[i].h_adjPeakLocations.size() && j < 20; j++) {

            if (h_diffs.size() < j) {
                h_diffs.push_back(new TH1F("h_diff", "h_diff", 30, 0, 30));
                leg->AddEntry(h_diffs[j-1], Form("%i and %i", j, j+1));
            }

            h_diffs[j-1]->Fill(spectra[i].h_adjPeakLocations[j] - spectra[i].h_adjPeakLocations[j-1]);
            h_diffs[j-1]->SetLineColor(1);
            h_diffs[j-1]->SetLineWidth(2);
            h_diffs[j-1]->SetFillColor(j);
        }
    }
    for(auto h : h_diffs) { h_diff_stack->Add(h); }

    h_diff_stack->Draw();
    leg->Draw();

    h_diff_stack->SetTitle("Distances Between Consecutive Peaks in Relative Minima");
    h_diff_stack->GetXaxis()->SetTitle("Difference in Peak Sample Cell");
    h_diff_stack->GetYaxis()->SetTitle("Events");
    h_diff_stack->GetYaxis()->SetLabelSize(0.025);

    gPad->GetCanvas()->Modified();

    return h_diff_stack;
}


THStack * drawStackedByNpeaksHistSimple(vector<Spectrum> spectra) {

    THStack  * h_diff_stack = new THStack("h_diff_stack", "");
    vector<TH1F*> h_diffs;
    h_diffs.push_back(new TH1F("h_diff_less", "h_diff_less", 30, 0, 30));
    h_diffs.push_back(new TH1F("h_diff_greater", "h_diff_greater", 30, 0, 30));

    TLegend * leg = new TLegend(0.1, 0.7, 0.5, 0.9);
    leg->SetHeader("Peaks used to define difference");

    for (int i = 0; i < spectra.size(); i++) {
        for(int j = 1; j < spectra[i].h_adjPeakLocations.size(); j++) {
            h_diffs[(j <= 7) ? 0 : 1]->Fill(spectra[i].h_adjPeakLocations[j] - spectra[i].h_adjPeakLocations[j-1]);
        }
    }
    h_diffs[0]->SetLineColor(1);
    h_diffs[0]->SetLineWidth(2);
    h_diffs[0]->SetFillColor(8);
    h_diffs[1]->SetLineColor(1);
    h_diffs[1]->SetLineWidth(2);
    h_diffs[1]->SetFillColor(9);
    for(auto h : h_diffs) { h_diff_stack->Add(h); }
    h_diff_stack->Draw();

    leg->AddEntry(h_diffs[0], "1 through 7");
    leg->AddEntry(h_diffs[1], "8 through 12");
    leg->Draw();

    h_diff_stack->SetTitle("Distances Between Consecutive Peaks in Relative Minima");
    h_diff_stack->GetXaxis()->SetTitle("Difference in Peak Sample Cell");
    h_diff_stack->GetYaxis()->SetTitle("Events");
    h_diff_stack->GetYaxis()->SetLabelSize(0.025);

    gPad->GetCanvas()->Modified();

    return h_diff_stack;
}


void fft(vector<Spectrum> spectra) {

    int nBins = 100;
    TH1 * fft_h = new TH1F("fft_h", "fft_h", nBins, 0, nBins); 

    spectra[2].adj_h->FFT(fft_h, "MAG");
    
    // Start debugging
    //TF1 *  f_sin = new TF1("f_sin", "sin(10*x)", 0, 2*TMath::Pi());
    //TH1F * h_sin = new TH1F("h_sin", "h_sin", nBins + 1, 0, 2*TMath::Pi());  
    //float x;
    //for (int i = 0; i <= nBins; i++) {
        //x = (float)i / nBins * 2*TMath::Pi();  
        //h_sin->SetBinContent(i+1, f_sin->Eval(x));
        //cout << i << " " <<  x << " " << f_sin->Eval(x) << endl;
    //}
    //h_sin->FFT(fft_h, "MAG");

    //h_sin->Draw();
    // End debugging

    fft_h->Draw();
}


void reflections() {

    vector<Spectrum> spectra;

    // Read spectra
    for (int i = 0; i < fileNames.size(); i++) {
        spectra.push_back(Spectrum(fileNames[i], significances[i]));
    }

    // Testing
    for (int i = 0; i < spectra[0].h_adjPeakLocations.size(); i++) { 
        cout << setw(6) << i*16 << setw(6) << i*11 << setw(6) << spectra[0].h_adjPeakLocations[i] << endl; }

    //drawHists(spectra);
    //drawAdjHists(spectra);
    //drawDiffHist(spectra);
    //drawStackedBySigHist(spectra);
    //drawStackedByNpeaksHist(spectra);
    //drawStackedByNpeaksHistSimple(spectra);
    //fft(spectra);
}







