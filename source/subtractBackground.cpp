
void subtractBackground() {
    
    TFile *f_noSource = TFile::Open("spectra/07_12_16_OutBox_NoSource.root");
    TFile *f_source   = TFile::Open("spectra/07_12_16_10cm.root");

    TTree *t_noSource = (TTree*)f_noSource->Get("Channel_7");
    TTree *t_source   = (TTree*)f_source->Get("Channel_7");
    
    t_noSource->Draw("Charge*(-1)>>h1(200, -50, 300)");
    TH1F * h_noSource = (TH1F*)gDirectory->Get("h1");
    t_source->Draw("Charge*(-1)>>h2(200, -50, 300)");
    TH1F * h_source = (TH1F*)gDirectory->Get("h2");
    
    h_source->SetTitle("Cesium Spectrum with Subtracted Background");
    h_source->GetXaxis()->SetTitle("Integrated Charge / (pC)");
    h_source->GetYaxis()->SetTitle("Events");

    TH1F * h_diff = (TH1F*)h_source->Clone();
    h_diff->Add(h_noSource, -1);

    h_noSource->SetLineWidth(2);
    h_noSource->SetLineColor(1);
    h_source->SetLineWidth(2);
    h_source->SetLineColor(2);
    h_diff->SetLineWidth(2);
    h_diff->SetLineColor(4);

    h_source->Draw();
    h_noSource->Draw("SAME");
    h_diff->Draw("SAME");
    
    TLegend * leg = new TLegend(0.60, 0.70, 0.88, 0.88);
    leg->AddEntry(h_noSource, "No Source");
    leg->AddEntry(h_source, "Source");
    leg->AddEntry(h_diff, "Subtracted background");
    leg->Draw();
}
