


void plotSpectraBlank() {
    
    TFile *f_charge = TFile::Open("spectra/06_28_16_ChargeMode_Blank.root");
    TFile *f_norm   = TFile::Open("spectra/06_28_16_NormalMode_Blank.root");
    TTree *t_charge = (TTree*)f_charge->Get("Channel_8");
    TTree *t_norm   = (TTree*)f_norm->Get("Channel_8");
    
    t_charge->Draw("Charge>>h1(100,-4, 3)");
    TH1F * h_charge = (TH1F*)gDirectory->Get("h1");
    t_norm->Draw("Charge>>h2(100,-4, 3)");
    TH1F * h_norm = (TH1F*)gDirectory->Get("h2");
    
    
    h_norm->SetTitle("Blank PMT Spectra in Normal Mode and Charge Mode");
    h_norm->GetXaxis()->SetTitle("Integrated Charge / (pC)");
    h_norm->GetYaxis()->SetTitle("Events");
    
    h_charge->SetLineWidth(2);
    h_norm->SetLineWidth(2);
    h_charge->SetLineColor(2);
    h_norm->Draw();
    h_charge->Draw("SAME");
    
    TLegend * leg = new TLegend(0.65, 0.75, 0.95, 0.92);
    leg->SetHeader("Blank PMT Spectra");
    leg->AddEntry(h_charge, "Charge Mode");
    leg->AddEntry(h_norm, "Normal Mode");
    leg->Draw();

}

