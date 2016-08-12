


void plotSpectraManyPe() {
    
    TFile *f_charge = TFile::Open("spectra/06_28_16_ChargeMode_10.0ns.root");
    TFile *f_norm   = TFile::Open("spectra/06_28_16_NormalMode_10.0ns.root");
    TTree *t_charge = (TTree*)f_charge->Get("Channel_8");
    TTree *t_norm   = (TTree*)f_norm->Get("Channel_8");
    
    t_charge->Draw("Charge>>h1(100,-50, 25)");
    TH1F * h_charge = (TH1F*)gDirectory->Get("h1");
    t_norm->Draw("Charge>>h2(100,-50, 25)");
    TH1F * h_norm = (TH1F*)gDirectory->Get("h2");
    
    
    h_charge->SetTitle("Pulsed-LED PMT Spectra in Normal Mode and Charge Mode");
    h_charge->GetXaxis()->SetTitle("Integrated Charge / (pC)");
    h_charge->GetYaxis()->SetTitle("Events");
    
    h_charge->SetLineWidth(2);
    h_norm->SetLineWidth(2);
    h_charge->SetLineColor(2);
    h_charge->Draw();
    h_norm->Draw("SAME");
    
    TLegend * leg = new TLegend(0.65, 0.75, 0.95, 0.92);
    leg->SetHeader("10.0ns LED Pulse PMT Spectra");
    leg->AddEntry(h_charge, "Charge Mode");
    leg->AddEntry(h_norm, "Normal Mode");
    leg->Draw();

}

