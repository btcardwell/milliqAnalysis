
void plotTimeOfDrop() {
    const vector<int> times = {0, 1650, 1500, 99999, 99999, 99999, 0, 99999, 930, 99999, 
                               1440, 0, 0, 0, 99999, 0, 0, 1380, 99999};
    
    TH1F * h = new TH1F("TimeOfDrop", "Time of Rate Drop", 15, 900, 1800);

    for(auto time : times) { h->Fill(time); }

    gStyle->SetOptStat(110011);
    h->GetXaxis()->SetTitle("Time since turning on PMT voltage / (s)");
    h->GetYaxis()->SetTitle("Runs");
    h->Draw();
}
