//#include "rate.C"

void SetRootPalette(Int_t pal = 0){
  // pal =  1: rainbow\n"                                                                                                                                                                              
  // pal =  2: reverse-rainbow\n"                                                                                                                                                                      
  // pal =  3: amber\n"                                                                                                                                                                                
  // pal =  4: reverse-amber\n"                                                                                                                                                                        
  // pal =  5: blue/white\n"                                                                                                                                                                           
  // pal =  6: white/blue\n"                                                                                                                                                                           
  // pal =  7: red temperature\n"                                                                                                                                                                      
  // pal =  8: reverse-red temperature\n"                                                                                                                                                              
  // pal =  9: green/white\n"                                                                                                                                                                          
  // pal = 10: white/green\n"                                                                                                                                                                          
  // pal = 11: orange/blue\n"                                                                                                                                                                          
  // pal = 12: blue/orange\n"                                                                                                                                                                          
  // pal = 13: white/black\n"                                                                                                                                                                          
  // pal = 14: black/white\n"                                                                                                                                                                           
  const Int_t NRGBs = 5;
  const Int_t NCont = 255;
  gStyle->SetNumberContours(NCont);

  if (pal < 1 && pal> 14) return;
  else pal--;

  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[14][NRGBs]   = {{ 0.00, 0.00, 0.87, 1.00, 0.51 },
                               { 0.51, 1.00, 0.87, 0.00, 0.00 },
                               { 0.17, 0.39, 0.62, 0.79, 1.00 },
                               { 1.00, 0.79, 0.62, 0.39, 0.17 },
                               { 0.00, 0.00, 0.00, 0.38, 1.00 },
                               { 1.00, 0.38, 0.00, 0.00, 0.00 },
                               { 0.00, 0.50, 0.89, 0.95, 1.00 },
                               { 1.00, 0.95, 0.89, 0.50, 0.00 },
                               { 0.00, 0.00, 0.38, 0.75, 1.00 },
                               { 0.00, 0.34, 0.61, 0.84, 1.00 },
                               { 0.75, 1.00, 0.24, 0.00, 0.00 },
                               { 0.00, 0.00, 0.24, 1.00, 0.75 },
                               { 0.00, 0.34, 0.61, 0.84, 1.00 },
                               { 1.00, 0.84, 0.61, 0.34, 0.00 }
  };
  Double_t green[14][NRGBs] = {{ 0.00, 0.81, 1.00, 0.20, 0.00 },
                               { 0.00, 0.20, 1.00, 0.81, 0.00 },
                               { 0.01, 0.02, 0.39, 0.68, 1.00 },
                               { 1.00, 0.68, 0.39, 0.02, 0.01 },
                               { 0.00, 0.00, 0.38, 0.76, 1.00 },
                               { 1.00, 0.76, 0.38, 0.00, 0.00 },
                               { 0.00, 0.00, 0.27, 0.71, 1.00 },
                               { 1.00, 0.71, 0.27, 0.00, 0.00 },
                               { 0.00, 0.35, 0.62, 0.85, 1.00 },
                               { 1.00, 0.75, 0.38, 0.00, 0.00 },
                               { 0.24, 1.00, 0.75, 0.18, 0.00 },
                               { 0.00, 0.18, 0.75, 1.00, 0.24 },
                               { 0.00, 0.34, 0.61, 0.84, 1.00 },
                               { 1.00, 0.84, 0.61, 0.34, 0.00 }
  };
  Double_t blue[14][NRGBs]  = {{ 0.51, 1.00, 0.12, 0.00, 0.00 },
                               { 0.00, 0.00, 0.12, 1.00, 0.51 },
                               { 0.00, 0.09, 0.18, 0.09, 0.00 },
                               { 0.00, 0.09, 0.18, 0.09, 0.00 },
                               { 0.00, 0.47, 0.83, 1.00, 1.00 },
                               { 1.00, 1.00, 0.83, 0.47, 0.00 },
                               { 0.00, 0.00, 0.00, 0.40, 1.00 },
                               { 1.00, 0.40, 0.00, 0.00, 0.00 },
                               { 0.00, 0.00, 0.00, 0.47, 1.00 },
                               { 1.00, 0.47, 0.00, 0.00, 0.00 },
                               { 0.00, 0.62, 1.00, 0.68, 0.12 },
                               { 0.12, 0.68, 1.00, 0.62, 0.00 },
                               { 0.00, 0.34, 0.61, 0.84, 1.00 },
                               { 1.00, 0.84, 0.61, 0.34, 0.00 }
  };

  TColor::CreateGradientColorTable(NRGBs, stops, red[pal], green[pal], blue[pal], NCont);
}


void Plots_new(){

  SetRootPalette(1);
    gStyle->SetOptStat(0);
    gStyle->SetLabelSize(0.05,"X");
    gStyle->SetLabelSize(0.05,"Y");
    gStyle->SetLabelSize(0.05,"Z");
    gStyle->SetTitleSize(0.05,"X");
    gStyle->SetTitleSize(0.05,"Y");
    gStyle->SetTitleSize(0.05,"Z");

    //    TTree* T = ((TTree*)(TFile::Open("gn_pimp_N10000_ng_tm1.root")->Get("T")));
    //       TTree* T = ((TTree*)(TFile::Open("gn_pimp_N100000_carbon_test_th40_140.root")->Get("T")));
         TTree* T = ((TTree*)(TFile::Open("gn_pimp_N100000_deu_test3_th40_140.root")->Get("T")));

  //nubmer of targets:
  TString nt;
  nt.Form("/5.");
  TString unit;
  unit.Form("1");

  TString w;
  w.Form("*weight_kk*0.8*0.8");// _kk   detection efficiency

  TString wst, wst_f;
  wst.Form("TMath::Abs(t)>2.&&TMath::Abs(u)>2.");
  wst_f = wst;
  wst_f.Prepend("(");
  wst_f.Append(")");
  wst_f.Append(w);
  cout<<"wst"<<endl;
  cout<<wst_f.Data()<<endl;

  TString wSRC90, wSRC, wMF90, wMF;
  TString wSRC90coh, wSRCcoh, wMF90coh, wMFcoh;
  TString rec, rec1, tcm90, fac, coh;
  rec.Form("(Pmiss<0.25&&");
  rec.Append(wst);
  wMF.Append(rec);
  wMFcoh.Append(wMF);
  wMF.Append(")");
  wMF.Append(w);
  fac.Form("*0.65"); // efficiency for recoil reconstruction
  cout<<"wMF"<<endl;
  cout<<  wMF.Data()<<endl;

  //    coh.Form("&&Original_E>=7.5&&Original_E<=11.7)"); // 9.1    7.8
  coh.Form(")");
  wMFcoh.Append(coh);
  wMFcoh.Append(w);
  cout<<"wMFcoh"<<endl;
  cout<<  wMFcoh.Data()<<endl;

  //  cout<<"weight_kk*(TMath::Abs(t)>3&&TMath::Abs(u)>3&&theta_recoil<160&&Precoil>0.3)*0.5*0.8"<<endl;
  rec1.Form("(theta_recoil<160.&&Pmiss>0.3");// Precoil > 0.3
  wSRC.Append(rec1);
  wSRC.Append("&&");
  wSRC.Append(wst);
  wSRCcoh.Append(wSRC);
  wSRC.Append(")");
  wSRC.Append(w);
  wSRC.Append(fac);
  cout<<"wSRC"<<endl;
  cout<<  wSRC.Data()<<endl;

  wSRCcoh.Append(coh);
  wSRCcoh.Append(w);
  wSRCcoh.Append(fac);
  cout<<"wSRCcoh"<<endl;
  cout<<  wSRCcoh.Data()<<endl;

  //Draw beam photons:
  TH1F* GammaBeamHist= new TH1F("GammaBeamHist","E_{beam} [GeV]; [#gamma/s]", 46, 2.8, 12.); // nbins = 46
  TCanvas* cb = new TCanvas("cb","beam",800,800);
  T->Draw("Original_E>>GammaBeamHist");

  // Draw s:
  TH1F* s1 = new TH1F("s1","; s [GeV^{2}]",120, -5.,55.);
  TH1F* s2 = new TH1F("s2","; s [GeV^{2}]",120, -5.,55.);
  TH1F* s3 = new TH1F("s3","; s [GeV^{2}]",120, -5.,55.);
  s1->SetLineColor(8);
  s2->SetLineColor(4);
  s3->SetLineColor(6);
  s1->SetLineWidth(3);
  s2->SetLineWidth(3);
  s3->SetLineWidth(3);
  //TCanvas* Cs = new TCanvas("Cs","s",800,800);
  //T->Draw("s>>s1", unit,"h");
  //  T->Draw("s_check>>s2", unit,"same");
  //  T->Draw("s_misak>>s3", unit,"same");

  TLegend* legs = new TLegend(0.5,0.5,0.75,0.75);
  legs->SetTextSize(0.04);
  legs->SetFillColor(0);
  legs->SetLineColor(0);
  legs->AddEntry(s1, "s for final state","l");
  legs->AddEntry(s2, "s for scattering on a SRC nucleon", "l");
  legs->AddEntry(s3, "s from Misak","l");
  legs->Draw();
  
    // Mean Field. Theta c.m. = All
    TH2F* theta_P_p_MF_All = new TH2F("theta_P_p_MF_All","; |P_{p}| [GeV/c]; #theta_{p} [#circ]",24,0,12 , 14,5,75);
    TH2F* theta_P_pi_MF_All = new TH2F("theta_P_pi_MF_All","; |P_{#pi}| [GeV/c]; #theta_{#pi} [#circ]",24,0,12 , 14,5,75);
    TH2F* theta_Ppi_MF_All = new TH2F("theta_Ppi_MF_All","; #theta_{p} [degrees]; #theta_{#pi} [#circ]",14,5,75 , 14,5,75);
    TH1F* phi_MF_All = new TH1F("phi_MF_All","; #Delta#phi [degrees]",40,160,200);
    TCanvas* C2 = new TCanvas("C2","MF, all",800,800);
    C2->Divide(2,2);
    C2->cd(1);
    gPad->SetFillColor(0);
    T->Draw("theta_P4:absP4>>theta_P_pi_MF_All"      , wMFcoh,"col");
    C2->cd(2);
    gPad->SetFillColor(0);
    T->Draw("theta_P3:absP3>>theta_P_p_MF_All"      , wMFcoh,"col");
    C2->cd(3);
    gPad->SetFillColor(0);
    T->Draw("theta_P4:theta_P3>>theta_Ppi_MF_All", wMFcoh,"col");
    C2->cd(4);
    gPad->SetFillColor(0);
    T->Draw("TMath::Abs(phi_P4-phi_P3)>>phi_MF_All",wMFcoh,"h");
    
    // SRC. Theta c.m. = All
    TH2F* theta_P_p_SRC_All = new TH2F("theta_P_p_SRC_All","; |P_{p}| [GeV/c]; #theta_{p} [#circ]",24,0,12 , 14,5,75);
    TH2F* theta_P_pi_SRC_All = new TH2F("theta_P_pi_SRC_All","; |P_{#pi}| [GeV/c]; #theta_{#pi} [#circ]",24,0,12 , 14,5,75);
    TH2F* theta_Ppi_SRC_All = new TH2F("theta_Ppi_SRC_All","; #theta_{p} [degrees]; #theta_{#pi} [#circ]",14,5,75 , 14,5,75);
    TH1F* phi_SRC_All = new TH1F("phi_SRC_All","; #Delta#phi [degrees]",80,140,220);
    TCanvas* C4 = new TCanvas("C4","SRC, All",800,800);
    C4->Divide(2,2);
    C4->cd(1);
    gPad->SetFillColor(0);
    T->Draw("theta_P4:absP4>>theta_P_pi_SRC_All",wSRCcoh,"col");
    C4->cd(2);
    gPad->SetFillColor(0);
    T->Draw("theta_P3:absP3>>theta_P_p_SRC_All",wSRCcoh,"col");
    C4->cd(3);
    gPad->SetFillColor(0);
    T->Draw("theta_P4:theta_P3>>theta_Ppi_SRC_All", wSRCcoh,"col");
    C4->cd(4);
    gPad->SetFillColor(0);
    T->Draw("TMath::Abs(phi_P4-phi_P3)>>phi_SRC_All", wSRCcoh,"h");
    
    // Print the number of events!!!!

             cout<<"MF, All, N events = "<<theta_P_p_MF_All->Integral()<<endl;
     cout<<"SRC, All, N events = "<<theta_P_p_SRC_All->Integral()<<endl;
     cout<<"For rho: MF events = "<< theta_P_p_MF_All->Integral()*10./0.8/0.8*0.9*0.3<<", SRC = "<< theta_P_p_SRC_All->Integral()*10./0.8/0.8*0.3*0.9 <<endl;     

    TH1F* Ebeam_MF_All = new TH1F("Ebeam_MF_All","; E_{beam} [GeV]", 46, 2.8, 12);    
    TH1F* Ebeam_SRC_All = new TH1F("Ebeam_SRC_All","; E_{beam} [GeV]", 46, 2.8, 12);
    TCanvas* Cbe = new TCanvas("Cbe","Cbe",800,800);
    gPad->SetLogy();
    T->Draw("Original_E>>Ebeam_MF_All",wMF,"h");
    Ebeam_MF_All->SetLineColor(8);
    Ebeam_MF_All->SetLineWidth(3);
    gROOT->ForceStyle();
    T->Draw("Original_E>>Ebeam_SRC_All",wSRC,"sameh");
    Ebeam_SRC_All->SetLineColor(2);
    Ebeam_SRC_All->SetLineWidth(3);
    gROOT->ForceStyle();
    TLegend* leg = new TLegend(0.5,0.5,0.75,0.75);
    leg->SetLineColor(0);
    leg->SetFillColor(0);
    leg->SetTextSize(0.04);
    leg->AddEntry(Ebeam_MF_All,"MF, #theta_{cm}=(50,130)#circ","l");
    leg->AddEntry(Ebeam_SRC_All,"SRC, #theta_{cm}=(50,130)#circ","l");
    leg->Draw("samel");

    //E_beam vs t
    gStyle->SetPalette(56);
    TH2F* Et_MF_All = new TH2F("Et_MF_All","MF, #theta_{cm} in (50,130)#circ; E_{beam} [GeV]; |t| [GeV^{2}]",46, 2.8,12, 17,3.,20.);
    TH2F* Et_SRC_All = new TH2F("Et_SRC_All","SRC, #theta_{cm} in (50,130)#circ; E_{beam} [GeV]; |t| [GeV^{2}]",46, 2.8,12, 17,3.,20.);
    TCanvas* C6 = new TCanvas("C6","Et",500,500);
    C6->Divide(2,1);
    C6->cd(1);
    T->Draw("TMath::Abs(t):Original_E>>Et_MF_All",wMF,"col2z");
    C6->cd(2);
    T->Draw("TMath::Abs(t):Original_E>>Et_SRC_All",wSRC,"col2z");
     
    // divide by 5 for 5 nuclear targets 'sharing' the beam time
    
    /*
    TCanvas TmpC;
    T->Draw("1","weight*(TMath::Abs(t)>3&&TMath::Abs(u)>3&&Precoil<0.1)/5.","h");
    T->Draw("1","weight*(TMath::Abs(t)>3&&TMath::Abs(u)>3&&Precoil<0.1&&TMath::Abs(theta_cm-90)<10)/5.","h same");

    TCanvas TmpC2;
    T->Draw("1"  ,"weight*(TMath::Abs(t)>3&&TMath::Abs(u)>3&&Precoil>0.1&&theta_recoil<160&&Precoil>0.3                            )*0.5*0.8/5.","h");
    T->Draw("1"  ,"weight*(TMath::Abs(t)>3&&TMath::Abs(u)>3&&Precoil>0.1&&theta_recoil<160&&Precoil>0.3&&TMath::Abs(theta_cm-90)<10)*0.5*0.8/5.","h same");
    */
  
    TCanvas* TmpC1 = new TCanvas("TmpC1","theta_recoil",500,500);
    TH1F* hth_rec = new TH1F("hth_rec","; #theta_{recoil} [#circ]; counts",18,0,180);
    TmpC1->cd();
    T->Draw("theta_recoil>>hth_rec"  ,wSRCcoh,"h");
    //    T->Draw("theta_recoil"  ,wSRCcoh,"samelh");

    TCanvas* TmpC3 = new TCanvas("TmpC3","t for all theta range",800,600);
    TH1F* ht = new TH1F("ht","; |t| [GeV^{2}]; counts",10,1.,21.);
    ht->SetTitle("");
    ht->SetLineColor(8);
    //    ht->GetYaxis()->SetRangeUser(0.1,500.);
    TH1F* ht2 = new TH1F("ht2","; |t| [GeV^{2}]; counts",10,1.,21.);
    ht2->SetLineColor(4);
    TmpC3->cd();
    gPad->SetLogy();
    T->Draw("TMath::Abs(t)>>ht"  ,wMFcoh,"h");
    //    T->Draw("TMath::Abs(t)>>ht2" ,wSRCcoh,"hsame");

    cout<<"# t [(GeV/c)^2]     counts MF"<<endl;
    for(int j=0; j<10; j++){
      cout<<ht->GetBinCenter(j+1)<<" "<<ht->GetBinContent(j+1)<<endl;
    }

    TLegend* leg1 = new TLegend(0.5,0.5,0.75,0.75);
    leg1->SetLineColor(0);
    leg1->SetFillColor(0);
    leg1->SetTextSize(0.04);
    leg1->AddEntry(ht,"MF, #theta_{cm} in (50, 130)#circ","l");
    leg1->AddEntry(ht2, "SRC, #theta_{cm} in (50, 130)#circ","l");
    leg1->Draw("same");
    
    TCanvas* Tcm = new TCanvas("Tcm","th_cm",500,500);
    TH1F* hth_cm = new TH1F("hth_cm","; #theta_{cm} [#circ]; counts",12,30.,150.);
    hth_cm->SetLineWidth(2);
    hth_cm->SetLineColor(8);
    Tcm->cd();
    T->Draw("theta_cm>>hth_cm",wMFcoh,"h");
    TH1F* hth_cm1 = new TH1F("hth_cm1","; #theta_{cm} [#circ]; counts",12,30.,150.);
    hth_cm1->SetLineWidth(2);
    hth_cm1->SetLineColor(4);
    T->Draw("theta_cm>>hth_cm1",wSRCcoh,"sameh");
    leg1->Draw("same");

    cout<<"# theta_cm [degrees]     counts MF      counts SRC"<<endl;
    for(int j=0; j<12; j++){
      cout<<hth_cm->GetBinCenter(j+1)<<" "<<hth_cm->GetBinContent(j+1)<<" "<<hth_cm1->GetBinContent(j+1)<<endl;
    }

    //sanity checks:
    TCanvas* sany = new TCanvas("sany","sanity checks",900,300);
    sany->Divide(3,2);
    sany->cd(1);
    T->Draw("Precoil:theta_recoil"  ,wst,"col2z");
    sany->cd(2);
    T->Draw("P_miss:theta_miss"  ,wst,"col2z");
    sany->cd(3);
    T->Draw("Precoil:P_miss"  ,wst,"col2z");
    sany->cd(4);
    T->Draw("theta_cm:t",wst,"col2z");


    sany->cd(5);
    T->Draw("P_miss",wst,"hist");
    sany->cd(6);
    T->Draw("P_miss");
    
    theta_P_p_MF_All->GetXaxis()->SetTitleSize(0.06);
    theta_P_p_MF_All->GetYaxis()->SetTitleSize(0.06);
    theta_P_pi_MF_All->GetXaxis()->SetTitleSize(0.06);
    theta_P_pi_MF_All->GetYaxis()->SetTitleSize(0.06);
    theta_Ppi_MF_All->GetXaxis()->SetTitleSize(0.06);
    theta_Ppi_MF_All->GetYaxis()->SetTitleSize(0.06);
    phi_MF_All->GetXaxis()->SetTitleSize(0.06);
    
    theta_P_p_MF_All->GetXaxis()->SetTitleOffset(0.8);
    theta_P_p_MF_All->GetYaxis()->SetTitleOffset(0.8);
    theta_P_pi_MF_All->GetXaxis()->SetTitleOffset(0.8);
    theta_P_pi_MF_All->GetYaxis()->SetTitleOffset(0.8);
    theta_Ppi_MF_All->GetXaxis()->SetTitleOffset(0.8);
    theta_Ppi_MF_All->GetYaxis()->SetTitleOffset(0.8);
    phi_MF_All->GetXaxis()->SetTitleOffset(0.8);
    
    theta_P_p_MF_All->GetXaxis()->SetLabelSize(0.06);
    theta_P_p_MF_All->GetYaxis()->SetLabelSize(0.06);
    theta_P_pi_MF_All->GetXaxis()->SetLabelSize(0.06);
    theta_P_pi_MF_All->GetYaxis()->SetLabelSize(0.06);
    theta_Ppi_MF_All->GetXaxis()->SetLabelSize(0.06);
    theta_Ppi_MF_All->GetYaxis()->SetLabelSize(0.06);
    phi_MF_All->GetXaxis()->SetLabelSize(0.06);
    

    theta_P_p_SRC_All->GetXaxis()->SetTitleSize(0.06);
    theta_P_p_SRC_All->GetYaxis()->SetTitleSize(0.06);
    theta_P_pi_SRC_All->GetXaxis()->SetTitleSize(0.06);
    theta_P_pi_SRC_All->GetYaxis()->SetTitleSize(0.06);
    theta_Ppi_SRC_All->GetXaxis()->SetTitleSize(0.06);
    theta_Ppi_SRC_All->GetYaxis()->SetTitleSize(0.06);
    phi_SRC_All->GetXaxis()->SetTitleSize(0.06);
    
    theta_P_p_SRC_All->GetXaxis()->SetTitleOffset(0.8);
    theta_P_p_SRC_All->GetYaxis()->SetTitleOffset(0.8);
    theta_P_pi_SRC_All->GetXaxis()->SetTitleOffset(0.8);
    theta_P_pi_SRC_All->GetYaxis()->SetTitleOffset(0.8);
    theta_Ppi_SRC_All->GetXaxis()->SetTitleOffset(0.8);
    theta_Ppi_SRC_All->GetYaxis()->SetTitleOffset(0.8);
    phi_SRC_All->GetXaxis()->SetTitleOffset(0.8);
    
    theta_P_p_SRC_All->GetXaxis()->SetLabelSize(0.06);
    theta_P_p_SRC_All->GetYaxis()->SetLabelSize(0.06);
    theta_P_pi_SRC_All->GetXaxis()->SetLabelSize(0.06);
    theta_P_pi_SRC_All->GetYaxis()->SetLabelSize(0.06);
    theta_Ppi_SRC_All->GetXaxis()->SetLabelSize(0.06);
    theta_Ppi_SRC_All->GetYaxis()->SetLabelSize(0.06);
    phi_SRC_All->GetXaxis()->SetLabelSize(0.06);   
  
}













