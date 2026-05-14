{
 
  TFile* data = new TFile("/P/stan/scratch/ksaldan/fcal_timing/test_p3pi/p3pi_FCALStudy_2018-01_fitData.root","read");
  TFile* dalt = new TFile("/P/stan/scratch/ksaldan/fcal_timing/test_p3pi/p3pi_FCALStudy_2018-01_oldSmear_preFiducial_fit.root","read"); 

  TGraphErrors* data_m1_E = data->Get<TGraphErrors>("gr_m1_effic_Ebins");  
  TGraphErrors* data_m2_E = data->Get<TGraphErrors>("gr_m2_effic_Ebins");  
  TGraphErrors* dalt_m1_E = dalt->Get<TGraphErrors>("gr_m1_effic_Ebins");  
  TGraphErrors* dalt_m2_E = dalt->Get<TGraphErrors>("gr_m2_effic_Ebins");  

  TGraphErrors* data_m1_Theta = data->Get<TGraphErrors>("gr_m1_effic_Thetabins");
  TGraphErrors* data_m2_Theta = data->Get<TGraphErrors>("gr_m2_effic_Thetabins");
  TGraphErrors* dalt_m1_Theta = dalt->Get<TGraphErrors>("gr_m1_effic_Thetabins");
  TGraphErrors* dalt_m2_Theta = dalt->Get<TGraphErrors>("gr_m2_effic_Thetabins");

  TGraphErrors* data_m1_Phi = data->Get<TGraphErrors>("gr_m1_effic_Phibins");
  TGraphErrors* data_m2_Phi = data->Get<TGraphErrors>("gr_m2_effic_Phibins");
  TGraphErrors* dalt_m1_Phi = dalt->Get<TGraphErrors>("gr_m1_effic_Phibins");
  TGraphErrors* dalt_m2_Phi = dalt->Get<TGraphErrors>("gr_m2_effic_Phibins");
  
  TMultiGraph* mgE1 = new TMultiGraph();
  TMultiGraph* mgT1 = new TMultiGraph();
  TMultiGraph* mgP1 = new TMultiGraph();

  mgE1->SetTitle(" ;E_{#gamma} GeV; Efficiency"); 
  mgT1->SetTitle(" ;Missing Photon Polar Angle (degrees); Efficiency");
  mgP1->SetTitle(" ; Missing Photon Azimunthal Angle (degrees)");
  mgE1->GetXaxis()->SetLimits(0.0,5.0);
  mgE1->GetYaxis()->SetRangeUser(0.25,1.1);
  mgT1->GetXaxis()->SetLimits(0.0,14.0);
  mgT1->GetYaxis()->SetRangeUser(0.25,1.1);
  mgP1->GetXaxis()->SetLimits(-180,180);
  mgP1->GetYaxis()->SetRangeUser(0.25,1.1);
  
  data_m1_E->RemovePoint(0);
  data_m2_E->RemovePoint(0);
  dalt_m1_E->RemovePoint(0);
  dalt_m2_E->RemovePoint(0);
  data_m1_E->SetLineColor(kRed); 
  dalt_m1_E->SetLineColor(kViolet); 
  data_m2_E->SetLineColor(kOrange-3); 
  dalt_m2_E->SetLineColor(kGreen+2); 
  data_m1_Theta->SetLineColor(kRed);
  dalt_m1_Theta->SetLineColor(kViolet);
  data_m2_Theta->SetLineColor(kOrange-3);
  dalt_m2_Theta->SetLineColor(kGreen+2); 
  data_m1_Phi->SetLineColor(kRed);
  dalt_m1_Phi->SetLineColor(kViolet);
  data_m2_Phi->SetLineColor(kOrange-3);
  dalt_m2_Phi->SetLineColor(kGreen+2);

  data_m1_E->SetMarkerColor(kRed);
  dalt_m1_E->SetMarkerColor(kViolet);
  data_m2_E->SetMarkerColor(kOrange-3);
  dalt_m2_E->SetMarkerColor(kGreen+2);
  data_m1_Theta->SetMarkerColor(kRed);
  dalt_m1_Theta->SetMarkerColor(kViolet);
  data_m2_Theta->SetMarkerColor(kOrange-3);
  dalt_m2_Theta->SetMarkerColor(kGreen+2);
  data_m1_Phi->SetMarkerColor(kRed);
  dalt_m1_Phi->SetMarkerColor(kViolet);
  data_m2_Phi->SetMarkerColor(kOrange-3);
  dalt_m2_Phi->SetMarkerColor(kGreen+2);

  data_m1_E->SetLineWidth(1);
  data_m2_E->SetLineWidth(1);
  dalt_m1_E->SetLineWidth(1);
  dalt_m2_E->SetLineWidth(1);
  data_m1_Theta->SetLineWidth(1);
  data_m2_Theta->SetLineWidth(1);
  dalt_m1_Theta->SetLineWidth(1);
  dalt_m2_Theta->SetLineWidth(1);
  data_m1_Phi->SetLineWidth(1);
  data_m2_Phi->SetLineWidth(1);
  dalt_m1_Phi->SetLineWidth(1);
  dalt_m2_Phi->SetLineWidth(1);
  data_m1_E->SetMarkerSize(1.2);
  data_m2_E->SetMarkerSize(1.2);
  data_m1_Theta->SetMarkerSize(1.2);
  data_m2_Theta->SetMarkerSize(1.2);
  data_m1_Phi->SetMarkerSize(1.2);
  data_m2_Phi->SetMarkerSize(1.2);
  dalt_m1_E->SetMarkerSize(1.2);
  dalt_m2_E->SetMarkerSize(1.2);
  dalt_m1_Theta->SetMarkerSize(1.2);
  dalt_m2_Theta->SetMarkerSize(1.2);
  dalt_m1_Phi->SetMarkerSize(1.2);
  dalt_m2_Phi->SetMarkerSize(1.2);

  data_m1_E->SetMarkerStyle(8);
  dalt_m1_E->SetMarkerStyle(25);
  data_m2_E->SetMarkerStyle(8);
  dalt_m2_E->SetMarkerStyle(25);
  data_m1_Theta->SetMarkerStyle(8);
  dalt_m1_Theta->SetMarkerStyle(25);
  data_m2_Theta->SetMarkerStyle(8);
  dalt_m2_Theta->SetMarkerStyle(25);
  data_m1_Phi->SetMarkerStyle(8);
  dalt_m1_Phi->SetMarkerStyle(25);
  data_m2_Phi->SetMarkerStyle(8);
  dalt_m2_Phi->SetMarkerStyle(25);


  TLegend* l1 = new TLegend(0.3,0.2,0.73,0.45);
  l1->AddEntry(data_m1_E,"Spring 2018 Data Method 1");
  l1->AddEntry(data_m2_E,"Spring 2018 Data Method 2");
  l1->AddEntry(dalt_m1_E,"MC Method 1");
  l1->AddEntry(dalt_m2_E,"MC Method 2");

  mgE1->Add(data_m1_E,"AP");
  mgE1->Add(dalt_m1_E,"AP");
  mgE1->Add(data_m2_E,"AP");
  mgE1->Add(dalt_m2_E,"AP");
  
  mgT1->Add(data_m1_Theta,"AP");
  mgT1->Add(dalt_m1_Theta,"AP");
  mgT1->Add(data_m2_Theta,"AP");
  mgT1->Add(dalt_m2_Theta,"AP");

  mgP1->Add(data_m1_Phi,"AP");
  mgP1->Add(dalt_m1_Phi,"AP");
  mgP1->Add(data_m2_Phi,"AP");
  mgP1->Add(dalt_m2_Phi,"AP");
  
  TCanvas* c1 = new TCanvas("c1","c1");
  TCanvas* c2 = new TCanvas("c2","c2");
  TCanvas* c3 = new TCanvas("c3","c3");
  
  c1->cd();  mgE1->Draw("AP");  l1->Draw("same");
  c2->cd();  mgT1->Draw("AP");  l1->Draw("same");
  c3->cd();  mgP1->Draw("AP");  l1->Draw("same");

}
 
