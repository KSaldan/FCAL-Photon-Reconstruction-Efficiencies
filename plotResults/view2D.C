{
 
  //TFile* data = new TFile("/P/stan/scratch/ksaldan/fcal_timing/test_p3pi/EBinFits/p3pi_FCALStudy_2017-01_42_fitEBins.root","read");
  //TFile* oldS = new TFile("/P/stan/scratch/ksaldan/fcal_timing/test_p3pi/EBinFits/p3pi_FCALStudy_2017-01_Fiducial_fitEBin.root","read"); 
  //TFile* data = new TFile("/P/stan/scratch/ksaldan/fcal_timing/test_p3pi/EBinFits/p3pi_FCALStudy_2018-01_Data_fitEBins.root","read");
  //TFile* oldS = new TFile("/P/stan/scratch/ksaldan/fcal_timing/test_p3pi/EBinFits/p3pi_FCALStudy_2018-01_oldSmear_Fiducial_fitEBins.root","read");
  //TFile* data = new TFile("/P/stan/scratch/ksaldan/fcal_timing/test_p3pi/EBinFits/p3pi_FCALStudy_2018-08_Data_fitEBins.root","read");
  //TFile* oldS = new TFile("/P/stan/scratch/ksaldan/fcal_timing/test_p3pi/EBinFits/p3pi_FCALStudy_2018-08_Fiducial_fitEBins.root","read");
  TFile* data = new TFile("/P/stan/scratch/ksaldan/fcal_timing/test_p3pi/EBinFits/p3pi_FCALStudy_2019-11_Data_fitEBins.root","read");
  TFile* oldS = new TFile("/P/stan/scratch/ksaldan/fcal_timing/test_p3pi/EBinFits/p3pi_FCALStudy_2019-11_oldSmear_Fiducial_fitEBins.root","read");

  TFile* fileOne = new TFile("FCAL2DPlots.root","update");

  TGraphErrors* data_m1_Bin1_Theta = data->Get<TGraphErrors>("gr_m1_Bin1_effic_Thetabins");
  TGraphErrors* data_m2_Bin1_Theta = data->Get<TGraphErrors>("gr_m2_Bin1_effic_Thetabins");
  TGraphErrors* oldS_m1_Bin1_Theta = oldS->Get<TGraphErrors>("gr_m1_Bin1_effic_Thetabins");
  TGraphErrors* oldS_m2_Bin1_Theta = oldS->Get<TGraphErrors>("gr_m2_Bin1_effic_Thetabins");
  
  Int_t TPoints = data_m1_Bin1_Theta->GetN();

  TGraphErrors* data_m1_Bin2_Theta = data->Get<TGraphErrors>("gr_m1_Bin2_effic_Thetabins");
  TGraphErrors* data_m2_Bin2_Theta = data->Get<TGraphErrors>("gr_m2_Bin2_effic_Thetabins");
  TGraphErrors* oldS_m1_Bin2_Theta = oldS->Get<TGraphErrors>("gr_m1_Bin2_effic_Thetabins");
  TGraphErrors* oldS_m2_Bin2_Theta = oldS->Get<TGraphErrors>("gr_m2_Bin2_effic_Thetabins");
    
  TGraphErrors* data_m1_Bin3_Theta = data->Get<TGraphErrors>("gr_m1_Bin3_effic_Thetabins");
  TGraphErrors* data_m2_Bin3_Theta = data->Get<TGraphErrors>("gr_m2_Bin3_effic_Thetabins");
  TGraphErrors* oldS_m1_Bin3_Theta = oldS->Get<TGraphErrors>("gr_m1_Bin3_effic_Thetabins");
  TGraphErrors* oldS_m2_Bin3_Theta = oldS->Get<TGraphErrors>("gr_m2_Bin3_effic_Thetabins");

  TGraphErrors* data_m1_Bin4_Theta = data->Get<TGraphErrors>("gr_m1_Bin4_effic_Thetabins");
  TGraphErrors* data_m2_Bin4_Theta = data->Get<TGraphErrors>("gr_m2_Bin4_effic_Thetabins");
  TGraphErrors* oldS_m1_Bin4_Theta = oldS->Get<TGraphErrors>("gr_m1_Bin4_effic_Thetabins");
  TGraphErrors* oldS_m2_Bin4_Theta = oldS->Get<TGraphErrors>("gr_m2_Bin4_effic_Thetabins");

  TGraphErrors* data_m1_Bin5_Theta = data->Get<TGraphErrors>("gr_m1_Bin5_effic_Thetabins");
  TGraphErrors* data_m2_Bin5_Theta = data->Get<TGraphErrors>("gr_m2_Bin5_effic_Thetabins");
  TGraphErrors* oldS_m1_Bin5_Theta = oldS->Get<TGraphErrors>("gr_m1_Bin5_effic_Thetabins");
  TGraphErrors* oldS_m2_Bin5_Theta = oldS->Get<TGraphErrors>("gr_m2_Bin5_effic_Thetabins");
    
  TGraphErrors* data_m1_Bin6_Theta = data->Get<TGraphErrors>("gr_m1_Bin6_effic_Thetabins");
  TGraphErrors* data_m2_Bin6_Theta = data->Get<TGraphErrors>("gr_m2_Bin6_effic_Thetabins");
  TGraphErrors* oldS_m1_Bin6_Theta = oldS->Get<TGraphErrors>("gr_m1_Bin6_effic_Thetabins");
  TGraphErrors* oldS_m2_Bin6_Theta = oldS->Get<TGraphErrors>("gr_m2_Bin6_effic_Thetabins");

  TGraphErrors* data_m1_Bin7_Theta = data->Get<TGraphErrors>("gr_m1_Bin7_effic_Thetabins");
  TGraphErrors* data_m2_Bin7_Theta = data->Get<TGraphErrors>("gr_m2_Bin7_effic_Thetabins");
  TGraphErrors* oldS_m1_Bin7_Theta = oldS->Get<TGraphErrors>("gr_m1_Bin7_effic_Thetabins");
  TGraphErrors* oldS_m2_Bin7_Theta = oldS->Get<TGraphErrors>("gr_m2_Bin7_effic_Thetabins");

  TGraphErrors* data_m1_Bin8_Theta = data->Get<TGraphErrors>("gr_m1_Bin8_effic_Thetabins");
  TGraphErrors* data_m2_Bin8_Theta = data->Get<TGraphErrors>("gr_m2_Bin8_effic_Thetabins");
  TGraphErrors* oldS_m1_Bin8_Theta = oldS->Get<TGraphErrors>("gr_m1_Bin8_effic_Thetabins");
  TGraphErrors* oldS_m2_Bin8_Theta = oldS->Get<TGraphErrors>("gr_m2_Bin8_effic_Thetabins");

  Double_t B1M2Y[40],B1M1Y[40],B2M2Y[40],B2M1Y[40],B3M2Y[40],B3M1Y[40],B4M2Y[40],B4M1Y[40],B5M2Y[40],B5M1Y[40],B6M2Y[40],B6M1Y[40],B7M2Y[40],B7M1Y[40],B8M2Y[40],B8M1Y[40];
  Double_t B1M2E[40],B1M1E[40],B2M2E[40],B2M1E[40],B3M2E[40],B3M1E[40],B4M2E[40],B4M1E[40],B5M2E[40],B5M1E[40],B6M2E[40],B6M1E[40],B7M2E[40],B7M1E[40],B8M2E[40],B8M1E[40];


  for(int j=0; j<TPoints; j++){
    
    if(oldS_m1_Bin1_Theta->GetY()[j]==0 || data_m1_Bin1_Theta->GetY()[j]==0){
      B1M1Y[j]=-100; B1M1E[j]=-100;}
    else{ 
      B1M1Y[j] = (oldS_m1_Bin1_Theta->GetY()[j]-data_m1_Bin1_Theta->GetY()[j])/data_m1_Bin1_Theta->GetY()[j];
      B1M1E[j] = sqrt(pow(oldS_m1_Bin1_Theta->GetEY()[j]/oldS_m1_Bin1_Theta->GetY()[j],2)+4.*pow(data_m1_Bin1_Theta->GetEY()[j]/data_m1_Bin1_Theta->GetY()[j],2))*B1M1Y[j];
    }
    if(oldS_m2_Bin1_Theta->GetY()[j]==0 || data_m2_Bin1_Theta->GetY()[j]==0){
      B1M2Y[j]=-100; B1M2E[j]=-100;}
    else{
      B1M2Y[j] = (oldS_m2_Bin1_Theta->GetY()[j]-data_m2_Bin1_Theta->GetY()[j])/data_m2_Bin1_Theta->GetY()[j];
      B1M2E[j] = sqrt(pow(oldS_m2_Bin1_Theta->GetEY()[j]/oldS_m2_Bin1_Theta->GetY()[j],2)+4.*pow(data_m2_Bin1_Theta->GetEY()[j]/data_m2_Bin1_Theta->GetY()[j],2))*B1M2Y[j]; 
    }
    if(oldS_m1_Bin2_Theta->GetY()[j]==0 || data_m1_Bin2_Theta->GetY()[j]==0){
      B2M1Y[j]=-100; B2M1E[j]=-100;}
    else{ 
      B2M1Y[j] = (oldS_m1_Bin2_Theta->GetY()[j]-data_m1_Bin2_Theta->GetY()[j])/data_m1_Bin2_Theta->GetY()[j];
      B2M1E[j] = sqrt(pow(oldS_m1_Bin2_Theta->GetEY()[j]/oldS_m1_Bin2_Theta->GetY()[j],2)+4.*pow(data_m1_Bin2_Theta->GetEY()[j]/data_m1_Bin2_Theta->GetY()[j],2))*B2M1Y[j];
    }
    if(oldS_m2_Bin2_Theta->GetY()[j]==0 || data_m2_Bin2_Theta->GetY()[j]==0){
      B2M2Y[j]=-100; B2M2E[j]=-100;}
    else{
      B2M2Y[j] = (oldS_m2_Bin2_Theta->GetY()[j]-data_m2_Bin2_Theta->GetY()[j])/data_m2_Bin2_Theta->GetY()[j];
      B2M2E[j] = sqrt(pow(oldS_m2_Bin2_Theta->GetEY()[j]/oldS_m2_Bin2_Theta->GetY()[j],2)+4.*pow(data_m2_Bin2_Theta->GetEY()[j]/data_m2_Bin2_Theta->GetY()[j],2))*B2M2Y[j];
    }
    if(oldS_m1_Bin3_Theta->GetY()[j]==0 || data_m1_Bin3_Theta->GetY()[j]==0){
      B3M1Y[j]=-100; B3M1E[j]=-100;}
    else{ 
      B3M1Y[j] = (oldS_m1_Bin3_Theta->GetY()[j]-data_m1_Bin3_Theta->GetY()[j])/data_m1_Bin3_Theta->GetY()[j];
      B3M1E[j] = sqrt(pow(oldS_m1_Bin3_Theta->GetEY()[j]/oldS_m1_Bin3_Theta->GetY()[j],2)+4.*pow(data_m1_Bin3_Theta->GetEY()[j]/data_m1_Bin3_Theta->GetY()[j],2))*B3M1Y[j];
    }
    if(oldS_m2_Bin3_Theta->GetY()[j]==0 || data_m2_Bin3_Theta->GetY()[j]==0){
      B3M2Y[j]=-100; B3M2E[j]=-100;}
    else{
      B3M2Y[j] = (oldS_m2_Bin3_Theta->GetY()[j]-data_m2_Bin3_Theta->GetY()[j])/data_m2_Bin3_Theta->GetY()[j];
      B3M2E[j] = sqrt(pow(oldS_m2_Bin3_Theta->GetEY()[j]/oldS_m2_Bin3_Theta->GetY()[j],2)+4.*pow(data_m2_Bin3_Theta->GetEY()[j]/data_m2_Bin3_Theta->GetY()[j],2))*B3M2Y[j];
    }
    if(oldS_m1_Bin4_Theta->GetY()[j]==0 || data_m1_Bin4_Theta->GetY()[j]==0){
      B4M1Y[j]=-100; B4M1E[j]=-100;}
    else{ 
      B4M1Y[j] = (oldS_m1_Bin4_Theta->GetY()[j]-data_m1_Bin4_Theta->GetY()[j])/data_m1_Bin4_Theta->GetY()[j];
      B4M1E[j] = sqrt(pow(oldS_m1_Bin4_Theta->GetEY()[j]/oldS_m1_Bin4_Theta->GetY()[j],2)+4.*pow(data_m1_Bin4_Theta->GetEY()[j]/data_m1_Bin4_Theta->GetY()[j],2))*B4M1Y[j];
    }
    if(oldS_m2_Bin4_Theta->GetY()[j]==0 || data_m2_Bin4_Theta->GetY()[j]==0){
      B4M2Y[j]=-100; B4M2E[j]=-100;}
    else{
      B4M2Y[j] = (oldS_m2_Bin4_Theta->GetY()[j]-data_m2_Bin4_Theta->GetY()[j])/data_m2_Bin4_Theta->GetY()[j];
      B4M2E[j] = sqrt(pow(oldS_m2_Bin4_Theta->GetEY()[j]/oldS_m2_Bin4_Theta->GetY()[j],2)+4.*pow(data_m2_Bin4_Theta->GetEY()[j]/data_m2_Bin4_Theta->GetY()[j],2))*B4M2Y[j];
    }
    if(oldS_m1_Bin5_Theta->GetY()[j]==0 || data_m1_Bin5_Theta->GetY()[j]==0){
      B5M1Y[j]=-100; B5M1E[j]=-100;}
    else{ 
      B5M1Y[j] = (oldS_m1_Bin5_Theta->GetY()[j]-data_m1_Bin5_Theta->GetY()[j])/data_m1_Bin5_Theta->GetY()[j];
      B5M1E[j] = sqrt(pow(oldS_m1_Bin5_Theta->GetEY()[j]/oldS_m1_Bin5_Theta->GetY()[j],2)+4.*pow(data_m1_Bin5_Theta->GetEY()[j]/data_m1_Bin5_Theta->GetY()[j],2))*B5M1Y[j];
    }
    if(oldS_m2_Bin5_Theta->GetY()[j]==0 || data_m2_Bin5_Theta->GetY()[j]==0){
      B5M2Y[j]=-100; B2M2E[j]=-100;}
    else{
      B5M2Y[j] = (oldS_m2_Bin5_Theta->GetY()[j]-data_m2_Bin5_Theta->GetY()[j])/data_m2_Bin5_Theta->GetY()[j];
      B5M2E[j] = sqrt(pow(oldS_m2_Bin5_Theta->GetEY()[j]/oldS_m2_Bin5_Theta->GetY()[j],2)+4.*pow(data_m2_Bin5_Theta->GetEY()[j]/data_m2_Bin5_Theta->GetY()[j],2))*B5M2Y[j];
    }
    if(oldS_m1_Bin6_Theta->GetY()[j]==0 || data_m1_Bin6_Theta->GetY()[j]==0){
      B6M1Y[j]=-100; B6M1E[j]=-100;}
    else{ 
      B6M1Y[j] = (oldS_m1_Bin6_Theta->GetY()[j]-data_m1_Bin6_Theta->GetY()[j])/data_m1_Bin6_Theta->GetY()[j];
      B6M1E[j] = sqrt(pow(oldS_m1_Bin6_Theta->GetEY()[j]/oldS_m1_Bin6_Theta->GetY()[j],2)+4.*pow(data_m1_Bin6_Theta->GetEY()[j]/data_m1_Bin6_Theta->GetY()[j],2))*B6M1Y[j];
    }
    if(oldS_m2_Bin6_Theta->GetY()[j]==0 || data_m2_Bin6_Theta->GetY()[j]==0){
      B6M2Y[j]=-100; B6M2E[j]=-100;}
    else{
      B6M2Y[j] = (oldS_m2_Bin6_Theta->GetY()[j]-data_m2_Bin6_Theta->GetY()[j])/data_m2_Bin6_Theta->GetY()[j];
      B6M2E[j] = sqrt(pow(oldS_m2_Bin6_Theta->GetEY()[j]/oldS_m2_Bin6_Theta->GetY()[j],2)+4.*pow(data_m2_Bin6_Theta->GetEY()[j]/data_m2_Bin6_Theta->GetY()[j],2))*B6M2Y[j];
    }
    if(oldS_m1_Bin7_Theta->GetY()[j]==0 || data_m1_Bin7_Theta->GetY()[j]==0){
      B7M1Y[j]=-100; B7M1E[j]=-100;}
    else{ 
      B7M1Y[j] = (oldS_m1_Bin7_Theta->GetY()[j]-data_m1_Bin7_Theta->GetY()[j])/data_m1_Bin7_Theta->GetY()[j];
      B7M1E[j] = sqrt(pow(oldS_m1_Bin7_Theta->GetEY()[j]/oldS_m1_Bin7_Theta->GetY()[j],2)+4.*pow(data_m1_Bin7_Theta->GetEY()[j]/data_m1_Bin7_Theta->GetY()[j],2))*B7M1Y[j];
    }
    if(oldS_m2_Bin7_Theta->GetY()[j]==0 || data_m2_Bin7_Theta->GetY()[j]==0){
      B7M2Y[j]=-100; B7M2E[j]=-100;}
    else{
      B7M2Y[j] = (oldS_m2_Bin7_Theta->GetY()[j]-data_m2_Bin7_Theta->GetY()[j])/data_m2_Bin7_Theta->GetY()[j];
      B7M2E[j] = sqrt(pow(oldS_m2_Bin7_Theta->GetEY()[j]/oldS_m2_Bin7_Theta->GetY()[j],2)+4.*pow(data_m2_Bin7_Theta->GetEY()[j]/data_m2_Bin7_Theta->GetY()[j],2))*B7M2Y[j];
    }
    if(oldS_m1_Bin8_Theta->GetY()[j]==0 || data_m1_Bin8_Theta->GetY()[j]==0){
      B8M1Y[j]=-100; B8M1E[j]=-100;}
    else{ 
      B8M1Y[j] = (oldS_m1_Bin8_Theta->GetY()[j]-data_m1_Bin8_Theta->GetY()[j])/data_m1_Bin8_Theta->GetY()[j];
      B8M1E[j] = sqrt(pow(oldS_m1_Bin8_Theta->GetEY()[j]/oldS_m1_Bin8_Theta->GetY()[j],2)+4.*pow(data_m1_Bin8_Theta->GetEY()[j]/data_m1_Bin8_Theta->GetY()[j],2))*B8M1Y[j]; 
    }
    if(oldS_m2_Bin8_Theta->GetY()[j]==0 || data_m2_Bin8_Theta->GetY()[j]==0){
      B8M2Y[j]=-100; B8M2E[j]=-100;}
    else{
      B8M2Y[j] = (oldS_m2_Bin8_Theta->GetY()[j]-data_m2_Bin8_Theta->GetY()[j])/data_m2_Bin8_Theta->GetY()[j];
      B8M2E[j] = sqrt(pow(oldS_m2_Bin8_Theta->GetEY()[j]/oldS_m2_Bin8_Theta->GetY()[j],2)+4.*pow(data_m2_Bin8_Theta->GetEY()[j]/data_m2_Bin8_Theta->GetY()[j],2))*B8M2Y[j];
    }
/*
    cout << B1M1Y[j] << endl;
    cout << B1M2Y[j] << endl;
    cout << B2M1Y[j] << endl;
    cout << B2M2Y[j] << endl;
    cout << B3M1Y[j] << endl;
    cout << B3M2Y[j] << endl;
    cout << B4M1Y[j] << endl;
    cout << B4M2Y[j] << endl;
    cout << B5M1Y[j] << endl;
    cout << B5M2Y[j] << endl;
    cout << B6M1Y[j] << endl;
    cout << B6M2Y[j] << endl;
    cout << B7M1Y[j] << endl;
    cout << B7M2Y[j] << endl;
    cout << B8M1Y[j] << endl;
    cout << B8M2Y[j] << "\n" << endl;
  */
  }
  
  TH2F* M1Plot = new TH2F("M1Plot","",18,0,9,7,0.6,4.8);
  TH2F* M2Plot = new TH2F("M2Plot","",18,0,9,7,0.6,4.8);

  for(int i=1; i<8;i++){
    for(int j=1; j<19; j++){
      if(i==1){
        M1Plot->SetBinContent(j,i,(B1M1Y[j-1]));
        M1Plot->SetBinError(j,i,(B1M1E[j-1]));
        M2Plot->SetBinContent(j,i,(B1M2Y[j-1]));
        M2Plot->SetBinError(j,i,(B1M2E[j-1]));
      }
      if(i==2){
        M1Plot->SetBinContent(j,i,(B2M1Y[j-1]));
        M1Plot->SetBinError(j,i,(B2M1E[j-1]));
        M2Plot->SetBinContent(j,i,(B2M2Y[j-1]));
        M2Plot->SetBinError(j,i,(B2M2E[j-1]));
      }
      if(i==3){
        M1Plot->SetBinContent(j,i,(B3M1Y[j-1]));
        M1Plot->SetBinError(j,i,(B3M1E[j-1]));
        M2Plot->SetBinContent(j,i,(B3M2Y[j-1]));
        M2Plot->SetBinError(j,i,(B3M2E[j-1]));
      }
      if(i==4){
        M1Plot->SetBinContent(j,i,(B4M1Y[j-1]));
        M1Plot->SetBinError(j,i,(B4M1E[j-1]));
        M2Plot->SetBinContent(j,i,(B4M2Y[j-1]));
        M2Plot->SetBinError(j,i,(B4M2E[j-1]));
      }
      if(i==5){
        M1Plot->SetBinContent(j,i,(B5M1Y[j-1]));
        M1Plot->SetBinError(j,i,(B5M1E[j-1]));
        M2Plot->SetBinContent(j,i,(B5M2Y[j-1]));
        M2Plot->SetBinError(j,i,(B5M2E[j-1]));
      }
      if(i==6){
        M1Plot->SetBinContent(j,i,(B6M1Y[j-1]));
        M1Plot->SetBinError(j,i,(B6M1E[j-1]));
        M2Plot->SetBinContent(j,i,(B6M2Y[j-1]));
        M2Plot->SetBinError(j,i,(B6M2E[j-1]));
      }
      if(i==7){
        M1Plot->SetBinContent(j,i,(B7M1Y[j-1]));
        M1Plot->SetBinError(j,i,(B7M1E[j-1]));
        M2Plot->SetBinContent(j,i,(B7M2Y[j-1]));
        M2Plot->SetBinError(j,i,(B7M2E[j-1]));
      }
      if(i==8){
        M1Plot->SetBinContent(j,i,(B8M1Y[j-1]));
        M1Plot->SetBinError(j,i,(B8M1E[j-1]));
        M2Plot->SetBinContent(j,i,(B8M2Y[j-1]));
        M2Plot->SetBinError(j,i,(B8M2E[j-1]));
      }
    
    }
  } 
  M1Plot->SetTitle("Method 1; #theta (deg); E_{#gamma} (GeV); ");
  M2Plot->SetTitle("Method 2; #theta (deg); E_{#gamma} (GeV); ");
  M1Plot->SetTitle("2019-11 M1");
  M2Plot->SetTitle("2019-11 M2");
  M1Plot->SetStats(0);
  M2Plot->SetStats(0);
  M1Plot->GetZaxis()->SetRangeUser(-0.05, 0.05);
  M2Plot->GetZaxis()->SetRangeUser(-0.05, 0.05);
  gStyle->SetPalette(kRainBow);


  gStyle->SetPaintTextFormat(".3f"); 

  M1Plot->Write();
  M2Plot->Write();
  fileOne->Close();
}
 
