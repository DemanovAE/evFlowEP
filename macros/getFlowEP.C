#include <string>
#include <iostream>
#include <vector>
#include <TGraph.h>
//#include <Math/SpecFuncMathMore.h>
#include <TMath.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TProfile3D.h>
#include "TFile.h"
#include <TH1.h>
#include <TH2.h>
#include <TPad.h>
#include <algorithm>
#include <TGraphErrors.h>
#include <TF1.h>


double RatioErrors(double x, double xe, double y, double ye){
   return TMath::Sqrt( TMath::Power( xe / y, 2.) + TMath::Power( (x*ye) / (y*y), 2.) );
}

// Resolution calculation
//S----------------------------------------------------------------
double GetRes(double _chi, double _harm)
{
    double con = TMath::Sqrt(TMath::Pi() / 2) / 2;
    double arg = _chi * _chi / 4.;
    double order1 = (_harm - 1) / 2.;
    double order2 = (_harm + 1) / 2.;
    double res = con * _chi * exp(-arg) * (ROOT::Math::cyl_bessel_i(order1, arg) + ROOT::Math::cyl_bessel_i(order2, arg));
    return res;
}

// Chi2 calculation
//S----------------------------------------------------------------
double GetChi(double _res, double _harm, int accuracy)
{
    double chi = 2.0;
    double delta = 1;
    for (int i = 0; i < accuracy; i++)
    {
        if (GetRes(chi, _harm) < _res)
            chi = chi + delta;
        else
            chi = chi - delta;
        delta = delta / 2.;
    }
    return chi;
}
//E-----

std::pair<double,double> GetResolution(std::string _method, double _Res2, double _Res2Err)
{

   double res2    = -999;
   double res     = -999;
   double res_err = -999;
   
   if(_method=="v2_TPC_NS" || _method=="v1_FHCal_NS" || _method=="v2_FHCal_NS")
   {
      res2 = TMath::Abs(_Res2);
      res  = TMath::Sqrt(res2);
      res_err = TMath::Abs(_Res2Err) / ( 2. * res );
   }

   if(_method=="v1_FHCal_RP" || _method=="v2_FHCal_RP")
   {
      res2 = 1.;
      res  =1.;
      res_err = 1.;
   }

   if(_method=="v1_FHCal_F")
   {
      res2 = TMath::Abs(_Res2);
      res  = TMath::Sqrt(res2);
      res_err = TMath::Abs(_Res2Err) / ( 2. * res );
      double _chi = TMath::Sqrt(2)*GetChi(res, 1., 50);
      res = GetRes(_chi,1.);
   }

   if(_method=="v2_FHCal_F")
   {
      res2 = TMath::Abs(_Res2);
      res  = TMath::Sqrt(res2);
      res_err = TMath::Abs(_Res2Err) / ( 2. * res );
      double _chi = TMath::Sqrt(2)*GetChi(res, 2., 50);
      res = GetRes(_chi,2.);
   }


   return {res,res_err};
}


TProfile *TProfileVnvsPt(TProfile2D *const &prVn, double eta_low, double eta_high)
{
  Int_t eta_bin_low = prVn->GetYaxis()->FindBin(eta_low);
  Int_t eta_bin_high = prVn->GetYaxis()->FindBin(eta_high-0.001);
  TProfile *prVnresult = (TProfile *)prVn->ProfileX(Form("%s_eta_%.1f_%.1f",prVn->GetName(),eta_low, eta_high), eta_bin_low, eta_bin_high);
  prVnresult->SetTitle(Form("%s, %.1f<#eta<%.1f;p_{T}, GeV/c;v_{n}", prVn->GetTitle(), eta_low, eta_high));
  return prVnresult;
}

TProfile *TProfileVnvsEta(TProfile2D *const &prVn, double pt_low, double pt_high)
{
  Int_t pt_bin_low = prVn->GetXaxis()->FindBin(pt_low);
  Int_t pt_bin_high = prVn->GetXaxis()->FindBin(pt_high-0.001);
  TProfile *prVnresult = (TProfile *)prVn->ProfileY(Form("%s_pt_%.1f_%.1f",prVn->GetName(),pt_low, pt_high), pt_bin_low, pt_bin_high);
  prVnresult->SetTitle(Form("%s, %.1f<p_{T}<%.1f;#eta;v_{n}", prVn->GetTitle(), pt_low, pt_high));
  return prVnresult;
}

TProfile2D *TProfile2D_VnPtEta(TFile *inFile, std::string _method, std::string _NameFLowProfile, std::string _NameResProfile, std::vector<std::pair<int,int>> _cent, std::string _pid)
{
   
   std::vector<TProfile2D*> fFLow;

   for(auto &CentBin : _cent)
   {
      TProfile2D *TPrCos = (TProfile2D*)inFile->Get(Form("%s_%i%i_%s", _NameFLowProfile.c_str(), CentBin.first,CentBin.second, _pid.c_str()))->Clone();
      if(TPrCos==nullptr){
         std::cout<<"No TProfile2D named "<< Form("%s_%i%i_%s", _NameFLowProfile.c_str(), CentBin.first,CentBin.second, _pid.c_str()) <<" was found in the ROOT file!"<<std::endl;
         return nullptr;
      }

      TProfile *TPrRes = (TProfile*)  inFile->Get(Form("%s_%i%i", _NameResProfile.c_str(),CentBin.first,CentBin.second))->Clone();
      if(TPrRes==nullptr){
         std::cout<<"No TProfile named "<< Form("%s_%i%i", _NameResProfile.c_str(),CentBin.first,CentBin.second) <<" was found in the ROOT file!"<<std::endl;
         return nullptr;
      }

      // Calculate resolution in one CentBin
      std::pair<double,double> fResolution = GetResolution( _method, TPrRes->GetBinContent(TPrRes->FindBin(0.)),TPrRes->GetBinError(TPrRes->FindBin(0.)));
      
      if(fResolution.first==0. || fResolution.first==-999){
         std::cout<<"Error!!! Resolution from "<<TPrRes->GetName()<<" in "<<CentBin.first<<"-"<<CentBin.second<<"%% equals "<<fResolution.first<<std::endl;
         return nullptr;
      }else{
         //Calculate Flow in One CentBin
         TProfile2D *FlowOneBin = (TProfile2D*)TPrCos->Clone();
         //std::cout<<fResolution.first<<std::endl;
         FlowOneBin->Scale(1./fResolution.first);
         fFLow.push_back(FlowOneBin);
      }
   }

   // Some centrality bins hadd into one
   for(int i=1; i<_cent.size(); i++){
      fFLow[0]->Add(fFLow[i]);
   }

   if(fFLow.size()!=0){
      fFLow[0]->SetName(Form("%s_%i%i_%s", _NameFLowProfile.c_str(), _cent[0].first,_cent.back().second, _pid.c_str()));
      fFLow[0]->SetTitle(Form("cent %i-%i",_cent[0].first,_cent.back().second));
      //fFLow[0]->Write();
      return fFLow[0];
   }

   return nullptr;

}


TGraphErrors *vn_diff_EP(TFile *inFile, std::string _method ,std::vector<std::pair<int,int>> _cent, std::pair<std::string, std::vector<double>> ProfileRange, std::string _pid, std::vector<double> AxisXBin){
   
   std::vector<double> x;
   std::vector<double> y;
   std::vector<double> xe;
   std::vector<double> ye;
   TProfile *fResult;
   TProfile2D *PrVn;
   
   if(_method=="v2_TPC_NS")   PrVn = TProfile2D_VnPtEta(inFile, "v2_TPC_NS",   "mhCos2_PhiPsiTPC_sNnS",   "mhCos2TpcNTpcS",     _cent, _pid);
   if(_method=="v2_FHCal_F")  PrVn = TProfile2D_VnPtEta(inFile, "v2_FHCal_F",  "mhCos2_PhiPsiFHCal_sFnF", "mhCos2FHCalSFHCalN", _cent, _pid);
   if(_method=="v2_FHCal_NS") PrVn = TProfile2D_VnPtEta(inFile, "v2_FHCal_NS", "mhCos2_PhiPsiFHCal_sNnS", "mhCos2FHCalSFHCalN", _cent, _pid);
   if(_method=="v2_FHCal_RP") PrVn = TProfile2D_VnPtEta(inFile, "v2_FHCal_RP", "mhCos2_PhiPsiRP_nFsF",    "mhCos2FHCalFPsiRP",  _cent, _pid);
   if(_method=="v1_FHCal_NS") PrVn = TProfile2D_VnPtEta(inFile, "v1_FHCal_NS", "mhCos1_PhiPsiFHCal_sNnS", "mhCos1FHCalSFHCalN", _cent, _pid);
   if(_method=="v1_FHCal_RP") PrVn = TProfile2D_VnPtEta(inFile, "v1_FHCal_RP", "mhCos1_PhiPsiRP_nFsF",    "mhCos1FHCalFPsiRP" , _cent, _pid);
   if(_method=="v1_FHCal_F")  PrVn = TProfile2D_VnPtEta(inFile, "v1_FHCal_F",  "mhCos1_PhiPsiFHCal_nFsF", "mhCos1FHCalSFHCalN" , _cent, _pid);

   if(PrVn==nullptr) return nullptr;
   
   if(ProfileRange.first=="eta"){
      if(ProfileRange.second.size()==2)
         fResult = TProfileVnvsPt(PrVn, ProfileRange.second[0], ProfileRange.second[1]);
      else
         fResult = TProfileVnvsPt(PrVn, 0.2, 1.0);
   }else if(ProfileRange.first=="pt"){
      if(ProfileRange.second.size()==2)
         fResult = TProfileVnvsEta(PrVn, ProfileRange.second[0], ProfileRange.second[1]);
      else
         fResult = TProfileVnvsEta(PrVn, 0.2, 3.2);
   }
   
   if((int)AxisXBin.size()>1){
      fResult = (TProfile*)fResult->Rebin(AxisXBin.size()-1,Form("%s_rebin",fResult->GetName()),&AxisXBin[0]);
   }

   for(int ipt=0; ipt<fResult->GetNbinsX();ipt++){
      x. push_back(fResult->GetBinCenter(ipt+1));
      xe.push_back(0.);
      y. push_back(fResult->GetBinContent(fResult->FindBin(x.back())));
      ye.push_back(fResult->GetBinError  (fResult->FindBin(x.back())));
   }

   if((int)x.size()==0){
      std::cout<<"Error! Empty graph!\n";
      return nullptr;
   }
   
   TGraphErrors *fGraphResult = new TGraphErrors(x.size(), &x[0], &y[0], &xe[0], &ye[0]);
   fGraphResult->SetName(Form("gr_%s",fResult->GetName()));
   return fGraphResult;
}

void SaveGraphInFile(TFile *inFile, TFile *outFile, std::string _method ,std::vector<std::pair<int,int>> _cent, std::pair<std::string, std::vector<double>> ProfileRange, std::string _pid, std::vector<double> AxisXBin){
   
   TGraphErrors *graph  = vn_diff_EP(inFile,_method,_cent,ProfileRange,_pid,AxisXBin);
   outFile->cd();
   if(ProfileRange.first=="eta"){
      graph->SetName(Form("%s_%s_%i-%i_pt",_method.c_str(),_pid.c_str(),_cent[0].first,_cent.back().second));
   }else{
      graph->SetName(Form("%s_%s_%i-%i_eta",_method.c_str(),_pid.c_str(),_cent[0].first,_cent.back().second));
   }
   if(_method=="v2_TPC_NS" || _method=="v2_FHCal_NS" || _method=="v2_FHCal_RP") graph->GetYaxis()->SetTitle("v_{2}");
   if(_method=="v1_FHCal_NS" || _method=="v1_FHCal_RP" ) graph->GetYaxis()->SetTitle("v_{1}");
   if(ProfileRange.first=="eta") graph->GetXaxis()->SetTitle("p_{T}, GeV/c");
   if(ProfileRange.first=="pt"){
      graph->GetXaxis()->SetTitle("y");
      if(_pid=="ch")graph->GetXaxis()->SetTitle("#eta");
   }
   graph->SetTitle("");
   graph->SetMarkerColor(kRed);
   graph->SetLineColor(kRed);
   graph->SetLineWidth(2);
   graph->Write();
}

void getFlowEP(std::string iFileName="pFlowEP.root", std::string oFileName="pFlowEP_graph.root") 
{
   gSystem->Load("libMathMore.so");
   // Get input information
   TFile *fout = new TFile(oFileName.c_str(),"RECREATE");
   TFile *fi = new TFile(iFileName.c_str(),"read");
   if (!fi) {
      cerr << "No input file was found!" << endl;
      return;
   }
   fout->cd();
   
   std::vector<double> PtBin     ={0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0};
   std::vector<double> EtaBin    ={-1.4,-1.2,-1.0,-0.8,-0.6,-0.4,-0.2,0.,0.2,0.4,0.6,0.8,1.0,1.2,1.4};
   std::vector<double> CentBin   ={0,10,20,30,40,50,60,70,80};

   std::string pid_txt = "PiP";

   SaveGraphInFile(fi,fout,"v1_FHCal_NS",{{0,10}},                  { "pt",{ 0.2 ,2.0}},pid_txt,EtaBin);
   SaveGraphInFile(fi,fout,"v1_FHCal_NS",{{10,20}},                 { "pt",{ 0.2 ,2.0}},pid_txt,EtaBin);
   SaveGraphInFile(fi,fout,"v1_FHCal_NS",{{10,20},{20,30},{30,40}}, { "pt",{ 0.2 ,2.0}},pid_txt,EtaBin);
   SaveGraphInFile(fi,fout,"v1_FHCal_NS",{{40,50},{50,60}},         { "pt",{ 0.2 ,2.0}},pid_txt,EtaBin);
   SaveGraphInFile(fi,fout,"v1_FHCal_NS",{{0,10}},                  {"eta",{ 0.1 ,1.4}},pid_txt,PtBin);
   SaveGraphInFile(fi,fout,"v1_FHCal_NS",{{10,20}},                 {"eta",{ 0.1 ,1.4}},pid_txt,PtBin);
   SaveGraphInFile(fi,fout,"v1_FHCal_NS",{{10,20},{20,30},{30,40}}, {"eta",{ 0.1 ,1.4}},pid_txt,PtBin);
   SaveGraphInFile(fi,fout,"v1_FHCal_NS",{{40,50},{50,60}},         {"eta",{ 0.1 ,1.4}},pid_txt,PtBin);
   
   SaveGraphInFile(fi,fout,"v1_FHCal_F", {{0,10}},                  { "pt",{ 0.2 ,2.0}},pid_txt,EtaBin);
   SaveGraphInFile(fi,fout,"v1_FHCal_F", {{10,20}},                 { "pt",{ 0.2 ,2.0}},pid_txt,EtaBin);
   SaveGraphInFile(fi,fout,"v1_FHCal_F", {{10,20},{20,30},{30,40}}, { "pt",{ 0.2 ,2.0}},pid_txt,EtaBin);
   SaveGraphInFile(fi,fout,"v1_FHCal_F", {{40,50},{50,60}},         { "pt",{ 0.2 ,2.0}},pid_txt,EtaBin);
   SaveGraphInFile(fi,fout,"v1_FHCal_F", {{0,10}},                  {"eta",{ 0.1 ,1.4}},pid_txt,PtBin);
   SaveGraphInFile(fi,fout,"v1_FHCal_F", {{10,20}},                 {"eta",{ 0.1 ,1.4}},pid_txt,PtBin);
   SaveGraphInFile(fi,fout,"v1_FHCal_F", {{10,20},{20,30},{30,40}}, {"eta",{ 0.1 ,1.4}},pid_txt,PtBin);
   SaveGraphInFile(fi,fout,"v1_FHCal_F", {{40,50},{50,60}},         {"eta",{ 0.1 ,1.4}},pid_txt,PtBin);

   SaveGraphInFile(fi,fout,"v1_FHCal_RP",{{0,10}},                  { "pt",{ 0.2 ,2.0}},pid_txt,EtaBin);
   SaveGraphInFile(fi,fout,"v1_FHCal_RP",{{10,20}},                 { "pt",{ 0.2 ,2.0}},pid_txt,EtaBin);
   SaveGraphInFile(fi,fout,"v1_FHCal_RP",{{10,20},{20,30},{30,40}}, { "pt",{ 0.2 ,2.0}},pid_txt,EtaBin);
   SaveGraphInFile(fi,fout,"v1_FHCal_RP",{{40,50},{50,60}},         { "pt",{ 0.2 ,2.0}},pid_txt,EtaBin);
   SaveGraphInFile(fi,fout,"v1_FHCal_RP",{{0,10}},                  {"eta",{ 0.1 ,1.4}},pid_txt,PtBin);
   SaveGraphInFile(fi,fout,"v1_FHCal_RP",{{10,20}},                 {"eta",{ 0.1 ,1.4}},pid_txt,PtBin);
   SaveGraphInFile(fi,fout,"v1_FHCal_RP",{{10,20},{20,30},{30,40}}, {"eta",{ 0.1 ,1.4}},pid_txt,PtBin);
   SaveGraphInFile(fi,fout,"v1_FHCal_RP",{{40,50},{50,60}},         {"eta",{ 0.1 ,1.4}},pid_txt,PtBin);

   SaveGraphInFile(fi,fout,"v2_TPC_NS",  {{0,10}},                  {"eta",{-1.4, 1.4}},pid_txt,PtBin);
   SaveGraphInFile(fi,fout,"v2_TPC_NS",  {{10,20},{20,30},{30,40}}, {"eta",{-1.4, 1.4}},pid_txt,PtBin);
   SaveGraphInFile(fi,fout,"v2_TPC_NS",  {{40,50},{50,60}},         {"eta",{-1.4, 1.4}},pid_txt,PtBin);
   SaveGraphInFile(fi,fout,"v2_TPC_NS",  {{0,10}},                  { "pt",{ 0.2 ,2.0}},pid_txt,EtaBin);
   SaveGraphInFile(fi,fout,"v2_TPC_NS",  {{10,20},{20,30},{30,40}}, { "pt",{ 0.2 ,2.0}},pid_txt,EtaBin);
   SaveGraphInFile(fi,fout,"v2_TPC_NS",  {{40,50},{50,60}},         { "pt",{ 0.2 ,2.0}},pid_txt,EtaBin);

   SaveGraphInFile(fi,fout,"v2_FHCal_NS",{{0,10}},                  {"eta",{-1.4, 1.4}},pid_txt,PtBin);
   SaveGraphInFile(fi,fout,"v2_FHCal_NS",{{10,20},{20,30},{30,40}}, {"eta",{-1.4, 1.4}},pid_txt,PtBin);
   SaveGraphInFile(fi,fout,"v2_FHCal_NS",{{40,50},{50,60}},         {"eta",{-1.4, 1.4}},pid_txt,PtBin);
   SaveGraphInFile(fi,fout,"v2_FHCal_NS",{{0,10}},                  { "pt",{ 0.2 ,2.0}},pid_txt,EtaBin);
   SaveGraphInFile(fi,fout,"v2_FHCal_NS",{{10,20},{20,30},{30,40}}, { "pt",{ 0.2 ,2.0}},pid_txt,EtaBin);
   SaveGraphInFile(fi,fout,"v2_FHCal_NS",{{40,50},{50,60}},         { "pt",{ 0.2 ,2.0}},pid_txt,EtaBin);

   SaveGraphInFile(fi,fout,"v2_FHCal_F", {{0,10}},                  {"eta",{-1.4 ,1.4}},pid_txt,PtBin);
   SaveGraphInFile(fi,fout,"v2_FHCal_F", {{10,20},{20,30},{30,40}}, {"eta",{-1.4 ,1.4}},pid_txt,PtBin);
   SaveGraphInFile(fi,fout,"v2_FHCal_F", {{40,50},{50,60}},         {"eta",{-1.4 ,1.4}},pid_txt,PtBin);
   SaveGraphInFile(fi,fout,"v2_FHCal_F", {{0,10}},                  { "pt",{ 0.2 ,2.0}},pid_txt,EtaBin);
   SaveGraphInFile(fi,fout,"v2_FHCal_F", {{10,20},{20,30},{30,40}}, { "pt",{ 0.2 ,2.0}},pid_txt,EtaBin);
   SaveGraphInFile(fi,fout,"v2_FHCal_F", {{40,50},{50,60}},         { "pt",{ 0.2 ,2.0}},pid_txt,EtaBin);

   SaveGraphInFile(fi,fout,"v2_FHCal_RP",{{0,10}},                  {"eta",{-1.4, 1.4}},pid_txt,PtBin);
   SaveGraphInFile(fi,fout,"v2_FHCal_RP",{{10,20},{20,30},{30,40}}, {"eta",{-1.4, 1.4}},pid_txt,PtBin);
   SaveGraphInFile(fi,fout,"v2_FHCal_RP",{{40,50},{50,60}},         {"eta",{-1.4, 1.4}},pid_txt,PtBin);
   SaveGraphInFile(fi,fout,"v2_FHCal_RP",{{0,10}},                  { "pt",{ 0.2 ,2.0}},pid_txt,EtaBin);
   SaveGraphInFile(fi,fout,"v2_FHCal_RP",{{10,20},{20,30},{30,40}}, { "pt",{ 0.2 ,2.0}},pid_txt,EtaBin);
   SaveGraphInFile(fi,fout,"v2_FHCal_RP",{{40,50},{50,60}},         { "pt",{ 0.2 ,2.0}},pid_txt,EtaBin);

   pid_txt = "PrP";

   SaveGraphInFile(fi,fout,"v1_FHCal_NS",{{0,10}},                  { "pt",{ 0.2 ,2.0}},pid_txt,EtaBin);
   SaveGraphInFile(fi,fout,"v1_FHCal_NS",{{10,20}},                 { "pt",{ 0.2 ,2.0}},pid_txt,EtaBin);
   SaveGraphInFile(fi,fout,"v1_FHCal_NS",{{10,20},{20,30},{30,40}}, { "pt",{ 0.2 ,2.0}},pid_txt,EtaBin);
   SaveGraphInFile(fi,fout,"v1_FHCal_NS",{{40,50},{50,60}},         { "pt",{ 0.2 ,2.0}},pid_txt,EtaBin);
   SaveGraphInFile(fi,fout,"v1_FHCal_NS",{{0,10}},                  {"eta",{ 0.1 ,1.4}},pid_txt,PtBin);
   SaveGraphInFile(fi,fout,"v1_FHCal_NS",{{10,20}},                 {"eta",{ 0.1 ,1.4}},pid_txt,PtBin);
   SaveGraphInFile(fi,fout,"v1_FHCal_NS",{{10,20},{20,30},{30,40}}, {"eta",{ 0.1 ,1.4}},pid_txt,PtBin);
   SaveGraphInFile(fi,fout,"v1_FHCal_NS",{{40,50},{50,60}},         {"eta",{ 0.1 ,1.4}},pid_txt,PtBin);
   
   SaveGraphInFile(fi,fout,"v1_FHCal_F", {{0,10}},                  { "pt",{ 0.2 ,2.0}},pid_txt,EtaBin);
   SaveGraphInFile(fi,fout,"v1_FHCal_F", {{10,20}},                 { "pt",{ 0.2 ,2.0}},pid_txt,EtaBin);
   SaveGraphInFile(fi,fout,"v1_FHCal_F", {{10,20},{20,30},{30,40}}, { "pt",{ 0.2 ,2.0}},pid_txt,EtaBin);
   SaveGraphInFile(fi,fout,"v1_FHCal_F", {{40,50},{50,60}},         { "pt",{ 0.2 ,2.0}},pid_txt,EtaBin);
   SaveGraphInFile(fi,fout,"v1_FHCal_F", {{0,10}},                  {"eta",{ 0.1 ,1.4}},pid_txt,PtBin);
   SaveGraphInFile(fi,fout,"v1_FHCal_F", {{10,20}},                 {"eta",{ 0.1 ,1.4}},pid_txt,PtBin);
   SaveGraphInFile(fi,fout,"v1_FHCal_F", {{10,20},{20,30},{30,40}}, {"eta",{ 0.1 ,1.4}},pid_txt,PtBin);
   SaveGraphInFile(fi,fout,"v1_FHCal_F", {{40,50},{50,60}},         {"eta",{ 0.1 ,1.4}},pid_txt,PtBin);

   SaveGraphInFile(fi,fout,"v1_FHCal_RP",{{0,10}},                  { "pt",{ 0.2 ,2.0}},pid_txt,EtaBin);
   SaveGraphInFile(fi,fout,"v1_FHCal_RP",{{10,20}},                 { "pt",{ 0.2 ,2.0}},pid_txt,EtaBin);
   SaveGraphInFile(fi,fout,"v1_FHCal_RP",{{10,20},{20,30},{30,40}}, { "pt",{ 0.2 ,2.0}},pid_txt,EtaBin);
   SaveGraphInFile(fi,fout,"v1_FHCal_RP",{{40,50},{50,60}},         { "pt",{ 0.2 ,2.0}},pid_txt,EtaBin);
   SaveGraphInFile(fi,fout,"v1_FHCal_RP",{{0,10}},                  {"eta",{ 0.1 ,1.4}},pid_txt,PtBin);
   SaveGraphInFile(fi,fout,"v1_FHCal_RP",{{10,20}},                 {"eta",{ 0.1 ,1.4}},pid_txt,PtBin);
   SaveGraphInFile(fi,fout,"v1_FHCal_RP",{{10,20},{20,30},{30,40}}, {"eta",{ 0.1 ,1.4}},pid_txt,PtBin);
   SaveGraphInFile(fi,fout,"v1_FHCal_RP",{{40,50},{50,60}},         {"eta",{ 0.1 ,1.4}},pid_txt,PtBin);

   SaveGraphInFile(fi,fout,"v2_TPC_NS",  {{0,10}},                  {"eta",{-1.4 ,1.4}},pid_txt,PtBin);
   SaveGraphInFile(fi,fout,"v2_TPC_NS",  {{10,20},{20,30},{30,40}}, {"eta",{-1.4 ,1.4}},pid_txt,PtBin);
   SaveGraphInFile(fi,fout,"v2_TPC_NS",  {{40,50},{50,60}},         {"eta",{-1.4 ,1.4}},pid_txt,PtBin);
   SaveGraphInFile(fi,fout,"v2_TPC_NS",  {{0,10}},                  { "pt",{ 0.2 ,2.0}},pid_txt,EtaBin);
   SaveGraphInFile(fi,fout,"v2_TPC_NS",  {{10,20},{20,30},{30,40}}, { "pt",{ 0.2 ,2.0}},pid_txt,EtaBin);
   SaveGraphInFile(fi,fout,"v2_TPC_NS",  {{40,50},{50,60}},         { "pt",{ 0.2 ,2.0}},pid_txt,EtaBin);

   SaveGraphInFile(fi,fout,"v2_FHCal_NS",{{0,10}},                  {"eta",{-1.4 ,1.4}},pid_txt,PtBin);
   SaveGraphInFile(fi,fout,"v2_FHCal_NS",{{10,20},{20,30},{30,40}}, {"eta",{-1.4 ,1.4}},pid_txt,PtBin);
   SaveGraphInFile(fi,fout,"v2_FHCal_NS",{{40,50},{50,60}},         {"eta",{-1.4 ,1.4}},pid_txt,PtBin);
   SaveGraphInFile(fi,fout,"v2_FHCal_NS",{{0,10}},                  { "pt",{ 0.2 ,2.0}},pid_txt,EtaBin);
   SaveGraphInFile(fi,fout,"v2_FHCal_NS",{{10,20},{20,30},{30,40}}, { "pt",{ 0.2 ,2.0}},pid_txt,EtaBin);
   SaveGraphInFile(fi,fout,"v2_FHCal_NS",{{40,50},{50,60}},         { "pt",{ 0.2 ,2.0}},pid_txt,EtaBin);

   SaveGraphInFile(fi,fout,"v2_FHCal_F", {{0,10}},                  {"eta",{-1.4 ,1.4}},pid_txt,PtBin);
   SaveGraphInFile(fi,fout,"v2_FHCal_F", {{10,20},{20,30},{30,40}}, {"eta",{-1.4 ,1.4}},pid_txt,PtBin);
   SaveGraphInFile(fi,fout,"v2_FHCal_F", {{40,50},{50,60}},         {"eta",{-1.4 ,1.4}},pid_txt,PtBin);
   SaveGraphInFile(fi,fout,"v2_FHCal_F", {{0,10}},                  { "pt",{ 0.2 ,2.0}},pid_txt,EtaBin);
   SaveGraphInFile(fi,fout,"v2_FHCal_F", {{10,20},{20,30},{30,40}}, { "pt",{ 0.2 ,2.0}},pid_txt,EtaBin);
   SaveGraphInFile(fi,fout,"v2_FHCal_F", {{40,50},{50,60}},         { "pt",{ 0.2 ,2.0}},pid_txt,EtaBin);

   SaveGraphInFile(fi,fout,"v2_FHCal_RP",{{0,10}},                  {"eta",{-1.4, 1.4}},pid_txt,PtBin);
   SaveGraphInFile(fi,fout,"v2_FHCal_RP",{{10,20},{20,30},{30,40}}, {"eta",{-1.4, 1.4}},pid_txt,PtBin);
   SaveGraphInFile(fi,fout,"v2_FHCal_RP",{{40,50},{50,60}},         {"eta",{-1.4, 1.4}},pid_txt,PtBin);
   SaveGraphInFile(fi,fout,"v2_FHCal_RP",{{0,10}},                  { "pt",{ 0.2 ,2.0}},pid_txt,EtaBin);
   SaveGraphInFile(fi,fout,"v2_FHCal_RP",{{10,20},{20,30},{30,40}}, { "pt",{ 0.2 ,2.0}},pid_txt,EtaBin);
   SaveGraphInFile(fi,fout,"v2_FHCal_RP",{{40,50},{50,60}},         { "pt",{ 0.2 ,2.0}},pid_txt,EtaBin);

   fi->Close();
   fout->Close();

}