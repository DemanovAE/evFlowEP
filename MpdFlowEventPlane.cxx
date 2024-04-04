#include <iostream>
#include <fstream> // std::ifstream

#include "MpdMCTrack.h"
#include "MpdKalmanFilter.h"
#include "MpdTpcKalmanTrack.h"
#include "MpdParticle.h"
#include "MpdVertex.h"
#include "MpdEvent.h"
#include "MpdZdcDigi.h"
#include "TRandom.h"

#include "MpdFlowEventPlane.h"
#include "TFile.h"

ClassImp(MpdFlowEventPlane);

MpdFlowEventPlane::MpdFlowEventPlane(const char *name, const char *outputName) : MpdAnalysisTask(name, outputName)
{
   mParamConfig = outputName;
}

MpdFlowEventPlane::~MpdFlowEventPlane() {}

void MpdFlowEventPlane::UserInit()
{

   mParams.ReadFromFile(mParamConfig);
   mParams.Print();

   // Prepare histograms etc.
   fOutputList = new TList();
   fOutputList->SetOwner(kTRUE);

   TH1::AddDirectory(kFALSE); // sets a global switch disabling the reference to histos in gROOT and their overwriting

   // Resize vector
   // Res
   mhCos1FHCalFPsiRP .resize(mInitCent.size());
   mhCos2FHCalFPsiRP .resize(mInitCent.size());
   mhCos1FHCalSFHCalN.resize(mInitCent.size());
   mhCos2TpcNTpcS    .resize(mInitCent.size());
   mhCos2FHCalSFHCalN.resize(mInitCent.size());

   // <cos()>
   mhCos1_PhiPsiRP_nFsF   .resize(mInitCent.size());
   mhCos2_PhiPsiRP_nFsF   .resize(mInitCent.size());
   mhCos2_PhiPsiTPC_sNnS  .resize(mInitCent.size());
   mhCos2_PhiPsiFHCal_sNnS.resize(mInitCent.size());
   mhCos2_PhiPsiFHCal_sFnF.resize(mInitCent.size());   
   mhCos1_PhiPsiFHCal_sNnS.resize(mInitCent.size()); 
   mhCos1_PhiPsiFHCal_nFsF.resize(mInitCent.size()); 
   
   for(auto &iCent : mInitCent)
   {
      mhCos1_PhiPsiRP_nFsF   [iCent.first].resize(mInitPID.size());
      mhCos2_PhiPsiRP_nFsF   [iCent.first].resize(mInitPID.size());
      mhCos2_PhiPsiTPC_sNnS  [iCent.first].resize(mInitPID.size());
      mhCos2_PhiPsiFHCal_sNnS[iCent.first].resize(mInitPID.size());
      mhCos2_PhiPsiFHCal_sFnF[iCent.first].resize(mInitPID.size());   
      mhCos1_PhiPsiFHCal_sNnS[iCent.first].resize(mInitPID.size()); 
      mhCos1_PhiPsiFHCal_nFsF[iCent.first].resize(mInitPID.size()); 
   }

   // General QA
   mhEvents = new TH1D("mhEvents", "Number of events", 10, 0., 10.);
   fOutputList->Add(mhEvents);
   mhTracks = new TH1D("mhTracks", "Number of tracks", 10, 0., 10.);
   fOutputList->Add(mhTracks);
   mhCent = new TH1D("mhCent", "Number of events", 100, 0., 100.);
   fOutputList->Add(mhCent);
   mhVertex = new TH1D("hVertex", "Event vertex distribution", 100, -200., 200.);
   fOutputList->Add(mhVertex);
   mhHits = new TH1D("mhHits", "Number of TPC hits", 100, -0.5, 99.5);
   fOutputList->Add(mhHits);
   mhEta = new TH1D("mhEta", "Eta", 100, -2., 2.);
   fOutputList->Add(mhEta);
   mhDca = new TH1D("mhDca", "DCA", 100, 0., 5.);
   fOutputList->Add(mhDca);
   mhPt = new TH1D("mhPt", "Pt", 300, 0., 3.);
   fOutputList->Add(mhPt);
   mhCh = new TH1D("mhCh", "Pt", 10, -5, 5.);
   fOutputList->Add(mhCh);
   mhBCent = new TH2D("mhBCent","b vs Cent;cent;b", 10,0,100,200,0,20);
   fOutputList->Add(mhBCent);

   mhRapidity.resize(mInitPID.size());
   mhRapidityPt.resize(mInitPID.size());
   for(auto &iPID : mInitPID){
      mhRapidity[iPID.first] = new TH1D(Form("mhRapidity_%s",iPID.second.c_str()), Form("%s;y (or #eta for CH)",iPID.second.c_str()), 100, -2., 2.);
      fOutputList->Add(mhRapidity[iPID.first]);
      mhRapidityPt[iPID.first] = new TH2D(Form("mhRapidityPt_%s",iPID.second.c_str()), Form("%s;y (or #eta for CH);p_{T}, GeV/c",iPID.second.c_str()), 100, -2., 2., 300, 0., 3.);
      fOutputList->Add(mhRapidityPt[iPID.first]);
   }

   if(mInitCent.size()==0)
      std::cout<<"Centrality not set!"<<std::endl;
   if(mInitPID.size()==0)
      std::cout<<"PID not set!"<<std::endl;

   mhCos2FHCalFPsiRP_All = new TProfile(Form("mhCos2FHCalFPsiRP_All"),"<cos(2(#Psi_{EP}^{FHCal F}-#Psi_{RP}))> vs. centrality",10,0.,100);
   fOutputList->Add(mhCos2FHCalFPsiRP_All);
   mhCos1FHCalFPsiRP_All = new TProfile(Form("mhCos1FHCalFPsiRP_All"),"<cos(#Psi_{EP}^{FHCal F}-#Psi_{RP})> vs. centrality",10,0.,100);
   fOutputList->Add(mhCos1FHCalFPsiRP_All);
   mhCos2TpcNTpcS_All = new TProfile(Form("mhCos2TpcNTpcS_All"),"<cos(2(#Psi_{EP}^{TPC N}-#Psi_{EP}^{TPC S}))> vs. centrality",10,0.,100);
   fOutputList->Add(mhCos2TpcNTpcS_All);
   mhCos1FHCalSFHCalN_All = new TProfile(Form("mhCos1FHCalSFHCalN_All"),"<cos(#Psi_{EP}^{FHCal N}-#Psi_{EP}^{FHCal S})> vs. centrality",10,0.,100);
   fOutputList->Add(mhCos1FHCalSFHCalN_All);
   mhCos2FHCalSFHCalN_All = new TProfile(Form("mhCos2FHCalSFHCalN_All"),"<cos(2(#Psi_{EP}^{FHCal N}-#Psi_{EP}^{FHCal S}))> vs. centrality",10,0.,100);
   fOutputList->Add(mhCos2FHCalSFHCalN_All);

   for(auto &iCent : mInitCent)
   {
      int ic = iCent.first;
      int iCentMin = iCent.second.first;
      int iCentMax = iCent.second.second;

      mhCos2FHCalFPsiRP[ic]  = new TProfile(Form("mhCos2FHCalFPsiRP_%i%i",iCentMin,iCentMax),
         Form("<cos(2(#Psi_{EP}^{FHCal F}-#Psi_{RP}))> %i-%i %%"    ,iCentMin,iCentMax),1,-0.5,0.5);
      fOutputList->Add(mhCos2FHCalFPsiRP[ic]);

      mhCos1FHCalFPsiRP[ic]  = new TProfile(Form("mhCos1FHCalFPsiRP_%i%i",iCentMin,iCentMax),
         Form("<cos(#Psi_{EP}^{FHCal F}-#Psi_{RP})> %i-%i %%"    ,iCentMin,iCentMax),1,-0.5,0.5);
      fOutputList->Add(mhCos1FHCalFPsiRP[ic]);

      mhCos1FHCalSFHCalN[ic]  = new TProfile(Form("mhCos1FHCalSFHCalN_%i%i",iCentMin,iCentMax),
         Form("<cos(#Psi_{EP}^{FHCal N}-#Psi_{EP}^{FHCal S})> %i-%i %%"    ,iCentMin,iCentMax),1,-0.5,0.5);
      fOutputList->Add(mhCos1FHCalSFHCalN[ic]);
      
      mhCos2TpcNTpcS[ic]      = new TProfile(Form("mhCos2TpcNTpcS_%i%i",iCentMin,iCentMax),
         Form("<cos(2(#Psi_{EP}^{TPC N}-#Psi_{EP}^{TPC S}))> %i-%i %%" ,iCentMin,iCentMax),1,-0.5,0.5);      
      fOutputList->Add(mhCos2TpcNTpcS[ic]);
      
      mhCos2FHCalSFHCalN[ic]  = new TProfile(Form("mhCos2FHCalSFHCalN_%i%i",iCentMin,iCentMax),
         Form("<cos(2(#Psi_{EP}^{FHCal N}-#Psi_{EP}^{FHCal S}))> %i-%i %%" ,iCentMin,iCentMax),1,-0.5,0.5);      
      fOutputList->Add(mhCos2FHCalSFHCalN[ic]);

      for(auto &iPID : mInitPID)
      {
         int ip = iPID.first;
         const char *ipidName = iPID.second.c_str();

         //v1(PsiRP)
         mhCos1_PhiPsiRP_nFsF[ic][ip] = new TProfile2D(Form("mhCos1_PhiPsiRP_nFsF_%i%i_%s"    ,iCentMin,iCentMax,ipidName),
            Form("<cos[1*(#phi-#Psi_{RP})]> cent %i-%i %%, %s; p_T, GeV/c; #eta" ,iCentMin,iCentMax,ipidName),100,0.,5.,40,-2.,2.);
         fOutputList->Add(mhCos1_PhiPsiRP_nFsF[ic][ip]);

         //v2(PsiRP)
         mhCos2_PhiPsiRP_nFsF[ic][ip] = new TProfile2D(Form("mhCos2_PhiPsiRP_nFsF_%i%i_%s"    ,iCentMin,iCentMax,ipidName),
            Form("<cos[2*(#phi-#Psi_{RP})]> cent %i-%i %%, %s; p_T, GeV/c; #eta" ,iCentMin,iCentMax,ipidName),100,0.,5.,40,-2.,2.);
         fOutputList->Add(mhCos2_PhiPsiRP_nFsF[ic][ip]);

         //v2(Psi2,TPC)
         mhCos2_PhiPsiTPC_sNnS[ic][ip] = new TProfile2D(Form("mhCos2_PhiPsiTPC_sNnS_%i%i_%s"    ,iCentMin,iCentMax,ipidName),
            Form("<cos[2*(#phi^{N,S}-#Psi_{2}^{TPC S,N})]> cent %i-%i %%, %s; p_T, GeV/c; #eta" ,iCentMin,iCentMax,ipidName),100,0.,5.,40,-2.,2.);
         fOutputList->Add(mhCos2_PhiPsiTPC_sNnS[ic][ip]);

         //v2(Psi1,FHCal)
         mhCos2_PhiPsiFHCal_sNnS[ic][ip] = new TProfile2D(Form("mhCos2_PhiPsiFHCal_sNnS_%i%i_%s" ,iCentMin,iCentMax,ipidName),
            Form("<cos[2*(#phi^{N,S}-#Psi_{1}^{FHCal S,N})]> cent %i-%i %%, %s; p_T, GeV/c; #eta",iCentMin,iCentMax,ipidName),100,0.,5.,40,-2.,2.);
         fOutputList->Add(mhCos2_PhiPsiFHCal_sNnS[ic][ip]);
         mhCos2_PhiPsiFHCal_sFnF[ic][ip] = new TProfile2D(Form("mhCos2_PhiPsiFHCal_sFnF_%i%i_%s" ,iCentMin,iCentMax,ipidName),
            Form("<cos[2*(#phi^{N,S}-#Psi_{1}^{FHCal F})]> cent %i-%i %%, %s; p_T, GeV/c; #eta"  ,iCentMin,iCentMax,ipidName),100,0.,5.,40,-2.,2.);
         fOutputList->Add(mhCos2_PhiPsiFHCal_sFnF[ic][ip]);
         
         //v1(Psi1, FHCal)
         mhCos1_PhiPsiFHCal_sNnS[ic][ip] = new TProfile2D(Form("mhCos1_PhiPsiFHCal_sNnS_%i%i_%s",iCentMin,iCentMax,ipidName),
            Form("<cos[(#phi^{N,S}-#Psi_{1}^{FHCal S,N})]> cent %i-%i %%, %s; p_T, GeV/c; #eta" ,iCentMin,iCentMax,ipidName),100,0.,5.,40,-2.,2.);
         fOutputList->Add(mhCos1_PhiPsiFHCal_sNnS[ic][ip]);
         mhCos1_PhiPsiFHCal_nFsF[ic][ip] = new TProfile2D(Form("mhCos1_PhiPsiFHCal_nFsF_%i%i_%s",iCentMin,iCentMax,ipidName),
            Form("<cos[(#phi^{N,S}-#Psi_{1}^{FHCal N})]> cent %i-%i %%, %s; p_T, GeV/c; #eta"   ,iCentMin,iCentMax,ipidName),100,0.,5.,40,-2.,2.);
         fOutputList->Add(mhCos1_PhiPsiFHCal_nFsF[ic][ip]);
      }
   }

}

//--------------------------------------
void MpdFlowEventPlane::ProcessEvent(MpdAnalysisEvent &event)
{
   if (!isInitialized) {
      mKF = MpdKalmanFilter::Instance();
      mKHit.SetType(MpdKalmanHit::kFixedR);
      isInitialized = true;
   }

   //float cent    = event.getCentrTPC();
   mCentValerii  = GetCentValerii(event);
   float cent    = mCentValerii;
   int   centBin = GetCentBin(cent);   
   if (centBin == -1) return;
   mCent10  = ( mInitCent.at(centBin).first + mInitCent.at(centBin).second ) / 2.;
   
   if (!selectEvent(event)) return;

   // Getting the event plane angle, then filling the Res2 histograms
   double Phi1EP_FHCal_F = event.fMpdEP.GetPhiEP_FHCal_F_all(); // event plane angle from FHCal N+S // Res full oproximation
   double Phi1EP_FHCal_N = event.fMpdEP.GetPhiEP_FHCal_N_all(); // event plane angle from FHCal N (eta<0)
   double Phi1EP_FHCal_S = event.fMpdEP.GetPhiEP_FHCal_S_all(); // event plane angle from FHCal S (eta>0)
   double Phi2EP_TPC_N   = event.fMpdEP.GetPhiEP_TPC_N_all();   // event plane angle from TPC N (eta<0)
   double Phi2EP_TPC_S   = event.fMpdEP.GetPhiEP_TPC_S_all();   // event plane angle from TPC S (eta>0)
   double PsiRP          = event.fMCEventHeader->GetRotZ();

   if(Phi1EP_FHCal_F!=-9999 && PsiRP!=-9999){
      mhCos1FHCalFPsiRP_All     -> Fill(mCent10, cos(1. * (Phi1EP_FHCal_F - PsiRP)));
      mhCos1FHCalFPsiRP[centBin]-> Fill(0.,      cos(1. * (Phi1EP_FHCal_F - PsiRP)));
      mhCos2FHCalFPsiRP_All     -> Fill(mCent10, cos(2. * (Phi1EP_FHCal_F - PsiRP)));
      mhCos2FHCalFPsiRP[centBin]-> Fill(0.,      cos(2. * (Phi1EP_FHCal_F - PsiRP)));
   }

   if(Phi1EP_FHCal_N!=-9999 && Phi1EP_FHCal_S!=-9999){
      mhCos1FHCalSFHCalN_All -> Fill(mCent10,  cos(1. * (Phi1EP_FHCal_N - Phi1EP_FHCal_S)));
      mhCos2FHCalSFHCalN_All -> Fill(mCent10,  cos(2. * (Phi1EP_FHCal_N - Phi1EP_FHCal_S)));
      mhCos1FHCalSFHCalN[centBin] -> Fill(0.,  cos(1. * (Phi1EP_FHCal_N - Phi1EP_FHCal_S)));
      mhCos2FHCalSFHCalN[centBin] -> Fill(0.,  cos(2. * (Phi1EP_FHCal_N - Phi1EP_FHCal_S)));
   }

   if(Phi2EP_TPC_N!=-9999 && Phi2EP_TPC_S!=-9999){
      mhCos2TpcNTpcS_All     -> Fill(mCent10, cos(2. * (Phi2EP_TPC_N   - Phi2EP_TPC_S)));
      mhCos2TpcNTpcS[centBin]-> Fill(0.    , cos(2. * (Phi2EP_TPC_N   - Phi2EP_TPC_S)));
   }

   mMpdGlobalTracks = event.fMPDEvent->GetGlobalTracks();
   mMCTracks        = event.fMCTrack;
   long int nTracks = 0;
   if(mParams.mUseMCTracks==true){
      nTracks = mMCTracks->GetEntriesFast();
   }else{
      nTracks = mMpdGlobalTracks->GetEntriesFast();
   }

   // Track circle
   for (long int i = 0; i < nTracks; i++)
   {
      MpdTrack   *mpdtrack; 
      MpdMCTrack *mctrack;
      TVector3    momentum;

      // Tracks parameters
      bool   fGoodTrack= true;
      double fPhi      = -999.;
      double fPt       = -999.;
      double fEta      = -999.;
      double fRapidity = -999.;
      std::vector<int> nPID;
      // Flow
      double v2Psi2    = -999.;
      double v2Psi1SN  = -999.;
      double v2Psi1F   = -999.;
      double v2Rp      = -999.;
      double v1SN      = -999.;
      double v1F       = -999.;
      double v1Rp      = -999.;

      if(mParams.mUseMCTracks==true)
      {
         mctrack = (static_cast<MpdMCTrack *>(mMCTracks->UncheckedAt(i)));
         fGoodTrack = selectTrack(event,mctrack);
         mctrack->GetMomentum(momentum);
         fPhi      = momentum.Phi();
         fEta      = momentum.Eta();
         fPt       = mctrack->GetPt();
         fRapidity = mctrack->GetRapidity();
         SetPIDvector(nPID,mctrack);
      }
      else
      {
         mpdtrack = (MpdTrack *)mMpdGlobalTracks->UncheckedAt(i);
         long int matchId = mpdtrack->GetID();
            //if (matchId == -1) fGoodTrack=false;
         mctrack  = (static_cast<MpdMCTrack *>(mMCTracks->UncheckedAt((int)mpdtrack->GetID())));
         fGoodTrack = selectTrack(mpdtrack, (int)mctrack->GetMotherId());
         fPhi      = mpdtrack->GetPhi();
         fEta      = mpdtrack->GetEta();
         fPt       = mpdtrack->GetPt();
         fRapidity = GetRapidity(mpdtrack, mctrack->GetPdgCode());
         SetPIDvector(nPID,mpdtrack,mctrack->GetPdgCode());
      }
      if(!fGoodTrack)continue;
      if(fPhi==-999. || fEta==-999. || fPt==-999.) continue;

      // v2(Psi2), 2-sub EP method
      if( Phi2EP_TPC_N != -9999 && Phi2EP_TPC_S != -9999 )
         v2Psi2 = cos(2. * (fPhi - ( fEta < 0. ? Phi2EP_TPC_S : Phi2EP_TPC_N) ));         

      // v2(Psi1), 2-sub EP method
      if( Phi1EP_FHCal_N != -9999 && Phi1EP_FHCal_S != -9999 )
         v2Psi1SN = cos(2. * (fPhi - ( fEta < 0. ? Phi1EP_FHCal_S : Phi1EP_FHCal_N) ));         
      
      // v2(Psi1), event plane angle from FHCal N+S
      if( Phi1EP_FHCal_F != -9999 )
         v2Psi1F = cos(2. * (fPhi - Phi1EP_FHCal_F ));

      // v2(PsiRP), event plane angle from FHCal N+S
      if( PsiRP != -9999 )
         v2Rp = cos(2. * (fPhi - PsiRP ));
      
      // v1(Psi1), 2-sub EP method
      if( Phi1EP_FHCal_N != -9999 && Phi1EP_FHCal_S != -9999 )
         v1SN = cos( fPhi - ( fEta < 0. ? Phi1EP_FHCal_S : Phi1EP_FHCal_N) );         
      
      // v1(Psi1), event plane angle from FHCal N+S
      if( Phi1EP_FHCal_F != -9999 )
         v1F  = cos( fPhi - Phi1EP_FHCal_F );       

      // v1(Rp), event plane angle from FHCal N+S
      if( PsiRP != -9999 )
         v1Rp  = cos( fPhi - PsiRP );


      for(auto &_pid : nPID)
      {   
         if(_pid==-999)continue;

         mhRapidity[_pid]  ->Fill(fRapidity);
         mhRapidityPt[_pid]->Fill(fRapidity,fPt);

         if(v2Psi2!=-999.)   mhCos2_PhiPsiTPC_sNnS[centBin][_pid]  ->Fill(fPt,fRapidity,v2Psi2);
         if(v2Psi1SN!=-999.) mhCos2_PhiPsiFHCal_sNnS[centBin][_pid]->Fill(fPt,fRapidity,v2Psi1SN);
         if(v2Psi1F!=-999.)  mhCos2_PhiPsiFHCal_sFnF[centBin][_pid]->Fill(fPt,fRapidity,v2Psi1F);
         if(v1SN!=-999.)     mhCos1_PhiPsiFHCal_sNnS[centBin][_pid]->Fill(fPt,fRapidity,v1SN);
         if(v1F!=-999.)      mhCos1_PhiPsiFHCal_nFsF[centBin][_pid]->Fill(fPt,fRapidity,v1F);
         if(v1Rp!=-999.)     mhCos1_PhiPsiRP_nFsF[centBin][_pid]   ->Fill(fPt,fRapidity,v1Rp);
         if(v2Rp!=-999.)     mhCos2_PhiPsiRP_nFsF[centBin][_pid]   ->Fill(fPt,fRapidity,v2Rp);

      }
      // Track circle
   }

}

void MpdFlowEventPlane::Finish()
{
   // Post-scan processing not needed
}

//--------------------------------------
bool MpdFlowEventPlane::selectEvent(MpdAnalysisEvent &event)
{
   mhEvents->Fill(0.5);

   // Reject empty events (UrQMD, PHSD)
   mMCTracks = event.fMCTrack;

   int nTrMc = 0;
   for (int i = 0; i < mMCTracks->GetEntriesFast(); i++) {
      MpdMCTrack *pr  = (static_cast<MpdMCTrack *>(mMCTracks->At(i)));
      float       eta = 0.5 * log((pr->GetP() + pr->GetPz()) / (pr->GetP() - pr->GetPz()));
      if (pr->GetMotherId() == -1 && fabs(eta) <= mParams.mEtaCut && fabs(eta) >= mParams.mEtaGapCut / 2.) {
         nTrMc++;
      }
   } // i

   if (nTrMc <= 2) { // check if there're at least 2 mc particles within eta cuts
      return false;
   }

   mhEvents->Fill(1.5);

   // Reject bad vertex
   if (!event.fVertex) {
      return false;
   }

   MpdVertex *vertex = (MpdVertex *)event.fVertex->First();
   vertex->Position(mPrimaryVertex);

   if (mPrimaryVertex.Z() == 0) { // not reconstructed (==0)
      return false;
   }

   if (fabs(mPrimaryVertex.Z()) > mParams.mZvtxCut) { // beyond the limits
      return false;
   }
   
   mhEvents->Fill(2.5);

   ///// Test 
   mhVertex->Fill(mPrimaryVertex.Z());
   mhCent->Fill(mCentValerii);
   //mhCent->Fill(event.getCentrTPC());
   mhBCent->Fill(mCent10, event.fMCEventHeader->GetB() );

   return true;
}

bool MpdFlowEventPlane::selectTrack(MpdTrack *mpdtrack, int fMotherId)
{
   
   mhTracks->Fill(0.5);

   if (mpdtrack->GetNofHits() < mParams.mNofHitsCut) return false; // nhits > 16
   mhTracks->Fill(1.5);
   
   if (fabs(mpdtrack->GetEta()) > mParams.mEtaCut) return false;         // |eta| < 1.5
   if (fabs(mpdtrack->GetEta()) < mParams.mEtaGapCut / 2.) return false; // |eta| > 0.05
   mhTracks->Fill(2.5);

   float pt = TMath::Abs(mpdtrack->GetPt());
   if (pt < mParams.mPtminCut) return false; // pT > 100 MeV/c
   if (pt > mParams.mPtmaxCut) return false; // pT < 2000 MeV/c
   mhTracks->Fill(3.5);

   if (fabs(mpdtrack->GetDCAX()) > mParams.mDcaCut) return false; // |DCAx| < 2.
   if (fabs(mpdtrack->GetDCAY()) > mParams.mDcaCut) return false; // |DCAy| < 2.
   if (fabs(mpdtrack->GetDCAZ()) > mParams.mDcaCut) return false; // |DCAz| < 2.
   mhTracks->Fill(4.5);

   if(fMotherId > -1) return false;
   mhTracks->Fill(5.5);

   mhCh->Fill(mpdtrack->GetCharge());
   mhHits->Fill(mpdtrack->GetNofHits());
   mhPt->Fill(pt);
   mhEta->Fill(mpdtrack->GetEta());
   mhDca->Fill(sqrt(mpdtrack->GetDCAX() * mpdtrack->GetDCAX() + mpdtrack->GetDCAY() * mpdtrack->GetDCAY() +
                    mpdtrack->GetDCAZ() * mpdtrack->GetDCAZ()));

   return true;
}

// Select MC tracks
bool MpdFlowEventPlane::selectTrack(MpdAnalysisEvent &event, MpdMCTrack *mctrack)
{
   
   TVector3 _track;
   mctrack->GetMomentum(_track);

   mhTracks->Fill(0.5);

   if (fabs(_track.Eta()) > mParams.mEtaCut) return false;         // |eta| < 1.5
   if (fabs(_track.Eta()) < mParams.mEtaGapCut / 2.) return false; // |eta| > 0.05
   mhTracks->Fill(2.5);

   float pt = TMath::Abs(mctrack->GetPt());
   if (pt < mParams.mPtminCut) return false; // pT > 100 MeV/c
   if (pt > mParams.mPtmaxCut) return false; // pT < 2000 MeV/c
   mhTracks->Fill(3.5);

   MpdVertex *vertex = (MpdVertex *)event.fVertex->First();
   vertex->Position(mPrimaryVertex);
   if(abs(mctrack->GetStartX() - mPrimaryVertex.X()) > mParams.mDcaCut) return false;
   if(abs(mctrack->GetStartY() - mPrimaryVertex.Y()) > mParams.mDcaCut) return false;
   if(abs(mctrack->GetStartZ() - mPrimaryVertex.Z()) > mParams.mDcaCut) return false;
   mhTracks->Fill(4.5);

   if((int)mctrack->GetMotherId() > -1) return false;
   mhTracks->Fill(5.5);

   mhPt->Fill(pt);
   mhEta->Fill(_track.Eta());
   mhDca->Fill(sqrt( pow(mctrack->GetStartX() - mPrimaryVertex.X(), 2) + pow(mctrack->GetStartY() - mPrimaryVertex.Y(), 2) +
                     pow(mctrack->GetStartZ() - mPrimaryVertex.Z(), 2)));

   return true;
}

int MpdFlowEventPlane::GetPidPDG(int pdg_code){

   if (pdg_code ==  211)  return 4;  // pion+
   if (pdg_code == -211)  return 5;  // pion-
   if (pdg_code ==  321)  return 7;  // kaon+
   if (pdg_code == -321)  return 8;  // kaon-
   if (pdg_code ==  2212) return 10;  // proton
   if (pdg_code == -2212) return 11;  // anti-proton

   return -999;
}

void MpdFlowEventPlane::SetPIDvector(std::vector<int> &_vect, MpdMCTrack *mctrack)
{
   int pdg_code   = mctrack->GetPdgCode();
   int n_pdg_code = MpdFlowEventPlane::GetPidPDG(pdg_code);
   _vect.clear();
   
   if(pdg_code < 0 || pdg_code > 0)_vect.push_back(0);// all charge hadrons
   if(pdg_code > 0)_vect.push_back(1);                // ch-
   if(pdg_code < 0)_vect.push_back(2);                // ch+
   if(n_pdg_code!=-999){
      _vect.push_back(n_pdg_code);                    // pion+,pion-,kaon+,kaon-,proton,antiproton
      _vect.push_back(n_pdg_code - n_pdg_code%3);     // all pions, all kaons, all protons
   }
}

void MpdFlowEventPlane::SetPIDvector(std::vector<int> &_vect, MpdTrack *mpdtrack, int pdg_code)
{
   _vect.clear();
   int n_pdg_code = MpdFlowEventPlane::GetPidPDG(pdg_code);

   if(mParams.mPidPDG==1)
   {
      if(pdg_code < 0 || pdg_code > 0)_vect.push_back(0);// all charge hadrons
      if(pdg_code > 0)_vect.push_back(1);                // ch-
      if(pdg_code < 0)_vect.push_back(2);                // ch+
      if(n_pdg_code!=-999){
         _vect.push_back(n_pdg_code);                    // pion+,pion-,kaon+,kaon-,proton,antiproton
         _vect.push_back(n_pdg_code - n_pdg_code%3);     // all pions, all kaons, all protons
      }
   }
   else if (mParams.mPidPDG==0)
   {
      if(mpdtrack->GetCharge() < 0 || mpdtrack->GetCharge() > 0)_vect.push_back(0); // all charge hadrons
      if(mpdtrack->GetCharge() > 0)_vect.push_back(1);    // ch-
      if(mpdtrack->GetCharge() < 0)_vect.push_back(2);    // ch+
      if(n_pdg_code!=-999){
         _vect.push_back(n_pdg_code);                    // pion+,pion-,kaon+,kaon-,proton,antiproton
         _vect.push_back(n_pdg_code - n_pdg_code%3);     // all pions, all kaons, all protons
      }
   }
}

double MpdFlowEventPlane::GetRapidity(MpdTrack *mpdtrack, int pdg_code){

   if( MpdFlowEventPlane::GetMass2(pdg_code) == -999) return -999.;
   double p  = sqrt( mpdtrack->GetPx()*mpdtrack->GetPx() + mpdtrack->GetPy()*mpdtrack->GetPy() + mpdtrack->GetPz()*mpdtrack->GetPz());
   double E  = sqrt( MpdFlowEventPlane::GetMass2(pdg_code) + p*p);
   double pz = mpdtrack->GetPz();
   return 0.5 * log( (E+pz) / (E-pz) );
}

double MpdFlowEventPlane::GetMass2(int pdg_code){
   if(abs(pdg_code) == 211)  return 0.13957061*0.13957061;
   if(abs(pdg_code) == 321)  return 0.493667*0.493667;
   if(abs(pdg_code) == 2212) return 0.938272*0.938272;
   return -999;
}

int MpdFlowEventPlane::GetCentBin(float cent)
{
   if (cent > 0. && cent <= 10.)
      return 0;
   else if (cent > 10. && cent <= 20.)
      return 1;
   else if (cent > 20. && cent <= 30.)
      return 2;
   else if (cent > 30. && cent <= 40.)
      return 3;
   else if (cent > 40. && cent <= 50.)
      return 4;
   else if (cent > 50. && cent <= 60.)
      return 5;
   else if (cent > 60. && cent <= 70.)
      return 6;
   else if (cent > 70. && cent <= 80.)
      return 7;
   else
      return -1;
}

float MpdFlowEventPlane::GetCentValerii(MpdAnalysisEvent &event){

   // Reject empty events (UrQMD, PHSD)
   mMCTracks = event.fMCTrack;
   mMpdGlobalTracks = event.fMPDEvent->GetGlobalTracks();

   int nTrMc = 0;
   for (int i = 0; i < mMCTracks->GetEntriesFast(); i++) {
      MpdMCTrack *pr  = (static_cast<MpdMCTrack *>(mMCTracks->At(i)));
      float       eta = 0.5 * log((pr->GetP() + pr->GetPz()) / (pr->GetP() - pr->GetPz()));
      if (pr->GetMotherId() == -1 && fabs(eta) <= mParams.mEtaCut && fabs(eta) >= mParams.mEtaGapCut / 2.) {
         nTrMc++;
      }
   } // i

   if (nTrMc <= 2) { // check if there're at least 2 mc particles within eta cuts
      return -1.;
   }

   // Reject bad vertex
   if (!event.fVertex) {
      return -1.;
   }

   MpdVertex *vertex = (MpdVertex *)event.fVertex->First();
   vertex->Position(mPrimaryVertex);

   //if (mPrimaryVertex.Z() == 0) return -1.;
   if (vertex->GetChi2() < 1) return -1.;
   if (fabs(mPrimaryVertex.Z()) > 50.) return -1.;

   
   float _Mult = 0.;
   
   for (long int i = 0; i < mMpdGlobalTracks->GetEntriesFast(); i++)
   {
      MpdTrack *mpdtrack = (MpdTrack *)mMpdGlobalTracks->UncheckedAt(i);
      long int matchId = mpdtrack->GetID();
      if(mpdtrack->GetNofHits()>16 && 
         fabs(mpdtrack->GetEta())<0.5 && 
         mpdtrack->GetPt()>0.15 && mpdtrack->GetCharge()!=0 && 
         fabs(mpdtrack->GetDCAX())<1. && 
         fabs(mpdtrack->GetDCAY())<1. && 
         fabs(mpdtrack->GetDCAZ())<1. )
      {
         _Mult=_Mult+1.;
      }
   }

   for(int i=0; i<14; i++){
      if (_Mult >= (float)fMinMult_Valerii[i] && _Mult < (float)fMaxMult_Valerii[i]){
         return ( fMinCentPercent_Valerii[i] + fMaxCentPercent_Valerii[i] )/2.;
      }
   }

   return -1.;

}