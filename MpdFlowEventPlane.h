#ifndef MPDFLOWEVENTPLANE_H
#define MPDFLOWEVENTPLANE_H

#include <deque>

#include "TChain.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "TProfile2D.h"

#include "MpdEvent.h"
#include "MpdPid.h"
#include "MpdHelix.h"
#include "FairMCEventHeader.h"
#include "MpdKalmanFilter.h"
#include "MpdAnalysisTask.h"
#include "MpdParticle.h"
#include "MpdFlowEventPlaneParams.h"
#include "TRandom.h"

class MpdTpcKalmanTrack;
class MpdMCTrack;

class MpdFlowEventPlane : public MpdAnalysisTask {

public:
   MpdFlowEventPlane() {}
   MpdFlowEventPlane(const char *name, const char *outputName = "taskName");
   ~MpdFlowEventPlane();

   void UserInit();
   void ProcessEvent(MpdAnalysisEvent &event);
   void Finish();

   void setOutFile(std::string filename = "histos.root") { mOutFile = filename; }

protected:
   bool selectEvent(MpdAnalysisEvent &event);
   bool selectTrack(MpdTrack *tr, int fMotherId);
   bool selectTrack(MpdAnalysisEvent &event, MpdMCTrack *mctrack);
   
   void processHistograms(MpdAnalysisEvent &event);
   int GetCentBin(float cent);
   int GetPidPDG(int pdg_code);
   void SetPIDvector(std::vector<int> &_vect, MpdMCTrack *mctrack);
   void SetPIDvector(std::vector<int> &_vect, MpdTrack *mpdtrack, int pdg_code);
   double GetRapidity(MpdTrack *mpdtrack, int pdg_code);
   double GetMass2(int pdg_code);
   float GetCentValerii(MpdAnalysisEvent &event);

private:
   // Event properties
   bool     isInitialized = false;
   TVector3 mPrimaryVertex;

   float mCentValerii = -1.;
   int mCent10 = -1;
   
   std::map<int,std::pair<int,int>> mInitCent={
      {0,{ 0,10}},
      {1,{10,20}},
      {2,{20,30}},
      {3,{30,40}},
      {4,{40,50}},
      {5,{50,60}},
      {6,{60,70}},
      {7,{70,80}}
   };
   
   std::map<int,std::string> mInitPID={
      {0,"Ch"}, {1,"ChP"},  {2,"ChM"},
      {3,"Pi"}, {4,"PiP"},  {5,"PiM"},
      {6,"K"},  {7,"KP"},   {8,"KM"},
      {9,"Pr"}, {10,"PrP"}, {11,"PrM"}
   };

   //UrQMD BiBi 9.2
   Int_t fMinMult_Valerii [14]          = { 204, 172, 143, 120, 101,  84, 70, 58, 39, 24, 14,  8,  4, 1};
   Int_t fMaxMult_Valerii [14]          = { 316, 204, 172, 143, 120, 101, 84, 70, 58, 39, 24, 14,  8, 4};
   Float_t fMinCentPercent_Valerii [14] = { 0,     5,  10,  15,  20,  25, 30, 35, 40, 50, 60, 70, 80, 90};
   Float_t fMaxCentPercent_Valerii [14] = { 5,    10,  15,  20,  25,  30, 35, 40, 50, 60, 70, 80, 90, 100};
   Float_t fMinB_Valerii [14] = {    0, 3.09, 4.40, 5.46, 6.33, 7.07, 7.73, 8.33, 8.90,  9.97, 10.95, 11.83, 12.61, 13.41};
   Float_t fMaxB_Valerii [14] = { 3.09, 4.40, 5.46, 6.33, 7.07, 7.73, 8.33, 8.90, 9.97, 10.95, 11.83, 12.61, 13.41, 14.49};

   std::string             mParamConfig;
   MpdFlowEventPlaneParams mParams;

   std::string mOutFile = "histos.root";

   MpdKalmanFilter      *mKF              = nullptr;
   TClonesArray         *mMCTracks        = nullptr;
   TObjArray            *mEMCClusters     = nullptr;
   TClonesArray         *mKalmanTracks    = nullptr;
   TClonesArray         *mMpdGlobalTracks = nullptr;
   vector<MpdParticle *> mPartK;
   MpdKalmanHit          mKHit;

   // Histograms
   TList mHistoList;

   // General QA
   TH1D *mhEvents            = nullptr;
   TH1D *mhTracks            = nullptr;
   TH1D *mhCent              = nullptr;
   TH1D *mhVertex            = nullptr;
   TH1D *mhHits              = nullptr;
   TH1D *mhEta               = nullptr;
   TH1D *mhDca               = nullptr;
   TH1D *mhPt                = nullptr;
   TH1D *mhCh                = nullptr;
   TH2D *mhBCent             = nullptr;
   std::vector<TH1D*> mhRapidity;
   std::vector<TH2D*> mhRapidityPt;

   // Resolution

   TProfile *mh_v1_1040_PiP_eta  = nullptr;
   TProfile *mh_v1_1040_PiP_pt   = nullptr;
   TProfile *mh_v1_1040_PrP_eta  = nullptr;
   TProfile *mh_v1_1040_PrP_pt   = nullptr;
   
   TProfile *mh_v2_1040_PiP_eta  = nullptr;
   TProfile *mh_v2_1040_PiP_pt   = nullptr;
   TProfile *mh_v2_1040_PrP_eta  = nullptr;
   TProfile *mh_v2_1040_PrP_pt   = nullptr;



   TProfile *mhCos1FHCalFPsiRP_All = nullptr;
   TProfile *mhCos1FHCalSFHCalN_All = nullptr;
   TProfile *mhCos2FHCalFPsiRP_All = nullptr;
   TProfile *mhCos2TpcNTpcS_All     = nullptr;
   TProfile *mhCos2FHCalSFHCalN_All = nullptr;

   std::vector<TProfile*> mhCos1FHCalFPsiRP;      // <cos(1(Psi_1^{N}-Psi_1^{S}))>
   std::vector<TProfile*> mhCos1FHCalSFHCalN;     // <cos(1(Psi_1^{N}-Psi_1^{S}))>
   std::vector<TProfile*> mhCos2FHCalFPsiRP;      // <cos(1(Psi_1^{N}-Psi_1^{S}))>
   std::vector<TProfile*> mhCos2TpcNTpcS;         // <cos(2(Psi_2^{N}-Psi_2^{S}))>
   std::vector<TProfile*> mhCos2FHCalSFHCalN;     // <cos(2(Psi_1^{N}-Psi_1^{S}))>
   
   std::vector<std::vector<TProfile2D*>> mhCos1_PhiPsiRP_nFsF;    // <cos(1(phi^{N,S}-Psi_1^{RP}))>
   std::vector<std::vector<TProfile2D*>> mhCos2_PhiPsiRP_nFsF;    // <cos(2(phi^{N,S}-Psi_1^{RP}))>
   std::vector<std::vector<TProfile2D*>> mhCos2_PhiPsiTPC_sNnS;      // <cos(2(phi^{S,N}-Psi_2^{N,S}))>
   std::vector<std::vector<TProfile2D*>> mhCos2_PhiPsiFHCal_sNnS;    // <cos(2(phi^{S,N}-Psi_1^{N,S}))>  
   std::vector<std::vector<TProfile2D*>> mhCos2_PhiPsiFHCal_sFnF;    // <cos(2(phi^{S,N}-Psi_1^{F}))> 
   std::vector<std::vector<TProfile2D*>> mhCos1_PhiPsiFHCal_sNnS;    // <cos(1(phi^{S,N}-Psi_1^{N,S}))>
   std::vector<std::vector<TProfile2D*>> mhCos1_PhiPsiFHCal_nFsF;    // <cos(1(phi^{N,S}-Psi_1^{F}))>


   ClassDef(MpdFlowEventPlane, 1);
};
#endif