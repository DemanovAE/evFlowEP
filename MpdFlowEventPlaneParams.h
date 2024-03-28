#ifndef MPDFLOWEVENTPLANEPARAMS_H
#define MPDFLOWEVENTPLANEPARAMS_H

#include <map>
#include <string>

#include "Rtypes.h"
#include "TObject.h"

class MpdFlowEventPlaneParams : public TObject {

public:
   //
   // Event selection cuts
   float mZvtxCut = 130.; //(V) event selection cut (cm)

   // Track selection cuts for EP from TPC
   bool  mUseMCTracks= false;//(V) Use MC or Reco tracks for analysis
   int   mNofHitsCut = 16;   //(V) minimal number of hits to accept track
   float mEtaCut     = 1.5;  //(V) maximal pseudorapidity accepted
   float mEtaGapCut  = 0.1;  //(V) pseudorapidity gap between 2 TPC sub events (deltaEtaGap). Default 0.1 -> from -0.05 to 0.05.
   float mPtminCut   = 0.1;  //(V) minimal pt used in analysis
   float mPtmaxCut   = 3.0;  //(V) maximal pt used in analysis
   float mDcaCut     = 2.0;  //(V) maximal DCA accepted
   // PID cuts
   int   mPidPDG     = 1;   // use pdg-code to PID (0 - use real PID, 1-use pdg-code)
   float mPIDsigTPC  = 3.0; // (V)
   float mPIDsigTOF  = 3.0; // (V)

   void ReadFromFile(std::string fname = "EpQa");
   void Print() const;

protected:
   void read(std::string name, bool &b);
   void read(std::string name, int &b);
   void read(std::string name, float &b);
   void read(std::string name, std::string &b);

   std::map<std::string, std::string> mMap;

   ClassDef(MpdFlowEventPlaneParams, 1);
};
#endif // MPDFLOWEVENTPLANEPARAMS_H