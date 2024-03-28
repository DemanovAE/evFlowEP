#include <iostream> // std::cout
#include <fstream>  // std::ifstream
#include <map>

#include "MpdFlowEventPlaneParams.h"

using namespace std;

ClassImp(MpdFlowEventPlaneParams);

void MpdFlowEventPlaneParams::ReadFromFile(std::string fname)
{
   // Read from file in format
   //  varname  value
   //  comments start from #

   std::string fInputFileTtx = fname;
   fInputFileTtx             = fInputFileTtx + ".txt";

   cout << endl << "Read input text file: " << fInputFileTtx << endl;

   std::ifstream ifs(fInputFileTtx);
   if (!ifs) {
      cout << "File " << fInputFileTtx << " can not be opened -> using input parameters from the header file instead"
           << endl;
      return;
   }

   std::string a, b;
   while (ifs.good()) {
      ifs >> a;
      // if comment, skip to the enf of line
      if (a.find_first_of('#') == 0) {
         ifs.ignore(999, '\n');
         continue;
      } else {
         ifs >> b;
         mMap.insert({a, b});
      }
   }
   ifs.close();

   // Parse prepared map
   read("mZvtxCut",     mZvtxCut);
   // Track cuts
   read("mUseMCTracks", mUseMCTracks);
   read("mNofHitsCut",  mNofHitsCut);
   read("mEtaCut",      mEtaCut);
   read("mEtaGapCut",   mEtaGapCut);
   read("mPtminCut",    mPtminCut);
   read("mPtmaxCut",    mPtmaxCut);
   read("mDcaCut",      mDcaCut);
   read("mPidPDG",      mPidPDG);
   read("mPIDsigTPC",   mPIDsigTPC);
   read("mPIDsigTOF",   mPIDsigTOF);

}

void MpdFlowEventPlaneParams::Print() const
{
   cout << "#-------Parameters used for eventPlane analysis------" << endl;
   cout << "# Event selection: " << endl;
   cout << "mZvtxCut       " << mZvtxCut << " // cut on vertex z coordinate" << endl;
   cout << "# Track cuts:  " << endl;
   cout << "mUseMCTracks   " << mUseMCTracks << "  // Use MC or Reco tracks for analysis" << endl;
   cout << "mNofHitsCut    " << mNofHitsCut << "  // minimal number of hits to accept track" << endl;
   cout << "mEtaGapCut     " << mEtaGapCut
        << "  // pseudorapidity gap between 2 TPC sub events (deltaEtaGap). Default 0.1 -> from -0.05 to 0.05." << endl;
   cout << "mEtaCut        " << mEtaCut << "  // maximal pseudorapidity accepted" << endl;
   cout << "mPtminCut      " << mPtminCut << "  // minimal pt used in analysis" << endl;
   cout << "mPtmaxCut      " << mPtmaxCut << "  // maximal pt used in analysis" << endl;
   cout << "mDcaCut        " << mDcaCut << "  // maximal DCA accepted" << endl;
   cout << "# PID cuts:   " << endl;
   cout << "mPidPDG:      " << mPidPDG << "  // use pdg-code to PID (0 - use real PID, 1-use pdg-code)" << endl;
   cout << "mPIDsigTPC:   " << mPIDsigTPC << "  // dEdx PID parameters" << endl;
   cout << "mPIDsigTOF:   " << mPIDsigTOF << "  // Beta PID parameters" << endl;


   
   cout << "------------------------------------------" << endl << endl;
}

void MpdFlowEventPlaneParams::read(std::string name, bool &b)
{
   auto search = mMap.find(name);
   if (search != mMap.end()) {
      if (search->second.compare("true") == 0 || search->second.compare("TRUE") == 0) {
         b = true;
      } else {
         b = false;
      }
   }
}
void MpdFlowEventPlaneParams::read(std::string name, int &b)
{
   auto search = mMap.find(name);
   if (search != mMap.end()) {
      b = atoi(search->second.data());
   }
}
void MpdFlowEventPlaneParams::read(std::string name, float &b)
{
   auto search = mMap.find(name);
   if (search != mMap.end()) {
      b = atof(search->second.data());
   }
}
void MpdFlowEventPlaneParams::read(std::string name, std::string &b)
{
   auto search = mMap.find(name);
   if (search != mMap.end()) {
      b = search->second;
   } else {
      cerr << "[WARNING] MpdFlowEventPlaneParams: could not find parameter " << name << ". Skipping." << endl;
   }
}
