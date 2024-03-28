# <b>evFlowEP wagon for flow measurements </b>

This wagon performs measurements of the directed (v1) and elliptic (v2) flow using the event-plane method.
Requires `evCentrality`, `evPlane` wagon.

## Basic structure
The main classes are in the evFlowEP directory:

* `MpdFlowEventPlaneParams`: stores input parameters from the config txt file, tools to parse input config txt files
* `MpdFlowEventPlane`: main class that performs flow measurements

## Input file parameters
Input txt file (see `evFlowEP/macros/pFLowEP.txt` for example) can be configured as follows:
```
#-------Parameters used for analysis------
# Event selection:
mZvtxCut 170 #  cut on vertex z coordinate

# Track selection:
mUseMCTracks false
mNofHitsCut   16  # minimal number of hits to accept track
mEtaCut      1.5  # maximal pseudorapidity accepted
mEtaGapCut   0.1  # minimal pseudorapidity accepted: abs(eta)>0.05 for mEtaGap=0.1
mPtminCut    0.1  # minimal pt used in analysis
mPtmaxCut    2.0  # maximal pt used in analysis
mDcaCut      1.0  # maximal DCA accepted
mPidPDG      1.0  # use PDG-code to PID
``` 

One can comment lines or parts of the line by using `#`.

## Usage

The main macro to run is in `evPlane/macros/RunAnalyses.C`.

The final flow measurement result can be calculated using `evPlane/macros/getFlowEP.C`, which has 2 arguments: an input file (pFlowEP.root by default) and an output file (pFlowEP_graph.root by default).