void RunAnalyses(int nEvents = -1){

   gSystem->Load("libZdc.so") ;
   gSystem->Load("libMpdPhysics.so") ;

   MpdAnalysisManager man("ManagerAnal", nEvents) ;
   man.InputFileList("list25.txt") ;
   man.ReadBranches("*") ; 
   man.SetOutput("histos.root") ;
   
   MpdCentralityAll pCentr("pCentr","pCentr") ;
   man.AddTask(&pCentr) ;
   
   MpdEventPlaneAll pEP("pEP","pEP") ;
   man.AddTask(&pEP) ;

   MpdFlowEventPlane pFlowEP("pFlowEP","pFlowEP") ;
   man.AddTask(&pFlowEP) ;

   man.Process() ;

}
