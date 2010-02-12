{
  if (gSystem->Getenv("CMSSW_VERSION")) {
    TString str = gSystem->GetMakeSharedLib();
    if (str.Contains("-m32")==0 && str.Contains("-m64")==0) {
      str.ReplaceAll("g++", "g++ -m32");
      gSystem->SetMakeSharedLib(str);
    }      
    
    gROOT->Macro("$CMSSW_BASE/src/MitAna/macros/setRootEnv.C+");
    loadLibraries("libZAnaMods.so");
  }  
       
  gROOT->Macro("CPlot.cc+");
  gROOT->Macro("MitStyleRemix.cc+");
       
  // Show which process needs debugging
  gInterpreter->ProcessLine(".! ps |grep root.exe");
}
