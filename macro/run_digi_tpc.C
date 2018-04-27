#if !defined(__CLING__) || defined(__ROOTCLING__)
#include <iostream>

#include "Rtypes.h"
#include "TString.h"
#include "TStopwatch.h"
#include "TGeoManager.h"
#include "TFile.h"
#include "TH3.h"
#include "TSystem.h"

#include "FairLogger.h"
#include "FairRunAna.h"
#include "FairFileSource.h"
#include "FairSystemInfo.h"
#include "FairRuntimeDb.h"
#include "FairParRootFileIo.h"

#include "Field/MagneticField.h"

#include "TPCSimulation/DigitizerTask.h"
#include "TPCSimulation/SpaceChargeInterface.h"
#endif

void run_digi_tpc(Int_t nEvents = 10, TString mcEngine = "TGeant3", Int_t isContinuous = 1, Int_t scDistortionsMode=0)
{
  // Initialize logger
  FairLogger* logger = FairLogger::GetLogger();
  logger->SetLogVerbosityLevel("LOW");
  logger->SetLogScreenLevel("DEBUG");

  // Input and output file name
  std::stringstream inputfile, outputfile, paramfile;
  inputfile << "AliceO2_" << mcEngine << ".tpc.mc_" << nEvents << "_event.root";
  paramfile << "AliceO2_" << mcEngine << ".tpc.params_" << nEvents << ".root";
  outputfile << "AliceO2_" << mcEngine << ".tpc.digi_" << nEvents << "_event.root";

  // Setup timer
  TStopwatch timer;

  // Setup FairRoot analysis manager
  FairRunAna* run = new FairRunAna();

  FairFileSource* fFileSource = new FairFileSource(inputfile.str().c_str());
  run->SetSource(fFileSource);
  run->SetOutputFile(outputfile.str().c_str());

  // Setup Runtime DB
  FairRuntimeDb* rtdb = run->GetRuntimeDb();
  FairParRootFileIo* parInput1 = new FairParRootFileIo();
  parInput1->open(paramfile.str().c_str());
  rtdb->setFirstInput(parInput1);

  fFileSource->SetEventMeanTime(20 * 1000); // is in us
  //  TGeoManager::Import("geofile_full.root");

  o2::field::MagneticField* magField =
    new o2::field::MagneticField("Maps", "Maps", -1., -1., o2::field::MagFieldParam::k5kG);
  run->SetField(magField);

  // Setup digitizer
  o2::TPC::DigitizerTask* digiTPC = new o2::TPC::DigitizerTask(0);
  digiTPC->setContinuousReadout(isContinuous);
  digiTPC->setDebugOutput("DigitMCDebug");

  // Switch on distortions and get initial space-charge density histogram if provided in environment variables
  if (scDistortionsMode>0) {
    o2::TPC::SpaceChargeInterface::SCDistortionType distortionType = scDistortionsMode==1 ? o2::TPC::SpaceChargeInterface::SCDistortionType::SCDistortionsRealistic : o2::TPC::SpaceChargeInterface::SCDistortionType::SCDistortionsConstant;
    TH3 *hisSCDensity = nullptr;
    if (gSystem->Getenv("O2TPCSCDensityHisFilePath") && gSystem->Getenv("O2TPCSCDensityHisName")){
      TFile *fileSCInput = TFile::Open(gSystem->Getenv("O2TPCSCDensityHisFilePath"));
      hisSCDensity = (TH3*)fileSCInput->Get(gSystem->Getenv("O2TPCSCDensityHisName"));
    }
    if (distortionType==o2::TPC::SpaceChargeInterface::SCDistortionType::SCDistortionsConstant && !hisSCDensity){
      std::cout << "Constant space-charge distortions require an initial space-charge density histogram. Please provide the path to the root file (O2TPCSCDensityHisFilePath) and the histogram name (O2TPCSCDensityHisName) in your environment variables." << std::endl;
      return;
    }
    digiTPC->enableSCDistortions(distortionType, hisSCDensity);
  }

  run->AddTask(digiTPC);

  run->Init();

  timer.Start();
  run->Run();

  std::cout << std::endl << std::endl;

  // Extract the maximal used memory an add is as Dart measurement
  // This line is filtered by CTest and the value send to CDash
  FairSystemInfo sysInfo;
  Float_t maxMemory = sysInfo.GetMaxMemory();
  std::cout << R"(<DartMeasurement name="MaxMemory" type="numeric/double">)";
  std::cout << maxMemory;
  std::cout << "</DartMeasurement>" << std::endl;

  timer.Stop();
  Double_t rtime = timer.RealTime();
  Double_t ctime = timer.CpuTime();

  Float_t cpuUsage = ctime / rtime;
  std::cout << R"(<DartMeasurement name="CpuLoad" type="numeric/double">)";
  std::cout << cpuUsage;
  std::cout << "</DartMeasurement>" << std::endl;
  std::cout << std::endl << std::endl;
  std::cout << "Macro finished succesfully." << std::endl;
  std::cout << std::endl << std::endl;
  std::cout << "Output file is " << outputfile.str() << std::endl;
  // std::cout << "Parameter file is " << parFile << std::endl;
  std::cout << "Real time " << rtime << " s, CPU time " << ctime << "s" << std::endl << std::endl;
}
