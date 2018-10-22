//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: exampleB1.cc 86065 2014-11-07 08:51:15Z gcosmo $
//
/// \file exampleB1.cc
/// \brief Main program of the B1 example

#include "B1DetectorConstruction.hh"
#include "B1ActionInitialization.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"
//#include "QBBC.hh"
#include "B1PhysicsList.hh"

#include "G4SystemOfUnits.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "Randomize.hh"
#include "G4StepLimiter.hh"
#include "G4UserLimits.hh"
#include "G4StepLimiterPhysics.hh"

#include <stdio.h>
#include <stdlib.h>
#include "SteppingVerbose.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
using namespace std;


int main(int argc,char** argv)
{
	
	G4bool VisFlag=false;
	G4bool QuickFlag=false;

	// Detect interactive mode (if no arguments) and define UI session
	G4UIExecutive* ui = 0;
	
	G4double x0Scan=0, ZValue=1.73*mm, AbsorberDiam=0*mm,AbsorberThickness=1*mm, TBRvalue=1;
	G4int FilterFlag=1, SourceChoice=1, SrSourceFlag=0, SensorChoice=1, AbsorberMaterial=1, QuickFlagCommandLine=0;
	
	G4String MacroName ="";
	G4String FileNameLabel="";

	G4double DTmis=600; //in s
	G4double AttSorg[7]={
		350, //PSr
		2.53e3, //ExtSr Rm
		26.5e3, //Y - according to rosa's data for Z scan (9/1/18)
		23.24e3, //Co60
		16.04e3, //Na22
		33.89e3, //Ba133
		37.84e3, //Cs137
	};
	
	G4bool NoOfPrimToGenChangeFlag=false;
	/*
	 2.53e3 - ExtSrRome
	 23.24e3 - Co60
	 16.04e3 - Na22
	 33.89e3 - Ba133
	 37.84e3 - Cs137
	 1 - Pointlike Sr
	 2 - Extended Sr
	 3 - ExtY
	 4 - Co60
	 5 - Na22
	 6 - Ba133
	 7 - Cs137
	 8 - FlatEle 0-3MeV
	 9 - FlatGamma 0-1 MeV
	 */
	
	G4double PixelThickness=2.5; //um
	G4int NoOfPrimToGen=99, Verbose=0;

	for(int i=1;i<argc;i++)
		if(argv[i][0] =='-')
		{
			G4String option(argv[i]);
			G4cout<<"option: "<<i<<" "<<option<<G4endl;
			if(option.compare("-AbsD")==0)
			{
				AbsorberDiam=strtod (argv[++i], NULL);;
			}
			else if(option.compare("-AbsT")==0)
			{
				AbsorberThickness=strtod (argv[++i], NULL);;
			}
			else if(option.compare("-AbsMat")==0)
			{
				AbsorberMaterial=strtod (argv[++i], NULL);;
			}
			else if(option.compare("-Z")==0)
			{
				ZValue=strtod (argv[++i], NULL);;
			}
			else if(option.compare("-Fil")==0)
			{
				FilterFlag=strtod (argv[++i], NULL);;
			}
			else if(option.compare("-TBR")==0)
			{
				TBRvalue=strtod (argv[++i], NULL);;
			}
			else if(option.compare("-Source")==0)
			{
				SourceChoice=strtod (argv[++i], NULL);;
			}
			else if(option.compare("-X")==0)
			{
				x0Scan=strtod (argv[++i], NULL);;
			}
			else if(option.compare("-Vis")==0)
			{
				VisFlag=stoi (argv[++i], NULL);;
			}
			else if(option.compare("-Verbose")==0)
			{
				Verbose=stoi (argv[++i], NULL);;
			}
			else if(option.compare("-Sensor")==0)
			{
				SensorChoice=strtod (argv[++i], NULL);;
			}
			else if(option.compare("-PixT")==0)
			{
				PixelThickness=strtod (argv[++i], NULL);;
			}
			else if(option.compare("-NPrim")==0)
			{
				NoOfPrimToGen=stoi (argv[++i], NULL);;
			}
			else if(option.compare("-Quick")==0)
			{
				QuickFlagCommandLine=strtod (argv[++i], NULL);;
			}else if(option.compare("-Label")==0)
			{
				FileNameLabel= argv[++i];;
			}
			
		}
		else
		{
			MacroName = argv[i]; //se ho trovato una macro (senza il "-" davanti) significa che NON voglio l'interattivo
			VisFlag=false;
			QuickFlag=false;
		}
	
	if (VisFlag) QuickFlag=true;
	
	if (NoOfPrimToGen!=99) NoOfPrimToGenChangeFlag=true;
	
	if (!NoOfPrimToGenChangeFlag && (SourceChoice==2 || (SourceChoice>=4 && SourceChoice<=7))) { // If still 99 it means I did not choose a precise value via command line, so let's compute it! -   To be fixed: what to do in case of PSr/Y
		NoOfPrimToGen=DTmis*AttSorg[(int)SourceChoice-1];
	}
	G4cout<<"\n############## \nI WILL GENERATE n= "<<NoOfPrimToGen<<" primaries \n##############"<<G4endl;
	if (QuickFlagCommandLine) QuickFlag=true;

	G4long seed = time(NULL);
	if (1||VisFlag) seed=12345; //If vis was requested same always the same seed to have reproducibility
	G4Random::setTheSeed(seed);

	G4int SourceSelect=SourceChoice;
	if (SourceSelect==1|| SourceSelect==2) SrSourceFlag=1; //if it is a Sr source... tell to DetCons
	
	G4String MaterialiAssorbitore[3]= {"Cu","Al","ABS"};
	
//	G4String FileNamePrim="Primaries";
	G4String OutFileName="CMOSmc";
	G4String FileNameCommonPart;
	
	FileNameCommonPart.append("_X"+ std::to_string((G4int)x0Scan));
	FileNameCommonPart.append("_Z"+ std::to_string((G4int)(100*ZValue)));
	
	if (AbsorberDiam>=0) FileNameCommonPart.append("_AbsDz" + std::to_string((G4int)(1000*AbsorberThickness))+"_AbsHole" + std::to_string((G4int)(100*AbsorberDiam)) +"_AbsMat" + MaterialiAssorbitore[AbsorberMaterial-1]);
	else FileNameCommonPart.append("_NoAbs");
	
	FileNameCommonPart.append("_Fil" + std::to_string((G4int)FilterFlag));
	
	if (SourceSelect==1) FileNameCommonPart.append("_PSr");
	if (SourceSelect==2) FileNameCommonPart.append("_ExtSr");
	if (SourceSelect==3) FileNameCommonPart.append("_ExtY_TBR"+ std::to_string((G4int)TBRvalue));
	if (SourceSelect==4) FileNameCommonPart.append("_PCo60");
	if (SourceSelect==5) FileNameCommonPart.append("_PNa22");
	if (SourceSelect==6) FileNameCommonPart.append("_PBa133");
	if (SourceSelect==7) FileNameCommonPart.append("_PCs137");
	if (SourceSelect==8) FileNameCommonPart.append("_FlatEle");
	if (SourceSelect==9) FileNameCommonPart.append("_FlatGamma");
	
	if (SensorChoice==1) FileNameCommonPart.append("_011");
	if (SensorChoice==2) FileNameCommonPart.append("_115");
	if (SensorChoice==3) FileNameCommonPart.append("_60035");
	
		FileNameCommonPart.append("_PxT" + std::to_string((G4int)(100*PixelThickness)));
	
	if (VisFlag) FileNameCommonPart.append("TEST"); //if it was a TEST run under vis
	if (QuickFlagCommandLine) FileNameCommonPart.append("_Quick"); //if it was a TEST run under vis

	if (FileNameLabel!="") FileNameCommonPart.append("_" + FileNameLabel);
	
	if (NoOfPrimToGenChangeFlag) FileNameCommonPart.append("_N"+std::to_string((G4int)NoOfPrimToGen)); //if it was a TEST run under vis
	
//	FileNamePrim.append(FileNameCommonPart+".dat");
	OutFileName.append(FileNameCommonPart);
	
//	std::ofstream primFile(FileNamePrim, std::ios::out);
	
	// Choose the Random engine
	G4Random::setTheEngine(new CLHEP::RanecuEngine);
	
	// Construct the default run manager
	//
	//#ifdef G4MULTITHREAD
	//  G4MTRunManager* runManager = new G4MTRunManager;
	//#else
//	G4VSteppingVerbose::SetInstance(new SteppingVerbose); //to use my SteppingVerbose
	G4RunManager* runManager = new G4RunManager;
	//#endif
	
	// Set mandatory initialization classes
	// Detector construction

	runManager->SetUserInitialization(new B1DetectorConstruction(x0Scan, ZValue, AbsorberDiam, AbsorberThickness, AbsorberMaterial, FilterFlag, SourceChoice, SensorChoice, QuickFlag, PixelThickness)); //DetectorConstruction needs to know if it is a SrSource to place the right geometry
	
	// Physics list
	//G4VModularPhysicsList* physicsList = new QBBC;
	//physicsList->SetVerboseLevel(1);
	
	//  runManager->SetUserInitialization(new B1PhysicsList);
	
	B1PhysicsList* physicsList=new B1PhysicsList;
	physicsList->RegisterPhysics(new G4StepLimiterPhysics());
	runManager->SetUserInitialization(physicsList);
	
	// User action initialization
	//	runManager->SetUserInitialization(new B1ActionInitialization(x0Scan, ZValue, CollHoleDiam, FilterFlag, primFile, TBRvalue,SourceSelect, SourceSelect));
	runManager->SetUserInitialization(new B1ActionInitialization(x0Scan, ZValue, AbsorberDiam, FilterFlag, TBRvalue, SourceSelect, SensorChoice, OutFileName));
	
	// Initialize visualization
	//
	G4VisManager* visManager = new G4VisExecutive;
	// G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
	// G4VisManager* visManager = new G4VisExecutive("Quiet");
	visManager->Initialize();
	
	// Get the pointer to the User Interface manager
	G4UImanager* UImanager = G4UImanager::GetUIpointer();
	
	runManager->Initialize();

	if ( VisFlag ) { //Prepare for vis
		ui = new G4UIExecutive(argc, argv);
	}

	if ( ! ui ) {
		// batch mode
		if (MacroName!="") {
			G4String command = "/control/execute ";
			UImanager->ApplyCommand(command+MacroName);
		} else {
			UImanager->ApplyCommand("/tracking/verbose " + std::to_string(Verbose));
			UImanager->ApplyCommand("/run/beamOn " + std::to_string(NoOfPrimToGen));
			//			UImanager->ApplyCommand("/run/beamOn 100");
		}
	}
	else {
		// interactive mode
		UImanager->ApplyCommand("/control/execute init_vis.mac");
		ui->SessionStart();
		delete ui;
	}

	delete visManager;
	delete runManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
