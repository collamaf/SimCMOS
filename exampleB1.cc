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

#include <stdio.h>      /* printf, NULL */
#include <stdlib.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
	
	
	// Detect interactive mode (if no arguments) and define UI session
	G4UIExecutive* ui = 0;
	if ( argc == 10 ) {  //was argc==1, 7 to see geom using input parameters, 8 once added sensorchoice
		ui = new G4UIExecutive(argc, argv);
	}
	
	G4double x0Scan=0, ZValue=2*mm, AbsorberDiam=5*mm,AbsorberThickness=1*mm, TBRvalue=1;
	G4int FilterFlag=1, SourceChoice=1, SrSourceFlag=0, SensorChoice=1, AbsorberMaterial=1;
	
	//arguments list: CuZ, Zval, Filter, TBR, Source, X0, Sensor
	
	if ( argc >1 ) {
		//    x0Scan=(*argv[2]-48)*mm;
		AbsorberDiam=strtod (argv[1], NULL);
		AbsorberThickness=strtod (argv[2], NULL);
		AbsorberMaterial=strtod (argv[3], NULL);
		ZValue=strtod (argv[4], NULL);
		FilterFlag=strtod (argv[5], NULL);
		TBRvalue=strtod (argv[6],NULL);
		SourceChoice=strtod (argv[7], NULL);
		x0Scan=strtod (argv[8], NULL);
		SensorChoice=strtod(argv[9],NULL);
		
		G4cout<<"DEBUG Initial parameter check x0= "<<x0Scan<<G4endl;
		G4cout<<"DEBUG Initial parameter check z= "<<ZValue<<G4endl;
		G4cout<<"DEBUG Initial parameter check AbsorberDiam= "<<AbsorberDiam<<G4endl;
		G4cout<<"DEBUG Initial parameter check AbsorberThickness= "<<AbsorberThickness<<G4endl;
		G4cout<<"DEBUG Initial parameter check CuMaterial= "<<AbsorberMaterial<<G4endl;
		G4cout<<"DEBUG Initial parameter check TBRvalue= "<<TBRvalue<<G4endl;
		G4cout<<"DEBUG Initial parameter check FilterFlag= "<<FilterFlag<<G4endl;
		G4cout<<"DEBUG Initial parameter check SourceChoice= "<<SourceChoice<<G4endl;
		G4cout<<"DEBUG Initial parameter check SensorChoice= "<<SensorChoice<<G4endl;
		
	}
	
	G4int SourceSelect=SourceChoice;
	if (SourceSelect==1|| SourceSelect==2) SrSourceFlag=1; //if it is a Sr source... tell to DetCons
	
	G4String MaterialiAssorbitore[3]= {"Cu","Al","PVC"};
	
	G4String FileNamePrim="Primaries";
	G4String OutFileName="CMOSmc";
	G4String FileNameCommonPart;
	
	FileNameCommonPart.append("_X"+ std::to_string((G4int)x0Scan));
	FileNameCommonPart.append("_Z"+ std::to_string((G4int)(100*ZValue)));
	
	if (AbsorberDiam>=0) FileNameCommonPart.append("_AbsDz" + std::to_string((G4int)(1000*AbsorberThickness))+"_AbsHole" + std::to_string((G4int)(100*AbsorberDiam)) +"_AbsMat" + MaterialiAssorbitore[AbsorberMaterial]);
	else FileNameCommonPart.append("_NoAbs");
	

	if (SourceSelect==1) FileNameCommonPart.append("_PSr");
	if (SourceSelect==2) FileNameCommonPart.append("_ExtSr");
	if (SourceSelect==3) FileNameCommonPart.append("_ExtY_TBR"+ std::to_string((G4int)TBRvalue));
//	if (SourceSelect==4) FileNameCommonPart.append("_ExtGa_Diam" + std::to_string((G4int)SourceDiameter) + "_Dz" + std::to_string((G4int)SourceThickness));
	
	if (SensorChoice==1) FileNameCommonPart.append("_011");
	if (SensorChoice==2) FileNameCommonPart.append("_115");
	if (SensorChoice==3) FileNameCommonPart.append("_60035");

	if (argc<11) FileNameCommonPart.append("TEST"); //if it was a TEST run under vis

	
	FileNamePrim.append(FileNameCommonPart);
	OutFileName.append(FileNameCommonPart);


	std::ofstream primFile(FileNamePrim, std::ios::out);
	
	// Choose the Random engine
	G4Random::setTheEngine(new CLHEP::RanecuEngine);
	
	// Construct the default run manager
	//
	//#ifdef G4MULTITHREAD
	//  G4MTRunManager* runManager = new G4MTRunManager;
	//#else
	G4RunManager* runManager = new G4RunManager;
	//#endif
	
	// Set mandatory initialization classes
	// Detector construction
	runManager->SetUserInitialization(new B1DetectorConstruction(x0Scan, ZValue, AbsorberDiam, AbsorberThickness, AbsorberMaterial, FilterFlag, SrSourceFlag, SensorChoice)); //DetectorConstruction needs to know if it is a SrSource to place the right geometry
	
	// Physics list
	//G4VModularPhysicsList* physicsList = new QBBC;
	//physicsList->SetVerboseLevel(1);
	
	//  runManager->SetUserInitialization(new B1PhysicsList);
	
	B1PhysicsList* physicsList=new B1PhysicsList;
	physicsList->RegisterPhysics(new G4StepLimiterPhysics());
	runManager->SetUserInitialization(physicsList);
	
	// User action initialization
	//	runManager->SetUserInitialization(new B1ActionInitialization(x0Scan, ZValue, CuDiam, FilterFlag, primFile, TBRvalue,SourceSelect, SourceSelect));
	runManager->SetUserInitialization(new B1ActionInitialization(x0Scan, ZValue, AbsorberDiam, FilterFlag, primFile, TBRvalue, SourceSelect, SensorChoice, OutFileName));
	
	// Initialize visualization
	//
	G4VisManager* visManager = new G4VisExecutive;
	// G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
	// G4VisManager* visManager = new G4VisExecutive("Quiet");
	visManager->Initialize();
	
	// Get the pointer to the User Interface manager
	G4UImanager* UImanager = G4UImanager::GetUIpointer();
	
	// Process macro or start UI session
	//
	
	if ( ! ui ) {
		// batch mode
		G4String command = "/control/execute ";
		G4String fileName = argv[10];
		UImanager->ApplyCommand(command+fileName);
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
