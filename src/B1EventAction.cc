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
// $Id: B1EventAction.cc 93886 2015-11-03 08:28:26Z gcosmo $
//
/// \file B1EventAction.cc
/// \brief Implementation of the B1EventAction class

#include "B1EventAction.hh"
#include "B1RunAction.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"

#include "B1Analysis.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1EventAction::B1EventAction(B1RunAction* runAction)
: G4UserEventAction(),
fRunAction(runAction),
fEdep(0.),
fEdkin(0.),
fno(0),
fPreNo(0),
fPreFilterNo(0),
fEdepEle(0.),
fEdepPos(0.),
fEdepFot(0.),
fEnteringParticle(0),
fPassCounterSource(0.),
fPassCounterCmos(0.),
fNSourceExit(0.),
fStoreTrackIDSource(0),
fStoreTrackIDCmos(0),
fStoreEventIDCmosPrim(-10),
fPrimDecayTime(0)
{}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1EventAction::~B1EventAction()
{}

// Outputfiles definition

//std::ofstream fileOut("cmosenergy.dat", std::ios::out);
//std::ofstream fileOut1("precmosenergy.dat", std::ios::out);

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1EventAction::BeginOfEventAction(const G4Event* )
{
	fEdep = 0.;
	fEdkin = 0.;
	// xstop=-1000;
	//ystop=-1000;
	//zstop=-1000;
	(fRunAction->GetRunEnPre()).clear();
	(fRunAction->GetRunEnPrePrim()).clear();
	(fRunAction->GetRunPart()).clear();
	(fRunAction->GetRunPreCmosX()).clear();
	(fRunAction->GetRunPreCmosY()).clear();
	(fRunAction->GetRunPreCmosZ()).clear();
	(fRunAction->GetRunPreCmosTrackID()).clear();
	(fRunAction->GetRunEnCmos()).clear();
	(fRunAction->GetRunEnCmosPrim()).clear();
	(fRunAction->GetRunEnCmosTime()).clear();
	(fRunAction->GetRunXCmos()).clear();
	(fRunAction->GetRunYCmos()).clear();
	(fRunAction->GetRunZCmos()).clear();
	(fRunAction->GetRunPartCmos()).clear();
	
	//	if(evento->GetEventID()%7000==0)
	//	{
	(fRunAction->GetRunPixNo()).clear();
//	(fRunAction->GetRunPixEneDep()).clear();
	(fRunAction->GetRunPixXpos()).clear();
	(fRunAction->GetRunPixYpos()).clear();
	
	(fRunAction->GetRunCosX()).clear();
	(fRunAction->GetRunCosY()).clear();
	(fRunAction->GetRunCosZ()).clear();

	(fRunAction->GetRunEnGen()).clear();
	(fRunAction->GetRunPartGen()).clear();
	(fRunAction->GetRunIsotopeGen()).clear();

	(fRunAction->SetMotherIsotope(-10));
	(fRunAction->SetMotherEnergy(-10));
	(fRunAction->SetMotherTime(0));

	(fRunAction->GetRunEnExit()).clear();
	(fRunAction->GetRunXExit()).clear();
	(fRunAction->GetRunYExit()).clear();
	(fRunAction->GetRunZExit()).clear();
	(fRunAction->GetRunCosXExit()).clear();
	(fRunAction->GetRunCosYExit()).clear();
	(fRunAction->GetRunCosZExit()).clear();
	(fRunAction->GetRunPartExit()).clear();
	(fRunAction->GetRunParentIDExit()).clear();
	
	(fRunAction->GetRunExitProcess()).clear();
	
	(fRunAction->GetRunEAbsComp()).clear();
	
	(fRunAction->GetRunPreFilterEn()).clear();
	(fRunAction->GetRunPreFilterPart()).clear();
	

	
	
	//	}
	fno=0;
	fPreNo=0;
	fPreFilterNo=0;
	
	fEdepEle=0.;
	fEdepPos=0;
	fEdepFot=0.;
	fNSourceExit=0;
	fPassCounterSource=0;
	fPassCounterCmos=0;
	fStoreTrackIDSource=0;
	fStoreTrackIDCmos=0;
	fStoreEventIDCmosPrim=-10;
	fPassCounterCmosPrim=0;
	fEnteringParticle=0;
	
	fPrimDecayTime=0;
	/*
	 fSourceX=0;
	 fSourceY=0;
	 fSourceZ=0;
	 fSourceEne=0;
	 */
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//std::ofstream primFile("Primaries.dat", std::ios::out);


void B1EventAction::EndOfEventAction(const G4Event* evento)
{
	// accumulate statistics in run action
	
	fRunAction->AddEdep(fEdep);
	fRunAction->AddEdkin(fEdkin);

	(fRunAction->GetRunEAbsComp()).push_back(fEdepEle/keV);
	(fRunAction->GetRunEAbsComp()).push_back(fEdepPos/keV);
	(fRunAction->GetRunEAbsComp()).push_back(fEdepFot/keV);

	G4int NevTot=fRunAction->GetEventNumber();
	
	if ((100*evento->GetEventID())%NevTot==0) G4cout<<"Progress status: "<<(evento->GetEventID()/(G4double)NevTot)*100<<" %, Nev= "<<evento->GetEventID()<<", NTotEv= "<<NevTot<<G4endl;
	
	// get analysis manager
	
	auto analysisManager = G4AnalysisManager::Instance();
	
	// fill ntuple
	
	if(1||fEdep>0)analysisManager->FillNtupleDColumn(0, 0, fEdep/keV);
	analysisManager->FillNtupleDColumn(0, 2, fPreFilterNo); //number of particles entering filter
	analysisManager->FillNtupleDColumn(0, 5, fPreNo); //number of particles entering CMOS
	analysisManager->FillNtupleDColumn(0, 13, fno); //number of hits into the detector
//	if(1/*fEdepSr>0*/)analysisManager->FillNtupleDColumn(0, 10, fEdepSr/keV);
//	if(fEdepEl>0)analysisManager->FillNtupleDColumn(0, 11, fEdepEl/keV);
//	if(1/*fEdepY>0*/)analysisManager->FillNtupleDColumn(0, 11, fEdepY/keV);
	analysisManager->FillNtupleDColumn(0,24, fSourceX/mm);
	analysisManager->FillNtupleDColumn(0,25, fSourceY/mm);
	analysisManager->FillNtupleDColumn(0,26, fSourceZ/mm);
	/*
	analysisManager->FillNtupleDColumn(0,19, fSourceCosX/mm);
	analysisManager->FillNtupleDColumn(0,20, fSourceCosY/mm);
	analysisManager->FillNtupleDColumn(0,21, fSourceCosZ/mm);
	analysisManager->FillNtupleDColumn(0,22, fSourceEne/keV);
	analysisManager->FillNtupleDColumn(0,23, fSourceIsotope);
	*/
	
	//	G4cout<<"CMOSDEBUG - SourceX scritto in root: "<<fSourceX<<G4endl<<G4endl;
	
	if(1||fEdep>0) analysisManager->AddNtupleRow(0);
	
	if(evento->GetEventID()<=1e5){ //to write to proper ntuple all the source particles info
		analysisManager->FillNtupleDColumn(1,0, fSourceX/mm);
		analysisManager->FillNtupleDColumn(1,1, fSourceY/mm);
		analysisManager->FillNtupleDColumn(1,2, fSourceZ/mm);
//		analysisManager->FillNtupleDColumn(1,3, fSourceCosX/mm);
//		analysisManager->FillNtupleDColumn(1,4, fSourceCosY/mm);
//		analysisManager->FillNtupleDColumn(1,5, fSourceCosZ/mm);
//		analysisManager->FillNtupleDColumn(1,6, fSourceEne/keV);
//		analysisManager->FillNtupleDColumn(1,7, fSourceIsotope);
// 		analysisManager->FillNtupleDColumn(1,16, fSourceIsotope);
		analysisManager->FillNtupleDColumn(1,19, fNSourceExit);
		
		analysisManager->AddNtupleRow(1);
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
