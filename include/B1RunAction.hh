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
// $Id: B1RunAction.hh 93886 2015-11-03 08:28:26Z gcosmo $
//
/// \file B1RunAction.hh
/// \brief Definition of the B1RunAction class

#ifndef B1RunAction_h
#define B1RunAction_h 1

#include "G4UserRunAction.hh"
#include "G4AccumulableManager.hh"
#include "globals.hh"
#include <vector>

#include <iostream>

class G4Run;

class B1RunAction : public G4UserRunAction
{
public:
	B1RunAction(G4String);
	virtual ~B1RunAction();
	
	// virtual G4Run* GenerateRun();
	virtual void BeginOfRunAction(const G4Run*);
	virtual void   EndOfRunAction(const G4Run*);
	
	void AddEdep (G4double edep);
	void AddEdkin (G4double edkin);
	
	// void AddEdepPhot (G4double edepPhot);
	// void AddEdepEl (G4double edepEl);
	
	
	std::vector<G4double>& GetRunEnCmos() {return fRunEnCmos; }
	std::vector<G4double>& GetRunEnCmosPrim() {return fRunEnCmosPrim; }
	std::vector<G4float>& GetRunEnCmosTime() {return fRunEnCmosTime; }
	std::vector<G4double>& GetRunXCmos() {return fRunXCmos; }
	std::vector<G4double>& GetRunYCmos() {return fRunYCmos; }
	std::vector<G4double>& GetRunZCmos() {return fRunZCmos; }
	std::vector<G4double>& GetRunPartCmos() {return fRunPartCmos; }
	
	std::vector<G4double>& GetRunEnExit() {return fRunEnExit; }
	std::vector<G4double>& GetRunXExit() {return fRunXExit; }
	std::vector<G4double>& GetRunYExit() {return fRunYExit; }
	std::vector<G4double>& GetRunZExit() {return fRunZExit; }
	std::vector<G4double>& GetRunCosXExit() {return fRunCosXExit; }
	std::vector<G4double>& GetRunCosYExit() {return fRunCosYExit; }
	std::vector<G4double>& GetRunCosZExit() {return fRunCosZExit; }
	std::vector<G4double>& GetRunPartExit() {return fRunPartExit; }
	std::vector<G4double>& GetRunParentIDExit() {return fRunParentIDExit; }
	
	std::vector<G4int>& GetRunExitProcess() {return fRunExitProcess; }
	
	std::vector<G4double>& GetRunCosX() {return fRunCosX; }
	std::vector<G4double>& GetRunCosY() {return fRunCosY; }
	std::vector<G4double>& GetRunCosZ() {return fRunCosZ; }
	std::vector<G4double>& GetRunEnGen() {return fRunEnGen; }
	std::vector<G4double>& GetRunPartGen() {return fRunPartGen; }
	std::vector<G4double>& GetRunIsotopeGen() {return fRunIsotopeGen; }
	
	std::vector<G4double>& GetRunEAbsComp() {return fRunEAbsComp; }

	std::vector<G4double>& GetRunEnPre() {return fRunEnPre; }
	std::vector<G4double>& GetRunPart() {return fRunPart; }
	std::vector<G4double>& GetRunPreCmosX() {return fRunPreCmosX;}
	std::vector<G4double>& GetRunPreCmosY() {return fRunPreCmosY;}
	std::vector<G4double>& GetRunPreCmosZ() {return fRunPreCmosZ;}
	std::vector<G4double>& GetRunEnPrePrim() {return fRunEnPrePrim; }
	std::vector<G4double>& GetRunPreCmosTrackID() {return fRunPreCmosTrackID;}
	
	std::vector<G4double>& GetRunPreFilterEn() {return fRunPreFilterEn; }
	std::vector<G4double>& GetRunPreFilterPart() {return fRunPreFilterPart; }
	
	std::vector<G4double>& GetRunPixNo() {return fRunPixNo; }
//	std::vector<G4double>& GetRunPixEneDep() {return fRunPixEneDep; }
	std::vector<G4double>& GetRunPixXpos() {return fRunPixXpos; }
	std::vector<G4double>& GetRunPixYpos() {return fRunPixYpos; }
	
	G4int GetEventNumber() {return nbEventInRun;}
	
	void SetMotherIsotope(G4double miso) {fMotherIsotope=miso;}
	G4double GetMotherIsotope() {return fMotherIsotope;}
	
	void SetMotherEnergy(G4double mene) {fMotherEnergy=mene;}
	G4double GetMotherEnergy() {return fMotherEnergy;}

	void SetMotherPart(G4double mpart) {fMotherPart=mpart;}
	G4double GetMotherPart() {return fMotherPart;}
	
	void SetMotherTime(G4double mtime) {fMotherTime=mtime;}
	G4float GetMotherTime() {return fMotherTime;}
	

private:
	G4Accumulable<G4double> fEdep;
	G4Accumulable<G4double> fEdep2;
	G4Accumulable <G4double> fEdkin;
	
	G4int nbEventInRun;

	G4int fMotherIsotope=-10;

	G4double fMotherEnergy=-10;
	G4double fMotherPart=0;
	G4float fMotherTime=0;
	
	/////////////////
	// Histogramming
	//
	void CreateHistogram();
	void WriteHistogram();
	
	std::vector<G4double> fRunEnCmos;
	std::vector<G4double>	fRunEnCmosPrim;
	std::vector<G4float>	fRunEnCmosTime;
	std::vector<G4double> fRunXCmos;
	std::vector<G4double> fRunYCmos;
	std::vector<G4double> fRunZCmos;
	std::vector<G4double> fRunPartCmos;
	
	std::vector<G4double> fRunEnPre;
	std::vector<G4double>	fRunEnPrePrim;
	std::vector<G4double> fRunPart;
	std::vector<G4double> fRunPreCmosX;
	std::vector<G4double> fRunPreCmosY;
	std::vector<G4double> fRunPreCmosZ;
	std::vector<G4double> fRunPreCmosTrackID;

	std::vector<G4double> fRunPreFilterEn;
	std::vector<G4double> fRunPreFilterPart;
	
	std::vector<G4double> fRunPixNo;
//	std::vector<G4double> fRunPixEneDep;
	std::vector<G4double> fRunPixXpos;
	std::vector<G4double> fRunPixYpos;
	
	std::vector<G4double> fRunCosX;
	std::vector<G4double> fRunCosY;
	std::vector<G4double> fRunCosZ;
	
	std::vector<G4double> fRunEnGen;
	std::vector<G4double> fRunPartGen;
	std::vector<G4double> fRunIsotopeGen;
	
	
	std::vector<G4double> fRunEnExit;
	std::vector<G4double> fRunXExit;
	std::vector<G4double> fRunYExit;
	std::vector<G4double> fRunZExit;
	std::vector<G4double> fRunCosXExit;
	std::vector<G4double> fRunCosYExit;
	std::vector<G4double> fRunCosZExit;
	
	std::vector<G4double> fRunPartExit;
	std::vector<G4double> fRunParentIDExit;
	
	std::vector<G4int> fRunExitProcess;
	
	std::vector<G4double> fRunEAbsComp;

	
	G4String fFileName;
	
};

#endif

