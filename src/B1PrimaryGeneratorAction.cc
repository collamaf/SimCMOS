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
// $Id: B1PrimaryGeneratorAction.cc 94307 2015-11-11 13:42:46Z gcosmo $
//
/// \file B1PrimaryGeneratorAction.cc
/// \brief Implementation of the B1PrimaryGeneratorAction class

#include "B1PrimaryGeneratorAction.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4PhysicalConstants.hh"

#include "G4IonTable.hh"
#include "G4ChargedGeantino.hh"

#include "B1RunAction.hh"
#include "B1Analysis.hh"


#include "G4Event.hh"

#include <iostream>
#include <fstream>


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using std::ofstream;
using std::ios;
using std::endl;


B1PrimaryGeneratorAction::B1PrimaryGeneratorAction(B1EventAction* eventAction, G4double TBR, G4int SourceSelect)
: G4VUserPrimaryGeneratorAction(),
fParticleGun(0) ,
evtPrimAction(eventAction), fTBR(TBR), fSourceSelect(SourceSelect)

{
	G4int n_particle = 1;
	fParticleGun  = new G4ParticleGun(n_particle);

	switch (fSourceSelect) {
		case 1: //PSr
			fRadiusInt=0*mm;
			fDZInt=0*mm;
			fRadiusExt=0*mm;
			fDZExt=0*mm;
			break;
			
		case 8: //FlatEle
		case 9: //FlatGamma
		case 2: //ExtSr
			fRadiusInt=8*mm;  //8 for RM, 10.5mm PG source
			fDZInt=0*mm;
			fRadiusExt=8*mm;
			fDZExt=0*mm;
			break;
			
		case 3: //ExtY
			fRadiusInt=3*mm;
			fDZInt=1*mm;
			fRadiusExt=10.48*mm; //10.48 per Rosa, 6.65 per PG
			fDZExt=4.57*mm;   //4.4 per Rosa, 5.5 per PG, metto 4.57 per rosa dato V = 1.58
			break;
			
		case 4: //Co60
		case 5: //Na22
		case 6: //Ba133
		case 7: //Cs137
			fRadiusInt=0.5*mm;
			fDZInt=0*mm;
			fRadiusExt=0.5*mm; 
			fDZExt=0*mm;
			break;
			
		case 10: //Na22 nuda
			fRadiusInt=1.5*mm;
			fDZInt=0*mm;
			fRadiusExt=1.5*mm;
			fDZExt=0*mm;
			break;
			
		default:
			break;
	}
	
	if (fSourceSelect==8) FlatEle = true;
	if (fSourceSelect==9) FlatGamma=true;
	
	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	G4ParticleDefinition* particle = particleTable->FindParticle("geantino");
	
	fParticleGun->SetParticleDefinition(particle);
	
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1PrimaryGeneratorAction::~B1PrimaryGeneratorAction()
{
	delete fParticleGun;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1PrimaryGeneratorAction::GeneratePrimaries (G4Event* anEvent)
{

	
	switch (fSourceSelect) {
		case 3:
			fSourceZ=39;
			fSourceA=90;
			break;

		case 4:
			fSourceZ=27;
			fSourceA=60;
			break;

		case 5:
		case 10:
			fSourceZ=11;
			fSourceA=22;
			break;
			
		case 6:
			fSourceZ=56;
			fSourceA=133;
			break;

		case 7:
			fSourceZ=55;
			fSourceA=137;
			break;
			
		default: //Sr
			fSourceZ=38;
			fSourceA=90;
			break;
	}
	
	G4double ionCharge   = 0.*eplus;
	G4double excitEnergy = 0.*keV;
	
	G4ParticleDefinition* sourceION
	= G4IonTable::GetIonTable()->GetIon(fSourceZ,fSourceA,excitEnergy);
	if (FlatEle) {
		fParticleGun->SetParticleDefinition(	G4ParticleTable::GetParticleTable()->FindParticle("e-"));
	} else if (FlatGamma) {
		fParticleGun->SetParticleDefinition(	G4ParticleTable::GetParticleTable()->FindParticle("gamma"));
		fParticleGun->SetParticleCharge(0);
	} else  {
		fParticleGun->SetParticleDefinition(sourceION);
		fParticleGun->SetParticleCharge(ionCharge);
	}
	
	
	//###################################################
	// Compute volumes for TBR source
	//##########################
	G4double VolA=CLHEP::pi*fDZExt*(fRadiusExt*fRadiusExt-fRadiusInt*fRadiusInt);
	G4double VolB=CLHEP::pi*fRadiusInt*fRadiusInt*fDZInt;
	G4double VolC=CLHEP::pi*fRadiusInt*fRadiusInt*(fDZExt-fDZInt);
	G4double denominatore=VolA+VolB*fTBR+VolC;
	G4double ProbA=VolA/denominatore;
	G4double ProbB=VolB*fTBR/denominatore;
	G4double ProbC=VolC/denominatore;
	//###################################################
	
	G4double zSource=0;
	G4double zSourceOffset=1e-6*mm; //to avoid generating particles at the very boundary of source!
	
	if (fRadiusExt==fRadiusInt) { //se ho un solo raggio ignoro il TBR e faccio la pasticca di sorgente
		fRadiusMax=fRadiusInt;
		fRadiusMin=0*mm;
		zSource = -zSourceOffset;
	} else { //se ho una sorgente con TBR (lo capisco dal fatto che ho due raggi diversi per Int e Ext)
		G4double random=G4UniformRand();
		if (random<=ProbA) {  //faccio il cilindretto cavo esterno al centro (VolA)
			fRadiusMax=fRadiusExt;
			fRadiusMin=fRadiusInt;
			fZ=fDZExt;
			zSource = -G4UniformRand()*fZ-zSourceOffset;
		} else if (random>ProbA && random<=ProbA+ProbB) {    //faccio il cilindretto attivo al centro (VolB) SEGNALE!!!!
			fRadiusMax=fRadiusInt;
			fRadiusMin=0*mm;
			fZ=fDZInt;
			zSource = -G4UniformRand()*fZ-zSourceOffset;
		} else if (random>ProbA+ProbB) {     //faccio il cilindretto dietro a quello attivo al centro (VolC)
			fRadiusMax=fRadiusInt;
			fRadiusMin=0*mm;
			fZ=fDZExt-fDZInt;
			zSource = -G4UniformRand()*fZ-fDZInt-zSourceOffset;
		}
	}
	
	if (fSourceSelect>=4 && fSourceSelect<=7 ) zSource=-1.5*mm; //Co60 Na Ba Cs sources are in the middle of the plastic thickness
	
	//###################################################
	// Sampling particle energy
	//##########################
	if (FlatEle) {
		G4double randomEne=G4UniformRand()*3;
		fParticleGun->SetParticleEnergy(randomEne*MeV); //SetParticleEnergy uses kinetic energy
		evtPrimAction->SetSourceEne(randomEne);
	} else if (FlatGamma) {
		G4double randomEne=G4UniformRand()*1;
		fParticleGun->SetParticleEnergy(randomEne*MeV);
		evtPrimAction->SetSourceEne(randomEne);
	} else {
		fParticleGun->SetParticleEnergy(0*MeV);
	}
	//###################################################

	//###################################################
	// Sampling particle initial position
	//##########################
	G4double rho = sqrt(fRadiusMin*fRadiusMin + G4UniformRand()*(fRadiusMax*fRadiusMax-fRadiusMin*fRadiusMin));   //fixed square problem by collamaf with internal radius!
	G4double alpha = G4UniformRand()*CLHEP::pi*2.;
	const G4ThreeVector position = G4ThreeVector(rho*cos(alpha), rho*sin(alpha), zSource);
	//###################################################

	//###################################################
	// Sampling particle initial direction
	//##########################
	if (FlatEle || FlatGamma) { //If FlatSource (for Eff) was requested, generate only towards up
		G4double phi = G4UniformRand()*CLHEP::pi*2.;
		G4double costheta = G4UniformRand();
		G4double theta = acos(costheta);
		G4double xDirection = sin(theta)*cos(phi);
		G4double yDirection = sin(theta)*sin(phi);
		G4double zDirection = costheta;
		const G4ThreeVector momentumDirection = G4ThreeVector(xDirection,yDirection,zDirection);
		fParticleGun->SetParticleMomentumDirection(momentumDirection);
	} else {
		G4ThreeVector momentumDirection = G4ThreeVector(0,0,0);
		fParticleGun->SetParticleMomentumDirection(momentumDirection);
	}
	fParticleGun->SetParticlePosition(position);
	//###################################################

	evtPrimAction->SetSourceX((position.x())/mm);
	evtPrimAction->SetSourceY((position.y())/mm);
	evtPrimAction->SetSourceZ((position.z())/mm);
	
	//###################################################
	// GENERATE PRIMARY VERTEX
	//##########################
	fParticleGun->GeneratePrimaryVertex(anEvent);
	//###################################################

	if(anEvent->GetEventID()==1) {  //stampo informazioni sorgente solo per il primo evento
		G4cout<<"Dimensioni sorgente: Raggio interno = "<<fRadiusInt<<", Raggio esterno = "<<fRadiusExt<<", H = "<<fZ<<G4endl;
		if (fSourceSelect==3) { //solo se è la sorgente ExtY..
			G4cout<<"TBR richiesto= "<<fTBR<<G4endl;
			G4cout<<"VolA= "<<VolA<<", ProbA= "<<ProbA<<G4endl;
			G4cout<<"VolB= "<<VolB<<", ProbB= "<<ProbB<<G4endl;
			G4cout<<"VolC= "<<VolC<<", ProbC= "<<ProbC<<G4endl;
			G4cout<<"Volume sorgente tot= "<<VolA+VolB+VolC<<G4endl;
		}
	}
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

