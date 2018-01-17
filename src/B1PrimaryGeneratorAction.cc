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
	G4bool fPointLike=true;
	G4bool fExtended=false;
	G4bool fSTB=false; 
	
	if (fSourceSelect==1) {  //pointlike Sr
		fPointLike=true;
		fExtended=false;
		fSTB=false;
	} else if (fSourceSelect==2) { //extended Sr
		fPointLike=false;
		fExtended=true;
		fSTB=false;
	} else if (fSourceSelect==3) { //DOTA Sr
		fPointLike=false;
		fExtended=false;
		fSTB=true;
	}
	
	if (fSTB) {
	fRadiusInt=3*mm;
	fDZInt=1*mm;
	fRadiusExt=8*mm;
	fDZExt=2*mm;
	} else if (fPointLike) {
		fRadiusInt=0*mm;
		fDZInt=0*mm;
		fRadiusExt=0*mm;
		fDZExt=0*mm;
	} else if (fExtended) {
		fRadiusInt=10.5*mm;  //8 for RM, 10.5mm PG source
		fDZInt=0*mm;
		fRadiusExt=10.5*mm;
		fDZExt=0*mm;
	}
	
	
	ofstream SourceFile;
	
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
	
	//Stronzium
	G4int Z = 38, A = 90;
	G4double ionCharge   = 0.*eplus;
	G4double excitEnergy = 0.*keV;
	
	G4ParticleDefinition* ion
	= G4IonTable::GetIonTable()->GetIon(Z,A,excitEnergy);
	fParticleGun->SetParticleDefinition(ion);
	fParticleGun->SetParticleCharge(ionCharge);
	
	
	G4double VolInt=CLHEP::pi*fRadiusInt*fRadiusInt*fDZInt;
	G4double VolExt=CLHEP::pi*(fRadiusExt*fRadiusExt*fDZExt-fRadiusMin*fRadiusMin*fDZInt);
	fRatio=VolInt/VolExt;
	

	G4double zSource=0;
	G4double zSourceOffset=1e-6; //to avoid generating particles at the very boundary of source!

	
	if (fRadiusExt==fRadiusInt) { //se ho un solo raggio ignoro il TBR
		fRadiusMax=fRadiusInt;
		fRadiusMin=0*mm;
		zSource = -zSourceOffset;
	} else {
		if (G4UniformRand()<1/(fRatio*fTBR+1) || fTBR<0) {  //parte esterna
			if(G4UniformRand()<.5) { //in metà dei casi faccio l'anello esterno
				fRadiusMax=fRadiusExt;
				fRadiusMin=fRadiusInt;
				fZ=fDZExt;
				zSource = -G4UniformRand()*fZ-zSourceOffset;
			} else { //nell'altra metà faccio il cilindretto sotto la parte centrale
				fRadiusMax=fRadiusInt;
				fRadiusMin=0*mm;
				fZ=fDZExt-fDZInt;
				zSource = -G4UniformRand()*fZ-zSourceOffset;
			}
		} else {                           //parte interna
			fRadiusMax=fRadiusInt;
			fRadiusMin=0*mm;
			fZ=fDZInt;
			zSource = -G4UniformRand()*fZ-fDZInt-zSourceOffset;
		}
	}

	fParticleGun->SetParticleEnergy(0*MeV); //SetParticleEnergy seems use kinetic energy
	evtPrimAction->SetSourceEne(fParticleGun->GetParticleEnergy());
	
	G4double rho = sqrt(fRadiusMin*fRadiusMin + G4UniformRand()*(fRadiusMax*fRadiusMax-fRadiusMin*fRadiusMin));   //fixed square problem by collamaf with internal radius!
	G4double alpha = G4UniformRand()*CLHEP::pi*2.;
//	zSource = -G4UniformRand()*fZ;

	const G4ThreeVector position = G4ThreeVector(rho*cos(alpha), rho*sin(alpha), zSource);

	evtPrimAction->SetSourceCosX(0);
	evtPrimAction->SetSourceCosY(0);
	evtPrimAction->SetSourceCosZ(0);
	
	G4ThreeVector momentumDirection = G4ThreeVector(0,0,0);
	
	fParticleGun->SetParticleMomentumDirection(momentumDirection);
	
	evtPrimAction->SetSourceX((position.x())/mm);
	evtPrimAction->SetSourceY((position.y())/mm);
	evtPrimAction->SetSourceZ((position.z())/mm);
	
	fParticleGun->SetParticlePosition(position);
	
	fParticleGun->GeneratePrimaryVertex(anEvent);
	
	
	if(anEvent->GetEventID()==1) {  //stampo informazioni sorgente
		G4cout<<"SIMULATA SORGENTE DI SR90 "<<G4endl;
		G4cout<<"Dimensioni sorgente: Raggio interno = "<<fRadiusInt<<", Raggio esterno = "<<fRadiusExt<<", H = "<<fZ<<G4endl;
		G4cout<<"TBR richiesto= "<<fTBR<<G4endl;
		G4cout<<"Volume sorgente interna= "<<VolInt<<G4endl;
		G4cout<<"Volume sorgente esterna= "<<VolExt<<G4endl;
		G4cout<<"Volume sorgente tot= "<<VolExt+VolInt<<G4endl;
		G4cout<<"Rapporto Volume int/ext = "<<fRatio<<G4endl;
		G4cout<<"Condizione 1/(fRatio*fTBR+1) = "<<1/(fRatio*fTBR+1)<<G4endl;
	}
	
	
	// ################################################################################
	
	
	
}



G4double  B1PrimaryGeneratorAction::BetaDecaySpectrum(G4double Ek, G4double EndPoint)
{
	G4double ElMassMev = 0.510998928*MeV;
	
	G4double res=0.;
	
	G4double omega= Ek /ElMassMev  +1.;
	G4double omegamax= EndPoint /ElMassMev +1.;
	if(omega>1 && omega<omegamax)
	{
		res=(omegamax-omega)*(omegamax-omega) *Ek*sqrt(omega*omega-1.);
	}
	return res;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

