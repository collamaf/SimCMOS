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
// $Id: B1SteppingAction.cc 74483 2013-10-09 13:37:06Z gcosmo $
//
/// \file B1SteppingAction.cc
/// \brief Implementation of the B1SteppingAction class

#include "B1SteppingAction.hh"
#include "B1EventAction.hh"
#include "B1RunAction.hh"
#include "B1DetectorConstruction.hh"


#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

#include "B1Analysis.hh"



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1SteppingAction::B1SteppingAction(B1EventAction* eventAction, B1RunAction* runAction, G4double CollHoleDiam)
: G4UserSteppingAction(),
fEventAction(eventAction),
fScoringVolume(0),
runStepAction(runAction),
fCollHoleDiam(CollHoleDiam)
{}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1SteppingAction::~B1SteppingAction()
{}

//std::ofstream pixelOut("PixelTest.dat", std::ios::out);

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1SteppingAction::UserSteppingAction(const G4Step* step)
{
	
	G4VPhysicalVolume* ThisVol = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
	G4VPhysicalVolume* NextVol = step->GetPostStepPoint()->GetTouchableHandle()->GetVolume();
	
	G4int debug=0;
	G4int EvtNumber = G4EventManager::GetEventManager()->GetNonconstCurrentEvent()->GetEventID();
	
	G4bool newDebug = false;
	
#pragma mark What Enters CMOS
	// ###########################################################################
	// ###################### ENTERING CMOS
	
	//	if((NextVol && ThisVol->GetName()=="DummyCMOS" && NextVol->GetName()=="CMOS")) { //what enters CMOS from front after 2018.05.29
	if((NextVol && ThisVol->GetName()!="CMOS" && NextVol->GetName()=="CMOS")) { //what enters CMOS from wherever - 2018.10.10 - "Good" one for Efficiencies
		if (debug) G4cout<<"\nSTEPDEBUG\n Particella entrata in CMOS da dummy - fEventAction->GetEnteringParticle() ERA = "<<fEventAction->GetEnteringParticle();
		
		fEventAction->SetEnteringParticle(step->GetTrack()->GetDynamicParticle() ->GetPDGcode());
		if (debug) G4cout<<" SETTO fEventAction->GetEnteringParticle()= "<<fEventAction->GetEnteringParticle()<<G4endl<<G4endl;
		
		if (fEventAction->GetStoreTrackIDCmos()==step->GetTrack()->GetTrackID()) { //if I already saw this track entering CMOS...
			fEventAction->AddPassCounterCmos(1);  //increase the counter
			if (debug) G4cout<<"STEPDEBUG!!! Già entrata!!! "<<step->GetTrack()->GetTrackID()<<G4endl;
		}else {
			fEventAction->SetStoreTrackIDCmos(step->GetTrack()->GetTrackID());
		}

		if (fEventAction->GetStoreEventIDCmosPrim() ==  EvtNumber) {
			fEventAction->AddPassCounterCmosPrim(1);
			if (newDebug) G4cout<<"CIAONE!!! evento gia visto dare segnale in COS!!! "<<EvtNumber<<G4endl;
			
		} else {
			fEventAction->SetStoreEventIDCmosPrim(EvtNumber);
			if (newDebug) G4cout<<"CIAONE!!! Nuovo evento che dà segnale in CMOS, Incremento contatore !!! "<<EvtNumber<<G4endl;
		}
		
		// Salvo le info solo della prima volta che una particella entra nel CMOS
		if (fEventAction->GetPassCounterCmos()==0) {
			G4double eKinPre = step->GetPostStepPoint()->GetKineticEnergy();
			//Fill vector
			(runStepAction->GetRunEnPre()).push_back(eKinPre/keV);
			fEventAction->AddNoPre(1); //update the counter of particles entering CMOS in the event
			(runStepAction->GetRunPart()).push_back(step->GetTrack()->GetDynamicParticle() ->GetPDGcode()); //add PID of particle enetering CMOS
			(runStepAction->GetRunPreCmosX()).push_back(step->GetPostStepPoint()->GetPosition().x()/mm);
			(runStepAction->GetRunPreCmosY()).push_back(step->GetPostStepPoint()->GetPosition().y()/mm);
			(runStepAction->GetRunPreCmosZ()).push_back(step->GetPostStepPoint()->GetPosition().z()/mm);
			
			(runStepAction->GetRunPreCmosTrackID()).push_back(fEventAction->GetPassCounterCmosPrim());
		}
		
		if (fEventAction->GetPassCounterCmos()==0 && fEventAction->GetPassCounterCmosPrim()==0) {
			(runStepAction->GetRunEnPrePrim()).push_back(fEventAction->GetSourceEne()/keV);
		}
	}
	// ###################### END ENTERING CMOS
	// ###########################################################################
	
#pragma mark What Enters Filter

	// ###########################################################################
	// ###################### ENTERING FILTER
	
	if((NextVol && ThisVol->GetName()=="World" && NextVol->GetName()=="Resin")) { //what enters Resin
		if (debug) G4cout<<"STEPDEBUG!!! Entrato in Resin!!! "<<G4endl;
		G4double eKinPre = step->GetPostStepPoint()->GetKineticEnergy();
		
		runStepAction->GetRunPreFilterEn().push_back(eKinPre/keV);
		runStepAction->GetRunPreFilterPart().push_back(step->GetTrack()->GetDynamicParticle() ->GetPDGcode());
		fEventAction->AddNoPreFilter(1); //update the counter of particles entering RESIN in the event
	}
	
	// ###################### END ENTERING FILTER
	// ###########################################################################
	
	// ###########################################################################
	// ###################### EXITING SOURCE i.e.
	//Modified on 2017-11-17 by collamaf: now the condition works for both cases: with or without Cu collimator.
	//If there is not collimator save what goes from source to dummy. If there is a collimator save what goes from world (the hole) into dummy
	if( NextVol && ( (fCollHoleDiam<0 &&  ( (ThisVol->GetName()=="Source" && NextVol->GetName()=="Dummy"))) || ( (fCollHoleDiam>=0 &&   (ThisVol->GetName()=="CuCollimator" && NextVol->GetName()=="Dummy") ) )) ) { //what actually exits the source
		
		if (debug) G4cout<<"STEPDEBUG!!! Uscito da source!!! "<<G4endl;
		
		//collamaf: to avoid double counting same track going back and forth, check if I already counted it
		if (fEventAction->GetStoreTrackIDSource()==step->GetTrack()->GetTrackID()) { //if I already saw this track exiting the source...
			fEventAction->AddPassCounterSource(1);  //increase the counter
			if (debug) G4cout<<"STEPDEBUG!!! Già contata!!! "<<G4endl;
		}else {
			fEventAction->SetStoreTrackIDSource(step->GetTrack()->GetTrackID());
		}
		
		// Salvo le info solo della prima volta che una particella esce dalla sorgente
		if (fEventAction->GetPassCounterSource()==0) {
			if (debug) G4cout<<"STEPDEBUG!!! Nuova traccia uscente la sorgente!!! "<<G4endl;
			//			if (step->GetPostStepPoint()->GetKineticEnergy()/keV<50) G4cout<<"STEPDEBUG!!! Nuova traccia uscente la sorgente con <50keV!!! "<<G4endl;
			
			fEventAction->AddNSourceExit(1);
			(runStepAction->GetRunEnExit()).push_back(step->GetPostStepPoint()->GetKineticEnergy()/keV);
			(runStepAction->GetRunXExit()).push_back(step->GetPostStepPoint()->GetPosition().x()/mm);
			(runStepAction->GetRunYExit()).push_back(step->GetPostStepPoint()->GetPosition().y()/mm);
			(runStepAction->GetRunZExit()).push_back(step->GetPostStepPoint()->GetPosition().z()/mm);
			(runStepAction->GetRunCosXExit()).push_back(step->GetPreStepPoint()->GetMomentumDirection().x());
			(runStepAction->GetRunCosYExit()).push_back(step->GetPreStepPoint()->GetMomentumDirection().y());
			(runStepAction->GetRunCosZExit()).push_back(step->GetPreStepPoint()->GetMomentumDirection().z());
			(runStepAction->GetRunPartExit()).push_back(step->GetTrack()->GetDynamicParticle() ->GetPDGcode());
			(runStepAction->GetRunParentIDExit()).push_back(step->GetTrack()->GetParentID());
			if (step->GetTrack()->GetCreatorProcess()) (runStepAction->GetRunExitProcess().push_back((step->GetTrack()->GetCreatorProcess()->GetProcessType())));
			if (debug && step->GetTrack()->GetParentID()>3) G4cout<<"CONTROLLA!!!!! TrackID= "<<step->GetTrack()->GetTrackID()<<" ParentID= "<<step->GetTrack()->GetParentID()<<G4endl;
		}
		/*
		 We have to use PreStepPoint to save the exit cosines, otherwise we already have particles flipped..
		 */
	}
	
	if (!fScoringVolume) {
		const B1DetectorConstruction* detectorConstruction
		= static_cast<const B1DetectorConstruction*>
		(G4RunManager::GetRunManager()->GetUserDetectorConstruction());
		fScoringVolume = detectorConstruction->GetScoringVolume();
	}

	
	// get volume of the current step
	G4LogicalVolume* volume
	= step->GetPreStepPoint()->GetTouchableHandle()
	->GetVolume()->GetLogicalVolume();
	
	// ###########################################################################
	// ###################### If I want to keep some kind of event for visualization
	
	G4Event* evt = G4EventManager::GetEventManager()->GetNonconstCurrentEvent();
	if (ThisVol->GetName()=="DummyCMOS" /*&& step->GetTrack()->GetDynamicParticle() ->GetPDGcode()==22*/ && (NextVol && NextVol->GetName()=="CMOS")) {
		evt->KeepTheEvent();
	}
	// ###########################################################################

	if (debug) { //Debug time stuff
		G4cout<<"CIAONE GlobalT [ns]= "<<std::setprecision(30) << step->GetTrack()->GetGlobalTime()/ns<<G4endl;
		G4cout<<"GlobalT [ns]= "<<step->GetTrack()->GetGlobalTime()/ns<<G4endl;
		G4cout<<" GetPrimDecayTime [ns]= " << fEventAction->GetPrimDecayTime()/ns<<G4endl;
		G4cout<<"LocalT= " <<  G4BestUnit(step->GetTrack()->GetLocalTime(),"Time") <<G4endl;
	}
	
	// ###########################################################################
	// ###################### INSIDE CMOS - Per each hit into sensitive detector
	// check if we are in scoring volume
	if (volume== fScoringVolume) { //fScoringVolume is "logicPix"
		
		//pixel information collection
		G4int CopyNB=step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber();
		fEventAction->AddNo(1);
		
		G4ThreeVector pixCenter;
		G4TouchableHandle touchHandle =step->GetPreStepPoint()->GetTouchableHandle();
		G4ThreeVector vec_origin(0.,0.,0.);
		G4ThreeVector globalPos = touchHandle->GetHistory()-> GetTopTransform().Inverse().TransformPoint(vec_origin);
		pixCenter = globalPos;
		
		if (CopyNB>0) {
			//fill vectors
			(runStepAction->GetRunPixNo()).push_back(CopyNB);
			(runStepAction->GetRunPixXpos()).push_back(pixCenter.getX()/mm);
			(runStepAction->GetRunPixYpos()).push_back(pixCenter.getY()/mm);
		}
		
		// collect energy deposited in this step
		G4StepPoint* postPoint = step->GetPostStepPoint();
		G4double edepStep = step->GetTotalEnergyDeposit();
		G4ThreeVector post=postPoint->GetPosition();
		
		//Fill vector
		(runStepAction->GetRunEnCmos()).push_back(step->GetTotalEnergyDeposit()/keV);
		(runStepAction->GetRunEnCmosPrim()).push_back(runStepAction->GetMotherEnergy());
		//		(runStepAction->GetRunEnCmosTime()).push_back(step->GetTrack()->GetLocalTime()/ns);
		//		(runStepAction->GetRunEnCmosTime()).push_back(step->GetTrack()->GetGlobalTime()/ns-runStepAction->GetMotherTime());
		(runStepAction->GetRunEnCmosTime()).push_back(step->GetTrack()->GetGlobalTime()/ns-fEventAction->GetPrimDecayTime());
		//		G4cout<<"CMOSDEBUG  MotherTime= "<< runStepAction->GetMotherTime()<<" PostDiff= "<<  step->GetTrack()->GetGlobalTime()/ns-runStepAction->GetMotherTime() <<G4endl;
		(runStepAction->GetRunXCmos()).push_back(step->GetPreStepPoint()->GetPosition().x()/mm);
		(runStepAction->GetRunYCmos()).push_back(step->GetPreStepPoint()->GetPosition().y()/mm);
		(runStepAction->GetRunZCmos()).push_back(step->GetPreStepPoint()->GetPosition().z()/mm);
		(runStepAction->GetRunPartCmos()).push_back(step->GetTrack()->GetDynamicParticle() ->GetPDGcode());
		if (debug)  G4cout<<"CIAODEBUG ! Evento n: "<< G4EventManager::GetEventManager()->GetNonconstCurrentEvent()->GetEventID() <<"Ho un rilascio di energia ("<< step->GetTotalEnergyDeposit()/keV<<" [KeV]) dovuto ad una particella entrata nel CMOS di tipo: "<<fEventAction->GetEnteringParticle()<<G4endl;
		if (debug && edepStep>0)  G4cout<<"CIAODEBUG MAGGIORE DI 0!! "<<G4endl;
		
		//Collect deposited energy in CMOS  due to Sr electons
		if (fEventAction->GetEnteringParticle()==11) {  //if son of electron
			fEventAction->AddEdepEle(step->GetTotalEnergyDeposit());
		}
		else if (fEventAction->GetEnteringParticle()==-11) {  //if son of positron
			fEventAction->AddEdepPos(step->GetTotalEnergyDeposit());
		} else if (fEventAction->GetEnteringParticle()==22) {  //if son of photon
			if (debug&& step->GetTotalEnergyDeposit()>0) G4cout<<"CONTROLLA RILASCIO DA FOTONE! Evento n "<<evt->GetEventID()<<G4endl;
			fEventAction->AddEdepFot(step->GetTotalEnergyDeposit());
			if (debug&&step->GetTotalEnergyDeposit()>0) G4cout<<"CONTROLLA"<<G4endl;
		}
		
		fEventAction->AddEdep(edepStep);
	}
	// ###########################################################################

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

