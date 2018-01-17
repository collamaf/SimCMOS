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
// $Id: B1ActionInitialization.hh 68058 2013-03-13 14:47:43Z gcosmo $
//
/// \file B1ActionInitialization.hh
/// \brief Definition of the B1ActionInitialization class

#ifndef B1ActionInitialization_h
#define B1ActionInitialization_h 1

#include "G4VUserActionInitialization.hh"
#include "globals.hh"
/// Action initialization class.

class B1ActionInitialization : public G4VUserActionInitialization
{
public:
	B1ActionInitialization(G4double, G4double, G4double, G4int, std::ofstream &, G4double/*, G4bool*/, G4int);
	virtual ~B1ActionInitialization();
	
	virtual void BuildForMaster() const;
	virtual void Build() const;
	
protected:
	G4double fX0Scan;
	G4double fZValue;
	G4double fCuDiam;
	G4int fFilterFlag;
	std::ofstream &FilePrimaries;
	G4double fTBR;
//	G4bool fSrSourceFlag;
	G4int fSourceSelect;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


