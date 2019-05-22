# Full simulation of CMOS beta- probe in Geant4

## HOW TO RUN:
```
cd build
cmake -DGeant4_DIR=$G4INSTALL ../
make
./exampleB1

./exampleB1 -AbsD -1 -Z 0.010 -Fil 0 -Source 8 -Sensor 2 -PixT 1.75 -NPrim 1000000
./exampleB1 -AbsD -1 -Z 0.445 -Fil 1 -Source 8 -Sensor 2 -PixT 1.75 -NPrim 1000000

```
### Source Choice:\
1 - Pointlike Sr\
2 - Extended Sr\
3 - ExtY\
4 - Co60\
5 - Na22\
6 - Ba133\
7 - Cs137\
8 - Flat Ele 0-3MeV\
9 - Flat Gamma 0-1MeV\
10 - Na22 "nude" (PG)

### Sensor Choice:
1 - MT9V011
2 - MT9V115
3 - bare 60035 SiPm



## GEOMETRY (maybe no more so accurate)
- Extended  Source ending at Z=0
- Cu collimator on top of Sr source (toggleble)
- CMOS Detector starting at Z offset (Z distance is from source surface to possible resin in case of sensor 1, or up to the sensor in case of sensor 2, even if with filter)
- Sensor resin filter in contact with CMOS (towards source)
- Dummy volume for scoring purposes between source (or if present CU Collimator) and World to score what exits the primary generation
- "infinite" carrier volume behind CMOS to simulate mechanical support


## PHYSICS
Process				Type		SubType
RadioActiveDecay	6			210
eIoni				2			2
eBrem				2			3

###CPU TIMES NEEDED FOR 1e5 PRIMARIES:

#### To send the sim to FARM:
from build/
```
rsync -avzh --exclude '*.root' --exclude '*.dat' ../* collamaf@farm-login.roma1.infn.it:NewLife/CMOS/

```


## OUTPUT:
A root file named CMOSmc_{XX}.root is created, reporting the several parameters used for the run, in which on an event (i.e. a primary particle) by event basis it is stored:

### SOURCE vector (one entry per primary particle, only first 100k events are written for disk space sake):
- AllX: X coordinate of primary particle [mm];
- AllY: Y coordinate of primary particle [mm];
- AllZ: Z coordinate of primary particle [mm];
- AllCosX[]: X directive cosine of produced electron;
- AllCosY[]: Y directive cosine of produced electron;
- AllCosZ[]: Z directive cosine of produced electron;
- AllEne[]: kinetic energy of produced particle  [keV];
- AllIsotope[]: parentID - 1 of track: so for example for Sr source is 0 if track is son of Sr, 1 if of Y;
- ExitX[]: X coordinate of primary particle exiting the source volume [mm]; ("Exiting" means going from source to dummy or from absorber to dummy if there is an absorber)
- ExitY[]: Y coordinate of primary particle exiting the source volume [mm];
- ExitZ[]: Z coordinate of primary particle exiting the source volume [mm];
- ExitCosX[]: X directive cosine of primary particle exiting the source volume;
- ExitCosY[]: Y directive cosine of primary particle exiting the source volume;
- ExitCosZ[]: Z directive cosine of primary particle exiting the source volume;
- ExitEne[]: kinetic energy of primary particle exiting the source volume [keV];
- ExitPart[]: kind of primary particle (11=e-, -11=e+, 22=gamma, 13=mu-...) exiting the source volume;
- ExitParentID[]: partent-id of particle exiting the source
- ExitProcess[]: process that created the particles that exits the source (see table above)
- ExitTrackN: number of different tracks exiting the source per event

### B1 vector (one entry per primary particle):
- Eabs: energy absorbed in CMOS [keV];
- EAbsComp[2]: vector containing energy absorbed in CMOS [keV] due to Sr (comp 1) and to Y (comp 2)
- PreFilterTrackN: number of tracks entering the Filter (if present, empty otherwise);
- PreFilterPart: kind of particle entering the Filter (if present, empty otherwise);
- PreFilterEn: energy of particle entering the Filter (if present, empty otherwise);
- PreCmosTrackN: number of tracks entering Cmos per primary (from front resin) (it's the length of the following vector);
- PreCmosPart[PreCmosTrackN]: kind of particle of each track entering Cmos (from front resin);
- PreCmosEn[PreCmosTrackN]: kinetic energy of particle of each tracks entering Cmos (from front resin) [keV];
- PreCmosX[PreCmosTrackN]: x position of each tracks entering Cmos (from front resin) [mm];
- PreCmosY[PreCmosTrackN]: y position of each tracks entering Cmos (from front resin) [mm];
- PreCmosZ[PreCmosTrackN]: z position of each tracks entering Cmos (from front resin) [mm];
- PreCmosEnPrim[PreCmosTrackN]: energy of the primary particle that originated the track that is now entering CMOS (from front resin) [keV];
- PreCmosEventNum[PreCmosTrackN]: number of times that the event of the primary particle that originated the track that is now entering CMOS (from front resin)  gave a track entering CMOS (useful to check for double countings) ;
- InCmosTrackN: number of hits inside Cmos (it's the length of the following vector);
- InCmosPart[InCmosTrackN]: kind of particle of hit inside Cmos;
- InCmosEn[InCmosTrackN]: energy deposit of single hit of particle inside Cmos;
- InCmosEnPrim[InCmosTrackN]: energy of the primary particle that origined the hit of particle inside Cmos [keV];
- InCmosTime[InCmosTrackN]: time of interaction of hit inside Cmos [ns] (To be really undersood);
- InCmosX[InCmosTrackN]: X position of hit inside Cmos [mm];
- InCmosY[InCmosTrackN]: Y position of hit inside Cmos [mm];
- InCmosZ[InCmosTrackN]: Z position of hit inside Cmos [mm];
- PixelID[InCmosTrackN]: number of pixel in which the hit occurred (from 1 to NpixMax);
- PixXPos[InCmosTrackN]: x position of the pixel in which the hit occurred [mm];
- PixYPos[InCmosTrackN]: y position of the pixel in which the hit occurred [mm];
- SourceX: X coordinate of primary particle (isotope) giving a signal in Cmos [mm];
- SourceY: Y coordinate of primary particle (isotope) giving a signal in Cmos [mm];
- SourceZ: Z coordinate of primary particle (isotope) giving a signal in Cmos [mm];
- SourceCosX[]: X directive cosine of decay electron(s) giving a signal in Cmos;
- SourceCosY[]: Y directive cosine of  decay electron(s) giving a signal in Cmos;
- SourceCosZ[]: Z directive cosine of decay electron(s) giving a signal in Cmos;
- SourceEne[]: kinetic energy of produced particle  [keV];
- SourceIsotope: parentID - 1 of track: so for example for Sr source is 0 if track is son of Sr, 1 if of Y;
- Nev: storing number of events generated


For example, to see energy spectrum of electrons created by Sr/Y that exit the source
```
Source->Draw("ExitEne","ExitPart==11&&ExitProcess==6")
```


## Pass The MC output to data-analysis algorithm

da CMOS/CodiceCMOS
```

Riduzione /Users/francesco/MonteCarlo/Sonda/SimCMOS/build/CMOSmcX0_Z173_NOCuD_Fil1_TBR10_ExtSr_115_Frame100.root -noise 2018-04-20_MT9V115_stronzioRM22gradi_0000_noise_100.root -frameSize 488x648 -t 7 -mc


Riduzione /Users/francesco/MonteCarlo/Sonda/SimCMOS/build/CMOSmc_X0_Z44_NoAbs_Fil1_FlatEle_115_PxT175_N1000000_Frame1800.root -noise 2018-04-20_MT9V115_stronzioRM22gradi_0000_noise_100.root -seedSize 9 -edge 4 -checkLocalMaximumSide 9

Riduzione /Users/francesco/MonteCarlo/Sonda/SimCMOS/build/CMOSmc_X0_Z44_NoAbs_Fil1_FlatEle_115_PxT175_N1000000_Frame1800.root -noise 2018-04-20_MT9V115_stronzioRM22gradi_0000_noise_100.root



DataAnalysis /Users/francesco/MonteCarlo/Sonda/SimCMOS/build/CMOSmc_X0_Z50_NoAbs_Fil1_ExtSr_115Eff_Frame800_Reduced.root -frameSize 488x648  -mc


```
## CREATE EFFICIENCY CURVES with macro CreaEff.C

cd /build/

### GammaFilter
```
 .X CreaEff.C("../build/CMOSmc_X0_Z1_NoAbs_Fil0_FlatGamma_115_PxT175_N10000000","../build/CMOSmc_X0_Z44_NoAbs_Fil1_FlatGamma_115_PxT175_N10000000_Frame1800_Analized",2)
 ```

### GammaNoFilter
 ```
.X CreaEff.C("../build/CMOSmc_X0_Z1_NoAbs_Fil0_FlatGamma_115_PxT175_N10000000","../build/CMOSmc_X0_Z1_NoAbs_Fil0_FlatGamma_115_PxT175_N10000000_Frame1800_Analized",2)
```

### LowEneGammaNoFilter
```
.X CreaEff.C("../build/CMOSmc_X0_Z1_NoAbs_Fil0_FlatGamma_115_PxT175_LowEne_N10000000","../build/CMOSmc_X0_Z1_NoAbs_Fil0_FlatGamma_115_PxT175_LowEne_N10000000_Frame1800_Analized",3)
```

### LowEneGammaFilter
```
.X CreaEff.C("../build/CMOSmc_X0_Z1_NoAbs_Fil0_FlatGamma_115_PxT175_LowEne_N10000000","../build/CMOSmc_X0_Z44_NoAbs_Fil1_FlatGamma_115_PxT175_LowEne_N10000000_Frame1800_Analized",3)
```

### EleFilter
```
.X CreaEff.C("../build/CMOSmc_X0_Z1_NoAbs_Fil0_FlatEle_115_PxT175_N10000000","../build/CMOSmc_X0_Z44_NoAbs_Fil1_FlatEle_115_PxT175_N10000000_Frame1800_Analized",1)
```

### EleNoFilter
```
.X CreaEff.C("../build/CMOSmc_X0_Z1_NoAbs_Fil0_FlatEle_115_PxT175_N10000000","../build/CMOSmc_X0_Z1_NoAbs_Fil0_FlatEle_115_PxT175_N10000000_Frame1800_Analized",1)
```





## CHANGELOG
2017.12.1 by collamaf
- Try to fix problem of primary particles double counting by putting a check in Stepping Action. Seems to reduce by about ~9% the number of exiting particles (NoCudZ2 test)

2017.12.4 by collamaf
- Added reset of pass counter in stacking action. Seems to increase back of about another 5% the number of exiting particles. now we have about 96% of the beginning. Verified this should be the right approach.

2017.12.12 by collamaf
- Corrected z_resin:  no need to add 1mm of copper since distance is always from source top

2018.01.17 by collamaf
- Deep reorganization of DetConstr, much clearer now
- Code extended to both sensors smoothly (new argument to be passed by terminal, default sensor 1).
- Fixed ExtY source problem

2018.01.30 by collamaf
- Quite deep reorganization of output root file, cleaned a little in ordering and name (and updated readme)
- Changed "DOTA" to "ExtY" in output file naming

2018.04.23 by collamaf
- Now resin is always present in front of CMOS, if flag not selected made of air (useful for scoring)
- Added double crossing check also for particles entering CMOS

2018.04.26 by collamaf
- Added InCmosEnPrim to bring primary particle info to Riduzione and DataAnalysis

2018.05.7 by collamaf
- first implementation of storage of time of interection. Still not clear which time to save...
- Introduced possibility to simulate bare SiPm (assumed to be a particular version of CMOS detector): SensorChoice=3

2018.05.29 by collamaf
- Added dummy volume in front of CMOS (towards source) to score what enters it after some problem in case of space between the CMOS and the filter (like in sensor 2-115)
- Now in SteppingAction to score what enters CMOS we use the new Dummy volume
- Modified size of resin in case of sensor 2 (is larger than the active area), according to datasheet
- Fixed error in positioning of filter in case of sensor 1 (was after cmos)
- Fixed why filter in Sensor2 was not working.. the unit for density=4 was missing, chissa come interpretava lui...

2018.05.30 by collamaf
- Fixed an error in the positioning of the carrier behind CMOS
- Added as argument to pass the thickness of the absorber
- Added new parts from SimDaVinci (StackingAction: kill neutrinos and change in primary particle storing, fileName generated in main)
- Added as argument to pass material of absorber (1: Cu, 2:Al, 3:PVC)
- Now if run in vis mode the file root ends with xxxTEST not to overwrite important ones
- Now arguments are taken with labels, not necessary to give them all! If a macro is provided no visualization is init.. cool!

2018.06.01 by collamaf
- Added "QuickFlag" in main to automatically scale number of pixel if I want to visualize geometry
- ".dat" was missing after Primaries file
- Fixed problem with positioning of dummy after source (or possible absorber)
- Added QuickFlag to be passed as an argument to program to call for scaled down geometry
- Added structure to classify energy release due to e+/e-/gamma

2018.06.06 by collamaf
- Added scoring of what enters the filter (PreFilterEn/TrackN/Part) - THERE IS SOME PROBLEM! If made by hand empty the spectrum of particles hitting the filter changes... TO BE UNDERSTOOD

2018.07.25 by collamaf
- Added Co60 source {SourceChoice= 4}
- Whole sources volumes reorganisation
- Added Na, Ba and Cs sources {SourceChoice= 5, 6, 7}

2018.08.10 by collamaf
- Gamma sources have now 1mm diameter (active spot)
- Added possibility to pass "label" to filename from command line
- No more need for macros, now number of primaries is calculated from source activity
- Pixel thickness is now passed as argument

2018.09.17 by collamaf
- Finally working for efficiency!!!!

2019.05.08 by collamaf
- Addedd Source 10, Na22 "nude", so that the beta+ exits

2019.05.21 by collamaf
- README cleaning



## TO DO's


