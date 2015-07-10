#include "DetectorConstruction.h"

//Material manager
#include "G4NistManager.hh"
//Basic Units.
#include "G4SystemOfUnits.hh"

//Basic Volume and Placement
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

//Types of Volumes
#include "G4Box.hh"
#include "G4Tubs.hh"

//layer visible attributes
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

//Include necessary scorers
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"

//Optical photon surfaces
#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"

//Scoring
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PSNofSecondary.hh"
#include "G4SDParticleFilter.hh"
#include "OpticalDetection.h"

DetectorConstruction::DetectorConstruction() :
	G4VUserDetectorConstruction()
{ }

DetectorConstruction::~DetectorConstruction()
{ }

/**This routine is used to construct the geometry to be simulated. This includes the
 * necessary materials to produce objects.
 */
G4VPhysicalVolume *DetectorConstruction::Construct()
{
	// Option to switch on/off checking of volumes overlaps
	G4bool checkOverlaps = false;

	G4bool fUseOptical = true;

	//-------------------------------------------------------
	//------------ Material Definitions ---------------------
	//-------------------------------------------------------

	// Get nist material manager
	G4NistManager* nist = G4NistManager::Instance();
	G4Material* aluminum = nist->FindOrBuildMaterial("G4_Al");
	//G4Material* stainlessSteel = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");

	//------------ Air --------------------------------------
	G4Material* air = nist->FindOrBuildMaterial("G4_AIR");

	if (fUseOptical) {
		//Create Material Properties Table to store optical properties
		G4MaterialPropertiesTable* air_MPT = new G4MaterialPropertiesTable();

		//Assume constant index of refraction
		const G4int nEntries = 2;
		G4double air_Energy[nEntries] = {3.82495*eV, 2.91425*eV};
		G4double air_rIndex[nEntries] = {1.0002853,1.0002785};
		air_MPT->AddProperty("RINDEX", air_Energy, air_rIndex, nEntries);

		air->SetMaterialPropertiesTable(air_MPT);
	}	

	//--------------- Teflon ---------------------------------
	G4Material* teflonMaterial = nist->FindOrBuildMaterial("G4_TEFLON");

	G4OpticalSurface *teflonSurface = NULL;
	if (fUseOptical) {
		G4MaterialPropertiesTable * teflon_MPT = new G4MaterialPropertiesTable();

		//\note Need to find a good source of teflon optical properties.
		//Assume constant index of refraction
		const G4int nEntries = 2;
		G4double teflon_Energy[nEntries] = {3.82495*eV, 2.91425*eV};

		G4double teflon_rIndex[nEntries] = {1.35,1.35};
		G4double teflon_absorption[nEntries] = {2.*mm,2.*mm};
		teflon_MPT->AddProperty("RINDEX", teflon_Energy, teflon_rIndex, nEntries);
		teflon_MPT->AddProperty("ABSLENGTH", teflon_Energy, teflon_absorption, nEntries);

		teflonMaterial->SetMaterialPropertiesTable(teflon_MPT);

		//Create surface properties
		G4MaterialPropertiesTable *teflonSurface_MPT = new G4MaterialPropertiesTable();
		G4double teflon_reflectivity[nEntries] = {0.98, 0.98};
		teflonSurface_MPT->AddProperty("REFLECTIVITY", teflon_Energy, teflon_reflectivity, nEntries);

		//Create the optical surface
		teflonSurface = new G4OpticalSurface("teflonSurface");
		teflonSurface->SetType(dielectric_dielectric);
		//GroundFrontPainted finish specifies only reflection or absorption. 100% Lambertian
		teflonSurface->SetFinish(groundfrontpainted);
		teflonSurface->SetMaterialPropertiesTable(teflonSurface_MPT);
	}


	//--------------- Glass Window ---------------------------
	G4Material* pyrex = nist->FindOrBuildMaterial("G4_Pyrex_Glass");

	if (fUseOptical) {
		//Create Material Properties Table to store optical properties
		G4MaterialPropertiesTable *glassMTP = new G4MaterialPropertiesTable();

		const G4int nEntries = 2;
		G4double LaBr_Energy[nEntries] = {3.82495*eV, 2.91425*eV};
		//Index of reflection
		G4double rIndex_Glass[nEntries]=
			{1.49,1.49};
		glassMTP->AddProperty("RINDEX",LaBr_Energy,rIndex_Glass,nEntries);
		//Absorbstion length
		G4double absLength_Glass[nEntries]={420.*cm,420.*cm};
		glassMTP->AddProperty("ABSLENGTH",LaBr_Energy,absLength_Glass,nEntries);

		pyrex->SetMaterialPropertiesTable(glassMTP);
	}

	//--------------- LaBr3(Ce) ------------------------------
	G4Element* La = nist->FindOrBuildElement("La");
	G4Element* Br = nist->FindOrBuildElement("Br");
	G4Element* Ce = nist->FindOrBuildElement("Ce"); 

	G4Material* LaBr3 = new G4Material("LaBr3",5.07*g/cm3,2);
	LaBr3->AddElement(La,1);
	LaBr3->AddElement(Br,3);

	//Create LaBr3(Ce) material by doping LaBr3
	G4Material* scintMaterial = new G4Material("LaBr3_Ce",5.08*g/cm3, 2);
	scintMaterial->AddMaterial(LaBr3, 95*perCent);
	scintMaterial->AddElement(Ce, 5*perCent);

	if (fUseOptical) {
		//Create Material Properties Table to store optical properties
		G4MaterialPropertiesTable* scintMaterial_MPT = new G4MaterialPropertiesTable();

		// Data from 10.1109/TNS.2012.2193597
		const G4int nEntries = 21;
		G4double LaBr_Energy[nEntries] =  
		{3.82495*eV, 3.76611*eV, 3.70904*eV, 3.65368*eV, 3.59995*eV, 
			3.54778*eV, 3.4971*eV, 3.44784*eV, 3.39996*eV, 3.35338*eV, 
			3.30807*eV, 3.26396*eV, 3.22101*eV, 3.17918*eV, 3.13842*eV, 
			3.09869*eV, 3.05996*eV, 3.02218*eV, 2.98533*eV, 2.94936*eV, 
			2.91425*eV};

		//Fast compoennt
		G4double LaBr_Fast[nEntries]  =  
		{0.00171744, 0.0169987, 0.104353, 0.413358, 1.09586, 2.01313, 2.65712, 
			2.64391, 2.16732, 1.73403, 1.60549, 1.65549, 1.63166, 1.41242, 
			1.05117, 0.674215, 0.375993, 0.184151, 0.0799886, 0.0310985, 
			0.0109151};     
		scintMaterial_MPT->AddProperty("FASTCOMPONENT",LaBr_Energy, LaBr_Fast, nEntries);

		//Index of reflection
		G4double LaBr_RIND[nEntries] =
			{2.53548, 2.50933, 2.48529, 2.4631, 2.44257, 2.42352, 2.40581, 
			2.38929, 2.37386, 2.35942, 2.34587, 2.33314, 2.32116, 2.30986, 
			2.2992, 2.28912, 2.27958, 2.27054, 2.26196, 2.25381, 2.24605} ;
		scintMaterial_MPT->AddProperty("RINDEX",LaBr_Energy, LaBr_RIND,    nEntries);

		//Absorbstion length
		G4double LaBr_ABSL[nEntries] = 
			{0.00214455*cm, 0.00485183*cm, 0.0130231*cm, 0.0408024*cm, 
			0.147022*cm, 0.601092*cm, 2.75427*cm, 13.9861*cm, 77.9037*cm, 
			471.533*cm, 3074.95*cm, 21435.2*cm, 158585.*cm, 1.23706E6*cm, 
			1.01134E7*cm, 8.61758E7*cm, 7.61483E8*cm, 6.94559E9*cm, 
			6.51156E10*cm, 6.25021E11*cm, 6.12049E12*cm};
		scintMaterial_MPT->AddProperty("ABSLENGTH", LaBr_Energy, LaBr_ABSL, nEntries);

		G4double LaBr_RayleighAbs[nEntries] = 
			{8.736*cm, 9.035*cm, 9.334*cm, 9.633*cm, 9.932*cm, 10.231*cm, 
			10.53*cm, 10.829*cm, 11.128*cm, 11.427*cm, 11.726*cm, 12.025*cm, 
			12.324*cm, 12.623*cm, 12.922*cm, 13.221*cm, 13.52*cm, 13.819*cm, 
			14.118*cm, 14.417*cm, 14.716*cm};
		scintMaterial_MPT->AddProperty("RAYLEIGH", LaBr_Energy, LaBr_RayleighAbs, nEntries);

		//Following is a calculated quantity that most likely underestiamtes the abs length.
		//G4double LaBr_WLSABS[nEntries] = 
		//	{0.0123851*cm, 0.00283096*cm, 0.0012378*cm, 0.000979046*cm, 
		//	 0.00133068*cm, 0.00296152*cm, 0.0102811*cm, 0.0524681*cm, 
		//	 0.356516*cm, 2.69712*cm, 18.9965*cm, 128.423*cm, 963.998*cm, 
		//	 8687.01*cm, 95426.1*cm, 1.26774E6*cm, 2.00875E7*cm, 
		//	 3.74093E8*cm, 8.07423E9*cm, 1.99343E11*cm, 5.56163E12*cm};
		//scintMaterial_MPT->AddProperty("WLSABSLENGTH", LaBr_Energy, LaBr_ABSL, nEntries);

		//63 photons / keV
		scintMaterial_MPT->AddConstProperty("SCINTILLATIONYIELD",6.3/keV);
		scintMaterial_MPT->AddConstProperty("RESOLUTIONSCALE", 0.71);
		scintMaterial_MPT->AddConstProperty("FASTTIMECONSTANT", 16.*ns);

		scintMaterial->SetMaterialPropertiesTable(scintMaterial_MPT);
	}

	//--------------- PMT ------------------------------------
	//PMT cathode is made of Super Bialkali from Hamatsu R6231-100
	//Just use bialkali for now
	G4Element* Sb = nist->FindOrBuildElement("Sb"); 
	G4Element* Rb = nist->FindOrBuildElement("Rb"); 
	G4Element* Cs = nist->FindOrBuildElement("Cs"); 

	G4Material* pmtMaterial = new G4Material("SuperBialkali",0,3);
	pmtMaterial->AddElement(Sb,1);
	pmtMaterial->AddElement(Rb,1);
	pmtMaterial->AddElement(Cs,1);

	G4OpticalSurface* opSurfPMT;

	if (fUseOptical) {
		//--------------- PMT Optical Properties -----------------
		opSurfPMT = new G4OpticalSurface("opSurfPMT");
		opSurfPMT->SetModel(unified);
		//Must be dielectric metal for detection
		opSurfPMT->SetType(dielectric_metal);
		opSurfPMT->SetFinish(polished);

		//Define optical surface optical properties 
		G4MaterialPropertiesTable* opSurfPMT_MPT = new G4MaterialPropertiesTable();
		const G4int numEntries = 22;
		G4double pmtEnergies[numEntries] = 
			{4.531417*eV, 4.214736*eV, 4.059972*eV, 3.870475*eV, 3.677380*eV, 3.502637*eV, 3.326978*eV,
			 3.168097*eV, 3.065617*eV, 2.891889*eV, 2.794437*eV, 2.692368*eV, 2.587361*eV, 2.518601*eV,
			 2.417636*eV, 2.365904*eV, 2.308279*eV, 2.193781*eV, 2.083541*eV, 1.977935*eV, 1.882518*eV,
			 1.748549*eV};

		G4double pmtReflectivity[numEntries] = {0}; 
		opSurfPMT_MPT->AddProperty("REFLECTIVITY",pmtEnergies,pmtReflectivity,numEntries);

		//Efficiencncy from http://www.hamamatsu.com/us/en/technology/innovation/photocathode/index.html
		//Quantum Efficiency 
		G4double pmtEfficiency[numEntries] = 
			{0.198864, 0.284091, 0.313920, 0.339489, 0.348011, 0.349432, 0.348011, 
			0.342330, 0.335227, 0.315341, 0.299716, 0.276989, 0.237216, 0.207386, 
			0.163352, 0.140625, 0.117898, 0.075284, 0.042614, 0.019886, 0.007102,
			0.000000};
			//{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
		//We can use the efficiecny as the reflectivity is 0
		opSurfPMT_MPT->AddProperty("EFFICIENCY",pmtEnergies,pmtEfficiency,numEntries);
		
		opSurfPMT->SetMaterialPropertiesTable(opSurfPMT_MPT);

	}

	G4Material* containerMaterial = aluminum;
	G4OpticalSurface *containerSurface = teflonSurface;

	// Construct the World
	G4double world_sizeXY = 1*m;
	G4double world_sizeZ  = 1*m;

	G4Box* solidWorld =    
		new G4Box("World",                       //its name
				0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);    //its size

	G4LogicalVolume* logicWorld =                         
		new G4LogicalVolume(solidWorld,          //its solid
				air,           //its material
				"World");      //its name

	G4VPhysicalVolume* physWorld = new G4PVPlacement(
			0, 					//no rotation
			G4ThreeVector(),	//at (0,0,0)
			logicWorld,			//contataing logical volume
			"World",				//name
			0,						//no mother volume
			false,				//no boolean operation
			0,						//copy number
			checkOverlaps);	//overlap checking

	//--------------Build the required geometry here-----------------//
	G4double in = 2.54 * cm;
	G4double ft = 12 * in;

	G4double detOffset = 2*cm;
	G4double barLength = 1*ft; 			//Active area of bar
	G4double barCrossSection = 2 * in; 	//Radius or length of side
	G4double containerThickness = 0.0065 * in; 		//Thickness of aluminum container hold active liquid
	G4double windowThickness = 1/8. * in;	//Thickness of window  
	G4double pmtDiameter = 51*mm;
	G4double pmtThickness = 5 * cm;

	//Computed dimensions
	G4double containerCrossSection = barCrossSection + 2 * containerThickness;
	G4double containerLength = barLength + 2 * windowThickness;
	G4double windowCrossSection = barCrossSection; //Crosssection of window

	//Build a volume of air to contain the detector.
	G4ThreeVector detectorPos = G4ThreeVector(0,0,detOffset + barCrossSection / 2.);
	G4RotationMatrix *detectorRot = new G4RotationMatrix();
	detectorRot->rotateX(90.*deg);

	G4Tubs *detector = new G4Tubs(
		"detector",
		0,
		containerCrossSection / 2., //Cylinder Radius
		(barLength + 2 * windowThickness + 2 * pmtThickness) / 2., //Half height
		0. * deg,
		360. * deg);

	G4LogicalVolume* detectorLogic = new G4LogicalVolume(
		detector,
		air,
		"detector");					

	new G4PVPlacement(
			detectorRot,
			detectorPos,
			detectorLogic,
			"detector",
			logicWorld,
			false,
			0,
			checkOverlaps);


	//Highlight the detector volume with a red wire frame.
	G4VisAttributes *detectorVis = new G4VisAttributes(G4Colour(1,0,0));
	detectorVis->SetForceWireframe(true);
	detectorVis->SetVisibility(false);
	detectorLogic->SetVisAttributes(detectorVis);


	//Build the active material container
	G4ThreeVector containerPos=G4ThreeVector(0,0,0);
	G4Tubs *container = 
		new G4Tubs("container",
				0,
				containerCrossSection / 2., //Cylinder radius
				containerLength / 2., //Half height
				0.*deg,
				360.*deg);		

	G4LogicalVolume* logicContainer =
		new G4LogicalVolume(container,
				containerMaterial,
				"alumExterior");						
	logicContainer->SetVisAttributes(new G4VisAttributes(G4Colour(0.75,0.75,0.75)));

	G4PVPlacement *physContainer = new G4PVPlacement(0,
			containerPos,
			logicContainer,
			"container",
			detectorLogic,
			false,
			0,
			checkOverlaps);					

	//Build the scintillation material.
	G4ThreeVector activeVolumePos = G4ThreeVector(0,0,0);

	G4Tubs *activeVolume = 
		new G4Tubs("activeVolume",
				0,
				barCrossSection / 2., //cylinder radius
				barLength / 2., //half height
				0.*deg,
				360.*deg);
	G4LogicalVolume* logicActiveVolume =
		new G4LogicalVolume(activeVolume,
				scintMaterial,
				"activeVolume");						

	G4PVPlacement *physActiveVolume = 
	new G4PVPlacement(0,
			activeVolumePos,
			logicActiveVolume,
			"activeVolume",
			logicContainer, 
			false,
			0,
			checkOverlaps);		

	logicActiveVolume->SetVisAttributes(new G4VisAttributes(G4Colour(1,0.33,1)));

	if (fUseOptical) {
		//Add optical border surface between LaBr3 and teflon.
		new G4LogicalBorderSurface(
				"ContainerActiveVolumeSurface",
				physActiveVolume,
				physContainer,
				containerSurface);
	}



	//Build the container windows
	G4ThreeVector windowPos = G4ThreeVector(0,0,containerLength/2. - windowThickness/2.);
	G4Tubs *window = new G4Tubs(
		"window",		//name
		0.,					//inner radius
		windowCrossSection/2.,		//outer radius
		windowThickness/2.,//half thickness
		0.*deg,				//starting angle
		360.*deg);			//ending angle

	G4LogicalVolume* logicWindow = new G4LogicalVolume(
		window,
		pyrex,
		"glassFace");						
	logicWindow->SetVisAttributes(new G4VisAttributes(G4Colour(0.33,1,1)));

	//Upper window placement
	G4PVPlacement *physGlassUpper = 
	new G4PVPlacement(
		0,
		windowPos,
		logicWindow,
		"window",
		logicContainer,
		false,
		0,
		checkOverlaps);

	//Lower window placement
	G4PVPlacement *physGlassLower = 
	new G4PVPlacement(
		0,
		-windowPos,
		logicWindow,
		"window",
		logicContainer,
		false,
		0,
		checkOverlaps);


	//--------------- PMT ------------------------------------
	G4ThreeVector pmtPos = G4ThreeVector(0,0,containerLength/ 2.+pmtThickness/2.);

	G4Tubs *pmt = new G4Tubs(
			"PMT",				//name
			0.,					//inner radius
			pmtDiameter/2,			//outer radius
			pmtThickness/2.,	//thickness
			0.*deg,				//starting angle
			360.*deg);			//ending angle

		/*G4Box *pmt = new G4Box(
			"PMT",				//name
			pmtRadius,			//outer radius
			pmtRadius,			//outer radius
			pmtThickness/2.);	//thickness
*/

		G4LogicalVolume* logicUpperPMT = new G4LogicalVolume(
			pmt,
			pmtMaterial,
			"UpperPMT");						

		G4LogicalVolume* logicLowerPMT = new G4LogicalVolume(
			pmt,
			pmtMaterial,
			"LowerPMT");						

	G4PVPlacement *physPMTUpper = 
			new G4PVPlacement(
				0,
				pmtPos,
				logicUpperPMT,
				"UpperPMT",
				detectorLogic,
				false,
				0,
				checkOverlaps);
	G4PVPlacement *physPMTLower = 
			new G4PVPlacement(
				0,
				-pmtPos,
				logicLowerPMT,
				"LowerPMT",
				detectorLogic,
				false,
				0,
				checkOverlaps);



	logicUpperPMT->SetVisAttributes(new G4VisAttributes(G4Colour(0.33,1,0.33)));
	logicLowerPMT->SetVisAttributes(new G4VisAttributes(G4Colour(0.33,1,0.33)));

	if (fUseOptical) {
		//Create the optical surface as the directional border from first
		//	physical volume to the second physical volume
		new G4LogicalBorderSurface("opSurfGlassPMT", physGlassUpper, physPMTUpper, opSurfPMT);
		new G4LogicalBorderSurface("opSurfGlassPMT", physGlassLower, physPMTLower, opSurfPMT);
	}

	//Return the world
	return physWorld;
}

/**In this routine we construct the Senditive Detector (SD). This is where the scorers are registered.
 */
void DetectorConstruction::ConstructSDandField() 
{
	G4MultiFunctionalDetector* activeVolume = new G4MultiFunctionalDetector("activeVolume");
	activeVolume->RegisterPrimitive(new G4PSEnergyDeposit("edep"));
	activeVolume->RegisterPrimitive(new G4PSNofSecondary("numSecondary"));
	SetSensitiveDetector("activeVolume",activeVolume);

	G4MultiFunctionalDetector* upperPMT = new G4MultiFunctionalDetector("UpperPMT");
	upperPMT->SetFilter(new G4SDParticleFilter("opPhotonFilter","opticalphoton"));
	upperPMT->RegisterPrimitive(new OpticalDetection("opDet"));
	SetSensitiveDetector("UpperPMT",upperPMT);

	G4MultiFunctionalDetector* lowerPMT = new G4MultiFunctionalDetector("LowerPMT");
	lowerPMT->SetFilter(new G4SDParticleFilter("opPhotonFilter","opticalphoton"));
	lowerPMT->RegisterPrimitive(new OpticalDetection("opDet"));
	SetSensitiveDetector("LowerPMT",lowerPMT);
}
