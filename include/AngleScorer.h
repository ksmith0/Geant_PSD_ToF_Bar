#ifndef ANGLESCORER_H
#define ANGLESCORER_H

#include "G4VPrimitiveScorer.hh"
#include "G4THitsMap.hh"

//Forward declarations
class G4HCofThisEvent;
class G4Step;
class G4TouchableHistory;

class AngleScorer : public G4VPrimitiveScorer
{
	private:
      G4int  fHCID;
      G4THitsMap<G4double>* fEvtMap;

	public:
		AngleScorer(G4String name, G4int depth=0);
		virtual ~AngleScorer() {};
      virtual void Initialize(G4HCofThisEvent*);
      virtual void EndOfEvent(G4HCofThisEvent*);

	protected:
		virtual G4bool ProcessHits(G4Step*, G4TouchableHistory*);
};

#endif
