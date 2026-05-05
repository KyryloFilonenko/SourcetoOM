#ifndef MICD_HH
#define MICD_HH

// MiHeaders
#include "MiCDCaloHit.h"
#include "MiCDTrackerHit.h"

// ROOT headers
#include "TObject.h"

// Standard headers
#include "vector"

using namespace std;

class MiCD: public TObject
{
	public:
		//! Constructor
		MiCD();

		//! Destructor
		~MiCD();
		
		int addcalohit(MiCDCaloHit& in_calo);

		int getnoofcaloh();

		MiCDCaloHit*  getcalohit(int in_no);

		void print();		

		//// From Filip
		int  getnooftrackerh();
		void addtrackerhit(MiCDTrackerHit& in_tracker);
		vector<MiCDTrackerHit>*  gettrackerhitv();

	private:

		int noOfCaloHits;

		vector<MiCDCaloHit> calohit;

		int noOfTrackerHits;
		vector<MiCDTrackerHit> trackerhit;

	ClassDef(MiCD,1);		
};

#endif // MICD_HH
