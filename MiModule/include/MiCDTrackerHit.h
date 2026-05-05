#ifndef MICDTRACKERHIT_HH
#define MICDTRACKERHIT_HH

// ROOT headers
#include "TObject.h"

using namespace std;

class MiCDTrackerHit: public TObject
{
	public:
		//! Constructor
		MiCDTrackerHit();

		//! Destructor
		~MiCDTrackerHit();

		double getX();
		double getY();
		double getZ();
		double getR();

		int setX(double in_X);
		int setY(double in_Y);
		int setZ(double in_Z);
		int setR(double in_R);
		
	private:

		double X;
		double Y;
		double Z;
		double R;

	ClassDef(MiCDTrackerHit,1);		
};

#endif 
