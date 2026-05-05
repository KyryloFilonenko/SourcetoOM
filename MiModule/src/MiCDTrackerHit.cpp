// MiHeaders
#include "MiCDTrackerHit.h"

ClassImp(MiCDTrackerHit);

MiCDTrackerHit::MiCDTrackerHit()
{
}

MiCDTrackerHit::~MiCDTrackerHit()
{
}
	
double MiCDTrackerHit::getX()
{
	return X;
}

double MiCDTrackerHit::getY()
{
	return Y;
}

double MiCDTrackerHit::getZ()
{
	return Z;
}

double MiCDTrackerHit::getR()
{
	return R;
}

int MiCDTrackerHit::setX(double in_X)
{
	X = in_X;
	return 0;
}

int MiCDTrackerHit::setY(double in_Y)
{
	Y = in_Y;
	return 0;
}

int MiCDTrackerHit::setZ(double in_Z)
{
	Z = in_Z;
	return 0;
}

int MiCDTrackerHit::setR(double in_R)
{
	R = in_R;
	return 0;
}
