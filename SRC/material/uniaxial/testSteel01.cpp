/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */

// $Revision: 1.6 $
// $Date: 2006-08-15 00:41:05 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/testSteel01.cpp,v $

// Written: MHS
// Created: Aug 2001
//
// Description: This file contains the class implementation for 
// testSteel01.

#include <testSteel01.h>
#include <Vector.h>
#include <Channel.h>
#include <math.h>
#include <float.h>

testSteel01::testSteel01(int tag,double Ep, double stressY, double bp)
  :UniaxialMaterial(tag,MAT_TAG_testSteel01),
   trialStrain(0.0), trialStress(0.0), trialTangent(0.0),
   commitStrain(0.0), commitStress(0.0), commitTangent(0.0)
{
	b = bp;
	stressYield = stressY;
	E = Ep;
}

testSteel01::testSteel01()
  :UniaxialMaterial(0,MAT_TAG_testSteel01),
   trialStrain(0.0), trialStress(0.0), trialTangent(0.0),
   commitStrain(0.0), commitStress(0.0), commitTangent(0.0)
{

}

testSteel01::~testSteel01()
{

}

int 
testSteel01::setTrialStrain(double strain, double strainRate)
{
	// elastic predictor
	trialStrain = strain;
	double dStrain = trialStrain - commitStrain;
	trialStress = commitStress + E * dStrain;
	trialTangent = E;


	// plastic corrector
	double ultimateElasticStrain = stressYield / E;
	
	if (trialStress > (stressYield + (trialStrain - ultimateElasticStrain)*b*E))
	{

		trialStress = stressYield + (trialStrain - ultimateElasticStrain)*b*E;
		trialTangent = b * E;
	}
	if (trialStress < -stressYield + (trialStrain + ultimateElasticStrain)*b*E)
	{
		trialStress = -stressYield + (trialStrain + ultimateElasticStrain)*b*E;
		trialTangent = b * E;
	}

  return 0;
}

double 
testSteel01::getStress(void)
{
  return trialStress;
}

double 
testSteel01::getTangent(void)
{
  return trialTangent;
}

double 
testSteel01::getInitialTangent(void)
{
  // return the initial tangent
  return 0.0;
}

double 
testSteel01::getStrain(void)
{
  return trialStrain;
}

int 
testSteel01::commitState(void)
{
  commitStrain  = trialStrain;
  commitStress  = trialStress;
  commitTangent = trialTangent;

  return 0;
}

int 
testSteel01::revertToLastCommit(void)
{
  trialStrain = commitStrain;
  trialStress = commitStress;
  trialTangent = commitTangent;

  return 0;
}

int 
testSteel01::revertToStart(void)
{
  trialStrain = 0.;
  trialStress = 0.0;
  trialTangent = 0.0;
  commitStrain = 0.;
  commitStress = 0.0;
  commitTangent = 0.0;

  return 0;
}

UniaxialMaterial *
testSteel01::getCopy(void)
{
  testSteel01 *theCopy = new testSteel01(this->getTag(),E,stressYield,b);

  return theCopy;
}

int 
testSteel01::sendSelf(int cTag, Channel &theChannel)
{
  return -1;
}

int 
testSteel01::recvSelf(int cTag, Channel &theChannel, 
			      FEM_ObjectBroker &theBroker)
{
  return -1;
}

void 
testSteel01::Print(OPS_Stream &s, int flag)
{
  s << "testSteel01 : " << this->getTag();

  return;
}


