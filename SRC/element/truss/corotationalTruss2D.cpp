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
                                                                        
// $Revision$
// $Date$
// $URL$
                                                                        
                                                                        
// Written: fmk 
// Created: 07/98
// Revision: A
//
// Description: This file contains the implementation for the corotationalTruss2D class.
//
// What: "@(#) corotationalTruss2D.C, revA"

#include <corotationalTruss2D.h>
#include <Information.h>
#include <Parameter.h>

#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <UniaxialMaterial.h>
#include <Renderer.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <ElementResponse.h>

//#include <fstream>

// initialise the class wide variables

Matrix corotationalTruss2D::corotationalTruss2DM4(4,4);
Vector corotationalTruss2D::corotationalTruss2DV4(4);
// constructor:
//  responsible for allocating the necessary space needed by each object
//  and storing the tags of the corotationalTruss2D end nodes.

#include <elementAPI.h>
#define OPS_Export 

OPS_Export void *
OPS_corotationalTruss2D()
{
  Element *theElement = 0;

  int numRemainingArgs = OPS_GetNumRemainingInputArgs();

  if (numRemainingArgs < 4) {
    opserr << "Invalid Args want: element corotationalTruss2D $tag $iNode $jNode $sectTag \n";
    opserr << " or: element corotationalTruss2D $tag $iNode $jNode $A $matTag\n";
    return 0;	
  }


  int iData[3];
  double A = 0.0;
  int matTag = 0;
  int ndm = OPS_GetNDM();

  int numData = 3;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING invalid integer (tag, iNode, jNode) in element corotationalTruss2D " << endln;
    return 0;
  }

  numData = 1;
  if (OPS_GetDouble(&numData, &A) != 0) {
    opserr << "WARNING: Invalid A: element corotationalTruss2D " << iData[0] << 
      " $iNode $jNode $A $matTag\n";
    return 0;	
  }

  numData = 1;
  if (OPS_GetInt(&numData, &matTag) != 0) {
    opserr << "WARNING: Invalid matTag: element corotationalTruss2D " << iData[0] << 
      " $iNode $jNode $A $matTag\n";
    return 0;
  }

  UniaxialMaterial *theUniaxialMaterial = OPS_GetUniaxialMaterial(matTag);
    
  if (theUniaxialMaterial == 0) {
    opserr << "WARNING: Invalid material not found element corotationalTruss2D " << iData[0] << " $iNode $jNode $A " << 
      matTag << "\n";
    return 0;
  }
  
  // now create the corotationalTruss2D
  theElement = new corotationalTruss2D(iData[0], ndm, iData[1], iData[2], *theUniaxialMaterial, A);

  if (theElement == 0) {
    opserr << "WARNING: out of memory: element corotationalTruss2D " << iData[0] << 
      " $iNode $jNode $A $matTag\n";
  }

  return theElement;
}

corotationalTruss2D::corotationalTruss2D(int tag, int dim,
         int Nd1, int Nd2, 
         UniaxialMaterial &theMat,
         double a)
 :Element(tag,ELE_TAG_corotationalTruss2D),
  theMaterial(0), connectedExternalNodes(2),
  dimension(dim), numDOF(0),
   theMatrix(0), theVector(0),
  Lo(0.0), initialDisp(0)
{
    // get a copy of the material and check we obtained a valid copy
    theMaterial = theMat.getCopy();
    if (theMaterial == 0) {
      opserr << "FATAL corotationalTruss2D::corotationalTruss2D - " << tag <<
	"failed to get a copy of material with tag " << theMat.getTag() << endln;
      exit(-1);
    }
    
    // ensure the connectedExternalNode ID is of correct size & set values
    if (connectedExternalNodes.Size() != 2) {
      opserr << "FATAL corotationalTruss2D::corotationalTruss2D - " <<  tag << "failed to create an ID of size 2\n";
      exit(-1);
    }

    connectedExternalNodes(0) = Nd1;
    connectedExternalNodes(1) = Nd2;        

    // set node pointers to NULL
    for (int i=0; i<2; i++)
      theNodes[i] = 0;

    cosX[0] = 0.0;
    cosX[1] = 0.0;
    cosX[2] = 0.0;

	cosBeta = 0;
	sinBeta = 0;
}

// constructor:
//   invoked by a FEM_ObjectBroker - blank object that recvSelf needs
//   to be invoked upon
corotationalTruss2D::corotationalTruss2D()
:Element(0,ELE_TAG_corotationalTruss2D),     
 theMaterial(0),connectedExternalNodes(2),
 dimension(0), numDOF(0),
  theMatrix(0), theVector(0),
 Lo(0.0), A(0.0), initialDisp(0)
{
    // ensure the connectedExternalNode ID is of correct size 
  if (connectedExternalNodes.Size() != 2) {
      opserr << "FATAL corotationalTruss2D::corotationalTruss2D - failed to create an ID of size 2\n";
      exit(-1);
  }

  for (int i=0; i<2; i++)
    theNodes[i] = 0;

  cosX[0] = 0.0;
  cosX[1] = 0.0;
  cosX[2] = 0.0;


}

//  destructor
//     delete must be invoked on any objects created by the object
//     and on the matertial object.
corotationalTruss2D::~corotationalTruss2D()
{
    // invoke the destructor on any objects created by the object
    // that the object still holds a pointer to
    if (theMaterial != 0)
	delete theMaterial;

    if (initialDisp != 0)
      delete [] initialDisp;
}


int
corotationalTruss2D::getNumExternalNodes(void) const
{
    return 2;
}

const ID &
corotationalTruss2D::getExternalNodes(void) 
{
    return connectedExternalNodes;
}

Node **
corotationalTruss2D::getNodePtrs(void) 
{
  return theNodes;
}

int
corotationalTruss2D::getNumDOF(void) 
{
    return numDOF;
}


// method: setDomain()
//    to set a link to the enclosing Domain and to set the node pointers.
//    also determines the number of dof associated
//    with the corotationalTruss2D element, we set matrix and vector pointers,
//    allocate space for t matrix, determine the length
//    and set the transformation matrix.
void
corotationalTruss2D::setDomain(Domain *theDomain)
{
    // check Domain is not null - invoked when object removed from a domain
    if (theDomain == 0) {
	theNodes[0] = 0;
	theNodes[1] = 0;
	Lo = 0;
	return;
    }

    // first set the node pointers
    int Nd1 = connectedExternalNodes(0);
    int Nd2 = connectedExternalNodes(1);
    theNodes[0] = theDomain->getNode(Nd1);
    theNodes[1] = theDomain->getNode(Nd2);	
    
    // if can't find both - send a warning message
    if ((theNodes[0] == 0) || (theNodes[1] == 0)) {
      if (theNodes[0] == 0)
	opserr <<"corotationalTruss2D::setDomain() - corotationalTruss2D" << this->getTag() << " node " << Nd1 <<
	  "does not exist in the model\n";
      else
	opserr <<"corotationalTruss2D::setDomain() - corotationalTruss2D" << this->getTag() << " node " << Nd2 <<
	  "does not exist in the model\n";
      return;
    }


    // call the base class method
    this->DomainComponent::setDomain(theDomain);
    // now set the number of dof for element and set matrix and vector pointer

	numDOF = 4;
	theMatrix = &corotationalTruss2DM4;
	theVector = &corotationalTruss2DV4;	

    // now determine the length, cosines and fill in the transformation
    // NOTE t = -t(every one else uses for residual calc)
    const Vector &end1Crd = theNodes[0]->getCrds();
    const Vector &end2Crd = theNodes[1]->getCrds();	
    const Vector &end1Disp = theNodes[0]->getDisp();
    const Vector &end2Disp = theNodes[1]->getDisp();

    double dx = end2Crd(0)-end1Crd(0);
    double dy = end2Crd(1)-end1Crd(1);	
    
      if (initialDisp == 0) {
	double iDispX = end2Disp(0)-end1Disp(0);
	double iDispY = end2Disp(1)-end1Disp(1);
	if (iDispX != 0 || iDispY != 0) {
	  initialDisp = new double[2];
	  initialDisp[0] = iDispX;
	  initialDisp[1] = iDispY;
	  dx += iDispX;
	  dy += iDispY;
	}
      }
      
      Lo = sqrt(dx*dx + dy*dy);
	  Ln = Lo;
      if (Lo == 0.0) {
	opserr <<"WARNING corotationalTruss2D::setDomain() - corotationalTruss2D " << this->getTag() << " has zero length\n";
	return;
      }
	
      cosX[0] = dx/Lo;
      cosX[1] = dy/Lo;
    }

 	 


int
corotationalTruss2D::commitState()
{
  int retVal = 0;
  // call element commitState to do any base class stuff
  if ((retVal = this->Element::commitState()) != 0) {
    opserr << "corotationalTruss2D::commitState () - failed in base class";
  }    
  retVal = theMaterial->commitState();
  return retVal;
}

int
corotationalTruss2D::revertToLastCommit()
{
    return theMaterial->revertToLastCommit();
}

int
corotationalTruss2D::revertToStart()
{
    return theMaterial->revertToStart();
}

int
corotationalTruss2D::update(void)
{
	// Nodal displacements
	const Vector &end1Disp = theNodes[0]->getTrialDisp();
	const Vector &end2Disp = theNodes[1]->getTrialDisp();
	const Vector &crdNode1 = theNodes[0]->getCrds();
	const Vector &crdNode2 = theNodes[1]->getCrds();
	Vector currentCrd1 = crdNode1 + end1Disp;
	Vector currentCrd2 = crdNode2 + end2Disp;
	Ln = pow((currentCrd2(0) - currentCrd1(0)), 2) + pow((currentCrd2(1) - currentCrd1(1)), 2);
	
	// Compute engineering strain and strain rate
	double strain = (Ln - Lo) / Lo;

	//it need to calculate the rigid body rotation! beta!
	//beta=deltaUx2/(Lo+deltaUx1)
	//line funciton: using crdNode1,crdNode2
	//disp from point to line:
	double A = (crdNode2(1) - crdNode1(1)) / (crdNode2(0) - crdNode1(0));
	double C = (crdNode1(1)*crdNode2(0) - crdNode2(1)*crdNode1(0)) / (crdNode2(0) - crdNode1(0));
	double u2 = (A* currentCrd1(0) + currentCrd1(1) + C) / (sqrt(pow(A, 2) + 1));
	double u5 = (A* currentCrd2(0) + currentCrd2(1) + C) / (sqrt(pow(A, 2) + 1));

	double u2X = (currentCrd1(0) - A * C - A * currentCrd1(1)) / (pow(A, 2) + 1);
	double u2Y = -A*u2X-C;

	double u5X = (currentCrd2(0) - A * C - A * currentCrd2(1)) / (pow(A, 2) + 1);
	double u5Y = -A * u5X - C;
	
	double deltaUx1 = sqrt(pow(u5X - u2X, 2) + pow(u5Y - u2Y, 2));
	double deltaUx2 = u5 - u2;
	
	cosBeta = Lo + deltaUx1/Ln;
	sinBeta= deltaUx2/ Ln;
	// Set material trial strain
	return theMaterial->setTrialStrain(strain);
}


const Matrix &
corotationalTruss2D::getTangentStiff(void)
{
	Matrix &K = *theMatrix;
	K.Zero();
	
	static Matrix kgm(4, 4);
	static Matrix kgt(4, 4);
	kgt.Zero();
	kgm.Zero();
	// Material stiffness
	//
	// Get material tangent
	double EA = A * theMaterial->getTangent();
	EA /= Ln ;
	Vector &rbm = *theVector;
	rbm.Zero();
	rbm(0) = -cosBeta; rbm(1) = -sinBeta;
	rbm(2) = cosBeta; rbm(3) = sinBeta;

	//to check Lei!!!
	int i, j;
	for (i = 0; i < 2; i++)
		for (j = 0; j < 2; j++)
			kl(i, j) = EA * cosX[i] * cosX[j];

	// Geometric stiffness
	//
	// Get material stress
	double q = A * theMaterial->getStress();
	//double SA = q / (Ln*Ln*Ln);
	double SL = q / Ln;
	Matrix tRBM = this->calculateTRBM();
	Matrix rRBMT = this->transT(tRBM);

	kgt.addMatrix(1.0, tRBM, SL);


	
	


	// Compute R'*kl*R
	static Matrix kg(3, 3);
	kg.addMatrixTripleProduct(0.0, R, kl, 1.0);

	

	// Copy stiffness into appropriate blocks in element stiffness
	int numDOF2 = numDOF / 2;
	for (i = 0; i < 2; i++) {
		for (j = 0; j < 2; j++) {
			K(i, j) = kg(i, j);
			K(i, j + numDOF2) = -kg(i, j);
			K(i + numDOF2, j) = -kg(i, j);
			K(i + numDOF2, j + numDOF2) = kg(i, j);
		}
	}
	
	return *theMatrix;
}

const Matrix &
corotationalTruss2D::calculateTRBM(void)
{
	Matrix &rbm = *theMatrix;
	rbm.Zero();
	double s2 = pow(sinBeta, 2);
	double sc = sinBeta * cosBeta;
	double c2 = pow(cosBeta, 2);
	rbm(0, 0) = s2;  rbm(1, 0) = -sc; rbm(2, 0) = -s2;  rbm(1, 0) = sc;
	rbm(0, 1) = -sc; rbm(1, 1) = c2;  rbm(2, 1) = sc; rbm(1, 1) = -c2;
	rbm(0, 2) = -s2; rbm(1, 2) = sc;  rbm(2, 2) = s2; rbm(1, 2) = -sc;
	rbm(0, 3) = sc;	 rbm(1, 3) = -c2; rbm(2, 3) = -sc;	 rbm(1, 3) = c2;

	return rbm;
}

const Matrix &
corotationalTruss2D::transT(Matrix matr)
{
	Matrix  &tMatr = *theMatrix;
	tMatr.Zero();
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			tMatr(j, i) = matr(i, j);
		}
	}

	return tMatr;
}

const Vector &
corotationalTruss2D::getResistingForce()
{	
    if (Lo == 0.0) { // - problem in setDomain() no further warnings
	theVector->Zero();
	return *theVector;
    }
    
    // R = Ku - Pext
    // Ku = F * transformation
    double force = A*theMaterial->getStress();
    int numDOF2 = numDOF/2;
    double temp;
    for (int i = 0; i < dimension; i++) {
      temp = cosX[i]*force;
      (*theVector)(i) = -temp;
      (*theVector)(i+numDOF2) = temp;
    }


    
  return *theVector;
}


int
corotationalTruss2D::sendSelf(int commitTag, Channel &theChannel)
{
  int res;

  return 0;
}

int
corotationalTruss2D::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int res;
  int dataTag = this->getDbTag();

  // corotationalTruss2D creates a Vector, receives the Vector and then sets the 
  // internal data with the data in the Vector

  static Vector data(12);
  res = theChannel.recvVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr <<"WARNING corotationalTruss2D::recvSelf() - failed to receive Vector\n";
    return -1;
  }	      

  this->setTag((int)data(0));
  dimension = (int)data(1);
  numDOF = (int)data(2);
  A = data(3);

  initialDisp = new double[dimension];
  for (int i=0; i<dimension; i++)
    initialDisp[i] = 0.0;

  int initial = 0;
  for (int i=0; i<dimension; i++) {
    if (data(9+i) != 0.0) {
      initial = 1;
    }
  }
  
  if (initial != 0) {
    for (int i=0; i<dimension; i++) {
      initialDisp[i] = data(9+i);
    }    
  }
  
  // corotationalTruss2D now receives the tags of it's two external nodes
  res = theChannel.recvID(dataTag, commitTag, connectedExternalNodes);
  if (res < 0) {
    opserr <<"WARNING corotationalTruss2D::recvSelf() - " << this->getTag() << " failed to receive ID\n";
    return -2;
  }

  // finally corotationalTruss2D creates a material object of the correct type,
  // sets its database tag and asks this new object to recveive itself.

  int matClass = (int)data(4);
  int matDb = (int)data(5);

  // check if we have a material object already & if we do if of right type
  if ((theMaterial == 0) || (theMaterial->getClassTag() != matClass)) {

    // if old one .. delete it
    if (theMaterial != 0)
      delete theMaterial;

    // create a new material object
    theMaterial = theBroker.getNewUniaxialMaterial(matClass);
    if (theMaterial == 0) {
      opserr <<"WARNING corotationalTruss2D::recvSelf() - " << this->getTag() 
	<< " failed to get a blank Material of type " << matClass << endln;
      return -3;
    }
  }

  theMaterial->setDbTag(matDb); // note: we set the dbTag before we receive the material
  res = theMaterial->recvSelf(commitTag, theChannel, theBroker);
  if (res < 0) {
    opserr <<"WARNING corotationalTruss2D::recvSelf() - "<< this->getTag() << "failed to receive its Material\n";
    return -3;    
  }

  return 0;
}

int
corotationalTruss2D::displaySelf(Renderer &theViewer, int displayMode, float fact, 
		   const char **displayModes, int numModes)
{
  int res = 0;
  
  return res;
}



void
corotationalTruss2D::Print(OPS_Stream &s, int flag)
{
  
}

double
corotationalTruss2D::computeCurrentStrain(void) const
{
	// NOTE method will not be called if Lo == 0

	// determine the strain
	const Vector &disp1 = theNodes[0]->getTrialDisp();
	const Vector &disp2 = theNodes[1]->getTrialDisp();

	double dLength = 0.0;
	for (int i = 0; i < dimension; i++) {
	    dLength += (disp2(i) - disp1(i))*cosX[i];
    }
  
    // this method should never be called with Lo == 0
    return dLength/Lo;
}

double
corotationalTruss2D::computeCurrentStrainRate(void) const
{
    // NOTE method will not be called if Lo == 0

    // determine the strain
    const Vector &vel1 = theNodes[0]->getTrialVel();
    const Vector &vel2 = theNodes[1]->getTrialVel();	

    double dLength = 0.0;
    for (int i = 0; i < dimension; i++)
      dLength += (vel2(i)-vel1(i))*cosX[i];

    // this method should never be called with Lo == 0
    return dLength/Lo;
}

Response*
corotationalTruss2D::setResponse(const char **argv, int argc, OPS_Stream &output)
{

    Response *theResponse = 0;

    output.tag("ElementOutput");
    output.attr("eleType","corotationalTruss2D");
    output.attr("eleTag",this->getTag());
    output.attr("node1",connectedExternalNodes[0]);
    output.attr("node2",connectedExternalNodes[1]);

    //
    // we compare argv[0] for known response types for the corotationalTruss2D
    //


    if ((strcmp(argv[0],"force") == 0) || (strcmp(argv[0],"forces") == 0) 
        || (strcmp(argv[0],"globalForce") == 0) || (strcmp(argv[0],"globalForces") == 0)){
            char outputData[10];
            int numDOFperNode = numDOF/2;
            for (int i=0; i<numDOFperNode; i++) {
                sprintf(outputData,"P1_%d", i+1);
                output.tag("ResponseType", outputData);
            }
            for (int j=0; j<numDOFperNode; j++) {
                sprintf(outputData,"P2_%d", j+1);
                output.tag("ResponseType", outputData);
            }
            theResponse =  new ElementResponse(this, 1, Vector(numDOF));

    } else if ((strcmp(argv[0],"axialForce") == 0) || 
	       (strcmp(argv[0],"basicForce") == 0) || 
	       (strcmp(argv[0],"localForce") == 0) || 
	       (strcmp(argv[0],"basicForces") == 0)) {
            output.tag("ResponseType", "N");
            theResponse =  new ElementResponse(this, 2, Vector(1));

    } else if (strcmp(argv[0],"defo") == 0 || strcmp(argv[0],"deformation") == 0 ||
        strcmp(argv[0],"deformations") == 0 || strcmp(argv[0],"basicDefo") == 0 ||
        strcmp(argv[0],"basicDeformation") == 0 || strcmp(argv[0],"basicDeformations") == 0) {

            output.tag("ResponseType", "U");
            theResponse = new ElementResponse(this, 3, Vector(1));

    } else if (strcmp(argv[0],"basicStiffness") == 0) {

      output.tag("ResponseType", "K");
      theResponse = new ElementResponse(this, 4, Matrix(1,1));
	    
    // a material quantity
    } else if (strcmp(argv[0],"material") == 0 || strcmp(argv[0],"-material") == 0) {

        theResponse =  theMaterial->setResponse(&argv[1], argc-1, output);
    }

    output.endTag();
    return theResponse;
}

int 
corotationalTruss2D::getResponse(int responseID, Information &eleInfo)
{
  double strain, force;
    static Vector fVec(1);
    static Matrix kVec(1,1);

    switch (responseID) {
    case 1:
        return eleInfo.setVector(this->getResistingForce());

    case 2:
      fVec(0) = A*theMaterial->getStress();
      return eleInfo.setVector(fVec);

    case 3:
        if (Lo == 0.0) {
            strain = 0.0;
        } else {
            strain = theMaterial->getStrain();
        }
	fVec(0) = Lo*strain;
        return eleInfo.setVector(fVec);

    case 4:
      force = 0.0;
      if (Lo > 0.0)
	force = theMaterial->getTangent();
      kVec(0,0) = A*force/Lo;
      return eleInfo.setMatrix(kVec);
      
    default:
      return 0;
    }
}

const Matrix &
corotationalTruss2D::getInitialStiff(void)
{
	if (Lo == 0.0) { // - problem in setDomain() no further warnings
		theMatrix->Zero();
		return *theMatrix;
	}

	double E = theMaterial->getInitialTangent();

	// come back later and redo this if too slow
	Matrix &stiff = *theMatrix;

	int numDOF2 = numDOF / 2;
	double temp;
	double EAoverL = E * A / Lo;
	for (int i = 0; i < dimension; i++) {
		for (int j = 0; j < dimension; j++) {
			temp = cosX[i] * cosX[j] * EAoverL;
			stiff(i, j) = temp;
			stiff(i + numDOF2, j) = -temp;
			stiff(i, j + numDOF2) = -temp;
			stiff(i + numDOF2, j + numDOF2) = temp;
		}
	}

	return *theMatrix;
}
