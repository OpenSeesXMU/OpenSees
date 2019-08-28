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
                                                                        
                                                                        
#ifndef corotationalTruss2D_h
#define corotationalTruss2D_h

// Written: fmk 
// Created: 07/98
// Revision: A
//
// Description: This file contains the class definition for corotationalTruss2D. A corotationalTruss2D object
// provides the abstraction of the small deformation bar element. Each corotationalTruss2D
// object is associated with a material object. This corotationalTruss2D element will work
// in 1d, 2d or 3d problems.
//
// What: "@(#) corotationalTruss2D.h, revA"

#include <Element.h>
#include <Matrix.h>

class Node;
class Channel;
class UniaxialMaterial;

class corotationalTruss2D : public Element
{
  public:
    corotationalTruss2D(int tag, int dimension,
	  int Nd1, int Nd2, 
	  UniaxialMaterial &theMaterial,
	  double A);
    
    corotationalTruss2D();    
    ~corotationalTruss2D();

    const char *getClassType(void) const {return "corotationalTruss2D";};

    // public methods to obtain information about dof & connectivity    
    int getNumExternalNodes(void) const;
    const ID &getExternalNodes(void);
    Node **getNodePtrs(void);

    int getNumDOF(void);	
    void setDomain(Domain *theDomain);

    // public methods to set the state of the element    
    int commitState(void);
    int revertToLastCommit(void);        
    int revertToStart(void);        
    int update(void);
	const Matrix &getInitialStiff(void);
    const Matrix &getTangentStiff(void);
    const Vector &getResistingForce(void);
	const Matrix &calculateTRBM(void);
	const Matrix &transT(Matrix matr);
    // public methods for element output
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    int displaySelf(Renderer &, int mode, float fact, const char **displayModes=0, int numModes=0);
    void Print(OPS_Stream &s, int flag =0);    

    Response *setResponse(const char **argv, int argc, OPS_Stream &s);
    int getResponse(int responseID, Information &eleInformation);


  protected:
    
  private:
    double computeCurrentStrain(void) const;
    double computeCurrentStrainRate(void) const;
    
    // private attributes - a copy for each object of the class
    UniaxialMaterial *theMaterial;  // pointer to a material
    ID  connectedExternalNodes;     // contains the tags of the end nodes
    int dimension;                  // corotationalTruss2D in 2 or 3d domain
    int numDOF;	                    // number of dof for corotationalTruss2D

    Matrix *theMatrix;  // pointer to objects matrix (a class wide Matrix)
    Vector *theVector;  // pointer to objects vector (a class wide Vector)

    double Lo;               // length of corotationalTruss2D based on undeformed configuration
	double Ln;
    double A;               // area of corotationalTruss2D

    double cosX[3];  // direction cosines

    Node *theNodes[2];
    double *initialDisp;

	double cosBeta;
	double sinBeta;
    // static data - single copy for all objects of the class	
    static Matrix corotationalTruss2DM4;   // class wide matrix for 4*4
    static Vector corotationalTruss2DV4;   // class wide Vector for size 4

};

#endif
