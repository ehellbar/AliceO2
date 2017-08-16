#include "TStopwatch.h"
#include "TPCSimulation/AliTPCPoissonSolver.h"
#include "TPCSimulation/AliTPC3DCylindricalInterpolatorFull.h"
#include "TVector.h"
#include "TVectorD.h"
#include "TMatrix.h"
#include "TMatrixD.h"
#include "TDecompSVD.h"



/// \cond CLASSIMP3
ClassImp(o2::TPC::AliTPC3DCylindricalInterpolatorFull)
/// \endcond

using namespace o2::TPC;

/// constructor
AliTPC3DCylindricalInterpolatorFull::AliTPC3DCylindricalInterpolatorFull
(
	Int_t nr, 
	Int_t nz, 
	Int_t nphi,
	Int_t stepr,
	Int_t stepz, 
	Int_t stepphi, 
	Int_t type
)
{	
	//printf("creating 3D Cylindrical Interpolator: %d,%d,%d,%d,%d,%d,%d\n",nr,nz,nphi,stepr,stepz,stepphi,type);
	fOrder = 1;
	fIsAllocatingLookUp = kFALSE;
	fIsInitCubic = kFALSE;
	fMinZIndex = 0;
	fNR = nr;
	fNZ = nz;
	fNPhi = nphi;
	
	fRBFWeightLookUp = new Int_t[nr * nz * nphi];

	Int_t nd = stepr*stepz * stepphi;
	fStepR = stepr;
	fStepZ = stepz;
	fStepPhi = stepphi;

	fType = type;
	fRBFWeight = new Double_t[nr * nz * nphi*nd];
	for (Int_t i=0;i<nr * nz * nphi;i++) 	fRBFWeightLookUp[i] = 0;

	SetKernelType(kRBFInverseMultiQuadratic);
}

/// constructor
AliTPC3DCylindricalInterpolatorFull::AliTPC3DCylindricalInterpolatorFull()
{	
	fOrder = 1;
	fIsAllocatingLookUp = kFALSE;
	fIsInitCubic = kFALSE;
	fMinZIndex = 0;
}

/// destructor
AliTPC3DCylindricalInterpolatorFull::~AliTPC3DCylindricalInterpolatorFull()
{
	
	if  (fIsAllocatingLookUp) {		
		delete fVals;		
		//		delete fValsNormalized;
		delete fRlist;
		delete fPhilist;
		delete fZlist;
	}
	if (fIsInitCubic) {
		delete fSecondDerZ;	}

	delete[]  fRBFWeightLookUp;
	delete[]  fRBFWeight;
}

/// main operation
/// Int_terpolation at the region 
Double_t AliTPC3DCylindricalInterpolatorFull::Interpolate3DTableCylIDW
( 
	Double_t r,   
	Double_t z,
	Double_t phi,
	Int_t rindex,
	Int_t zindex,
	Int_t phiindex,
	Int_t stepR,
	Int_t stepPhi,
	Int_t stepZ
) 
{
  //  Int_t stepPhi = 3;
  // Int_t stepR    = 3;
  //Int_t stepZ    = 3;

  Double_t r0,z0,phi0,d;
  Double_t MIN_DIST = 1e-3;
  Double_t val = 0.0;
  
  // guess the index of r,z,phi
  


//  Search( fNZ,&fZlist[phiindex  * (fNZ * fNR) + (index_r * fNZ)], 1, z, jlowlist[index_phi * fNR + index_r] ) ;            
//  Search( fNPhi,&fPhilist[rindex * (fNZ) + zindex], (fNR * fNZ), phi, phiindex) ;            

  Int_t startPhi = phiindex - stepPhi/2;
  Int_t index_phi;

  Int_t startR =  rindex  - stepR/2;  
  Int_t startZ =  zindex  - stepZ/2;
  
  if (startPhi < 0) startPhi = fNPhi + startPhi;
//  if (startPhi + stepPhi >= fNR ) startPhi = startPhi fNR - stepR ;


  if (startR    < 0) startR = 0;
  if (startR + stepR >= fNR ) startR = fNR - stepR ;

  if (startZ    < fMinZIndex) startZ = fMinZIndex;
  if (startZ + stepZ  >= fNZ  ) startZ = fNZ - stepZ;

  Int_t index;
  Double_t sum_w = 0.0;
  Double_t sum_d  = 0.0;
  Double_t shortest_d = 10000.0; 
  Int_t new_rindex;
  Int_t new_zindex;
  Int_t new_phiindex;



  for (Int_t iphi = startPhi;iphi < startPhi + stepPhi;iphi++) 
    {
      index_phi = iphi % fNPhi;
      for (Int_t index_r= startR ;index_r< startR + stepR;index_r++)
	{
	  for (Int_t index_z=startZ; index_z < startZ + stepZ;index_z++) 
	    {
	      // check for the closest poInt_t
	      index = index_phi * (fNZ * fNR) + index_r * fNZ + index_z;
	
	      r0 =  fRlist[index];
	      z0 =  fZlist[index];
	      phi0 =  fPhilist[index];
	

	      d = Distance(r0,phi0,z0,r,phi,z);
	      if (d < shortest_d) 
		{
			shortest_d = d;
			new_rindex = index_r;
			new_phiindex = index_phi;
			new_zindex = index_z;	
	      	}
	      
	      
	    }
	}
    }

//  if ((new_phiindex != phiindex) ||(new_rindex != rindex) ||(new_zindex != zindex))
//	printf("(%d,%d,%d) == (%d,%d,%d)\n",rindex,phiindex,zindex,new_rindex,new_phiindex,new_zindex);
  stepPhi = 3;
  stepR = 3;
	
  
  startPhi = new_phiindex - stepPhi/2;
  
  startR =  new_rindex  - stepR/2;  
  startZ =  new_zindex  - stepZ/2;
  
  if (startPhi < 0) startPhi = fNPhi + startPhi;
//  if (startPhi + stepPhi >= fNR ) startPhi = startPhi fNR - stepR ;


  if (startR    < 0) startR = 0;
  if (startR + stepR >= fNR ) startR = fNR - stepR ;

  if (startZ    < fMinZIndex) startZ = fMinZIndex;
  if (startZ + stepZ  >= fNZ  ) startZ = fNZ - stepZ;


  for (Int_t iphi = startPhi;iphi < startPhi + stepPhi;iphi++) 
    {
      index_phi = iphi % fNPhi;
      for (Int_t index_r= startR ;index_r< startR + stepR;index_r++)
	{
	  for (Int_t index_z=startZ; index_z < startZ + stepZ;index_z++) 
	    {
	      // check for the closest poInt_t
	      index = index_phi * (fNZ * fNR) + index_r * fNZ + index_z;
	
	      r0 =  fRlist[index_phi * (fNR * fNZ) + index_r * fNZ + index_z];
	      z0 =  fZlist[index_phi * (fNR * fNZ) + index_r * fNZ + index_z];
	      phi0 =  fPhilist[index_phi * (fNR * fNZ) + index_r * fNZ + index_z];
	
//	      if (phi0 < 0) phi0 = TMath::TwoPi() + phi0;
////	      if (phi < 0) phi = TMath::TwoPi() + phi;
//	      if (phi0 > TMath::TwoPi()) phi0 = phi0 - TMath::TwoPi();
//	      if (phi > TMath::TwoPi()) phi = phi - TMath::TwoPi();

	      d = Distance(r0,phi0,z0,r,phi,z);
//	if (index_phi == 0) printf("d=%f\n",d);
	      if (d < MIN_DIST) {
		return fVals[index];
	      }
	      d = 1.0/d;
	      sum_w +=  (fVals[index] * d *d*d *d);
	      sum_d += d*d*d*d;

	      
	      
	    }
	}
  }



  return (sum_w / sum_d);

}

// distance in Cyl coordinate 
Double_t AliTPC3DCylindricalInterpolatorFull::Distance
(
 Double_t r0,
 Double_t phi0,
 Double_t z0,
 Double_t r,
 Double_t phi,
 Double_t z
)
{
  Double_t x0  = r0 * TMath::Cos(phi0);
  Double_t y0  = r0 * TMath::Sin(phi0);  
  Double_t x  = r * TMath::Cos(phi);
  Double_t y  = r * TMath::Sin(phi);



//  Double_t ret = ((x-x0) * (x-x0) +  (y-y0) * (y-y0) + (z - z0) * (z - z0));

/**  Double_t phi0Corrected = phi0 - phi;
  // normalize dphi

  while ( dphi < -TMath::TwoPi() ) dphi += TMath::TwoPi() ;                   // Table uses phi from 0 to 2*Pi
  while ( dphi > TMath::TwoPi() ) dphi = dphi - TMath::TwoPi() ;                   // Table uses phi from 0 to 2*Pi
  
**/

  if (phi < 0) phi = TMath::TwoPi() + phi;
  if (phi > TMath::TwoPi()) phi = phi - TMath::TwoPi();

  if (phi0 < 0) phi0 = TMath::TwoPi() + phi0;
  if (phi0 > TMath::TwoPi()) phi0 = phi0 - TMath::TwoPi();

  Double_t dphi = phi  - phi0;
  if (dphi > TMath::Pi())
		dphi = TMath::TwoPi() - dphi;
  if (dphi < -TMath::Pi())
		dphi = TMath::TwoPi() + dphi;
	
  Double_t ret = (r-r0) * (r-r0) + (dphi   * dphi) * ((r +r0)/2.0) + (z - z0) * (z - z0);
  return TMath::Sqrt(ret);
  //return ret;
}



/// main operation
/// Int_terpolation at the region 
Double_t AliTPC3DCylindricalInterpolatorFull::Interpolate3DTableCylRBFCartesian
( 
	Double_t r,   
	Double_t z,
	Double_t phi,
	Int_t rindex,
	Int_t zindex,
	Int_t phiindex,
	Int_t stepR,
	Int_t stepPhi,
	Int_t stepZ,
	Double_t radiusRBF0
) 
{
  

  const Float_t gridSizeR   =  (AliTPCPoissonSolver::fgkOFCRadius-AliTPCPoissonSolver::fgkIFCRadius) / (fNR - 1) ;
  const Float_t gridSizeZ   =  AliTPCPoissonSolver::fgkTPCZ0 / (fNZ -1) ;
  const Float_t gridSizePhi =  TMath::TwoPi() / fNPhi;
	
  //  Int_t stepPhi = 3;
  // Int_t stepR    = 3;
  //Int_t stepZ    = 3;

  Double_t r0,z0,phi0,d;
  Double_t MIN_DIST = 1e-3;
  Double_t val = 0.0;
  
  // guess the index of r,z,phi
  
  
  //
  

  //printf("radiusRBF0 = %f\n",radiusRBF0);
//  radiusRBF0 = 0.0; 
//  Search( fNZ,&fZlist[phiindex  * (fNZ * fNR) + (index_r * fNZ)], 1, z, jlowlist[index_phi * fNR + index_r] ) ;            
//  Search( fNPhi,&fPhilist[rindex * (fNZ) + zindex], (fNR * fNZ), phi, phiindex) ;            

  Int_t startPhi = phiindex - stepPhi/2;
  Int_t index_phi;

  Int_t startR =  rindex  - stepR/2;  
  Int_t startZ =  zindex  - stepZ/2;
  
  if (startPhi < 0) startPhi = fNPhi + startPhi;
//  if (startPhi + stepPhi >= fNR ) startPhi = startPhi fNR - stepR ;


  if (startR    < 0) startR = 0;
  if (startR + stepR >= fNR ) startR = fNR - stepR ;

  if (startZ    < fMinZIndex) startZ = fMinZIndex;
  if (startZ + stepZ  >= fNZ  ) startZ = fNZ - stepZ;

  Int_t index;
  Double_t sum_w = 0.0;
  Double_t sum_d  = 0.0;
  Double_t shortest_d = 10000.0; 
  Int_t new_rindex;
  Int_t new_zindex;
  Int_t new_phiindex;




  Double_t x = r * TMath::Cos(phi);
  Double_t y = r * TMath::Sin(phi);

  printf("(%f,%f,%f) => (%f,%f,%f)\n",r,phi,z,x,y,z);   	

  for (Int_t iphi = startPhi;iphi < startPhi + stepPhi;iphi++) 
  {
      index_phi = iphi % fNPhi;
      for (Int_t index_r= startR ;index_r< startR + stepR;index_r++)
	{
	  for (Int_t index_z=startZ; index_z < startZ + stepZ;index_z++) 
	    {
	      // check for the closest poInt_t
	      index = index_phi * (fNZ * fNR) + index_r * fNZ + index_z;
	
	      r0 =  fRlist[index];
	      z0 =  fZlist[index];
	      phi0 =  fPhilist[index];
	

	      d = Distance(r0,phi0,z0,x,y,z);
	      if (d < shortest_d) 
		{
			shortest_d = d;
			new_rindex = index_r;
			new_phiindex = index_phi;
			new_zindex = index_z;	
	      	}
	      
	      
	    }
	}
  }

  printf("(%d,%d,%d) => (%d,%d,%d)\n",rindex,phiindex,zindex,new_rindex,new_phiindex,new_zindex);
// do Int_terpolant


  stepPhi = fStepPhi;
  stepR = fStepR;
  stepZ = fStepZ;
//  if (stepZ > 3) stepZ =3;


  startPhi = new_phiindex - stepPhi/2;
  
  startR =  new_rindex  - stepR/2;  
  startZ =  new_zindex  - stepZ/2;
  
  if (startPhi < 0) startPhi = fNPhi + startPhi;
 

  if (startR    < 0) startR = 0;
  if (startR + stepR >= fNR ) startR = fNR - stepR ;

  if (startZ    < fMinZIndex) startZ = fMinZIndex;
  if (startZ + stepZ  >= fNZ  ) startZ = fNZ - stepZ;


  Double_t * w;	

  //Int_t nd = (stepPhi-1) + (stepR-1) + (stepZ-1) + 1;
  Int_t nd = stepPhi  *stepR  * stepZ ;

  w = new Double_t[nd];


	r0   = AliTPCPoissonSolver::fgkIFCRadius + (startR  * gridSizeR);

//  radiusRBF0 = Distance(r0,0.0,0.0,r0+gridSizeR,gridSizePhi,gridSizeZ);
  	//radiusRBF0 = 			fRadiusRBF0;
	//Distance(r0 * TMath::Cos(phi),r0 * TMath::Sin(phi),0.0,(r0 + gridSizeR) * TMath::Cos(phi + gridSizePhi),(r0 + gridSizeR) * TMath::Sin(phi + gridSizePhi),gridSizeZ);  

//	radiusRBF0 = iDistance(r0 * TMath::Cos(0.0),r0 * TMath::Sin(0.0),0.0,(r0+gridSizeR) * TMath::Cos(gridSizePhi),(r0 + gridSizeR) * TMath::Sin(gridSizePhi),gridSizeZ);
	radiusRBF0 = GetRadius0RBF(new_rindex,new_phiindex,new_zindex);
  if (fType == 1) {	

	for (Int_t i=0;i<nd;i++)  w[i]=0.0;


//	for (Int_t i=0;i<nd;i++) printf("w[%d]=%f\t",i,w[i]);
	
//	printf("get rbf weight\n");
	GetRBFWeightCartesian ( 
		new_rindex,
		new_zindex,
 		new_phiindex,
		stepR,
		stepPhi,
		stepZ,
		radiusRBF0,
		0,
		w
  	);


	
	//for (Int_t i=0;i<nd;i++) printf("w[%d]=%f\t",i,w[i]);
	
	val = InterpRBFCartesian
  	( 
		r,
		phi,
		z,
		startR,
		startPhi,
		startZ,
		stepR,
		stepPhi,
		stepZ,
		radiusRBF0,
		0,
		w
	);

//	printf("val = %f\n",val);
//	printf("result val= %f\n",val);

  } else {
	GetRBFWeightHalf ( 
		new_rindex,
		new_zindex,
 		new_phiindex,
		stepR,
		stepPhi,
		stepZ,
		radiusRBF0,
		0,
		w
  	);


	
//	  printf("computing InterpRBFhalf %d,%d,%d\n",new_rindex,new_zindex,new_phiindex);
	
	//for(Int_t i=0;i<nd;i++) printf("w[%d]=%f\t",i,w[i]);
	//printf("\n");

	val = InterpRBFHalf
  	( 
		r,
		phi,
		z,
		startR,
		startPhi,
		startZ,
		stepR,
		stepPhi,
		stepZ,
		radiusRBF0,
		0,
		w
	);
  }
  //printf("nd=%d, value = %f\n",nd,val);  
  delete w;

  return val;
}


/// main operation
/// Int_terpolation at the region 
Double_t AliTPC3DCylindricalInterpolatorFull::Interpolate3DTableCylRBF
( 
	Double_t r,   
	Double_t z,
	Double_t phi,
	Int_t rindex,
	Int_t zindex,
	Int_t phiindex,
	Int_t stepR,
	Int_t stepPhi,
	Int_t stepZ,
	Double_t radiusRBF0
) 
{
  

  const Float_t gridSizeR   =  (AliTPCPoissonSolver::fgkOFCRadius-AliTPCPoissonSolver::fgkIFCRadius) / (fNR - 1) ;
  const Float_t gridSizeZ   =  AliTPCPoissonSolver::fgkTPCZ0 / (fNZ -1) ;
  const Float_t gridSizePhi =  TMath::TwoPi() / fNPhi;
	
  //  Int_t stepPhi = 3;
  // Int_t stepR    = 3;
  //Int_t stepZ    = 3;

  Double_t r0,z0,phi0,d;
  Double_t MIN_DIST = 1e-3;
  Double_t val = 0.0;
  
  // guess the index of r,z,phi
  
  
  //
  

  //printf("radiusRBF0 = %f\n",radiusRBF0);
//  radiusRBF0 = 0.0; 
//  Search( fNZ,&fZlist[phiindex  * (fNZ * fNR) + (index_r * fNZ)], 1, z, jlowlist[index_phi * fNR + index_r] ) ;            
//  Search( fNPhi,&fPhilist[rindex * (fNZ) + zindex], (fNR * fNZ), phi, phiindex) ;            

  Int_t startPhi = phiindex - stepPhi/2;
  Int_t index_phi;

  Int_t startR =  rindex  - stepR/2;  
  Int_t startZ =  zindex  - stepZ/2;
  
  if (startPhi < 0) startPhi = fNPhi + startPhi;
//  if (startPhi + stepPhi >= fNR ) startPhi = startPhi fNR - stepR ;


  if (startR    < 0) startR = 0;
  if (startR + stepR >= fNR ) startR = fNR - stepR ;

  if (startZ    < fMinZIndex) startZ = fMinZIndex;
  if (startZ + stepZ  >= fNZ  ) startZ = fNZ - stepZ;

  Int_t index;
  Double_t sum_w = 0.0;
  Double_t sum_d  = 0.0;
  Double_t shortest_d = 10000.0; 
  Int_t new_rindex;
  Int_t new_zindex;
  Int_t new_phiindex;


   	

  for (Int_t iphi = startPhi;iphi < startPhi + stepPhi;iphi++) 
  {
      index_phi = iphi % fNPhi;
      for (Int_t index_r= startR ;index_r< startR + stepR;index_r++)
	{
	  for (Int_t index_z=startZ; index_z < startZ + stepZ;index_z++) 
	    {
	      // check for the closest poInt_t
	      index = index_phi * (fNZ * fNR) + index_r * fNZ + index_z;
	
	      r0 =  fRlist[index];
	      z0 =  fZlist[index];
	      phi0 =  fPhilist[index];
	

	      d = Distance(r0,phi0,z0,r,phi,z);
	      if (d < shortest_d) 
		{
			shortest_d = d;
			new_rindex = index_r;
			new_phiindex = index_phi;
			new_zindex = index_z;	
	      	}
	      
	      
	    }
	}
  }

//  printf("(%d,%d,%d)\n",new_rindex,new_phiindex,new_zindex);
// do Int_terpolant


  stepPhi = fStepPhi;
  stepR = fStepR;
  stepZ = fStepZ;
//  if (stepZ > 3) stepZ =3;


  startPhi = new_phiindex - stepPhi/2;
  
  startR =  new_rindex  - stepR/2;  
  startZ =  new_zindex  - stepZ/2;
  
  if (startPhi < 0) startPhi = fNPhi + startPhi;
 

  if (startR    < 0) startR = 0;
  if (startR + stepR >= fNR ) startR = fNR - stepR ;

  if (startZ    < fMinZIndex) startZ = fMinZIndex;
  if (startZ + stepZ  >= fNZ  ) startZ = fNZ - stepZ;


  Double_t * w;	

  //Int_t nd = (stepPhi-1) + (stepR-1) + (stepZ-1) + 1;
  Int_t nd = stepPhi  *stepR  * stepZ ;

  w = new Double_t[nd];

	Float_t minTemp,minTemp2;

/**
	if (fKernelType == kRBFThinPlateSpline)
		r0   = AliTPCPoissonSolver::fgkIFCRadius + ((startR - 1)  * gridSizeR);	
	else
		r0   = AliTPCPoissonSolver::fgkIFCRadius + ((startR + fStepR/2)  * gridSizeR);


	radiusRBF0 = Distance(r0,0.0,0.0,r0 + gridSizeR, gridSizePhi, gridSizeZ);
**/
	radiusRBF0 = GetRadius0RBF(new_rindex,new_phiindex,new_zindex);


//	radiusRBF0 = GetRadius0RBF(new_rindex,new_phiindex,new_zindex);
//	radiusRBF0 = fRadiusRBF0;
//  radiusRBF0 = Distance(0.0,0.0,0.0,gridSizeR,gridSizePhi,gridSizeZ);  
//	minTemp = r0 * gridSizePhi;
//	minTemp2 = TMath::Min(minTemp,gridSizeR);
	//adiusRBF0 = TMath::Min(minTemp2,gridSizeZ);


  if (fType == 1) {	

	for (Int_t i=0;i<nd;i++)  w[i]=0.0;


//	for (Int_t i=0;i<nd;i++) printf("w[%d]=%f\t",i,w[i]);
	
//	printf("get rbf weight\n");
	GetRBFWeight ( 
		new_rindex,
		new_zindex,
 		new_phiindex,
		stepR,
		stepPhi,
		stepZ,
		radiusRBF0,
		0,
		w
  	);


	
	//for (Int_t i=0;i<nd;i++) printf("w[%d]=%f\t",i,w[i]);
	
	val = InterpRBF
  	( 
		r,
		phi,
		z,
		startR,
		startPhi,
		startZ,
		stepR,
		stepPhi,
		stepZ,
		radiusRBF0,
		0,
		w
	);

//	printf("val = %f\n",val);
//	printf("result val= %f\n",val);

  } else {
	GetRBFWeightHalf ( 
		new_rindex,
		new_zindex,
 		new_phiindex,
		stepR,
		stepPhi,
		stepZ,
		radiusRBF0,
		0,
		w
  	);


	
//	  printf("computing InterpRBFhalf %d,%d,%d\n",new_rindex,new_zindex,new_phiindex);
	
	//for(Int_t i=0;i<nd;i++) printf("w[%d]=%f\t",i,w[i]);
	//printf("\n");

	val = InterpRBFHalf
  	( 
		r,
		phi,
		z,
		startR,
		startPhi,
		startZ,
		stepR,
		stepPhi,
		stepZ,
		radiusRBF0,
		0,
		w
	);
  }
  //printf("nd=%d, value = %f\n",nd,val);  
  delete w;

  return val;
}



/// main operation
/// Int_terpolation at the region 
Double_t AliTPC3DCylindricalInterpolatorFull::Interpolate3DTableCyl
( 
	Double_t r,   
	Double_t z,
	Double_t phi,
	Int_t rindex,
	Int_t zindex,
	Int_t phiindex,
	Int_t stepR,
	Int_t stepPhi,
	Int_t stepZ

) 
{
  Int_t ilow = 0, jlow = 0, klow = 0, m=0, myilow=0, myjlow =0 ;
  Double_t r0,z0,phi0;

  Int_t indeks;

  
  // list all low in all direction
  Int_t jlowlist[fNPhi * fNR];
  Int_t klowlist[fNR * fNZ];

  // tricubic poInt_ts
  Double_t saveArray[fOrder + 1];
  Double_t savedArray[fOrder + 1];

  while (phi < 0.0) phi = TMath::TwoPi() + phi;
  while (phi > TMath::TwoPi()) phi = phi - TMath::TwoPi() ;
  

  

  Int_t startPhi = phiindex - stepPhi/2;
  Int_t index_phi;

  Int_t startR =  rindex  - stepR/2;  
  Int_t startZ =  zindex  - stepZ/2;
  

  if (startPhi < 0) startPhi = fNPhi + startPhi;

  if (startR    < 0) startR = 0;
  if (startR + stepR >= fNR ) startR = fNR - stepR ;

  if (startZ    < 0) startZ = 0;
  if (startZ + stepZ  >= fNZ  ) startZ = fNZ - stepZ;

  
  Double_t shortest_d = 10000.0; 
  Double_t d;
  Int_t new_rindex;
  Int_t new_zindex;
  Int_t new_phiindex;



  for (Int_t iphi = startPhi;iphi < startPhi + stepPhi;iphi++) 
    {
      index_phi = iphi % fNPhi;
      for (Int_t index_r= startR ;index_r< startR + stepR;index_r++)
	{
	  for (Int_t index_z=startZ; index_z < startZ + stepZ;index_z++) 
	    {
	      // check for the closest poInt_t
	      //index = index_phi * (fNZ * fNR) + index_r * fNZ + index_z;
	
	      r0 =  fRlist[index_phi * (fNR * fNZ) + index_r * fNZ + index_z];
	      z0 =  fZlist[index_phi * (fNR * fNZ) + index_r * fNZ + index_z];
	      phi0 =  fPhilist[index_phi * (fNR * fNZ) + index_r * fNZ + index_z];
	

	      d = Distance(r0,phi0,z0,r,phi,z);
	      if (d < shortest_d) 
		{
			shortest_d = d;
			new_rindex = index_r;
			new_phiindex = index_phi;
			new_zindex = index_z;	
	      	}
	      
	      
	    }
	}
    }

  // check phi 
  //while (phi < 0.0) phi = TMath::TwoPi() + phi;
  // while (phi > TMath::TwoPi()) phi = phi - TMath::TwoPi() ;

  ilow   	= new_rindex;
  klow       = new_phiindex;
  jlow        = new_zindex;  
  
	// order >= 3
  //	if (fOrder > 2) {		
 	klow -= (fOrder/2);
   ilow -= (fOrder/2);
   jlow -= (fOrder/2);
		//	} 

	
	// check if out of range
  if ( ilow < 0 ) ilow = 0 ;   
   if ( jlow < 0 ) jlow = 0;
   
  
  
  if ( klow < 0 ) klow = fNPhi + klow ;
  
  // check if out of range
   if ( ilow + fOrder  >=    fNR - 1 ) ilow =   fNR- 1 - fOrder ;

   if ( jlow + fOrder  >=    fNZ - 1 ) jlow =   fNZ - 1 - fOrder ;



	// do for each 
  for ( Int_t k = 0 ; k <  fOrder + 1 ; k++ ) {		
		m = (klow + k) % fNPhi;
		// Int_terpolate	
		for ( Int_t i = ilow ; i < ilow + fOrder + 1 ; i++ ) {

			if (fOrder <= 2) {
				  // test
				  indeks =  m* (fNZ * fNR) +  i *(fNZ)  + jlow;
				  saveArray[i-ilow] = Interpolate( &fZlist[indeks], &fVals[indeks], z );
			 	  
			}	else {
				indeks =  m* (fNZ * fNR) +  i *(fNZ);
				saveArray[i-ilow] = InterpolateCubicSpline(fZlist,&fVals[indeks],&fSecondDerZ[indeks],fNZ,fNZ,fNZ,z,1);
			}
		}
		//  printf("savear (%f,%f,%f)\n",saveArray[0],saveArray[1],saveArray[2]);

	 	indeks =  m* (fNZ * fNR) +  ilow *(fNZ)  + jlow;
	    savedArray[k] = Interpolate (&fRlist[indeks], fNZ, saveArray,	1,	r ); 
    }


  //  printf("(%f,%f,%f)\n",savedArray[0],savedArray[1],savedArray[2]);
 	indeks =   ilow *(fNZ)  + jlow;
	return( InterpolatePhi( &fPhilist[indeks], (fNZ * fNR), klow, fNPhi,  savedArray,  phi ) )   ;



  	


}

/// main operation
Double_t AliTPC3DCylindricalInterpolatorFull::Interpolate3DTableCyl
( 
	Double_t r,   
	Double_t z,
	Double_t phi
) 
{
  Int_t ilow = 0, jlow = 0, klow = 0, m=0;
  Int_t indeks;
  
  // tricubic poInt_ts
  Double_t saveArray[fOrder + 1];
  Double_t savedArray[fOrder + 1];
  Double_t zlistM1[fOrder + 1];
  Double_t valsM1[fOrder + 1];

  Bool_t neg = kFALSE;
  // check phi 
  while (phi < 0.0) phi = TMath::TwoPi() + phi;
  while (phi > TMath::TwoPi()) phi = phi - TMath::TwoPi() ;

   	
 // search lowest index related to r,z and phi
  Search( fNR, fRlistNormalized, r, ilow   ) ;
  Search( fNZ, fZlistNormalized, z, jlow   ) ; 
  Search( fNPhi, fPhilistNormalized, phi, klow   ) ;
  
  
	// order >= 3
  //	if (fOrder > 2) {		
   klow -= (fOrder/2);
   ilow -= (fOrder/2);
   jlow -= (fOrder/2);
		//	} 

	
	// check if out of range
  if ( ilow < 0 ) ilow = 0 ;   
  if ( jlow < 0 ) jlow = 0;
   
  //ilow -=  (fOrder/2);
  //jlow -=  (fOrder/2);
  
  
  
  if ( klow < 0 ) klow = fNPhi + klow ;
  
  // check if out of range
   if ( ilow + fOrder  >=    fNR - 1 ) ilow =   fNR- 1 - fOrder ;
   if ( jlow + fOrder  >=    fNZ - 1 ) jlow =   fNZ - 1 - fOrder ;

   
  
	// do for each 
  for ( Int_t k = 0 ; k <  fOrder + 1 ; k++ ) {		
		m = (klow + k) % fNPhi;
		// Int_terpolate	
		for ( Int_t i = ilow ; i < ilow + fOrder + 1 ; i++ ) {

			if (fOrder <= 2) {

				//				printf("fval (%f,%f,%f)\n",fVals[indeks],fVals[indeks + 1],fVals[indeks + 2]);
			  if (jlow >= 0) {
				  indeks =  m* (fNZ * fNR) +  i *(fNZ)  + jlow;
				  saveArray[i-ilow] = Interpolate( &fZlistNormalized[jlow], &fValsNormalized[indeks], z );
			  }
			  else  
			    {
				  indeks =  m* (fNZ * fNR) +  i *(fNZ);
				  zlistM1[0] = fZlistNormalized[0] - (fZlistNormalized[1] - fZlistNormalized[0]);
				  zlistM1[1] = fZlistNormalized[0];
				  zlistM1[2] = fZlistNormalized[1];
				  valsM1[0] = fValsNormalized[indeks] - (fValsNormalized[indeks + 1] - fValsNormalized[indeks]);
				  valsM1[1] = fValsNormalized[indeks];
				  valsM1[2] = fValsNormalized[indeks + 1];
				  saveArray[i-ilow] = Interpolate( &zlistM1[0], &valsM1[0], z );
				}
				  
			}	else {
				indeks =  m* (fNZ * fNR) +  i *(fNZ);
				saveArray[i-ilow] = InterpolateCubicSpline(fZlistNormalized,&fValsNormalized[indeks],&fSecondDerZ[indeks],fNZ,fNZ,fNZ,z,1);
			}
		}
		//  printf("savear (%f,%f,%f)\n",saveArray[0],saveArray[1],saveArray[2]);

    savedArray[k] = Interpolate( &fRlistNormalized[ilow], saveArray, r )  ;
    //table.PrInt_t();
  }

  //  printf("(%f,%f,%f)\n",savedArray[0],savedArray[1],savedArray[2]);
	return( InterpolatePhi( &fPhilistNormalized[0], klow, fNPhi,  savedArray,  phi ) )   ;
	//	return 0.0;
}


// search
void AliTPC3DCylindricalInterpolatorFull::Search
( 
	Int_t n, 
	const Double_t xArray[], 
	Double_t x, 
	Int_t &low 
) 
{
  /// Search an ordered table by starting at the most recently used poInt_t

  Long_t middle, high ;
  Int_t  ascend = 0, increment = 1 ;

  if ( xArray[n-1] >= xArray[0] ) ascend = 1 ;  // Ascending ordered table if true

  if ( low < 0 || low > n-1 ) {
    low = -1 ; high = n ;
  } else {                                            // Ordered Search phase
    if ( (Int_t)( x >= xArray[low] ) == ascend )  {
      if ( low == n-1 ) return ;
      high = low + 1 ;
      while ( (Int_t)( x >= xArray[high] ) == ascend ) {
				low = high ;
				increment *= 2 ;
				high = low + increment ;
				if ( high > n-1 )  {  high = n ; break ;  }
      }
    } else {
      if ( low == 0 )  {  low = -1 ;  return ;  }
      high = low - 1 ;
      while ( (Int_t)( x < xArray[low] ) == ascend ) {
				high = low ;
				increment *= 2 ;
				if ( increment >= high )  {  low = -1 ;  break ;  }
				else  low = high - increment ;
      }
    }
  }

  while ( (high-low) != 1 ) {                     // Binary Search Phase
    middle = ( high + low ) / 2 ;
    if ( (Int_t)( x >= xArray[middle] ) == ascend )
      low = middle ;
    else
      high = middle ;
  }

  if ( x > xArray[n-1] ) low = n-1 ;
  if ( x <  xArray[0]   ) low = 0 ;

}


// search
void AliTPC3DCylindricalInterpolatorFull::Search
( 
	Int_t n, 
	Double_t * xArray , 
	Int_t offset,
	Double_t x, 
	Int_t &low 
) 
{
  /// Search an ordered table by starting at the most recently used poInt_t

  Long_t middle, high ;
  Int_t  ascend = 0, increment = 1 ;

  if ( xArray[(n-1) * offset] >= xArray[0 * offset] ) ascend = 1 ;  // Ascending ordered table if true

  if ( low < 0 || low > n-1 ) {
    low = -1 ; high = n ;
  } else {                                            // Ordered Search phase
    if ( (Int_t)( x >= xArray[low * offset] ) == ascend )  {
      if ( low == n-1 ) return ;
      high = low + 1 ;
      while ( (Int_t)( x >= xArray[high * offset] ) == ascend ) {
				low = high ;
				increment *= 2 ;
				high = low + increment ;
				if ( high > n-1 )  {  high = n ; break ;  }
      }
    } else {
      if ( low == 0 )  {  low = -1 ;  return ;  }
      high = low - 1 ;
      while ( (Int_t)( x < xArray[low * offset] ) == ascend ) {
				high = low ;
				increment *= 2 ;
				if ( increment >= high )  {  low = -1 ;  break ;  }
				else  low = high - increment ;
      }
    }
  }

  while ( (high-low) != 1 ) {                     // Binary Search Phase
    middle = ( high + low ) / 2 ;
    if ( (Int_t)( x >= xArray[middle * offset] ) == ascend )
      low = middle ;
    else
      high = middle ;
  }

  if ( x >  xArray[n-1] ) low = n ;
  if ( x < xArray[0]   ) low = -1 ;
 

}



// Interpolate 1D
Double_t AliTPC3DCylindricalInterpolatorFull::Interpolate
( 
	Double_t xArray[], 
	Double_t yArray[],
	Double_t x 
) 
{
  // Interpolate function Y(x) using linear (order=1) or quadratic (order=2) Int_terpolation.
 
  Double_t y ;
  //printf("order: %d x: (%f,%f,%f)\n",fOrder,xArray[0],xArray[1],x);
  
  //printf("y: (%f,%f)\n",yArray[0],yArray[1]);
  if (fOrder > 2) {
		Double_t y2Array[fOrder + 1];	
		//printf("x: (%f,%f,%f,%f,%f)\n",xArray[0],xArray[1],xArray[2],xArray[3],x);
		InitCubicSpline(xArray,yArray,fOrder + 1,y2Array,1);
		y = InterpolateCubicSpline(xArray,yArray,y2Array,fOrder + 1,fOrder + 1,fOrder + 1,x,1);
		//printf("y: (%f,%f,%f,%f,%f)\n",yArray[0],yArray[1],yArray[2],yArray[3],y);
		
	} else if ( fOrder == 2 ) {                
		// Quadratic Interpolation = 2
    y  = (x-xArray[1]) * (x-xArray[2]) * yArray[0] / ( (xArray[0]-xArray[1]) * (xArray[0]-xArray[2]) ) ;
    y += (x-xArray[2]) * (x-xArray[0]) * yArray[1] / ( (xArray[1]-xArray[2]) * (xArray[1]-xArray[0]) ) ;
    y += (x-xArray[0]) * (x-xArray[1]) * yArray[2] / ( (xArray[2]-xArray[0]) * (xArray[2]-xArray[1]) ) ;
  } else {                           
		// Linear nterpolation = 1
    y  = yArray[0] + ( yArray[1]-yArray[0] ) * ( x-xArray[0] ) / ( xArray[1] - xArray[0] ) ;
  }

  return (y);

}


// Interpolate 1D
Double_t AliTPC3DCylindricalInterpolatorFull::Interpolate
( 
	Double_t xArray[], 
	Int_t offsetX,
	Double_t yArray[],
	Int_t offsetY,
	Double_t x 
) 
{
  // Interpolate function Y(x) using linear (order=1) or quadratic (order=2) Int_terpolation.
 
  Double_t y ;
  //printf("order: %d x: (%f,%f,%f)\n",fOrder,xArray[0],xArray[1],x);
  
  //printf("y: (%f,%f)\n",yArray[0],yArray[1]);
  if (fOrder > 2) {
		Double_t y2Array[fOrder + 1];	
		//printf("x: (%f,%f,%f,%f,%f)\n",xArray[0],xArray[1],xArray[2],xArray[3],x);
		InitCubicSpline(xArray,yArray,fOrder + 1,y2Array,1);
		y = InterpolateCubicSpline(xArray,yArray,y2Array,fOrder + 1,fOrder + 1,fOrder + 1,x,1);
		//printf("y: (%f,%f,%f,%f,%f)\n",yArray[0],yArray[1],yArray[2],yArray[3],y);
		
	} else if ( fOrder == 2 ) {                
		// Quadratic Interpolation = 2
    y  = (x-xArray[1 * offsetX]) * (x-xArray[2 * offsetX]) * yArray[0 * offsetY] / ( (xArray[0 * offsetX]-xArray[1 * offsetX]) * (xArray[0 * offsetX]-xArray[2 * offsetX]) ) ;
    y += (x-xArray[2 * offsetX]) * (x-xArray[0 * offsetX]) * yArray[1 * offsetY] / ( (xArray[1 * offsetX]-xArray[2 * offsetX]) * (xArray[1 * offsetX]-xArray[0 * offsetX]) ) ;
    y += (x-xArray[0 * offsetX]) * (x-xArray[1 * offsetX]) * yArray[2 * offsetY] / ( (xArray[2 * offsetX]-xArray[0 * offsetX]) * (xArray[2 * offsetX]-xArray[1 * offsetX]) ) ;
  } else {                           
		// Linear Interpolation = 1
    y  = yArray[0 * offsetY] + ( yArray[1 * offsetY]-yArray[0 * offsetY] ) * ( x-xArray[0 * offsetX] ) / ( xArray[1 * offsetX] - xArray[0 * offsetX] ) ;
  }

  return (y);

}


// Int_terpolate phi
Double_t AliTPC3DCylindricalInterpolatorFull::InterpolatePhi
( 
	Double_t xArray[], 
	const Int_t ilow,
	const Int_t nx,
	Double_t yArray[],
	Double_t x 
)
{
  /// Interpolate function Y(x) using linear (order=1) or quadratic (order=2) Int_terpolation.

	Int_t i0 = ilow;
	Double_t xi0 = xArray[ilow];	
	Int_t i1 = (ilow + 1) % nx;
	Double_t xi1 = xArray[i1];
	Int_t i2 = (ilow + 2) % nx;
	Double_t xi2 = xArray[i2];
	if ((ilow + 1) >= nx) {
		xi1 += TMath::TwoPi();
	}
	if ((ilow + 2) >= nx) {
		xi2 += TMath::TwoPi();
		
	}
  Double_t y ;
  if (fOrder > 2) {
		Double_t y2Array[fOrder + 1];	
		Double_t xArrayTemp[fOrder + 1];	
		Double_t dPhi = xArray[1] - xArray[0];
		for (Int_t i=0;i< fOrder + 1;i++) {
			xArrayTemp[i] = xArray[ilow] + (dPhi * i);
			
		
		}
		if ((xArrayTemp[0] - x) > 1e-10) x = x + TMath::TwoPi();
		
		InitCubicSpline(xArrayTemp,yArray,fOrder + 1,y2Array,1);
		y = InterpolateCubicSpline(xArrayTemp,yArray,y2Array,fOrder + 1,fOrder + 1,fOrder + 1,x,1);
	}
	else if ( fOrder == 2 ) {                // Quadratic Interpolation = 2		
    y  = (x-xi1) * (x-xi2) * yArray[0] / ( (xi0-xi1) * (xi0-xi2) ) ;
    y += (x-xi2) * (x-xi0) * yArray[1] / ( (xi1-xi2) * (xi1-xi0) ) ;
    y += (x-xi0) * (x-xi1) * yArray[2] / ( (xi2-xi0) * (xi2-xi1) ) ;
  } else {                           // Linear Interpolation = 1
    y  = yArray[0] + ( yArray[1]-yArray[0] ) * ( x-xArray[i0] ) / (xi1 - xArray[i0] ) ;
  }

  return (y);

}
// Int_terpolate phi
Double_t AliTPC3DCylindricalInterpolatorFull::InterpolatePhi
( 
	Double_t xArray[], 
	Int_t offsetX,
	const Int_t ilow,
	const Int_t nx,
	Double_t yArray[],
	Double_t x 
)
{
  /// Interpolate function Y(x) using linear (order=1) or quadratic (order=2) Int_terpolation.

	Int_t i0 = ilow;
	Double_t xi0 = xArray[ilow * offsetX];	
	Int_t i1 = (ilow + 1) % nx;

	Double_t xi1 = xArray[i1 * offsetX];
	Int_t i2 = (ilow + 2) % nx;
	Double_t xi2 = xArray[i2 *offsetX];
	Double_t y ;
	if (xi1 < xi0) xi1 = TMath::TwoPi() + xi1;
	if (xi2 < xi1) xi2 = TMath::TwoPi() + xi2;
	if (x   < xi0) x   = TMath::TwoPi() + x;


	if (fOrder > 2) {
		Double_t y2Array[fOrder + 1];	
		Double_t xArrayTemp[fOrder + 1];

	
		Double_t dPhi = xArray[1 * offsetX] - xArray[0 * offsetX];

		for (Int_t i=0;i< fOrder + 1;i++) {
			xArrayTemp[i] = xArray[ilow * offsetX] + (dPhi * i);
			
		
		}
		if ((xArrayTemp[0] - x) > 1e-10) x = x + TMath::TwoPi();
		
		InitCubicSpline(xArrayTemp,yArray,fOrder + 1,y2Array,1);

		y = InterpolateCubicSpline(xArrayTemp,yArray,y2Array,fOrder + 1,fOrder + 1,fOrder + 1,x,1);
	}
	else if ( fOrder == 2 ) {                // Quadratic Interpolation = 2		






	  y  = (x-xi1) * (x-xi2) * yArray[0] / ( (xi0-xi1) * (xi0-xi2) ) ;
	  y += (x-xi2) * (x-xi0) * yArray[1] / ( (xi1-xi2) * (xi1-xi0) ) ;
	  y += (x-xi0) * (x-xi1) * yArray[2] / ( (xi2-xi0) * (xi2-xi1) ) ;
	} else {                           // Linear Interpolation = 1
	  
   // make non negative
   
   

    y  = yArray[0] + ( yArray[1]-yArray[0] ) * ( x-xi0) / (xi1 - xi0 ) ;
  }

  return (y);

}

// get value
Double_t AliTPC3DCylindricalInterpolatorFull::GetValue
(
	Double_t r, 
	Double_t phi, 
	Double_t z,
	Int_t rindex,
	Int_t phiindex,
	Int_t zindex,
	Int_t stepR,
	Int_t stepPhi,
	Int_t stepZ
) 
{  

  // return 0.0;
//  return Interpolate3DTableCylIDW(r,z,phi,rindex,zindex,phiindex,stepR,stepPhi,stepZ);
	fMinZIndex = 0;
 	return Interpolate3DTableCylRBF(r,z,phi,rindex,zindex,phiindex,stepR,stepPhi,stepZ,0.0);
//  return Interpolate3DTableCylRBF(r,z,phi,rindex,zindex,phiindex,stepR,stepPhi,stepZ,0.0);
	
}



// get value
Double_t AliTPC3DCylindricalInterpolatorFull::GetValueCartesian
(
	Double_t r, 
	Double_t phi, 
	Double_t z,
	Int_t rindex,
	Int_t phiindex,
	Int_t zindex,
	Int_t stepR,
	Int_t stepPhi,
	Int_t stepZ
) 
{  

  // return 0.0;
//  return Interpolate3DTableCylIDW(r,z,phi,rindex,zindex,phiindex,stepR,stepPhi,stepZ);
	fMinZIndex = 0;
 	return Interpolate3DTableCylRBFCartesian(r,z,phi,rindex,zindex,phiindex,stepR,stepPhi,stepZ,0.0);
//  return Interpolate3DTableCylRBF(r,z,phi,rindex,zindex,phiindex,stepR,stepPhi,stepZ,0.0);
	
}

// get value
Double_t AliTPC3DCylindricalInterpolatorFull::GetValue
(
	Double_t r, 
	Double_t phi, 
	Double_t z,
	Int_t rindex,
	Int_t phiindex,
	Int_t zindex,
	Int_t stepR,
	Int_t stepPhi,
	Int_t stepZ,
	Int_t minzindex
) 
{  
  fMinZIndex = minzindex;
  
//  printf("GetValue (%f,%f,%f)\n",r,phi,z);
  return Interpolate3DTableCylRBF(r,z,phi,rindex,zindex,phiindex,stepR,stepPhi,stepZ,0.0);
	
}


// get value
Double_t AliTPC3DCylindricalInterpolatorFull::GetValueCartesian
(
	Double_t r, 
	Double_t phi, 
	Double_t z,
	Int_t rindex,
	Int_t phiindex,
	Int_t zindex,
	Int_t stepR,
	Int_t stepPhi,
	Int_t stepZ,
	Int_t minzindex
) 
{  
  fMinZIndex = minzindex;
  
//  printf("GetValue (%f,%f,%f)\n",r,phi,z);
  return Interpolate3DTableCylRBFCartesian(r,z,phi,rindex,zindex,phiindex,stepR,stepPhi,stepZ,0.0);
	
}


// get value
Double_t AliTPC3DCylindricalInterpolatorFull::GetValue
(
	Double_t r, 
	Double_t phi, 
	Double_t z
) 
{  

  // return 0.0;
  return Interpolate3DTableCyl(r,z,phi);
	
}



// spline from recepi
void AliTPC3DCylindricalInterpolatorFull::InitCubicSpline
(	
	Double_t * xArray,
	Double_t * yArray,
	const Int_t n,
	Double_t * y2Array,
	const Int_t skip
) 
{
	Double_t u[n];
	Double_t sig,p,qn,un;
	
	y2Array[0] = 0.0;
	u[0] = 0.0; //natural condition
		
		
	for (Int_t i=1;i<=n-2;i++) {
		sig = (xArray[i] - xArray[i-1])/(xArray[i+1] - xArray[i-1]);
		p = sig * y2Array[(i-1) * skip]+ 2.0;
		y2Array[i * skip] = (sig - 1.0)/p;
		u[i] = (yArray[(i+1) * skip] - yArray[i * skip])/(xArray[i+1] - xArray[i]) - (yArray[i * skip] - yArray[(i-1) * skip])/(xArray[i] - xArray[i-1]);
		u[i] = (6.0 * u[i]/(xArray[i+1] - xArray[i-1]) - sig * u[i-1])/p;	
	}
	
	
	
	
	qn = un = 0.0;
	
	y2Array[(n-1) * skip] = (un - qn * u[n-2])/(qn*y2Array[(n-2) * skip]+1.0);
	for (Int_t k=n-2;k>=1;k--)
		y2Array[k * skip] = y2Array[k * skip] * y2Array[(k+1) * skip] + u[k];	
}



// spline from recepi
//void AliTPC3DCylindricalInterpolator::InitCubicSplinePhi
//(	
//	Double_t * xArray,
//	Double_t * yArray,
//	const Int_t n,
//	Double_t * y2Array,
//	const Int_t skip
//) 
//{
//	/// find z
//	Double_t z[n];
//	Double_t u[n];
//	
//	Double_t alpha = (TMath::Sqrt(12.0) - 4.0)/2.0;
//	Double_t delta =  (1 + alpha*alpha);
//	delta = delta / (4.0 * (1 - (alpha * alpha)) *(1 - TMath::Power(alpha,n)));
	
	
	
//	for (Int_t i=0;i<n;i++) {
//		z[i] = delta * (TMath::Power(alpha,i) + TMath::Power(alpha,n-i));
//		
//		if (i == 0) {
//			u[i] = 3 * (yArray[(i+1) * skip] - yArray[(n-1) * skip]);	
//			continue;
//		}
//		if (i == n-1) {
//		u[i] = 3 * (yArray[0] - yArray[(n-1) * skip]);	
//			continue;
//		}
//		u[i] = 3 * (yArray[(i+1) * skip] - yArray[(i-1) * skip]);				
//	}
//	
//	y2Array[0] = z[0] * u[0];
//	Int_t m = (n+1)/2;
	
//	for (Int_t i=1;i<m;i++) 
//		y2Array[0] += z[i] * (u[i] + u[n-i]);
	
//	if ((n%2) == 0)
//		y2Array[0] += z[m] * u[m];
		
	
//	u[1] = u[1] - y2Array[0];
//	u[n-1] = u[n-1] - y2Array[0];
	
	// solve tridiagonalfor solving y2Array[1..n]
	
	// forward elimination
//	for (Int_t i = 2;i<n;i++)
//	{
//		u[i] = u[i] - 0.25 * u[i-1];
//	}
	
	// backward subs
//	y2Array[(n-1) * skip] = u[n-1] / (3.75);
//	for (Int_t i = n-2;i >= 2;i--)
//	{
//		y2Array[i * skip] = (u[i] - y2Array[(i+1) * skip])/(3.75);
//	}
//	y2Array[1 * skip] = (u[1] - y2Array[2 * skip])/(4.0);
//}


// spline from recepi
Double_t AliTPC3DCylindricalInterpolatorFull::InterpolateCubicSpline
(	
	Double_t * xArray, 
	Double_t * yArray,
	Double_t * y2Array,	
	const Int_t nxArray,
	const Int_t nyArray,
	const Int_t ny2Array,
	Double_t x,
	Int_t skip
) 
{
	Int_t klo,khi,k;
	Float_t h,b,a;
	
	
	klo = 0;
	khi = nxArray-1;
	while (khi-klo > 1) {
		k = (khi + klo) >> 1;
		if (xArray[k] > x) khi=k;
		else klo = k;
	}
	
	h = xArray[khi] - xArray[klo];
	
	
	if (h < 1e-20)  {
		printf("not found\n");
		return 0.0;
	}
	
	a = (xArray[khi] - x) / h;
	b = (x - xArray[klo]) / h;
	
	Double_t y = a *yArray[klo ] + b * yArray[khi ] + ((a*a*a -a) * y2Array[klo * skip] + (b*b*b -b) * y2Array[khi * skip]) * (h*h)/6.0;
	
	return y;
	
}



// spline from recepi
Double_t AliTPC3DCylindricalInterpolatorFull::InterpolateCubicSplinePhi
(	
	Double_t * xArray, 
	Double_t * yArray,
	Double_t * y2Array,	
	const Int_t n,
	Double_t x,
	const Int_t skip
) 
{
	Int_t klo,khi,k;
	Float_t h,b,a,t;
	
	
	klo = 0;
	khi = n-1;
	while (khi-klo > 1) {
		k = (khi + klo) >> 1;
		if (xArray[k] > x) khi=k;
		else klo = k;
	}
	
	h = xArray[1] - xArray[0]; // delta phi	
	t = (x - xArray[klo]) / h;
	
	
	Double_t ai = yArray[0];
	Double_t bi = y2Array[klo * skip];
	Double_t ci = 3 * (yArray[1] - yArray[0]) - (2 * bi) - y2Array[((klo + 1)%n) * skip];
	Double_t di = 2 * (yArray[0] - yArray[1]) + bi + y2Array[((klo + 1)%n) * skip];
	Double_t y = ai + bi * t + ci * t * t + di * t * t *t;
	return y;
	
}

void AliTPC3DCylindricalInterpolatorFull::InitCubicSpline() 
{
	if (fIsInitCubic != kTRUE) {
		
		fSecondDerZ = new Double_t[fNR * fNZ * fNPhi]; 
		
		TStopwatch w;
		
		w.Start();
		// Init at Z direction
		for (Int_t m=0;m < fNPhi;m++)
			for (Int_t i=0; i< fNR;i++)		
				InitCubicSpline
				(		
					fZlist, 
					&fVals[m * (fNZ * fNR) + i * fNZ],
					fNZ,
					&fSecondDerZ[m * (fNZ * fNR) + i * fNZ],
					1
				);

				
		// Init at R direction
		//for (Int_t m=0;m < fNPhi;m++)
		//	for (Int_t j=0; j< fNZ;j++)		
		//		InitCubicSpline
		//		(		
		//			fRlist, 
		//			&fVals[m * (fNZ * fNR) + j ],
		//			fNR,
		//			&fSecondDerR[m * (fNZ * fNR) + j ],
		//			fNZ
		//		);
		
		
		
		// Init at  Phi direction
		//for (Int_t i=0;i < fNR;i++)
		//	for (Int_t j=0; j< fNZ;j++)		
		//		InitCubicSplinePhi
		//		(		
		//			fPhilist, 
		//			&fVals[i * fNZ + j ],
		//			fNPhi,
		//			&fSecondDerPhi[i * fNZ + j ],
		//			fNZ * fNR
		//		);
				
				
		w.Stop();
		
		//for (Int_t m=0;m < fNPhi;m++)		
		//	for (Int_t i=0;i < fNR;i++)
		//		for (Int_t j=0; j< fNZ;j++)		
		//			printf("%f,%f,%f\n",fSecondDerR[m * fNR * fNZ + i * fNZ + j],fSecondDerZ[m * fNR * fNZ + i * fNZ + j],fSecondDerPhi[m * fNR * fNZ + i * fNZ + j]);
		
		
		printf("%f s\n",w.CpuTime());
		
		fIsInitCubic = kTRUE;
		
	}
}


void AliTPC3DCylindricalInterpolatorFull::SetVals
(
	TMatrixD ** mvals, 
	TMatrixD ** mrlist, 
	TMatrixD ** mphilist, 
	TMatrixD ** mzlist
) 
{
	Int_t indeksm;
	Int_t indeks;
	TMatrixD * mat;
	TMatrixD * matr;
	TMatrixD * matphi;
	TMatrixD * matz;

	if (!fIsAllocatingLookUp) {
	  fVals = new Double_t[fNPhi * fNR * fNZ];
	  fRlist = new Double_t[fNPhi * fNR * fNZ];
	  fPhilist = new Double_t[fNPhi * fNR * fNZ];
	  fZlist  = new Double_t[fNPhi * fNR * fNZ];
	  fIsAllocatingLookUp = kTRUE;		
	} 
	  
		for (Int_t m=0;m < fNPhi;m++)	{
			indeksm = m * fNR * fNZ;
			mat = mvals[m];
			matr = mrlist[m];
			matphi = mphilist[m];
			matz = mzlist[m];

			for (Int_t i=0;i < fNR;i++)	{
				indeks =  indeksm  + i * fNZ;
				for (Int_t j=0; j< fNZ;j++)	{
					fVals[indeks + j] = (*mat)(i,j);

					fRlist[indeks + j] = (*matr)(i,j);
					fPhilist[indeks + j] = (*matphi)(i,j);
					fZlist[indeks + j] = (*matz)(i,j);

	//				printf("(%d,%d,%d)= (%f,%f,%f)\n",i,j,m,(*matr)(i,j),(*matz)(i,j),(*matphi)(i,j));
					// init RBF weights here
				}	
			}	
		}	

		InitRBFWeight();

}

// Kartesian
void AliTPC3DCylindricalInterpolatorFull::SetValsCartesian
(
	TMatrixD ** mvals, 
	TMatrixD ** mrlist, 
	TMatrixD ** mphilist, 
	TMatrixD ** mzlist
) 
{
	Int_t indeksm;
	Int_t indeks;
	TMatrixD * mat;
	TMatrixD * matr;
	TMatrixD * matphi;
	TMatrixD * matz;

	if (!fIsAllocatingLookUp) {
	  fVals = new Double_t[fNPhi * fNR * fNZ];
	  fRlist = new Double_t[fNPhi * fNR * fNZ];
	  fPhilist = new Double_t[fNPhi * fNR * fNZ];
	  fZlist  = new Double_t[fNPhi * fNR * fNZ];
	  fIsAllocatingLookUp = kTRUE;		
	} 
	  
		for (Int_t m=0;m < fNPhi;m++)	{
			indeksm = m * fNR * fNZ;
			mat = mvals[m];
			matr = mrlist[m];
			matphi = mphilist[m];
			matz = mzlist[m];

			for (Int_t i=0;i < fNR;i++)	{
				indeks =  indeksm  + i * fNZ;
				for (Int_t j=0; j< fNZ;j++)	{
					fVals[indeks + j] = (*mat)(i,j);

					fRlist[indeks + j] = (*matr)(i,j) * TMath::Cos((*matphi)(i,j));
					fPhilist[indeks + j] = (*matr)(i,j) * TMath::Sin((*matphi)(i,j));
					fZlist[indeks + j] = (*matz)(i,j);

					printf("fRlist[%d] = %f\n",indeks + j,(*mat)(i,j));
					// init RBF weights here
				}	
			}	
		}	

		InitRBFWeightCartesian();

}

// init RBF Weights assume, vals already been set
void AliTPC3DCylindricalInterpolatorFull::InitRBFWeight()
{

	Int_t indeksm;
	Int_t indeksr;
	Int_t indeks;
	Int_t startR;
	Int_t nd;


	const Double_t gridSizeR   =  (AliTPCPoissonSolver::fgkOFCRadius-AliTPCPoissonSolver::fgkIFCRadius) / (fNR - 1) ;
	const Double_t gridSizeZ   =  AliTPCPoissonSolver::fgkTPCZ0 / (fNZ -1) ;
    const Double_t gridSizePhi =  TMath::TwoPi() / fNPhi;

	Float_t r0;
	Double_t radiusRBF0,minTemp,minTemp2;

	nd = fStepR * fStepPhi * fStepZ;
	// init RBF weights
	for (Int_t m=0;m < fNPhi;m++)	{
		printf("RBF weights init %d\n",m);
		indeksm = m * fNR * fNZ;
		for (Int_t i=0;i < fNR;i++)	{
			indeksr =  indeksm  + i * fNZ;

			startR =  i  - fStepR/2;  

			if (startR    < 0) startR = 0;
			if (startR + fStepR >= fNR ) startR = fNR - fStepR ;


			//r0   = AliTPCPoissonSolver::fgkIFCRadius +  i * gridSizeR;	
			//maxTemp = TMath::Min(r0 * gridSizePhi,gridSizeR);
			//radiusRBF0 = TMath::Min(maxTemp,gridSizeZ); 
//			minTemp = Distance(r0,0.0,0.0,r0,gridSizePhi,0.0);
/**			minTemp = r0 * gridSizePhi;
			minTemp2 = TMath::Min(minTemp,gridSizeR);
			radiusRBF0 = TMath::Min(minTemp2,gridSizeZ);
			fRadiusRBF0 = radiusRBF0;
**//**
			radiusRBF0 = Distance(r0,0.0,0.0,r0 + gridSizeR, 0.0, 0.0) ;
			radiusRBF0 += Distance(r0,0.0,0.0,r0, gridSizePhi, 0.0) ;
			radiusRBF0 += Distance(r0,0.0,0.0,r0, 0.0, gridSizeZ) ;

			radiusRBF0 += Distance(r0,0.0,0.0,r0 + gridSizeR, 0.0, gridSizeZ);
			radiusRBF0 += Distance(r0,0.0,0.0,r0 + gridSizeR, gridSizePhi, 0.0);
			radiusRBF0 += Distance(r0,0.0,0.0,r0 , gridSizePhi, gridSizeZ);


**//**
			if (fKernelType == kRBFThinPlateSpline)
				r0   = AliTPCPoissonSolver::fgkIFCRadius + ((startR-1)   * gridSizeR);
			else
				r0   = AliTPCPoissonSolver::fgkIFCRadius + ((startR + fStepR/2)  * gridSizeR);


			radiusRBF0 = Distance(r0,0.0,0.0,r0 + gridSizeR , gridSizePhi,  gridSizeZ);
**/
	//		radiusRBF0 /= 7.0;


			for (Int_t j=0; j< fNZ;j++)	{
				indeks = indeksr + j;


				radiusRBF0 = GetRadius0RBF(i,j,m);

  				RBFWeight ( 
					i,
					j,
			 		m,
					fStepR,
					fStepPhi,
					fStepZ,
					radiusRBF0,
					fKernelType,
					&fRBFWeight[indeks * nd]
				);

				fRBFWeightLookUp[indeks] = 1;

			}	

		}
	}	

}


// init RBF Weights assume, vals already been set
void AliTPC3DCylindricalInterpolatorFull::InitRBFWeightCartesian()
{

	Int_t indeksm;
	Int_t indeksr;
	Int_t indeks;
	Int_t startR;
	Int_t nd;


	const Double_t gridSizeR   =  (AliTPCPoissonSolver::fgkOFCRadius-AliTPCPoissonSolver::fgkIFCRadius) / (fNR - 1) ;
	const Double_t gridSizeZ   =  AliTPCPoissonSolver::fgkTPCZ0 / (fNZ -1) ;
    const Double_t gridSizePhi =  TMath::TwoPi() / fNPhi;

	Float_t r0;
	Double_t radiusRBF0, rdistance, rphidistance, zdistance;

	nd = fStepR * fStepPhi * fStepZ;
	// init RBF weights
	for (Int_t m=0;m < fNPhi;m++)	{
		printf("RBF weights init %d\n",m);
		indeksm = m * fNR * fNZ;
		for (Int_t i=0;i < fNR;i++)	{
			indeksr =  indeksm  + i * fNZ;

			startR =  i  - fStepR/2;  

			if (startR    < 0) startR = 0;
			if (startR + fStepR >= fNR ) startR = fNR - fStepR ;


			r0   = AliTPCPoissonSolver::fgkIFCRadius + (startR * gridSizeR);	
			//maxTemp = TMath::Min(r0 * gridSizePhi,gridSizeR);
			//radiusRBF0 = TMath::Min(maxTemp,gridSizeZ); 

//			rdistance = Distance(r0 * TMath::Cos(0.0),r0 * TMath::Sin(0.0),0.0,(r0 + gridSizeR) * TMath::Cos(0),(r0 + gridSizeR) * TMath::Sin(0),0);  
//			rphidistance = Distance(r0 * TMath::Cos(0.0),r0 * TMath::Sin(0.0),0.0,r0* TMath::Cos(gridSizePhi),r0  * TMath::Sin(gridSizePhi),0);  

//			zdistance = Distance(r0 * TMath::Cos(0.0),r0 * TMath::Sin(0.0),0.0,r0* TMath::Cos(0.0),r0  * TMath::Sin(0.0),gridSizeZ);  
	
//			radiusRBF0 = TMath::Min(rdistance,rphidistance);
//			rdistance  = TMath::Min(radiusRBF0,zdistance);
//			radiusRBF0  = rdistance;
//			fRadiusRBF0 = radiusRBF0;
			
		

//			radiusRBF0 = Distance(r0 * TMath::Cos(0.0),r0 * TMath::Sin(0.0),0.0,(r0+gridSizeR) * TMath::Cos(gridSizePhi),(r0 + gridSizeR) * TMath::Sin(gridSizePhi),gridSizeZ);

			//radiusRBF0 = GetRadius0RBF(i,j,m);
	//		printf("Radius0 = %f\n",radiusRBF0);

			fRadiusRBF0 = radiusRBF0;

			for (Int_t j=0; j< fNZ;j++)	{
				indeks = indeksr + j;

				radiusRBF0 = GetRadius0RBF(i,j,m);


  				RBFWeightCartesian ( 
					i,
					j,
			 		m,
					fStepR,
					fStepPhi,
					fStepZ,
					radiusRBF0,
					fKernelType,
					&fRBFWeight[indeks * nd]
				);

				fRBFWeightLookUp[indeks] = 1;

			}	

		}
	}	

}

// set vals 
// 
void AliTPC3DCylindricalInterpolatorFull::SetVals
(
	TMatrixD ** mvals, 
	TMatrixD ** mrlist, 
	TMatrixD ** mphilist, 
	TMatrixD ** mzlist, 
	Int_t jy
) 
{
	Int_t indeksm;
	Int_t indeks;
	TMatrixD * mat;
	TMatrixD * matr;
	TMatrixD * matphi;
	TMatrixD * matz;

	if (!fIsAllocatingLookUp) {
	  fVals = new Double_t[fNPhi * fNR * fNZ];
	  fValsNormalized = new Double_t[fNPhi * fNR * fNZ];
	  fRlist = new Double_t[fNPhi * fNR * fNZ];
	  fPhilist = new Double_t[fNPhi * fNR * fNZ];
	  fZlist  = new Double_t[fNPhi * fNR * fNZ];


	  fIsAllocatingLookUp = kTRUE;		
	} 
	  
	for (Int_t m=0;m < fNPhi;m++)	{
		indeksm = m * fNR * fNZ;
		mat = mvals[m];
		matr = mrlist[m];
		matphi = mphilist[m];
		matz = mzlist[m];

		for (Int_t i=0;i < fNR;i++)	{
			indeks =  indeksm  + i * fNZ;
//			for (Int_t j=0; j< fNZ;j++)	{
				fVals[indeks + jy] = (*mat)(i,jy);

				fRlist[indeks + jy] = (*matr)(i,jy);
				fPhilist[indeks + jy] = (*matphi)(i,jy);
				fZlist[indeks + jy] = (*matz)(i,jy);
//				printf("fRlist[%d] = %f\n",indeks + jy,(*mat)(i,jy));
//			}	
		}
	}	



	       
}


// calculate
// RBFWeight for all points in the interpolation
void AliTPC3DCylindricalInterpolatorFull::RBFWeight 
( 
	Int_t rindex,
	Int_t zindex,
	Int_t phiindex,
	Int_t stepR,
	Int_t stepPhi,
	Int_t stepZ,
	Double_t radius0,
	Int_t kernelType,
	Double_t * w
)

{


   	
  Double_t *a;
  Int_t i;
  Int_t j;
  Int_t k;
  Int_t ii;
  Int_t jj;
  Int_t kk;

  Int_t index0,index1;
  Int_t indexCyl0,indexCyl1;
  Double_t *r;
  Double_t *v;
  
  Double_t phi0;
  Double_t z0;
  Double_t r0;


  Double_t phi1;
  Double_t z1;
  Double_t r1;

  Int_t nd = stepR * stepPhi * stepZ;

  a = new Double_t[nd*nd];
  r = new Double_t[nd];
  v = new Double_t[nd];


  Int_t startPhi = phiindex - stepPhi/2;
  Int_t index_phi;
  Int_t index_phi1;

  Int_t startR =  rindex  - stepR/2;  
  Int_t startZ =  zindex  - stepZ/2;
  
  if (startPhi < 0) startPhi = fNPhi + startPhi;

  if (startR    < 0) startR = 0;
  if (startR + stepR >= fNR ) startR = fNR - stepR ;

  if (startZ    < fMinZIndex) startZ = fMinZIndex;
  if (startZ + stepZ  >= fNZ  ) startZ = fNZ - stepZ;


/// calculate distance between each poInt_t

  index0 = 0;

	for ( i=startPhi; i < startPhi + stepPhi;i++) 
	{
	    index_phi = i % fNPhi;	

		for (j=startR; j < startR + stepR;j++) 
		{
			for(k=startZ; k < startZ + stepZ;k++) 
			{
				indexCyl0 = index_phi * fNR * fNZ + j * fNZ + k;

			

			    r0 =  fRlist[indexCyl0];
		      	z0 =  fZlist[indexCyl0];
		      	phi0 =  fPhilist[indexCyl0];

				index1 = 0;	
			//	printf("init rbf for (%f,%f,%f)\n",r0,z0,phi0);
				for ( ii=startPhi; ii < startPhi + stepPhi;ii++) {
					index_phi1 = ii % fNPhi;
					for (jj=startR; jj < startR + stepR;jj++) {
						for(kk=startZ; kk < startZ + stepZ;kk++) {
							indexCyl1 = index_phi1 * fNR * fNZ + jj * fNZ + kk;
			    			r1 =  fRlist[indexCyl1];
		      				z1 =  fZlist[indexCyl1];
		      				phi1 =  fPhilist[indexCyl1];

						
				//			printf("(%f,%f)\n",phi0,phi1);
							r[index1] = Distance(r0,phi0,z0,r1,phi1,z1);
	
						//	if (zindex < 14)
						//		printf("[%d,%d](%f,%f,%f)(%f,%f,%f) = r[%d] = %f\n",indexCyl0,indexCyl1,r0,phi0,z0,r1,phi1,z1,index1,r[index1]);
							index1++;
						}
					}
				}			
		
/**
		  	TVectorD vec_r;

			vec_r.Use(nd,r);

			printf("r: \n");
			vec_r.Print();



		  	TVectorD vec_v;

			vec_v.Use(nd,v);

			printf("v: \n");
			vec_v.Print();
**/
	
			    Phi ( nd, r, radius0, v );		

				index1 = 0;	
				for ( ii=startPhi; ii < startPhi + stepPhi;ii++) 
				{
					index_phi1 = ii % fNPhi;
					for (jj=startR; jj < startR + stepR;jj++) 
					{
						for(kk=startZ; kk < startZ + stepZ;kk++) 
						{
				      		a[index0 * nd+ index1] = v[index1];
							index1++;
						}
					}
				}

				w[index0] = fVals[indexCyl0];
				index0++;
				

			}
		}

		
	 }
	
//  Solve for the weights.

  	TMatrixD mat_a;
  	mat_a.Use(nd,nd,a);
	TVectorD vec_w;

//	mat_a.Print();


  	vec_w.Use(nd,w);
  	TDecompSVD svd(mat_a);

   	svd.Solve(vec_w);


 	

	delete [] a;
 	delete [] r;
  	delete [] v;
	

	
	//weight = w;
 // return  w;
}




// calculate
// RBFWeight for all points in the interpolation
void AliTPC3DCylindricalInterpolatorFull::RBFWeightCartesian
( 
	Int_t rindex,
	Int_t zindex,
	Int_t phiindex,
	Int_t stepR,
	Int_t stepPhi,
	Int_t stepZ,
	Double_t radius0,
	Int_t kernelType,
	Double_t * w
)

{


   	
  Double_t *a;
  Int_t i;
  Int_t j;
  Int_t k;
  Int_t ii;
  Int_t jj;
  Int_t kk;

  Int_t index0,index1;
  Int_t indexCyl0,indexCyl1;
  Double_t *r;
  Double_t *v;
  
  Double_t phi0;
  Double_t z0;
  Double_t r0;


  Double_t phi1;
  Double_t z1;
  Double_t r1;

  Int_t nd = stepR * stepPhi * stepZ;

  a = new Double_t[nd*nd];
  r = new Double_t[nd];
  v = new Double_t[nd];


  Int_t startPhi = phiindex - stepPhi/2;
  Int_t index_phi;
  Int_t index_phi1;

  Int_t startR =  rindex  - stepR/2;  
  Int_t startZ =  zindex  - stepZ/2;
  
  if (startPhi < 0) startPhi = fNPhi + startPhi;

  if (startR    < 0) startR = 0;
  if (startR + stepR >= fNR ) startR = fNR - stepR ;

  if (startZ    < fMinZIndex) startZ = fMinZIndex;
  if (startZ + stepZ  >= fNZ  ) startZ = fNZ - stepZ;


/// calculate distance between each poInt_t

  index0 = 0;

	for ( i=startPhi; i < startPhi + stepPhi;i++) 
	{
	    index_phi = i % fNPhi;	

		for (j=startR; j < startR + stepR;j++) 
		{
			for(k=startZ; k < startZ + stepZ;k++) 
			{
				indexCyl0 = index_phi * fNR * fNZ + j * fNZ + k;

			

			    r0 =  fRlist[indexCyl0];
			    z0 =  fZlist[indexCyl0];
		      	    phi0 =  fPhilist[indexCyl0];

				index1 = 0;	
				for ( ii=startPhi; ii < startPhi + stepPhi;ii++) {
					index_phi1 = ii % fNPhi;
					for (jj=startR; jj < startR + stepR;jj++) {
						for(kk=startZ; kk < startZ + stepZ;kk++) {
							indexCyl1 = index_phi1 * fNR * fNZ + jj * fNZ + kk;
			    			r1 =  fRlist[indexCyl1];
		      				z1 =  fZlist[indexCyl1];
		      				phi1 =  fPhilist[indexCyl1];

						
							r[index1] = Distance(r0 ,phi0,z0,r1,phi1,z1);
	
						//	if (zindex < 14)
						//		printf("[%d,%d](%f,%f,%f)(%f,%f,%f) = r[%d] = %f\n",indexCyl0,indexCyl1,r0,phi0,z0,r1,phi1,z1,index1,r[index1]);
							index1++;
						}
					}
				}			
		
/**
		  	TVectorD vec_r;

			vec_r.Use(nd,r);

			printf("r: \n");
			vec_r.Print();



		  	TVectorD vec_v;

			vec_v.Use(nd,v);

			printf("v: \n");
			vec_v.Print();
**/
	
			    Phi ( nd, r, radius0, v );		

				index1 = 0;	
				for ( ii=startPhi; ii < startPhi + stepPhi;ii++) 
				{
					index_phi1 = ii % fNPhi;
					for (jj=startR; jj < startR + stepR;jj++) 
					{
						for(kk=startZ; kk < startZ + stepZ;kk++) 
						{
				      		a[index0 * nd+ index1] = v[index1];
							index1++;
						}
					}
				}

				w[index0] = fVals[indexCyl0];
				index0++;
				

			}
		}

		
	 }
	
//  Solve for the weights.

  	TMatrixD mat_a;
  	mat_a.Use(nd,nd,a);
	TVectorD vec_w;

//	mat_a.Print();


  	vec_w.Use(nd,w);
  	TDecompSVD svd(mat_a);

   	svd.Solve(vec_w);


 	

	delete [] a;
 	delete [] r;
  	delete [] v;
	

	
	//weight = w;
 // return  w;
}




void AliTPC3DCylindricalInterpolatorFull::rbf1 ( Int_t n, Double_t r[], Double_t r0, Double_t v[] )

{
  Int_t i;

  for ( i = 0; i < n; i++ )
  {
    v[i] = sqrt ( 1 + (r[i] * r[i] + r0 * r0 ));
  }
  return;
}
//****************************************************************************80

void AliTPC3DCylindricalInterpolatorFull::rbf2 ( Int_t n, Double_t r[], Double_t r0, Double_t v[] )

//****************************************************************************80
//
//  Purpose:
//
//    PHI2 evaluates the inverse multiquadric radial basis function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 June 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
//    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
//    Third Edition,
//    Cambridge University Press, 2007,
//    ISBN13: 978-0-521-88068-8,
//    LC: QA297.N866.
//
//  Parameters:
//
//    Input, Int_t N, the number of points.
//
//    Input, Double_t R[N], the radial separation.
//    0 < R.
//
//    Input, Double_t R0, a scale factor.
//
//    Output, Double_t V[N], the value of the radial basis function.
//
{
  Int_t i;

  for ( i = 0; i < n; i++ )
  {
    v[i] = 1.0 / sqrt (1 + ( r[i] * r[i] + r0 * r0 ));
  }
  return;
}
//****************************************************************************80

void AliTPC3DCylindricalInterpolatorFull::rbf3 ( Int_t n, Double_t r[], Double_t r0, Double_t v[] )

//****************************************************************************80
//
//  Purpose:
//
//    PHI3 evaluates the thin-plate spline radial basis function.
//
//  Discussion:
//
//    Note that PHI3(R,R0) is negative if R < R0.  Thus, for this basis function,
//    it may be desirable to choose a value of R0 smaller than any possible R.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 June 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
//    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
//    Third Edition,
//    Cambridge University Press, 2007,
//    ISBN13: 978-0-521-88068-8,
//    LC: QA297.N866.
//
//  Parameters:
////    Input, Int_t N, the number of points.
//
//    Input, Double_t R[N], the radial separation.
//    0 < R.
//
//    Input, Double_t R0, a scale factor.
//
//    Output, Double_t V[N], the value of the radial basis function.
//
{
  Int_t i;

  for ( i = 0; i < n; i++ )
  {
    if ( r[i] <= 0.0 )
    {
      v[i] = 0.0;
    }
    else
    {
      v[i] = r[i] * r[i] * log ( r[i] / r0 );
    }
  }
  return;
}
//****************************************************************************80

void AliTPC3DCylindricalInterpolatorFull::rbf4 ( Int_t n, Double_t r[], Double_t r0, Double_t v[] )

//****************************************************************************80
//
//  Purpose:
//
//    PHI4 evaluates the gaussian radial basis function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 June 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
//    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
//    Third Edition,
//    Cambridge University Press, 2007,
//    ISBN13: 978-0-521-88068-8,
//    LC: QA297.N866.
//
//  Parameters:
//
//    Input, Int_t N, the number of points.
//
//    Input, Double_t R[N], the radial separation.
//    0 < R.
//
//    Input, Double_t R0, a scale factor.
//
//    Output, Double_t V[N], the value of the radial basis function.
//
{
  Int_t i;

  for ( i = 0; i < n; i++ )
  {
    v[i] = TMath::Exp ( - 0.5 * r[i] * r[i] / (r0 * r0) );
  }
  return;
}




// RBF based interpolation
// return interpolated value
Double_t AliTPC3DCylindricalInterpolatorFull::InterpRBF
( 
	Double_t r,
	Double_t phi,
	Double_t z,
	Int_t startR,
	Int_t startPhi,
	Int_t startZ,
	Int_t stepR,
	Int_t stepPhi,
	Int_t stepZ,
	Double_t radius0,
	Int_t kernelType,
	Double_t * weight
)
{
	Double_t interpVal = 0.0;
	Double_t r0,z0,phi0;
	Double_t *dlist;
	Double_t *v;

	Int_t indexCyl0, index0, index_phi;

	Int_t nd = stepR * stepPhi * stepZ;

  	dlist = new Double_t[nd];
 	v = new Double_t[nd];


	index0 = 0;
	for (Int_t i=startPhi; i < startPhi + stepPhi;i++) 
	{
	    index_phi = i % fNPhi;	

		for (Int_t j=startR; j < startR + stepR;j++) 
		{
			for(Int_t k=startZ; k < startZ + stepZ;k++) 
			{

				indexCyl0 = index_phi * fNR * fNZ + j * fNZ + k;

				r0 =  fRlist[indexCyl0];
		      	z0 =  fZlist[indexCyl0];
		      	phi0 =  fPhilist[indexCyl0];

										
				dlist[index0] = Distance(r,phi,z,r0,phi0,z0);
				index0++;
				
			}
		}
	}

	Phi ( nd, dlist, radius0, v );		


  	TVectorD vec_v;
	vec_v.Use(nd,v);

  	TVectorD vec_w;
	vec_w.Use(nd,weight);
	
	interpVal = vec_v * vec_w;
	delete[] v;	
	delete[] dlist;
	return interpVal;
}


// RBF based interpolation
// return interpolated value
Double_t AliTPC3DCylindricalInterpolatorFull::InterpRBFCartesian
( 
	Double_t r,
	Double_t phi,
	Double_t z,
	Int_t startR,
	Int_t startPhi,
	Int_t startZ,
	Int_t stepR,
	Int_t stepPhi,
	Int_t stepZ,
	Double_t radius0,
	Int_t kernelType,
	Double_t * weight
)
{
	Double_t interpVal = 0.0;
	Double_t r0,z0,phi0;
	Double_t *dlist;
	Double_t *v;

	Int_t indexCyl0, index0, index_phi;

	Int_t nd = stepR * stepPhi * stepZ;

  	dlist = new Double_t[nd];
 	v = new Double_t[nd];

	Double_t x = r * TMath::Cos(phi);
	Double_t y = r * TMath::Sin(phi);

	index0 = 0;
	for (Int_t i=startPhi; i < startPhi + stepPhi;i++) 
	{
	    index_phi = i % fNPhi;	

		for (Int_t j=startR; j < startR + stepR;j++) 
		{
			for(Int_t k=startZ; k < startZ + stepZ;k++) 
			{

				indexCyl0 = index_phi * fNR * fNZ + j * fNZ + k;

				r0 =  fRlist[indexCyl0];
		      		z0 =  fZlist[indexCyl0];
		      		phi0 =  fPhilist[indexCyl0];

										
				dlist[index0] = Distance(x,y,z,r0,phi0,z0);
				index0++;
				
			}
		}
	}

	Phi ( nd, dlist, radius0, v );		


  	TVectorD vec_v;
	vec_v.Use(nd,v);

  	TVectorD vec_w;
	vec_w.Use(nd,weight);
	
	interpVal = vec_v * vec_w;
	delete[] v;	
	delete[] dlist;
	return interpVal;
}

// InitWeightRBF
// RBF init until z == jy
void AliTPC3DCylindricalInterpolatorFull::InitRBFWeight
(
	Int_t jy
) 
{

}





// calculate
// RBFWeight for all points in the interpolation
void AliTPC3DCylindricalInterpolatorFull::GetRBFWeight 
( 
	Int_t rindex,
	Int_t zindex,
	Int_t phiindex,
	Int_t stepR,
	Int_t stepPhi,
	Int_t stepZ,
	Double_t radius0,
	Int_t kernelType,
	Double_t * w
)
{

  Int_t index = phiindex * fNR * fNZ + rindex * fNZ + zindex;
  if (fRBFWeightLookUp[index] == 0)
  {
	//printf("compute RBF weights\n");
  	RBFWeight ( 
		rindex,
		zindex,
 		phiindex,
		stepR,
		stepPhi,
		stepZ,
		radius0,
		kernelType,
		w
  	);

//	if ((stepR == fStepR) && (stepZ == fStepZ) && (stepPhi == fStepPhi) && (zindex > fMinZIndex+fStepZ))
//	if ((stepR == fStepR) && (stepZ == fStepZ) && (stepPhi == fStepPhi) && (zindex > fMinZIndex+fStepZ))  
//	{
	fRBFWeightLookUp[index] = 1;
		// copy to lookup
	Int_t nd = stepR*stepZ*stepPhi;
		

	for (Int_t i=0;i<nd;i++) fRBFWeight[index * nd + i] = w[i];
//	}
  } else {

	//printf("loo. RBF weights\n");
	Int_t ndw = stepR*stepZ*stepPhi;

	Int_t nd  = fStepR*fStepZ*fStepPhi;
	Int_t indexw = phiindex * fNR * fNZ *nd + rindex * fNZ *nd + zindex *nd;


//	if (ndw == nd)
//	printf("lookup weight %d\n",zindex);
	for (Int_t i=0;i<nd;i++)  w[i] = fRBFWeight[indexw + i];
//	else
/**	RBFWeight ( 
		rindex,
		zindex,
 		phiindex,
		stepR,
		stepPhi,
		stepZ,
		radius0,
		kernelType,
		w
  	);
**/
  }
}


// calculate
// RBFWeight for all points in the interpolation
void AliTPC3DCylindricalInterpolatorFull::GetRBFWeightCartesian
( 
	Int_t rindex,
	Int_t zindex,
	Int_t phiindex,
	Int_t stepR,
	Int_t stepPhi,
	Int_t stepZ,
	Double_t radius0,
	Int_t kernelType,
	Double_t * w
)
{

  Int_t index = phiindex * fNR * fNZ + rindex * fNZ + zindex;
  if (fRBFWeightLookUp[index] == 0)
  {
	//printf("compute RBF weights\n");
  	RBFWeightCartesian ( 
		rindex,
		zindex,
 		phiindex,
		stepR,
		stepPhi,
		stepZ,
		radius0,
		kernelType,
		w
  	);

//	if ((stepR == fStepR) && (stepZ == fStepZ) && (stepPhi == fStepPhi) && (zindex > fMinZIndex+fStepZ))
//	if ((stepR == fStepR) && (stepZ == fStepZ) && (stepPhi == fStepPhi) && (zindex > fMinZIndex+fStepZ))  
//	{
	fRBFWeightLookUp[index] = 1;
		// copy to lookup
	Int_t nd = stepR*stepZ*stepPhi;
		

	for (Int_t i=0;i<nd;i++) fRBFWeight[index * nd + i] = w[i];
//	}
  } else {

	//printf("loo. RBF weights\n");
	Int_t ndw = stepR*stepZ*stepPhi;

	Int_t nd  = fStepR*fStepZ*fStepPhi;
	Int_t indexw = phiindex * fNR * fNZ *nd + rindex * fNZ *nd + zindex *nd;


//	if (ndw == nd)
//	printf("lookup weight %d\n",zindex);
	for (Int_t i=0;i<nd;i++)  w[i] = fRBFWeight[indexw + i];
//	else
/**	RBFWeight ( 
		rindex,
		zindex,
 		phiindex,
		stepR,
		stepPhi,
		stepZ,
		radius0,
		kernelType,
		w
  	);
**/
  }
}



// calculate
// RBFWeight for all points in the interpolation
void AliTPC3DCylindricalInterpolatorFull::GetRBFWeightHalf
( 
	Int_t rindex,
	Int_t zindex,
	Int_t phiindex,
	Int_t stepR,
	Int_t stepPhi,
	Int_t stepZ,
	Double_t radius0,
	Int_t kernelType,
	Double_t * w
)
{

  Int_t index = phiindex * fNR * fNZ + rindex * fNZ + zindex;

  if (fRBFWeightLookUp[index] == 0)
  {
  	RBFWeightHalf ( 
		rindex,
		zindex,
 		phiindex,
		stepR,
		stepPhi,
		stepZ,
		radius0,
		kernelType,
		w
  	);

	if ((stepR == fStepR) && (stepZ == fStepZ) && (stepPhi == fStepPhi) && (zindex > fMinZIndex+fStepZ)) {
		fRBFWeightLookUp[index] = 1;
		// copy to lookup
		Int_t nd = stepR+stepZ+stepPhi-2;
		

		for (Int_t i=0;i<nd;i++) fRBFWeight[index * nd + i] = w[i];
	}
  } else {


	//Int_t ndw = stepR*stepZ*stepPhi;
	Int_t nd = stepR+stepZ+stepPhi-2;
	Int_t indexw = phiindex * fNR * fNZ *nd + rindex * fNZ *nd + zindex *nd;


//	if (ndw == nd)
//	printf("lookup weight %d\n",zindex);
	for (Int_t i=0;i<nd;i++)  w[i] = fRBFWeight[indexw + i];
//	else
/**	RBFWeight ( 
		rindex,
		zindex,
 		phiindex,
		stepR,
		stepPhi,
		stepZ,
		radius0,
		kernelType,
		w
  	);
**/
  }
}

// calculate
// RBFWeight for all points in the interpolation
// half cubes (not included 
void AliTPC3DCylindricalInterpolatorFull::RBFWeightHalf
( 
	Int_t rindex,
	Int_t zindex,
	Int_t phiindex,
	Int_t stepR,
	Int_t stepPhi,
	Int_t stepZ,
	Double_t radius0,
	Int_t kernelType,
	Double_t * w
)

{
  Double_t *a;
  Int_t i;
  Int_t j;
  Int_t k;
  Int_t ii;
  Int_t jj;
  Int_t kk;

  Int_t index0,index1;
  Int_t indexCyl0,indexCyl1;
  Double_t *r;
  Double_t *v;
  
  Double_t phi0;
  Double_t z0;
  Double_t r0;


  Double_t phi1;
  Double_t z1;
  Double_t r1;

  Int_t nd = (stepR - 1) + (stepPhi - 1) + (stepZ - 1)  + 1;

  
  a = new Double_t[nd*nd];
  r = new Double_t[nd];
  v = new Double_t[nd];


  Int_t startPhi = phiindex - stepPhi/2;
  Int_t index_phi;
  Int_t index_phi1;

  Int_t startR =  rindex  - stepR/2;  
  Int_t startZ =  zindex  - stepZ/2;
  
  if (startPhi < 0) startPhi = fNPhi + startPhi;

  if (startR    < 0) startR = 0;
  if (startR + stepR >= fNR ) startR = fNR - stepR ;

  if (startZ    < fMinZIndex) startZ = fMinZIndex;
  if (startZ + stepZ  >= fNZ  ) startZ = fNZ - stepZ;


/// calculate distance between each poInt_t

  index0 = 0;

	
	



	for ( i=startPhi; i < startPhi + stepPhi;i++) 
	{
	    index_phi = i % fNPhi;	

		for (j=startR; j < startR + stepR;j++) 
		{
			for(k=startZ; k < startZ + stepZ;k++) 
			{

				if (
					(i ==  (startPhi + stepPhi/2)  && 	j ==  (startR + stepR/2)) ||
					(i ==  (startPhi + stepPhi/2)  && 	k ==  (startZ + stepZ/2)) ||
					(j ==  (startR + stepR/2)  && 	k ==  (startZ + stepZ/2)) 
				   )
				{
					indexCyl0 = index_phi * fNR * fNZ + j * fNZ + k;

			

				    r0 =  fRlist[indexCyl0];
			      	z0 =  fZlist[indexCyl0];
			      	phi0 =  fPhilist[indexCyl0];

					index1 = 0;	
					for ( ii=startPhi; ii < startPhi + stepPhi;ii++) {
						index_phi1 = ii % fNPhi;
						for (jj=startR; jj < startR + stepR;jj++) {
							for(kk=startZ; kk < startZ + stepZ;kk++) {
								if (
									(ii ==  (startPhi + stepPhi/2)  && 	jj ==  (startR + stepR/2)) ||
									(ii ==  (startPhi + stepPhi/2)  && 	kk ==  (startZ + stepZ/2)) ||
									(jj ==  (startR + stepR/2)  && 	kk ==  (startZ + stepZ/2)) 
								   )
								{

			
									indexCyl1 = index_phi1 * fNR * fNZ + jj * fNZ + kk;
					    			r1 =  fRlist[indexCyl1];
				      				z1 =  fZlist[indexCyl1];
				      				phi1 =  fPhilist[indexCyl1];

						
									r[index1] = Distance(r0,phi0,z0,r1,phi1,z1);

				//					printf("[%d,%d](%f,%f,%f)(%f,%f,%f) = r[%d] = %f\n",indexCyl0,indexCyl1,r0,phi0,z0,r1,phi1,z1,index1,r[index1]);
									index1++;
								}
							}
						}
					}			
		
/**
		  	TVectorD vec_r;

			vec_r.Use(nd,r);

			printf("r: \n");
			vec_r.Print();



		  	TVectorD vec_v;

			vec_v.Use(nd,v);

			printf("v: \n");
			vec_v.Print();
**/
	
			    	Phi ( nd, r, radius0, v );		

					index1 = 0;	
					for ( ii=startPhi; ii < startPhi + stepPhi;ii++) 
					{
						index_phi1 = ii % fNPhi;
						for (jj=startR; jj < startR + stepR;jj++) 
						{
							for(kk=startZ; kk < startZ + stepZ;kk++) 
							{
								if (
									(ii ==  (startPhi + stepPhi/2)  && 	jj ==  (startR + stepR/2)) ||
									(ii ==  (startPhi + stepPhi/2)  && 	kk ==  (startZ + stepZ/2)) ||
									(jj ==  (startR + stepR/2)  && 	kk ==  (startZ + stepZ/2)) 
								   )
								{
					      			a[index0 * nd+ index1] = v[index1];
									index1++;
								}
							}
						}
					}

					w[index0] = fVals[indexCyl0];
					index0++;
				
				}
			}
		}

		
	 }
	
//  Solve for the weights.

  	TMatrixD mat_a;
  	mat_a.Use(nd,nd,a);
	TVectorD vec_w;

  	vec_w.Use(nd,w);
  	TDecompSVD svd(mat_a);

   	svd.Solve(vec_w);


//  	vec_w.Print();


	delete [] a;
 	delete [] r;
  	delete [] v;
	

	
	//weight = w;
 // return  w;
}


// RBF based interpolation
// return interpolated value
// half points
Double_t AliTPC3DCylindricalInterpolatorFull::InterpRBFHalf
( 
	Double_t r,
	Double_t phi,
	Double_t z,
	Int_t startR,
	Int_t startPhi,
	Int_t startZ,
	Int_t stepR,
	Int_t stepPhi,
	Int_t stepZ,
	Double_t radius0,
	Int_t kernelType,
	Double_t * weight
)
{
	Double_t interpVal = 0.0;
	Double_t r0,z0,phi0;
	Double_t *dlist;
	Double_t *v;

	Int_t indexCyl0, index0, index_phi;

//	Int_t nd = stepR * stepPhi * stepZ;
  	Int_t nd = (stepR - 1) + (stepPhi - 1) + (stepZ - 1)  + 1;


  	dlist = new Double_t[nd];
 	v = new Double_t[nd];


	index0 = 0;
	for (Int_t i=startPhi; i < startPhi + stepPhi;i++) 
	{
	    index_phi = i % fNPhi;	

		for (Int_t j=startR; j < startR + stepR;j++) 
		{
			for(Int_t k=startZ; k < startZ + stepZ;k++) 
			{
				if (
					(i ==  (startPhi + stepPhi/2)  && 	j ==  (startR + stepR/2)) ||
					(i ==  (startPhi + stepPhi/2)  && 	k ==  (startZ + stepZ/2)) ||
					(j ==  (startR + stepR/2)  && 	k ==  (startZ + stepZ/2)) 
				) 
				{
	
					indexCyl0 = index_phi * fNR * fNZ + j * fNZ + k;

					r0 =  fRlist[indexCyl0];
			      	z0 =  fZlist[indexCyl0];
			      	phi0 =  fPhilist[indexCyl0];

										
					dlist[index0] = Distance(r,phi,z,r0,phi0,z0);
					index0++;
				}			
			}
		}
	}

	Phi ( nd, dlist, radius0, v );		


  	TVectorD vec_v;
	vec_v.Use(nd,v);

  	TVectorD vec_w;
	vec_w.Use(nd,weight);
	
	interpVal = vec_v * vec_w;
	delete[] v;	
	delete[] dlist;
	return interpVal;
}




// set Radius0
void AliTPC3DCylindricalInterpolatorFull::Phi
( 
	Int_t n, 
	Double_t r[], 
	Double_t r0, 
	Double_t v[] 
)
{


	switch (fKernelType) {
		case kRBFMultiQuadratic:		  				
			rbf1(n,r,r0,v);
		  	break;
		case kRBFInverseMultiQuadratic:
			rbf2(n,r,r0,v);
		  	break;
		case kRBFThinPlateSpline:
		  	rbf3(n,r,r0,v);
		  	break;
		case kRBFGaussian:
		  	rbf4(n,r,r0,v);
		  	break;

		default:
			rbf1(n,r,r0,v);
			break;
	}
}

//
Double_t AliTPC3DCylindricalInterpolatorFull::GetRadius0RBF
(
	const Int_t indexr,
	const Int_t indexphi,
	const Int_t indexz
)
{
  	const Float_t gridSizeR   =  (AliTPCPoissonSolver::fgkOFCRadius-AliTPCPoissonSolver::fgkIFCRadius) / (fNR - 1) ;
	const Float_t gridSizeZ   =  AliTPCPoissonSolver::fgkTPCZ0 / (fNZ -1) ;
  	const Float_t gridSizePhi =  TMath::TwoPi() / fNPhi;
 	Int_t startPhi 	= indexphi - fStepPhi/2;
  
	Int_t startR 	=  indexr  - fStepR/2;  
  	Int_t startZ 	=  indexz  - fStepZ/2;
  
  	if (startPhi < 0) startPhi = fNPhi + startPhi;
 

  	if (startR    < 0) startR = 0;
  	if (startR +fStepR >= fNR ) startR = fNR - fStepR ;

 	if (startZ    < 0) startZ = 0;
  	if (startZ + fStepZ  >= fNZ  ) startZ = fNZ - fStepZ;

	Double_t r0   = AliTPCPoissonSolver::fgkIFCRadius + (startR  * gridSizeR);
	Double_t phi0 = startPhi * gridSizePhi;
	Double_t z0 =	startZ * gridSizeZ;

	Double_t r1   = AliTPCPoissonSolver::fgkIFCRadius + (startR  * gridSizeR);
	Double_t phi1 = (startPhi + 1) * gridSizePhi;
	Double_t z1 =	(startZ+1) * gridSizeZ;

	if (fKernelType == kRBFThinPlateSpline)
		r0   = AliTPCPoissonSolver::fgkIFCRadius + ((startR - 1)  * gridSizeR);	
	else
		r0   = AliTPCPoissonSolver::fgkIFCRadius + (startR  * gridSizeR);

	return Distance(r0,0.0,0.0,r0 +  gridSizeR  ,gridSizePhi,   gridSizeR);
	
/**
	Int_t indexCyl0 = startPhi * (fNR * fNZ) + startR * fNZ + startZ;

	startPhi = (startPhi + 1) % fNPhi;

	Int_t indexCyl1 = startPhi * (fNR * fNZ) + (startR + 1) * fNZ + (startZ + 1);

	if (fKernelType != kRBFThinPlateSpline)
		return Distance(fRlist[indexCyl0],fPhilist[indexCyl0],fZlist[indexCyl0],fRlist[indexCyl1],fPhilist[indexCyl1],fZlist[indexCyl1]);
	else
		return Distance(r0,0.0,0.0,r0 + gridSizeR , gridSizePhi,  gridSizeZ);
**/
}
//return Distance(r0,0.0,0.0,r0 + gridSizeR , gridSizePhi,  gridSizeZ);
