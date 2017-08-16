#include <TStopwatch.h>
#include "TPCSimulation/AliTPC3DCylindricalInterpolator.h"


/// \cond CLASSIMP3
ClassImp(o2::TPC::AliTPC3DCylindricalInterpolator)
/// \endcond

using namespace o2::TPC;

/// constructor
AliTPC3DCylindricalInterpolator::AliTPC3DCylindricalInterpolator()
{	
  fOrder = 1;
  fIsAllocatingLookUp = kFALSE;
  fIsInitCubic = kFALSE;
}

/// destructor
AliTPC3DCylindricalInterpolator::~AliTPC3DCylindricalInterpolator()
{
	
  if  (fIsAllocatingLookUp) {		
    delete fVals;		
    //delete fRlist;
    //delete fPhilist;
    //delete fZlist;
  }
  if (fIsInitCubic) {
    delete fSecondDerZ;
  }
	
}

/// main operation
Double_t AliTPC3DCylindricalInterpolator::Interpolate3DTableCyl
( 
 Double_t r,   
 Double_t z,
 Double_t phi
  ) 
{
  Int_t ilow = 0, jlow = 0, klow = 0, m=0;
  Int_t indeks;
  
  // tricubic points
  Double_t saveArray[fOrder + 1];
  Double_t savedArray[fOrder + 1];
  Double_t zlistM1[fOrder + 1];
  Double_t valsM1[fOrder + 1];

  Bool_t neg = kFALSE;
  // check phi 
  while (phi < 0.0) phi = TMath::TwoPi() + phi;
  while (phi > TMath::TwoPi()) phi = phi - TMath::TwoPi() ;

   	
  // search lowest index related to r,z and phi
  Search( fNR, fRlist, r, ilow   ) ;
  Search( fNZ, fZlist, z, jlow   ) ; 
  Search( fNPhi, fPhilist, phi, klow   ) ;
  
  
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
    // interpolate	
    for ( Int_t i = ilow ; i < ilow + fOrder + 1 ; i++ ) {

      if (fOrder <= 2) {

	//				printf("fval (%f,%f,%f)\n",fVals[indeks],fVals[indeks + 1],fVals[indeks + 2]);
	if (jlow >= 0) {
	  indeks =  m* (fNZ * fNR) +  i *(fNZ)  + jlow;
	  saveArray[i-ilow] = Interpolate( &fZlist[jlow], &fVals[indeks], z );
	}
	else  
	  {
	    indeks =  m* (fNZ * fNR) +  i *(fNZ);
	    zlistM1[0] = fZlist[0] - (fZlist[1] - fZlist[0]);
	    zlistM1[1] = fZlist[0];
	    zlistM1[2] = fZlist[1];
	    valsM1[0] = fVals[indeks] - (fVals[indeks + 1] - fVals[indeks]);
	    valsM1[1] = fVals[indeks];
	    valsM1[2] = fVals[indeks + 1];
	    saveArray[i-ilow] = Interpolate( &zlistM1[0], &valsM1[0], z );
	  }
				  
      }	else {
	indeks =  m* (fNZ * fNR) +  i *(fNZ);
	saveArray[i-ilow] = InterpolateCubicSpline(fZlist,&fVals[indeks],&fSecondDerZ[indeks],fNZ,fNZ,fNZ,z,1);
      }
    }
    //  printf("savear (%f,%f,%f)\n",saveArray[0],saveArray[1],saveArray[2]);

    savedArray[k] = Interpolate( &fRlist[ilow], saveArray, r )  ;
    //table.Print();
  }

  //  printf("(%f,%f,%f)\n",savedArray[0],savedArray[1],savedArray[2]);
  return( InterpolatePhi( &fPhilist[0], klow, fNPhi,  savedArray,  phi ) )   ;
  //	return 0.0;
}



// search
void AliTPC3DCylindricalInterpolator::Search
( 
 Int_t n, 
 const Double_t xArray[], 
 Double_t x, 
 Int_t &low 
  ) 
{
  /// Search an ordered table by starting at the most recently used point

  Long_t middle, high ;
  Int_t  ascend = 0, increment = 1 ;

  if ( xArray[n-1] > xArray[0] ) ascend = 1 ;  // Ascending ordered table if true

  if ( low < 0 || low > n-1 ) {
    low = -1 ; high = n ;
  } else {                                            // Ordered Search phase
    if ( (Int_t)( x > xArray[low] ) == ascend )  {
      if ( low == n-1 ) return ;
      high = low + 1 ;
      while ( (Int_t)( x > xArray[high] ) == ascend ) {
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
    if ( (Int_t)( x > xArray[middle] ) == ascend )
      low = middle ;
    else
      high = middle ;
  }

  if ( x >  xArray[n-1] ) low = n;
  if ( x <  xArray[0]   ) low = -1;

}



// Interpolate 1D
Double_t AliTPC3DCylindricalInterpolator::Interpolate
( 
 Double_t xArray[], 
 Double_t yArray[],
 Double_t x 
  ) 
{
  // Interpolate function Y(x) using linear (order=1) or quadratic (order=2) interpolation.
 
  Double_t y ;
  //  printf("order: %d x: (%f,%f,%f)\n",fOrder,xArray[0],xArray[1],x);
  
  // printf("y: (%f,%f)\n",yArray[0],yArray[1]);
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
    // Linear Interpolation = 1
    /**
       if ( TMath::Abs(x - xArray[0]) < 1e-3 )
       {
       printf("%E == %E)\n",x,xArray[0]);
	
       y = yArray[0];
       }
       else if (TMath::Abs(x-xArray[1]) < 1e-3 ) 
       {      
       printf("%E == %E)\n",x,xArray[1]);
       y = yArray[1];
       }
       else
    **/
    y  = yArray[0] + ( yArray[1]-yArray[0] ) * ( x-xArray[0] ) / ( xArray[1] - xArray[0] ) ;
  }

  return (y);

}


// interpolate phi
Double_t AliTPC3DCylindricalInterpolator::InterpolatePhi
( 
 Double_t xArray[], 
 const Int_t ilow,
 const Int_t nx,
 Double_t yArray[],
 Double_t x 
  )
{
  /// Interpolate function Y(x) using linear (order=1) or quadratic (order=2) interpolation.

  Int_t i0 = ilow;
  Double_t xi0 = xArray[ilow];
	
  Int_t i1 = (ilow + 1) % nx;
  Double_t xi1 = xArray[i1];
  Int_t i2 = (ilow + 2) % nx;
  Double_t xi2 = xArray[i2];
  /**
     if ((ilow + 1) >= nx) {
     xi1 += TMath::TwoPi();
     }
     if ((ilow + 2) >= nx) {
     xi2 += TMath::TwoPi();
		
     }
     if (x < 0.0)
     printf("(%f,%f,%f,%f)\n",xi0,xi1,xi2,x);
  **/
  if (xi1 < xi0) xi1 = TMath::TwoPi() + xi1;
  if (xi2 < xi1) xi2 = TMath::TwoPi() + xi2;
  if (x   < xi0) x   = TMath::TwoPi() + x;

		
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
    //	printf("y (%f,%f,%f,%f)\n",xi0,xi1,xi2,y);

    y += (x-xi2) * (x-xi0) * yArray[1] / ( (xi1-xi2) * (xi1-xi0) ) ;
    //	printf("y (%f,%f,%f,%f)\n",xi0,xi1,xi2,y);
    y += (x-xi0) * (x-xi1) * yArray[2] / ( (xi2-xi0) * (xi2-xi1) ) ;
    //	printf("y (%f,%f,%f,%f)\n",xi0,xi1,xi2,y);



  } else {                           // Linear Interpolation = 1
    y  = yArray[0] + ( yArray[1]-yArray[0] ) * ( x-xi0 ) / (xi1 - xi0 ) ;
  }

  return (y);

}

// get value
Double_t AliTPC3DCylindricalInterpolator::GetValue
(
 Double_t r, 
 Double_t phi, 
 Double_t z
 ) 
{  
	
  return Interpolate3DTableCyl(r,z,phi);
	
}



// natural cubic spline
void AliTPC3DCylindricalInterpolator::InitCubicSpline
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
  for (Int_t k=n-2;k>=0;k--)
    y2Array[k * skip] = y2Array[k * skip] * y2Array[(k+1) * skip] + u[k];	
}




// clamped cubic spline
// yp0 =>
// ypn1 =>
void AliTPC3DCylindricalInterpolator::InitCubicSpline
(	
 Double_t * xArray,
 Double_t * yArray,
 const Int_t n,
 Double_t * y2Array,
 const Int_t skip,
 Double_t yp0,
 Double_t ypn1
	) 
{
  Double_t u[n];
  Double_t sig,p,qn,un;

		
  y2Array[0] = 0.0;
  u[0] = 0.0; //natural condition

  //	y2Array[0] = -0.5;
  //	u[0]=(3.0/(xArray[1]-xArray[0]))*((yArray[skip]-yArray[0])/(xArray[1]-xArray[0])-yp0);
		
  for (Int_t i=1;i<=n-2;i++) {
    sig = (xArray[i] - xArray[i-1])/(xArray[i+1] - xArray[i-1]);
    p = sig * y2Array[(i-1) * skip]+ 2.0;
    y2Array[i * skip] = (sig - 1.0)/p;
    u[i] = (yArray[(i+1) * skip] - yArray[i * skip])/(xArray[i+1] - xArray[i]) - (yArray[i * skip] - yArray[(i-1) * skip])/(xArray[i] - xArray[i-1]);
    u[i] = (6.0 * u[i]/(xArray[i+1] - xArray[i-1]) - sig * u[i-1])/p;	
  }
	
	
	
	
  qn = un = 0.0;
  //	qn=0.5;
  //s	un=(3.0/(xArray[n-1]-xArray[n-2]))*(ypn1-(yArray[(n-1) * skip]-yArray[(n-2) * skip])/(xArray[n-1]-xArray[n-2]));
	
  y2Array[(n-1) * skip] = (un - qn * u[n-2])/(qn*y2Array[(n-2) * skip]+1.0);
  for (Int_t k=n-2;k>=0;k--)
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
Double_t AliTPC3DCylindricalInterpolator::InterpolateCubicSpline
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
	
	
  if (h < 1e-10)  {
    printf("not found\n");
    return 0.0;
  }
	
  a = (xArray[khi] - x) / h;
  b = (x - xArray[klo]) / h;
	
  Double_t y = a *yArray[klo ] + b * yArray[khi ] + ((a*a*a -a) * y2Array[klo * skip] + (b*b*b -b) * y2Array[khi * skip]) * (h*h)/6.0;
	
  return y;
	
}



// spline from recepi
Double_t AliTPC3DCylindricalInterpolator::InterpolateCubicSplinePhi
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

void AliTPC3DCylindricalInterpolator::InitCubicSpline() 
{

  Double_t yp0, ypn1;
	
  if (fIsInitCubic != kTRUE) {
		
    fSecondDerZ = new Double_t[fNR * fNZ * fNPhi]; 
		
    TStopwatch w;
		
    w.Start();
    // Init at Z direction
    for (Int_t m=0;m < fNPhi;m++)
      {			
	for (Int_t i=0; i< fNR;i++)		
	  {
				
	    /// first derivative at endpoint

	    // forward derivative
	    //yp0 =  (fVals[(m * (fNZ * fNR) + i * fNZ) + 1] - fVals[(m * (fNZ * fNR) + i * fNZ) ]) / (fZlist[1] - fZlist[0]);


	    yp0 = ( -(11.0/6.0)* fVals[(m * (fNZ * fNR) + i * fNZ) ] + (3.0*fVals[(m * (fNZ * fNR) + i * fNZ) + 1]) - (1.5*fVals[(m * (fNZ * fNR) + i * fNZ) +2]) + ((1.0/3.0)* fVals[(m * (fNZ * fNR) + i * fNZ) + 4]) )  / (fZlist[1] - fZlist[0]);
	    ypn1 = ( -(11.0/6.0)*fVals[(m * (fNZ * fNR) + i * fNZ) + (fNZ - 1)]  + (3.0*fVals[(m * (fNZ * fNR) + i * fNZ) + (fNZ - 2)] ) - (1.5*fVals[(m * (fNZ * fNR) + i * fNZ) + (fNZ - 3)] ) + ((1.0/3.0)*fVals[(m * (fNZ * fNR) + i * fNZ) + (fNZ - 4)] ) )  /  (fZlist[0] - fZlist[1]);
	    InitCubicSpline
	      (		
	       fZlist, 
	       &fVals[m * (fNZ * fNR) + i * fNZ],
	       fNZ,
	       &fSecondDerZ[m * (fNZ * fNR) + i * fNZ],
	       1,
	       yp0,
	       ypn1
			);
	  }

      }				
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


void AliTPC3DCylindricalInterpolator::SetVals(TMatrixD ** mvals) 
{
  Int_t indeksm;
  Int_t indeks;
  TMatrixD * mat;
  if (!fIsAllocatingLookUp) {
    fVals = new Double_t[fNPhi * fNR * fNZ];
    for (Int_t m=0;m < fNPhi;m++)	{
      indeksm = m * fNR * fNZ;
      mat = mvals[m];
      for (Int_t i=0;i < fNR;i++)	{
	indeks =  indeksm  + i * fNZ;
	for (Int_t j=0; j< fNZ;j++)	{
	  fVals[indeks + j] = (*mat)(i,j);
	  //printf("fVaks[%d] = %f\n",indeks + j,(*mat)(i,j));
	}	
      }
    }
    fIsAllocatingLookUp = kTRUE;		
  }
}
