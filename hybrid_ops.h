/* Various Hybrid operators coded up in GRID */

///////////////////////////////////////////////////////////////
//		Staggered Sign Functions
///////////////////////////////////////////////////////////////

/*  LatticeInteger coor[4] = {&Grid, &Grid, &Grid, &Grid}; 
  LatticeInteger n[4] = {&Grid, &Grid, &Grid, &Grid, &Grid};
  LatticeComplex signs[4] =  {&Grid, &Grid, &Grid, &Grid, &Grid}; 
  LatticeComplex One(&Grid), minusOne(&Grid); One = 1; minusOne = -1;

  for(int m=0; m<4; m++) {
    LatticeCoordinate(coor[m], m);  
  }

  n[0] = coor[0]; n[1] = n[0] + coor[1]; n[2] = n[1] + coor[2]; n[3] = n[2] + coor[3];
  
  for(int i=0; i<4; i++) {
    signs[i] = where((mod(n[i],2)== (Integer) 1), minusOne, One) ;
  } 

// Decided to split the phases up a bit, here are the standard staggered sign functions, corresponding
// to the gamma matrices. the phase from the inversion of M is signs[3].  
//////////////////////////////////////////////////////////////
*/
LatticeStaggeredFermion symm_shift(LatticeGaugeField &Umu, LatticeStaggeredFermion q, int dir, int shift)
{
  GridBase *grid = Umu._grid;
  LatticeColourMatrix U(grid);
  U = peekLorentz(Umu, dir);
  LatticeStaggeredFermion tmp(grid);
  tmp = 0.5 * (U * Cshift(q, dir, shift) + adj(Cshift(U, dir, -shift))*Cshift(q, dir, -shift)) ;
  return tmp;
}

// remove the trace from a colour matrix
void make_traceless(LatticeColourMatrix &B)
{
  GridBase *grid = B._grid;
  LatticeColourMatrix con(grid) ;
  con = 1/3.0 ;

  LatticeComplex tt(grid);
  tt = trace(B);
  
  B -= tt * con;

}

//###########################################################################################################


// local operator taste gamma_i
LatticeStaggeredFermion onemm_local(LatticeGaugeField Umu, LatticeStaggeredFermion q, int dir)
{
  GridBase *grid = Umu._grid;
  LatticeColourMatrix Bi(grid);
  int i, j;
  if(dir==0) { i=2; j=1;} if(dir==1) {i=2; j=0;} if(dir==2) {i=1; j=0;}
  WilsonLoops<PeriodicGimplR>::FieldStrength(Bi, Umu, i, j); 
  //LatticeStaggeredFermion tmp(grid);
  //tmp = Bi * q;  
  return Bi * q;
}

// local operator taste gamma_i
LatticeStaggeredFermion zeromp_local(LatticeGaugeField Umu, LatticeStaggeredFermion q, LatticeComplex *signs)
{   
  GridBase *grid = Umu._grid;
  LatticeColourMatrix Bx(grid), By(grid), Bz(grid);
  WilsonLoops<PeriodicGimplR>::FieldStrength(Bx, Umu, 2, 1);
  WilsonLoops<PeriodicGimplR>::FieldStrength(By, Umu, 2, 0);
  WilsonLoops<PeriodicGimplR>::FieldStrength(Bz, Umu, 1, 0); 

  make_traceless(Bx);
  make_traceless(By);
  make_traceless(Bz);
  double milc_nrm = 8.0;
 
  return (milc_nrm * (signs[0]*Bx + signs[0]*signs[1]*By + signs[1]*signs[2]*Bz) * q);
  
}

// local operator from taste gamma_i
LatticeStaggeredFermion twomp_local(LatticeGaugeField Umu, LatticeStaggeredFermion q, LatticeComplex *signs, int dir)
{
  GridBase *grid = Umu._grid;
  //LatticeStaggeredFermion tmp(grid);
  LatticeColourMatrix Bi(grid); LatticeColourMatrix Bj(grid);
  
  int i, j;
  if(dir==0) { i=1; j=2;} if(dir==1) {i=2; j=0;} if(dir==2) {i=0; j=1;}
  // will need another conditional to ensure dir < i, dir < j.
  WilsonLoops<PeriodicGimplR>::FieldStrength(Bi, Umu, i, dir);  
  WilsonLoops<PeriodicGimplR>::FieldStrength(Bj, Umu, j, dir);
  
  LatticeStaggeredFermion tmp(grid);
  tmp = (signs[0]*signs[1]*Bi + signs[1]*signs[2]*Bj) * q;  // only works for x-dir
  return tmp;
}

// local operator taste gamma_i
LatticeStaggeredFermion onemp_local(LatticeGaugeField Umu, LatticeStaggeredFermion q, LatticeComplex *signs, int dir)
{
  GridBase *grid = Umu._grid;
  LatticeStaggeredFermion tmp(grid);
  LatticeColourMatrix Bi(grid); LatticeColourMatrix Bj(grid);
  
  int i, j;
  if(dir==2) {i=0; j=1;}
  else { cout << "ERROR! Other dirs not coded up " << endl; }
  // will need another conditional to ensure dir > i, dir > j.
  WilsonLoops<PeriodicGimplR>::FieldStrength(Bi, Umu, dir, i);  
  WilsonLoops<PeriodicGimplR>::FieldStrength(Bj, Umu, dir, j);

  make_traceless(Bi) ;
  make_traceless(Bj) ;
  double milc_nrm = 8.0 ;
  tmp = (signs[0]*Bi - signs[0]*signs[1]*Bj) * q * milc_nrm;  // only works for z-dir
  return tmp;
}

// ########################################################################################################
// ########################################################################################################
// taste singlet , shifted ops
LatticeStaggeredFermion onemp(LatticeGaugeField Umu, LatticeStaggeredFermion q, LatticeComplex *signs,
                                  int dir, int shift)
{
  GridBase *grid = Umu._grid;
  LatticeColourMatrix Bi(grid); LatticeColourMatrix Bj(grid);
  int i, j;
  if(dir==0) { i=1; j=2;} if(dir==1) {i=2; j=0;} if(dir==2) {i=0; j=1;}
  // will need another conditional to ensure dir < i, dir < j.
  WilsonLoops<PeriodicGimplR>::FieldStrength(Bi, Umu, j, dir);  
  WilsonLoops<PeriodicGimplR>::FieldStrength(Bj, Umu, i, dir);   

  LatticeStaggeredFermion tmp(grid);

  tmp = signs[i]*symm_shift(Umu, Bj*q , i, shift) + Bj*signs[i]*symm_shift(Umu, q, i, shift)
      - signs[j]*symm_shift(Umu, Bi*q, j, shift)  - Bi*signs[j]*symm_shift(Umu, q, j, shift); 
  return tmp;
}

LatticeStaggeredFermion zeromp(LatticeGaugeField Umu, LatticeStaggeredFermion q, LatticeComplex *signs, int shift)
{
  GridBase *grid = Umu._grid;
  LatticeColourMatrix Bx(grid), By(grid), Bz(grid);
  WilsonLoops<PeriodicGimplR>::FieldStrength(Bx, Umu, 2, 1);
  WilsonLoops<PeriodicGimplR>::FieldStrength(By, Umu, 2, 0);
  WilsonLoops<PeriodicGimplR>::FieldStrength(Bz, Umu, 1, 0);  

  LatticeStaggeredFermion tmp(grid);

  tmp = signs[0]*symm_shift(Umu, Bx*q , 0, shift) + Bx*signs[0]*symm_shift(Umu, q, 0, shift)
      + signs[1]*symm_shift(Umu, By*q, 1, shift)  + By*signs[1]*symm_shift(Umu, q, 1, shift) 
      + signs[2]*symm_shift(Umu, Bz*q, 2, shift)  + Bz*signs[2]*symm_shift(Umu, q, 2, shift); 
  return tmp;
}

LatticeStaggeredFermion twomp(LatticeGaugeField Umu, LatticeStaggeredFermion q, LatticeComplex *signs,
                                  int dir, int shift)
{
  GridBase *grid = Umu._grid;
  LatticeColourMatrix Bi(grid); LatticeColourMatrix Bj(grid);
  int i, j;
  if(dir==0) { i=1; j=2;} if(dir==1) {i=2; j=0;} if(dir==2) {i=0; j=1;}
  // will need another conditional to ensure dir < i, dir < j.
  WilsonLoops<PeriodicGimplR>::FieldStrength(Bi, Umu, j, dir);  
  WilsonLoops<PeriodicGimplR>::FieldStrength(Bj, Umu, i, dir);   

  LatticeStaggeredFermion tmp(grid);
  tmp = signs[i]*symm_shift(Umu, Bj*q , i, shift) + Bj*signs[i]*symm_shift(Umu, q, i, shift)
      + signs[j]*symm_shift(Umu, Bi*q, j, shift)  + Bi*signs[j]*symm_shift(Umu, q, j, shift); 
  return tmp;
}

