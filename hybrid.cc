<<<<<<< HEAD
    /*************************************************************************************


    Perhaps a better starting point:
    typedef typename ImprovedStaggeredFermionR::FermionField FermionField; 


    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./tests/Test_wilson_cg_unprec.cc
    Starting point:    Test_staggered_cg_unprec.cc 


    Copyright (C) 2015

Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
Author: Peter Boyle <paboyle@ph.ed.ac.uk>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

    See the full license in the file "LICENSE" in the top level distribution directory
    *************************************************************************************/
    /*  END LEGAL */
=======

/* This is a code to compute correlators
*/


>>>>>>> 3caa3e70f5bebd69be7c37565cb27800671a849f
#include <Grid/Grid.h>
#include <Grid/algorithms/iterative/BlockConjugateGradient.h>

using namespace std;
using namespace Grid;
using namespace Grid::QCD;


//////////////////////////////////////////////
// Load Fermion into propagator 
//////////////////////////////////////////////

void FermToProp_s(LatticeStaggeredPropagator &Qprop, LatticeStaggeredFermion &psi , const int c)
{
  const int nc = 3 ;  // later use the dimension 

  for(int i = 0; i < nc ; ++i)
    {
      pokeColour(Qprop, peekColour(psi, i), i, c);
    }

}

//////////////////////////////////////////////////
//  Apply anti-periodic boundary conditions in time
////////////////////////////////////////////////////

void anti_periodic( LatticeGaugeField & Umu , int nt)
{
  int mu = 3 ;  // time directiom
  GridBase *grid = Umu._grid;
  // conformable(grid,g._grid);
  LatticeColourMatrix U(grid);

  U= PeekIndex<LorentzIndex>(Umu,mu);

  // code hacked from Grid/Grid/qcd/action/fermion/FermionOperatorImpl.h`
  Lattice<iScalar<vInteger> > t(grid); LatticeCoordinate(t,3);
  LatticeComplex phases(grid); phases=1.0;

  --nt ;
  phases = where( t ==(Integer)nt, phases,-phases);
  U = phases * U ;

  PokeIndex<LorentzIndex>(Umu,U,mu);

  cout << "Anti-periodic boundary conditions in time applied" << endl ; 
}

<<<<<<< HEAD
void symm_shift(LatticeGaugeField &Umu, LatticeStaggeredFermion &q, int dir, int shift)
=======

LatticeStaggeredFermion symm_shift(LatticeGaugeField &Umu, LatticeStaggeredFermion q, int dir, int shift)
>>>>>>> 3caa3e70f5bebd69be7c37565cb27800671a849f
{
  GridBase *grid = Umu._grid;
  LatticeColourMatrix U(grid);
  U = peekLorentz(Umu, dir);
<<<<<<< HEAD
  q = 0.5 * (U * Cshift(q, dir, -shift) + adj(Cshift(U, dir, shift))*Cshift(q, dir, shift)) ;

}

=======
  LatticeStaggeredFermion tmp(grid);
  tmp = 0.5 * (U * Cshift(q, dir, shift) + adj(Cshift(U, dir, -shift))*Cshift(q, dir, -shift)) ;
  return tmp;
}

LatticeStaggeredFermion hybrid_op(LatticeGaugeField Umu, LatticeStaggeredFermion q, LatticeComplex *signs,
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
 // generalise this and/or make it more compact 
  return tmp;
}
>>>>>>> 3caa3e70f5bebd69be7c37565cb27800671a849f

int main (int argc, char ** argv)
{
  typedef typename ImprovedStaggeredFermionR::FermionField FermionField; 
  typedef typename ImprovedStaggeredFermionR::ComplexField ComplexField; 
  typename ImprovedStaggeredFermionR::ImplParams params; 


  Grid_init(&argc,&argv);

  std::vector<int> latt_size   = GridDefaultLatt();
  std::vector<int> simd_layout = GridDefaultSimd(Nd,vComplex::Nsimd());
  std::vector<int> mpi_layout  = GridDefaultMpi();
  GridCartesian               Grid(latt_size,simd_layout,mpi_layout);
  GridRedBlackCartesian     RBGrid(&Grid);

  for(int mu=0;mu<Nd;mu++) 
    {
      cout << "Lattice[ " << mu << " ]= " << latt_size[mu] << endl ;
    }


  std::vector<int> seeds({1,2,3,4});
  GridParallelRNG          pRNG(&Grid);  pRNG.SeedFixedIntegers(seeds);

  FermionField result(&Grid); result=zero;

///////////////////////////////////////////////////////////////
<<<<<<< HEAD
//		Staggered Phases
///////////////////////////////////////////////////////////////

  LatticeComplex phases[4] = {&Grid, &Grid, &Grid, &Grid}; 
  LatticeComplex one(&Grid), minusOne(&Grid); one = 1; minusOne = -1;
  LatticeInteger coor[4] = {&Grid, &Grid, &Grid, &Grid}; LatticeInteger r[4] = {&Grid, &Grid, &Grid, &Grid};
=======
//		Staggered Sign Functions
///////////////////////////////////////////////////////////////

  LatticeInteger coor[4] = {&Grid, &Grid, &Grid, &Grid}; 
  LatticeInteger n[5] = {&Grid, &Grid, &Grid, &Grid, &Grid};
  LatticeComplex signs[5] =  {&Grid, &Grid, &Grid, &Grid, &Grid}; 
  LatticeComplex One(&Grid), minusOne(&Grid); One = 1; minusOne = -1;
>>>>>>> 3caa3e70f5bebd69be7c37565cb27800671a849f

  for(int m=0; m<4; m++) {
    LatticeCoordinate(coor[m], m);  
  }
<<<<<<< HEAD
  r[0] = coor[0] + coor[1] + coor[2]; // + coor[3]?
  r[1] = r[0] - coor[0]; r[2] = r[1] - coor[1]; r[3] = coor[3];

  for(int i=0; i<4; i++) {
    phases[i] = where((mod(r[i],2)==1), minusOne, one);
  }

=======

  n[0] = 1; n[1] = coor[0]; n[2] = n[1] + coor[1]; n[3] = n[2] + coor[2];
  n[4] = n[3] + coor[3];
  
  for(int i=0; i<5; i++) {
    signs[i] = where((mod(n[i],2)== (Integer) 1), minusOne, One) ;
  } 

// Decided to split the phases up a bit, here are the standard staggered sign functions, corresponding
// to the gamma matrices. the phase from the inversion of M is signs[4].
  
>>>>>>> 3caa3e70f5bebd69be7c37565cb27800671a849f
//////////////////////////////////////////////////////////////

  LatticeGaugeField Umu(&Grid); 
  SU3::HotConfiguration(pRNG,Umu); 
  //SU3::ColdConfiguration(Umu); // Umu = 1

  int t_dir = 3;
  int nt = latt_size[t_dir];

  anti_periodic(Umu, nt);

<<<<<<< HEAD
  const int g_trans = 1;
=======
  const int g_trans = 0;
>>>>>>> 3caa3e70f5bebd69be7c37565cb27800671a849f
  if( g_trans == 1)
    {
      LatticeColourMatrix   g(&Grid); // Gauge xform
      SU3::RandomGaugeTransform(pRNG,Umu,g); // Unit gauge
      cout << "Random Gauge Transform applied "  << endl ; 
    }
  else
    {
   cout << "NO Gauge Transform applied "  << endl ; 
    }

  // This uses the Milc conventions. See gridStaggInvert.cc in the MILC code.
  // Naive staggered action
  RealD u0=1.0;
  RealD c1= 2.0 ;
  RealD c2= 0.0 ;

  // Naik coefficients (See reference in Grid library)
  // RealD c1=9.0/8.0; RealD c2=-1.0/24.0;

  RealD mass=0.1;

  ImprovedStaggeredFermionR Ds(Umu,Umu,Grid,RBGrid,2.0*mass,c1,c2,u0);

  MdagMLinearOperator<ImprovedStaggeredFermionR,FermionField> HermOp(Ds);
  ConjugateGradient<FermionField> CG(1.0e-8,10000);
  //LinearOperatorBase<FermionField> LinOp(Ds);

  //  ./Grid/qcd/QCD.h
  LatticeStaggeredFermion local_src(&Grid) ;
  LatticeStaggeredFermion out(&Grid) ;
  LatticeStaggeredPropagator Qprop[2] = {&Grid, &Grid}  ;

  Qprop[1] = Qprop[0] = zero ;
  cout << "Inversion for mass "   << mass << endl ; 
<<<<<<< HEAD

=======
  
  enum dirs {x=0, y=1, z=2, t=3};
  dirs dir = x; // just a bit more readable
>>>>>>> 3caa3e70f5bebd69be7c37565cb27800671a849f
  // Compute the staggered quark propagators
  for(int k=0; k<2; k++) {

    for(int ic = 0 ; ic < 3 ; ++ic)
      {
        cout << "---------------------------------------------" << endl ; 
        cout << "Inversion for colour " << ic << endl ; 
        cout << "---------------------------------------------" << endl ; 

        // create point source
<<<<<<< HEAD
        // tests/core/Test_staggered5Dvec.cc
      
=======
     
>>>>>>> 3caa3e70f5bebd69be7c37565cb27800671a849f
        std::vector<int> site({0,0,0,0});
        ColourVector cv = zero;
        cv()()(ic)=1.0;  
        local_src = zero;
        pokeSite(cv,local_src,site);
      
        // shift fermion field
<<<<<<< HEAD
        if(k) symm_shift(Umu, local_src, 0, 1); // do the unshifted qprop first
=======
        if(k) local_src = hybrid_op(Umu, local_src, signs, x, 1); // do the unshifted qprop first
>>>>>>> 3caa3e70f5bebd69be7c37565cb27800671a849f
//	cout << local_src;
//	break;
        Ds.Mdag(local_src, out) ; // apply Mdagger
        local_src = out;

        // invert 
        out = zero ;  // intial guess
<<<<<<< HEAD
//	if(k) symm_shift(Umu, out, 0, 1);
        CG(HermOp,local_src,out);

=======
        CG(HermOp,local_src,out);

	if(k) out = hybrid_op(Umu, out, signs, x, 1);
>>>>>>> 3caa3e70f5bebd69be7c37565cb27800671a849f
        // add solution to propagator structure
        FermToProp_s(Qprop[k], out , ic  ) ; 
      }
  }
  //  -----  Use the quark propagator to compute the correlator

<<<<<<< HEAD
  // 1-link rho correlator
  std::vector<vector<TComplex>> corr(3, vector<TComplex>(nt)) ;

  // contract the quark propagators
  LatticeComplex  c[3] = {&Grid, &Grid, &Grid};

  for(int j=0; j<3; j++) {
    c[j] = trace(adj(Qprop[0]) * Qprop[1]) ; 
    c[j] = c[j] * phases[j];	
=======
 
  std::vector<TComplex> corr(nt) ;

  // contract the quark propagators
  LatticeComplex  c(&Grid);

  for(int j=0; j<3; j++) {
    c = trace(adj(Qprop[0]) * Qprop[1]) ; 
    c = c * signs[4];	// phase from inversion of M (aka epsilon)
>>>>>>> 3caa3e70f5bebd69be7c37565cb27800671a849f
  }

  //  this correlator over the lattice is summed over the spatial
  //   lattice at each timeslice t.
  cout << "\nTp = " << Tp  << "\n"; 
<<<<<<< HEAD
  for(int k=0; k<3; k++) {
    sliceSum(c[k], corr[k], Tp);
  }

  // output the correlators
  cout << "\n\nSHIFTED rho meson \n\n";

  for(int m=0; m<3; m++) {
    cout << "\nCorrelator in spin component " << m << endl; 
    for(int tt = 0 ; tt < nt ; ++tt) {
      
        double ttt = real(corr[m][tt]) ;
        cout << tt << " "  <<  ttt  << endl ;
    }
  }

=======
  sliceSum(c, corr, Tp);
  

  // output the correlators
  cout << "\n\nHYBRID meson \n\n";

  cout << "\nCorrelator in spin component " << "x" << endl; 
  for(int tt = 0 ; tt < nt ; ++tt) {  
    double ttt = real(corr[tt]) ;
    cout << tt << " "  <<  ttt  << endl ;
  }
  
>>>>>>> 3caa3e70f5bebd69be7c37565cb27800671a849f
Grid_finalize();

}
