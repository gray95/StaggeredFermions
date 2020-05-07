/******** spectrum_hybrids5.c *************/
/* MIMD version 7*/
/* DT 11/21/95 started */
/* DT 2/21/96 version2.  Make it easy to add operators.  Keep collection
   of "mult_by_1mps" routines around.  Select only operators that look
   interesting in previous tests. 
   DT 7/10/96 multiple sink operators for each source.
   New output format:  "SINK:" line lists sink operators, 
    "operator_SRC:" lines give propagators from each source.
   DT 10/4/96 add 0+- "baryon number" operator
   DT 1/21/97 use defines to select desired source operators, add
	1-+ 4 quark operator QQQQ
   DT 4/4/97  convert to KS quarks
   CD 11/09/97 added parentheses to A & B == C -> changed results
   DT 9/01 fixed many bugs
   DT 1/02 need flavor structures that are singlets under spatial rotations
   DT 1/24/02 symmetrize - average over field at quark and antiquark
   CD 8/02, moved all nonhybrid q qbar code to spectrum_nlpi2.c

   CMN started from spectrum_hybrids5.c 

*/

/* Spectrum for Kogut-Susskind hybrid mesons. 
   Any antiperiodic boundary conditions (time direction usually)
   are absorbed in the link matrices (ie phases_in=ON )

   Also does spectrum for a few conventional mesons, as reference.

    Gauge fixing should be done before calling this function
    (Coulomb gauge, probably).
*/


/* Symbolic names for propagators.  prop_SOURCE_SINK */
/* operators:
   pion:	conventional 0-+:  qbar gamma_5 q
		(1)
   pion2:	conventional 0-+:  qbar gamma_0 gamma_5 q
		(-1)^(x+y+z+t)
   0mp:		0-+ hybrid
   rho:		conventional 1--: qbar gamma_i gamma_0 q
		(1)^(dir) (VT)
   rho2:	conventional 1--: qbar gamma_i q
		(1)^(x+y+z+t+dir) (PV)
   1mm:		1-- hybrid
   a1:		conventional 1++: qbar gamma_5 gamma_i q
   a1P:		P wave a_1: qbar eps_ijk gamma_j partial_k q
   1pp:		1++ hybrid
   0pm:		0+- hybrid
   0pmP		0+- P-wave hybrid
   0pmB		0+- bilinear, qbar gamma_0 q
   0mm		0-- hybrid
   0mmP:	0-- P-wave hybrid
   1mp0:	1-+ hybrid, qbar gamma_i q eps_ijk B_j, flavor gamma_0
   1mps:	1-+ hybrid, qbar gamma_i q eps_ijk B_j, flavor 1
   1mps:	1-+ hybrid, qbar gamma_0 q E_i
   qqqq	1-+, a_1 plus pion
*/
enum prop_name { 

    prop_1mp0_1mp0,
    prop_1mps_1mps,
    prop_qqqq_1mp0,
    prop_qqqq_1mps,

    prop_2pm_2pm,

    nprops		/* nprops = number of propagators */
};

#include "generic_ks_includes.h"
#include "../include/field_strength.h"


void mult_rhoi_field( int pdir,  int r0[], su3_vector *src, su3_vector *dest );

/* Multiply by the field strength.  Arguments are two indices on
  field strength tensor, source and dest vectors */
void mult_by_field_strength( int dir1, int dir2,
    su3_vector *src, su3_vector *dest ){
    su3_vector tvec;

    register site *s;
    register int i;
    int fs_dir=-99;	/* index in field strength tensor, FS_XY to FS_ZT */
    Real sign;	/* minus one if indices in "wrong" order */

    /* figure out what fs_dir is */
    if(dir1==dir2){
	if(this_node==0)printf("Bad call to mult_by_field_strength\n");
	terminate(0);
    }
    if(dir1>dir2){ 
	sign= -1.0;
	i=dir1; dir1=dir2; dir2=i;	/* switch dirs */
    }
    else sign= 1.0;
    switch( dir1+10*dir2 ){
	case XUP+10*YUP: fs_dir = FS_XY; break;
	case XUP+10*ZUP: fs_dir = FS_XZ; break;
	case XUP+10*TUP: fs_dir = FS_XT; break;
	case YUP+10*ZUP: fs_dir = FS_YZ; break;
	case YUP+10*TUP: fs_dir = FS_YT; break;
	case ZUP+10*TUP: fs_dir = FS_ZT; break;
	default: 
            if(this_node==0)printf("Bad call to mult_by_field_strength\n");
            terminate(0); break;
    }

    /* multiply by matrix */
    FORALLSITES(i,s){
	mult_su3_mat_vec( &(s->field_strength[fs_dir]), 
	    src+i, &tvec );
	scalar_mult_su3_vector( &tvec, sign, dest + i  );
    }
}


/***
DEBUG
 ***/

#include "../include/complex.h"
#include "../include/su3.h"

/* Copy a su3 matrix:  b <- a   */
void su3mat_dump(int dir, int ii,  su3_matrix *a, su3_matrix *b ){
register int i,j;
    for(i=0;i<3;i++)
      for(j=0;j<3;j++)
	{
	  printf("dir = %d, ii = %d colour%d %d ", dir, ii, i , j) ;
	  printf("a_re - b_re = %f \n", a->e[i][j].real - b->e[i][j].real) ;
	
	  //	b->e[i][j].imag = a->e[i][j].imag;
    }
}



int init_hybrids()
{ /* return the C.G. iteration number */

  int cgn;
  register int i,j;
  register site* s;
  register complex cc;
  register int t_source;
  int dir;	/* direction in lattice */
  int color;	/* color for source */

  node0_printf("init_hybrids: Creating the field strength tensor\n");
  /* uses field_strength[] for temporary space */
  const int ape_iter = 1 ;

/* Temporary - running in old version6 */
#define APE_PROJECT

  /**
   OLD ape parameters
   NHIT 10
    SMEAR 32
   staple_weight 0.25
   ape_smear( F_OFFSET(link[0]), F_OFFSET(tmplink[0]), staple_weight, u0, 1, 3*NHIT, 0.);
   **/


#ifdef APE_PROJECT
   if( ape_iter  != 0 )
     {
       if(this_node==0)printf("init_hybrid::creating field strength tensor WITH APE smearing)\n");
       // rephase(OFF);  //  hack
       // ape_links
       for(dir=XUP;dir<=TUP;dir++)
	 FORALLSITES(i,s){
	   //	   s->tmplink[dir] = *(ape_links + i+dir*sites_on_node) ;
	   // s->tmplink[dir] = *(ape_links + i+dir*sites_on_node) ;
	   //	   su3mat_copy( (ape_links + i+dir*sites_on_node) , &(s->tmplink[dir]) ) ;
   su3mat_copy( (ape_links + 4*i+ dir) , &(s->tmplink[dir]) ) ;
   //   su3mat_dump(dir, i,   &(s->link[dir])  , &(s->tmplink[dir]) ) ;
	   //	   s->tmplink[dir] = s->link[dir]  ;
	 }
       make_field_strength(F_OFFSET(tmplink[0]), F_OFFSET(field_strength[0]));

       //rephase(ON); // hack
     }
   else
     {
       if(this_node==0)printf("init_hybrid::creating field strength tensor WITH NO APE smearing)\n");
       make_field_strength(F_OFFSET(link[0]), F_OFFSET(field_strength[0]));
     }

#else
  if(this_node==0)printf("init_hybrid::creating field strength tensor (no APE smearing)\n");
  rephase(OFF);
  make_field_strength(F_OFFSET(link[0]), F_OFFSET(field_strength[0]));
  rephase(ON);
#endif

  return 0 ;
}



/* Apply the symmetric shift opperator in direction "dir" *
 * This is the explicit version                           *
 * Covariant shifts are used                              */
void sym_shift(int dir, field_offset src,field_offset dest)
{
  register int i ;
  register site *s ;
  msg_tag *tag[2];
  su3_vector *tvec;
  
  tvec = (su3_vector *)malloc( sites_on_node*sizeof(su3_vector) );

  tag[0] = start_gather_site( src, sizeof(su3_vector), dir, EVENANDODD ,gen_pt[0] );
  FORALLSITES(i,s)
    {
      mult_adj_su3_mat_vec( &(s->link[dir]), (su3_vector *)F_PT(s,src), 
			    &(tvec[i]) ) ;
    }
  tag[1] = start_gather_field(tvec, sizeof(su3_vector), OPP_DIR(dir), 
				  EVENANDODD ,gen_pt[1] );
  wait_gather(tag[0]);
  FORALLSITES(i,s)
    {
    mult_su3_mat_vec( &(s->link[dir]), (su3_vector *)gen_pt[0][i], 
		      (su3_vector *)F_PT(s,dest) ) ;    
    }
  wait_gather(tag[1]);
  FORALLSITES(i,s)
    {
      add_su3_vector( (su3_vector *)F_PT(s,dest), (su3_vector *)gen_pt[1][i], 
		      (su3_vector *)F_PT(s,dest) ) ;    
    }
  /* Now devide by 2 eq. (4.2b) of Golderman's Meson paper*/
 FORALLSITES(i,s)
   {
     scalar_mult_su3_vector( (su3_vector *)F_PT(s,dest), .5,
			     (su3_vector *)F_PT(s,dest) );
   }
  for(i=0;i<2;i++) cleanup_gather(tag[i]) ;
  free(tvec);
}


/* Apply the symmetric shift operator in direction "dir" *
 * This is the explicit version                           *
 * Covariant shifts are used                              *
 * The KS phases MUST BE in the links                     */
static void 
sym_shift_field(int dir, su3_vector *src, su3_vector *dest)
{
  register int i ;
  register site *s ;
  msg_tag *tag[2];
  su3_vector *tvec = create_v_field();


  tag[0] = start_gather_field( src, sizeof(su3_vector), dir, EVENANDODD, gen_pt[0] );
  /* With ONE_SIDED_SHIFT defined, the shift is asymmetric */
#ifndef ONE_SIDED_SHIFT
  FORALLSITES(i,s)
    {
      mult_adj_su3_mat_vec( &(s->link[dir]), src+i, tvec+i ) ;
    }
  tag[1] = start_gather_field(tvec, sizeof(su3_vector), OPP_DIR(dir), 
			      EVENANDODD ,gen_pt[1] );
#endif
  wait_gather(tag[0]);
  FORALLSITES(i,s)
    {
      mult_su3_mat_vec( &(s->link[dir]), (su3_vector *)gen_pt[0][i], dest+i );
    }
#ifndef ONE_SIDED_SHIFT
  wait_gather(tag[1]);
  FORALLSITES(i,s)
    {
      add_su3_vector( dest+i, (su3_vector *)gen_pt[1][i], dest+i ) ;    
    }
  /* Now divide by 2 eq. (4.2b) of Golterman's Meson paper*/
  FORALLSITES(i,s)
    {
      scalar_mult_su3_vector( dest+i, .5, dest+i );
    }
  for(i=0;i<2;i++) cleanup_gather(tag[i]) ;
#else
cleanup_gather(tag[0]) ;
#endif

  
  destroy_v_field(tvec);
}



/* "Multiply by" the second quark-antiquark one link rho operator */
void mult_rho0_field_all( int fdir, su3_vector *src,  su3_vector *dest )
{
  register int i;
  register site *s;  

  /* apply the symmetric shift opperator */
  //  rephase(OFF);
  sym_shift_field(fdir, src, dest);
  //x rephase(ON);

  //    copy_v_field(dest, src); // hack-craig

#if 0 
  FORALLSITES(i,s){
    /* \eta_k already in the phases                                     * 
     * only the \epsilon for the anti-quark is needed times (-1)^(x+y+z)*
     * this is equal to (-1)^t                                          */
    if((s->t)%2==1)
      scalar_mult_su3_vector( dest + i , -1.0, dest + i );
  }
#endif

}

/***
  Multiply by epsilon = (-1)**(x+y+z+yt)

 **/

void mult_epsilon(su3_vector *src,  su3_vector *dest )
{
  register int i;
  register site *s;  

  copy_v_field(dest, src); 

  FORALLSITES(i,s){
    if((s->t + s->x + s->y + s->z)%2==1)
      scalar_mult_su3_vector( dest + i , -1.0, dest + i );
  }


}






/* "Multiply by" a one-minus-plus operator 
    quark operator is "rho0", gluon operator is magnetic field B.
    Cross product is spin one.
    Z component = "rho_X*F_XZ + rho_Y*F_YZ.rho0
   "pdir" is the polarization direction of the meson */


void mult_1mp0_field( int pdir, su3_vector *src_in, su3_vector *dest )
{
  char debug[][4] = {"XUP", "YUP", "ZUP", "TUP"  };

#if 0
  node0_printf("DEBUG tup %d %s\n", TUP, debug[TUP] );
  node0_printf("DEBUG xup %d %s\n", XUP, debug[XUP] );
  node0_printf("DEBUG yup %d %s\n", YUP, debug[YUP] );
  node0_printf("DEBUG zup %d %s\n", ZUP, debug[ZUP] );

  node0_printf("DEBUG pdir = %s\n", debug[pdir] );
#endif

    /* use cg_p and ttt as temporary storage */
    register int dir,i;
    register site *s;

    static int do_init = 1 ;

    if( do_init )
      {
	int dummy = init_hybrids() ;
	do_init = 0 ;
      }


  su3_vector *ttt = create_v_field();
  su3_vector *cg_p = create_v_field();

  su3_vector *src = create_v_field();

#if 1 
  /** copy source so src_in and dfest can be the same ***/
  copy_v_field(src, src_in);

  /* set destination to zero */
  FORALLSITES(i,s){
    clearvec( dest + i );
  }
#endif

#if 0 
  rephase(OFF);
  mult_rho0_field_all( ZUP, src, dest );  // craig hack
   rephase(ON);
#endif

  /* Loop over spatial directions orthogonal to pdir */
  for(dir=XUP; dir<=ZUP ; dir++)
    {
      if(dir != pdir){
	/* Multiply by dir,pdir component of magnetic field, 
	   multiply by flavor gamma_0 rho (includes gamma_5 for antiquark
	   propagator) */

	/**     **/
	//	node0_printf("Field_strength [%3s,%3s] rho0[%3s]\n", debug[dir],debug[pdir], debug[dir] );


	 mult_by_field_strength( dir, pdir, src, cg_p );
	 mult_rho0_field_all( dir, cg_p, ttt );

	FORALLSITES(i,s){
	  add_su3_vector( dest + i,  ttt + i , dest + i ) ;
	}

        /**    **/
	mult_rho0_field_all( dir, src, cg_p );
	 mult_by_field_strength( dir, pdir,cg_p, ttt  );


	FORALLSITES(i,s){
	  add_su3_vector( dest + i, ttt + i , dest + i ) ;
	}

      } /* end loop on dir */
    }


  destroy_v_field(ttt);
  destroy_v_field(cg_p);
  destroy_v_field(src);


} /* end mult_1mp0_field */



/* Use local rho meson operator

    "Multiply by" a one-minus-plus operator 
    quark operator is "rho0", gluon operator is magnetic field B.
    Cross product is spin one.
    Z component = "rho_X*F_XZ + rho_Y*F_YZ.rho0
   "pdir" is the polarization direction of the meson */


void mult_1mp0_lrho_field( int pdir, su3_vector *src_in, su3_vector *dest )
{
  char debug[][4] = {"XUP", "YUP", "ZUP", "TUP"  };

#if 0
  node0_printf("DEBUG tup %d %s\n", TUP, debug[TUP] );
  node0_printf("DEBUG xup %d %s\n", XUP, debug[XUP] );
  node0_printf("DEBUG yup %d %s\n", YUP, debug[YUP] );
  node0_printf("DEBUG zup %d %s\n", ZUP, debug[ZUP] );

  node0_printf("DEBUG pdir = %s\n", debug[pdir] );
#endif

    /* use cg_p and ttt as temporary storage */
    register int dir,i;
    register site *s;
    int dir_a , dir_b ; 

    static int do_init = 1 ;

    if( do_init )
      {
	int dummy = init_hybrids() ;
	do_init = 0 ;
      }

    if( pdir == ZUP )
      {
	dir_a = XUP ;
	dir_b = YUP ;
      }
    else
      {
	if(this_node==0)printf("pdir = %d not coded up \n" , pdir);
	terminate(0);	
      }


  su3_vector *ttt = create_v_field();
  su3_vector *cg_p = create_v_field();

  su3_vector *src = create_v_field();

  int r0[4] = { 0 , 0 ,0 , 0 }  ;

  /** copy source so src_in and dfest can be the same ***/
  copy_v_field(src, src_in);

  /* set destination to zero */
  FORALLSITES(i,s){
    clearvec( dest + i );
  }

#if 0
  /* Loop over spatial directions orthogonal to pdir */
  for(dir=XUP; dir<=ZUP ; dir++)
    {
      if(dir != pdir)
	{
	  /* Multiply by dir,pdir component of magnetic field, 
	     multiply by flavor gamma_0 rho (includes gamma_5 for antiquark
	     propagator) */

	 mult_by_field_strength( dir, pdir, src, cg_p );
	// Apply the local rho operator
	 mult_rhoi_field(dir,r0, cg_p, ttt ) ;

	FORALLSITES(i,s){
	  add_su3_vector( dest + i,  ttt + i , dest + i ) ;
	}

      } 

    } /* end loop on dir */
#endif

  mult_by_field_strength( dir_a, pdir, src, cg_p );
  // Apply the local rho operator
  mult_rhoi_field(dir_a, r0, cg_p, ttt ) ;
  FORALLSITES(i,s){
    add_su3_vector( dest + i,  ttt + i , dest + i ) ;
  }

  mult_by_field_strength( dir_b, pdir, src, cg_p );
  // Apply the local rho operator
  mult_rhoi_field(dir_b, r0, cg_p, ttt ) ;
  FORALLSITES(i,s){
    sub_su3_vector( dest + i,  ttt + i , dest + i ) ;
  }



  destroy_v_field(ttt);
  destroy_v_field(cg_p);
  destroy_v_field(src);


} /* end mult_1mp0_lrho_field */
