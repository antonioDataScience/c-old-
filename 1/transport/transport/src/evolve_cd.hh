#ifndef EVOLVE_HH
#define EVOLVE_HH

#include <cassert>

/**
 * 
 * @param gv  = LeafGridView
 * @param current = rješenje na trenutnom vremenskom sloju
 * @param next = rješenje na sljedećem vremenskom sloju
 * @param bc = vrijednost rubnog uvjeta (na sljedećem vremenskom sloju)
 * @param lambda = v*dt/dx.
 */
 

template <typename GV, typename Vec, typename VMap>
void evolve_cd(GV const & gv, Vec const & current, Vec & next, double bc, double vel, double dt, double coeff_a, VMap const & vmapper){
    /* Invarijante:
     *  -- na ulazu je next = current
     */
    int dim = GV::dimension;
    assert(dim == 1);
    int last = 0;
    
    // po svim elementima
    auto it = gv.template begin<0>();
    for (; it != gv.template end<0>(); ++it){
	
	double dx = it->geometry().volume();
	int left  = vmapper.subIndex(*it,0,dim);
        int right = vmapper.subIndex(*it,1,dim);
	double lambda = vel*dt/dx;

  	//brzina>0 , gledamo diferencije u lijevo
	if(vel >= 0){
        
         //u(i,n+1)=u(i,n)-lambda*(u(i,n)-u(i-1,n))-(dt)*a*u(i,n)
	 //npr.:
         //next[ right ] +=-lambda * (current[ right ] - current[ left ]) - dt*coeff_a * next[ right ];
	   next[ right ] +=-lambda * ( next[ right ] - current[ left ] ) - dt*coeff_a * next[ right ];
	  
 	}
	
	//brzina<0 , gledamo diferencije u desno
	else if(vel < 0){
        		
	//u(i,n+1)=u(i,n)-lambda*(u(i+1,n)-u(i,n))-(dt)*a*u(i,n)
        //npr.:
        //next[ left ] += -lambda * ( current[ right ] - current[ left ] ) - dt*coeff_a * next[ left ];
	  next[ left ] += -lambda * ( next[ right ] - current[ left ] ) - dt*coeff_a * next[ left ];
	 
	}
	
	++last;
       
    }
	 
	if(vel >= 0) next[ 0 ] = bc;  // ispravi grešku napravljenu u x = left !
	
	else next[ last ] = bc; // ispravi grešku napravljenu u x = right !
	   
}

#endif /* EVOLVE_HH */
