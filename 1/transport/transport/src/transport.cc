#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <iostream>
#include <string>
#include <cstdlib>

#include <dune/common/parallel/mpihelper.hh>  
#include <dune/common/exceptions.hh>          
 
#include <dune/grid/onedgrid.hh>
#include <dune/common/parametertree.hh>
#include <dune/common/parametertreeparser.hh>

#include "driver_cd_simple.hh"

 

int main(int argc, char** argv)
{
    Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);


    std::string filename("1D.input");
    
     if (argc > 1) filename = argv[1];
    
    Dune::ParameterTree input_data;   
    try {
        Dune::ParameterTreeParser::readINITree(filename, input_data);
     }
     catch (...) {
           std::cerr << " Nemogu ucitati " << filename<<" ! "<<std::endl;
	       std::exit(1);
     }
    
    int N     = input_data.get<int>("N");
    double left  = input_data.get<double>("L");
    double right = input_data.get<double>("R");
    double T     = input_data.get<double>("T");
    double dt    = input_data.get<double>("dt");
    double vel   = input_data.get<double>("vel");
    double coeff_a   = input_data.get<double>("coeff_a");
    double x_right_or_left; //podatak s lijevog ili desnog ruba domene

    // Konstrukcija mreže
    using Grid = Dune::OneDGrid; 

    //prostorni korak, koristim ga za CFL ispod i eventualno radi promijene vrijednosti lijevog ruba
    double dx = ( right + abs(left) )/N ;  
 
    if(abs(vel)*dt/dx <= 1) std::cout<<abs(vel)*dt/dx<<" = CFL uvjet zadovoljen => STABILNOST !"<<std::endl;
 	 	 
    else {  
             std::cout<<"OPREZ ! "<<abs(vel)*dt/dx<<" = CFL uvjet nije zadovoljen => NESTABILNOST !"<<std::endl;
	     std::cout<<"Lijevi rub domene = "<<left<<" , desni rub domene = "<<right<<"."<<std::endl;
	     std::cout<<"Korekcija domene..."<<std::endl; 
	     
             while (abs(vel)*dt/dx > 1) dx = dx * 2; // povecavam prostorni korak, zadrzavam isti N
 	
	     //korekcija lijevog i desnog ruba, tj. prosirivanje u lijevu i desnu stranu prostora
             //imamo novi dx*N, za h idemo lijevo i desno
     	     double h = (N * dx - (right + abs(left)))/2;
             right = right + h; 
	     left = left - h;
	     
	     std::cout<<abs(vel)*dt/dx<<" = CFL uvjet zadovoljen => STABILNOST !"<<std::endl;
	     std::cout<<"Lijevi rub domene = "<<left<<" , desni rub domene = "<<right<<"."<<std::endl;
	     std::cout<<"Konstrukcija grida..."<<std::endl;
	     std::cout<<"Paraview...!"<<std::endl;
              
	      
        } 
 
    /*spremam podatak lijevog ili desnog ruba u ovisnosti o predznaku brzine koji ce biti proslijedjen u 
    funkciju boundary_left_or_right() zajedno s vremenom */
    if (vel>=0) x_right_or_left = left;
    else x_right_or_left = right; 

    Grid grid(N, left, right);
    auto gv = grid.leafGridView();  
 
  // simulacija se vrši u driver rutini
    driver(gv, T, dt, vel, coeff_a, x_right_or_left); 
    return 0;
}
