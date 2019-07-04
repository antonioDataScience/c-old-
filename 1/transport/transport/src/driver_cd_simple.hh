#ifndef DRIVER_CD_SIMPLE_HH
#define DRIVER_CD_SIMPLE_HH

#include <vector>
#include <cassert>

#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include "evolve_cd.hh"

//-------------------------FUNKCIJE----------------------------//  
 double initial(double x){ 
 return cos(x);
 }


double exact(double x, double t){
 return cos(x+2*t)*exp(-t);
 }
 

double boundary_left_or_right(double t, double x_right_or_left){
 return exact(x_right_or_left,t);
}
//-------------------------------------------------------------//


template <typename GV>
void driver(GV & gv, double T, double dt, double vel , double coeff_a,double x_right_or_left)
{
  const int dim = GV::dimension;
  Dune::LeafMultipleCodimMultipleGeomTypeMapper<typename GV::Grid, Dune::MCMGElementLayout> el_mapper(gv.grid());
  Dune::LeafMultipleCodimMultipleGeomTypeMapper<typename GV::Grid, Dune::MCMGVertexLayout>  ve_mapper(gv.grid());
  
  assert( el_mapper.size() + 1 == ve_mapper.size() );
  int N = ve_mapper.size(); // broj vrhova

  //prostorni korak, vezano konkretno za L2 normu greske
  double dx = gv.template begin<0>()->geometry().volume();
   
  // current = u^n, next = u^{n+1} 
  std::vector<double> current, next, egzaktno, razlika;
  current.resize(N);
  next.resize(N);
  egzaktno.resize(N);
  razlika.resize(N);
   
 
  
  // Inicijalizacija početnim uvjetom
  int index = 0;
  for(auto it = gv.template begin<dim>(); it != gv.template end<dim>(); ++it){ 

      current[index] = initial(it->geometry().center());
      
      egzaktno[index] = initial(it->geometry().center());

      razlika[index] = 0.0 ;
      
      index++;
  }
 
   
  // Ispis početnog podatka
   Dune::VTKWriter<GV> vtkwriter(gv);
   std::string fname("aproksimirano_vs_egzaktno_vs_greska-0");
   vtkwriter.addVertexData(current, "aproksimacija");
   vtkwriter.addVertexData(egzaktno, "egzaktno");
   vtkwriter.addVertexData(razlika, "razlika");
   vtkwriter.write(fname, Dune::VTK::ascii );
  
  next = current;
  double time = 0.0;
  int k_out = 0;
  double greska = 0.0; 

   // Vremenska petlja
  while( time < T){

      // rubni uvjet na sljedećem vremenskom sloju
      double bc = boundary_left_or_right(time + dt, x_right_or_left);  
       
     // evolucija sheme u novi vremenski korak
      evolve_cd(gv, current, next, bc, vel, dt, coeff_a, ve_mapper);
    
      //tocno rjesenje   
      int i=0;
      for(auto it = gv.template begin<dim>(); it != gv.template end<dim>(); ++it) {

      	egzaktno[i]=exact(it->geometry().center(),time+dt);

      	razlika[i]=egzaktno[i]-current[i];

        greska=greska + dt*dx*(razlika[i]*razlika[i]); // sqrt(greska) racunam na kraju 
     	
	i++;
       }
      
      // priprema za sljedeću iteraciju
      time += dt;
      current = next;
      k_out++;
    
      // Ispis trenutne iteracije
      Dune::VTKWriter<GV> vtkwriter(gv);
      std::string fname("aproksimirano_vs_egzaktno_vs_greska-"); 
      fname += std::to_string(k_out);   
      vtkwriter.addVertexData(current, "aproksimacija");
      vtkwriter.addVertexData(egzaktno, "egzaktno");
      vtkwriter.addVertexData(razlika, "razlika");
      vtkwriter.write(fname, Dune::VTK::ascii );
		 
      }
      std::cout<<"Greska u L2 normi = "<<sqrt(greska) <<" ! "<<std::endl;
  
}

#endif /* DRIVER_CD_SIMPLE_HH */

