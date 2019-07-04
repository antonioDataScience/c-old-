#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <array>
#include <bitset>

#include <dune/common/fvector.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/alugrid/grid.hh>
#include <dune/grid/io/file/gmshreader.hh>  
#include <dune/grid/common/gridinfo.hh>  
#include <dune/grid/io/file/dgfparser/gridptr.hh>
#if HAVE_ALBERTA
#include <dune/grid/albertagrid.hh>
#endif
#if HAVE_ALUGRID
#include <dune/grid/alugrid.hh>
#endif
#if HAVE_UG
#include <dune/grid/uggrid.hh>
#endif

#include "driver_simple.hh"

//===============================================================
// Main program with grid setup
//===============================================================

int main(int argc, char** argv)
{
	//Maybe initialize Mpi
	Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);

	if (argc!=2)
	{
		if(helper.rank()==0)
			std::cout << "usage: " << argv[0] << " <level>" << std::endl;
		return 1;
	}
 
  int nref = std::stoi(argv[1]);
  typedef Dune::ALUGrid<2,2,Dune::simplex,Dune::conforming> GridType;
  typedef typename GridType::LeafGridView GridView;
  typedef typename GridType::ctype ctype;
  GridType* gridptr = Dune::GmshReader<GridType>::read("torus.msh");
  gridptr->globalRefine(nref);
  auto gv = gridptr->leafGridView();
  std::cout << " Učitavanje mreže je gotovo.\n";
  driver(gv); 
//Dune::VTKWriter<GridView> vtkwriter(gv);
//vtkwriter.write("torus");
    return 0;
}
