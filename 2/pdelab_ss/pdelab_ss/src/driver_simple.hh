#ifndef _DRIVER_SIMPLE_HH_
#define _DRIVER_SIMPLE_HH_

#include <dune/istl/bvector.hh>
#include <dune/istl/io.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/superlu.hh>

#include <dune/pdelab/backend/istl/bcrsmatrixbackend.hh>
#include <dune/pdelab/backend/istl.hh>
#include <dune/pdelab/common/function.hh>
#include <dune/pdelab/common/vtkexport.hh>
#include <dune/pdelab/constraints/conforming.hh>
#include <dune/pdelab/finiteelementmap/pkfem.hh>
#include <dune/pdelab/gridfunctionspace/genericdatahandle.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include <dune/pdelab/gridfunctionspace/interpolate.hh>
#include <dune/pdelab/gridoperator/gridoperator.hh>
#include <dune/pdelab/stationary/linearproblem.hh>

#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

#include "bctype_simple.hh"
#include "operator_simple.hh"

#include <memory>
#include<cmath>
#include<dune/geometry/quadraturerules.hh>

//norma razlike tocnog i aproksimativnog rjesenja, to smo isto naknadno dodali
//verifikacija
 template <typename DGF, typename Exact>
double l2norma(DGF const & dgf, Exact const & exact)
{
	double norma2=0.0, norma2_1=1;

//tu dolazi račun
	auto const & gv=exact.getGridView(); 
	const int dim=Exact::Traits::GridViewType::Grid::dimension;
	
	for (auto it=gv.template begin<0>(); it !=gv.template end<0>(); ++it)
	{
	  // integracijska formula na stranici
    auto gt = it->geometry().type();
    const auto & rule = Dune::QuadratureRules<double,dim>::rule(gt,2); //za k=2 stavi 4, za k=1 stavi 2


	for (auto qit=rule.begin(); qit != rule.end(); ++qit)
		{	
			Dune::FieldVector<double,1> aprox = 0.0;		
			Dune::FieldVector<double,1> y_exact = 0.0;

			dgf.evaluate(*it, qit->position(), aprox);			
			exact.evaluate(*it, qit->position(),y_exact); //vrijednost rijesenja u integracijskoj tocki 
			aprox = aprox - y_exact;
			norma2 = norma2 + aprox.two_norm2()*qit->weight()*it->geometry().integrationElement(qit->position());
			
		}
		 
	}

return std::sqrt(norma2);
} 

template<class GV>
void driver(const GV& gv)
{
   using namespace Dune::PDELab;  // Da skratimo imena

  // <<<1>>> Tipovi u domeni i slici
  typedef typename GV::Grid::ctype Coord;
  typedef double Real;
//  const int dim = GV::dimension;

  // <<<2>>> Prostor konačnih elemenata (grid function space) i GridOperator.
  const int k = 1;
  typedef PkLocalFiniteElementMap<GV,Coord,Real,k>                FEM;
  typedef ConformingDirichletConstraints                          CONSTRAINTS;
  typedef istl::VectorBackend<>                                   VBE;
  typedef GridFunctionSpace<GV,FEM,CONSTRAINTS,VBE>               GFS;
  typedef typename GFS::template ConstraintsContainer<Real>::Type CC;
  typedef DiffusionLocalOperator<DirichletBdry>                   LOP;
  typedef istl::BCRSMatrixBackend<>                               MBE;
  typedef GridOperator<
    GFS,GFS,        /* prostor KE rješenja i test funkcije */
    LOP,            /* lokalni operator */
    MBE,            /* matrix backend */
    Real,Real,Real, /* tipovi u domeni, slici i jakobijanu */
    CC,CC           /* ograničenja za prostor rješenja i test funkcija. */
    > GO;

  FEM fem(gv);
  GFS gfs(gv,fem);
  CC cc;
  DirichletBdry bctype;
  constraints(bctype, gfs, cc); // asembliranje ograničenja Dirichletovog tipa
  std::cout << "constrained dofs=" << cc.size() << " of " << gfs.globalSize() << std::endl;
  LOP lop( bctype );
  MBE mbe(9);  // traži prosječan broj ne-nul elemenata u redu (=9)
  GO go(gfs,cc,gfs,cc,lop,mbe);

  // <<<3>>> Konstrukcija rješavača.
  typedef typename GO::Traits::Domain            U;
  typedef BCExtension<GV,Real>                   G;
  typedef ISTLBackend_SEQ_BCGS_SSOR              LS;
  typedef StationaryLinearProblemSolver<GO,LS,U> SLP;

  U u(gfs,0.0);
  G g(gv);
  interpolate(g,gfs,u);
  LS ls(5000,true);        // max 5000 iteracija, verbosity = true
  SLP slp(go,ls,u,1e-10);  // redukcija = 1e-10
  slp.apply();

  // <<<7>>> grafički izlaz (VTK)
  typedef DiscreteGridFunction<GFS,U> DGF;

  ExactSolution<GV> exact(gv);
  DGF udgf(gfs,u);
  Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,k-1);
  vtkwriter.addVertexData(std::make_shared<VTKGridFunctionAdapter<DGF>>(udgf,"solution"));
  //verifikacija
  //vtkwriter.addVertexData(std::make_shared<VTKGridFunctionAdapter<ExactSolution<GV>>>(exact,"exact")); //ODKOMENTIRATI ZA VERIFIKACIJU
  vtkwriter.write("rjesenje",Dune::VTK::ascii); //Dune::VTK::appendedraw);
	 

//verifikacija
// std::cout<< "L2 norma greske ="<< l2norma(udgf, exact)<<std::endl; //ODKOMENTIRATI ZA VERIFIKACIJU
}



#endif

