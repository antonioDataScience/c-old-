#ifndef _BCTYPE_SIMPLE_HH_
#define _BCTYPE_SIMPLE_HH_

#include <dune/common/fvector.hh>
#include <dune/pdelab/common/function.hh>
#include <dune/pdelab/constraints/common/constraintsparameters.hh>

// Klasa mora proširivati klasu PDELab::DirichletConstraintsParameters
// i u njoj prerađuje metodu isDirichlet() odlučuje je li neka točka
// na Dirichletovoj granici ili nije.
class DirichletBdry : public Dune::PDELab::DirichletConstraintsParameters
{
public:
  //  intersection = stranica elementa (u 3D) ili brid elementa (u 2D)
  //  coord        = lokalne koordinate točke na "intersectionu" koja se ispituje
  //  povratna vrijednost: true ako je točka na Dirichletovoj granici
  //                       false ako nije.
  template<typename I>
  bool isDirichlet(const I & intersection,
                   const Dune::FieldVector<typename I::ctype, I::dimension-1> & coord
                   ) const
  {
    // Globalne koordinate točke (uočite da su dimenzije lokalne i globalne točke različite )
    Dune::FieldVector<typename I::ctype, I::dimension> xg = intersection.geometry().global( coord );
     if((xg[0]*xg[0]+xg[1]*xg[1]) < 1+(1E-5))  return true;//x^2+y^2=r^2
     if((xg[0]*xg[0]+xg[1]*xg[1]) > 8.9+(1E-2))  return true;//x^2+y^2=r^2

	return false;

  }

};
/*  Klasa koja određuje vrijednost Dirichletovog rubnog uvjeta i
    njegovo proširenje na čitavu domenu.
    Template parametri:
       GV = GridView
       RF = Range Field Type (tip kojim su predstavljeni elementi slike funkcije)

       Treći parametar u GridFunctionTraits je dimenzija slike funkcije (1 jer su
       naše funkcije skalarne). Ta se dimenzija ponavlja u zadnjem parametru
       GridFunctionTraits klase.

    */
template<typename GV, typename RF>
class BCExtension
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::
                                          GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> >,
                                          BCExtension<GV,RF> > {
 // Klasa čuva referencu na GridView objekt.
  const GV& gv;
public:
  typedef Dune::PDELab::GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> > Traits;

  // Konstruktor samo uzima referencu na  GridView objekt.
  BCExtension (const GV& gv_) : gv(gv_) {}

  // Izračunaj Dirichletovu vrijednost na elementu. Ako točka nije na
  // Dirichletovoj granici, nda funkcija daje proširenje Dirichletovog rubnog
  // uvjeta na čitavu domenu. To je proširenje u osnovi proizvoljno.
  // e     = element
  // xlocal = lokalne koordinate točke u kojoj se računa Dirichletova vrijednost
  // y      = izračunata Dirichletova vrijednost
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& xlocal,
                        typename Traits::RangeType& y) const
  {
    const int dim = Traits::GridViewType::Grid::dimension;
    typedef typename Traits::GridViewType::Grid::ctype ctype;

    // Pretvori lokalne koordinate u globalne
    Dune::FieldVector<ctype,dim> x = e.geometry().global(xlocal);

	 //vezano za zadatak
	 if((x[0]*x[0]+x[1]*x[1]) < 1+(1E-5))  y=-70;
     if((x[0]*x[0]+x[1]*x[1]) > 8.9+(1E-2)) y=120;

	/*
	//dodano na zadatak, razlicite vrijednosti Dirichlea po ko kvadrantima za obje kruznice
    if((x[0]*x[0]+x[1]*x[1]) < 1+(1E-5) && x[0]>0 && x[1]>0) y =-70.0;
    if((x[0]*x[0]+x[1]*x[1]) < 1+(1E-5) && x[0]<0 && x[1]>0) y =10.0;
    if((x[0]*x[0]+x[1]*x[1]) < 1+(1E-5) && x[0]<0 && x[1]<0) y =-20.0;
    if((x[0]*x[0]+x[1]*x[1]) < 1+(1E-5) && x[0]>0 && x[1]<0) y =-70.0;
    if((x[0]*x[0]+x[1]*x[1]) > 8.9+(1E-2) && x[0]>0 && x[1]>0) y=120.0;
    if((x[0]*x[0]+x[1]*x[1]) > 8.9+(1E-2) && x[0]<0 && x[1]>0) y=70.0;
	if((x[0]*x[0]+x[1]*x[1]) > 8.9+(1E-2) && x[0]>0 && x[1]<0) y=20.0;
	if((x[0]*x[0]+x[1]*x[1]) > 8.9+(1E-2) && x[0]<0 && x[1]<0) y=0.0; */


   //verifikacija
   //if(((x[0]*x[0]+x[1]*x[1]) < 1+(1E-5)) || ((x[0]*x[0]+x[1]*x[1]) > 8.9+(1E-2)) ) y = x[0]*x[0]-x[1]*x[1] ; //ODKOMENTIRATI ZA VERIFIKACIJU

   //verifikacija
   //if(((x[0]*x[0]+x[1]*x[1]) < 1+(1E-5)) || ((x[0]*x[0]+x[1]*x[1]) > 8.9+(1E-2)) ) y = x[0]*x[0]*x[0]-x[1]*x[1]*x[1] ; //ODKOMENTIRATI ZA VERIFIKACIJU, drugi primjer


   //dodano na zadatak, Neumann na vasnjskoj kruznici
   //if((x[0]*x[0]+x[1]*x[1]) < 1+(1E-5)) y=-70;

    return;
  }

  // Vrati referencu na GridView
  inline const GV& getGridView () {return gv;}
};

//ova posebno napravljena klasa ce nam pomoci prikazivati tocno rjesenje
/* u driver datoteci moramo dodat na kraju kako je i stavljeno

ExactSolution<GV> exact(gv);
vtkwriter.addVertexData(std::make_shared<VTKGridFunctionAdapter<ExactSolution<GV>>>(exact,"exact"));*/

 template <typename GV>
class ExactSolution : public
					  Dune::PDELab::AnalyticGridFunctionBase<
						Dune::PDELab::AnalyticGridFunctionTraits< GV, double, 1>,
						 ExactSolution<GV>>
	{
public:
typedef typename Dune::PDELab::AnalyticGridFunctionTraits< GV, double, 1> Traits;
typedef typename Dune::PDELab::AnalyticGridFunctionBase<
						Dune::PDELab::AnalyticGridFunctionTraits< GV, double, 1>,
						 ExactSolution<GV>> Base; //čitava klasa nazvana Base
public:
ExactSolution(GV const & gv) : Base (gv){}

//tu bi trebali svoju formulu ubacit
void evaluate (const typename Traits::ElementType &e,
			 const typename Traits::DomainType &x,
			typename Traits::RangeType &y) const
{

const int dim=Traits::GridViewType::Grid::dimension;
 Dune::FieldVector<double,dim> xglob=e.geometry().global(x);
   y = xglob[0]*xglob[0]-xglob[1]*xglob[1];
 //y = xglob[0]*xglob[0]*xglob[0]-xglob[1]*xglob[1]*xglob[1];verifikacija drugi primjer
}
};
#endif
