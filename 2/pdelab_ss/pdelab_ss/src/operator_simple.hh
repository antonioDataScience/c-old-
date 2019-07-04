#ifndef _OPERATOR_HH_
#define _OPERATOR_HH_

#include<dune/geometry/quadraturerules.hh>
#include<dune/geometry/referenceelements.hh>

#include<dune/pdelab/localoperator/defaultimp.hh>
#include<dune/pdelab/localoperator/flags.hh>
#include<dune/pdelab/localoperator/pattern.hh>

template<class RF, int dim>
double k1(Dune::FieldVector<RF,dim> x)  // k(x)
{	
 /* if (x[0]>=0.0 && x[1]>0.0)  return 0.1; 
 	if (x[0]<0.0 && x[1]>0.0)  return 1.0; 
	if (x[0]<0.0 && x[1]<=0.0)  return 0.01; 
    return 10; */

 // std::cout << " return 1\n";   
    return 1;
}

template<class BCType>
class DiffusionLocalOperator : // derivacijska lista -- jakobijan i pattern računa PDELab
  public Dune::PDELab::NumericalJacobianApplyVolume  <DiffusionLocalOperator<BCType> >,
  public Dune::PDELab::NumericalJacobianVolume       <DiffusionLocalOperator<BCType> >,
  public Dune::PDELab::NumericalJacobianApplyBoundary<DiffusionLocalOperator<BCType> >,
  public Dune::PDELab::NumericalJacobianBoundary     <DiffusionLocalOperator<BCType> >,
  public Dune::PDELab::FullVolumePattern,
  public Dune::PDELab::LocalOperatorDefaultFlags
{
public:
  // Zastavice koje signaliziraju da na svakom elementu treba zvati: 
  enum { doPatternVolume = true };  // metodu za računanje patterna (iz volumnih doprinosa)
  enum { doAlphaVolume = true };    // alpha_volume
  enum { doAlphaBoundary = true };  // alpha_boundary         

  DiffusionLocalOperator(const BCType& bctype_, // boundary cond.type
                         unsigned int intorder_=2) :
    bctype( bctype_ ), intorder( intorder_ )
  {}

  // Računanje volumnog integrala
  // eg   = element (geometry)
  // lfsu = lokalni prostor funkcija za rješenje
  // lfsv = lokalni prostor funkcija za test funkciju
  // x    = vektor koeficijenata rješenja 
  // r    = lokalni rezidual
  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
  {
    // dimenzije
    const int dim = EG::Geometry::dimension;
    const int dimw = EG::Geometry::dimensionworld;

    // tipovi
    typedef typename LFSU::Traits::FiniteElementType::Traits::LocalBasisType::Traits::DomainFieldType DF;
    typedef typename LFSU::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType RF;
    typedef typename LFSU::Traits::FiniteElementType::Traits::LocalBasisType::Traits::JacobianType Jacobian;
    typedef typename LFSU::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeType Range;
    typedef Dune::FieldVector<RF,dimw> Gradient;
    typedef typename LFSU::Traits::SizeType size_type;

    // integracijska formula
    auto gt = eg.geometry().type();
    auto& rule = Dune::QuadratureRules<DF,dim>::rule(gt,intorder);

    // petlja po svim integracijskim točkama
    for (auto it=rule.begin(); it!=rule.end(); ++it)
      {
        // računanje baznih funckcija na referentnom elementu
        std::vector<Range> phi(lfsu.size());
        lfsu.finiteElement().localBasis().evaluateFunction(it->position(),phi);

        // rješenje u integracijskoj točki
        RF u=0.0;
        for (size_type i=0; i<lfsu.size(); ++i) u += x(lfsu,i)*phi[i];

        // gradijent baznih funkcija
        std::vector<Jacobian> js(lfsu.size());
        lfsu.finiteElement().localBasis().evaluateJacobian(it->position(),js);

        // transformacija gradijenata s referentnog na fizički element
        const Dune::FieldMatrix<DF,dimw,dim> &jac = eg.geometry().jacobianInverseTransposed(it->position());
        std::vector<Gradient> gradphi(lfsu.size()); //vrijednost rijesenja u integracijskoj tocki 
        for (size_type i=0; i<lfsu.size(); i++)
          jac.mv(js[i][0],gradphi[i]);

        //gradijent rješenja u integracijskoj točki
        Gradient gradu(0.0);
        for (size_type i=0; i<lfsu.size(); ++i)
          gradu.axpy(x(lfsu,i),gradphi[i]);
		 
	   
     	Dune::FieldVector<RF,dim> x_1 = eg.geometry().global(it->position());
        RF kk = k1(x_1);
	    RF f = 0.0; //x^2-y^2, za verifikaciju
        //RF f= -6*(x_1[0]-x_1[1]); // -x^3+y^3 za verifikaciju drugi primjer
	  
        // integriramo : k * grad u * grad phi_i - f phi_i
        RF factor = it->weight()*eg.geometry().integrationElement(it->position());

        for (size_type i=0; i<lfsu.size(); ++i)
          r.accumulate(lfsu, i, (k1(x_1)*(gradu*gradphi[i]) - f*phi[i]) * factor);
      }
  }

  // integral po rubu
  // ig     = intersection 😊 stranica elementa)
  // lfsu_s = lokalni prostor funkcija na stranici za rješenje
  // lfsu_v = lokalni prostor funkcija na stranici za test funkciju 
  // x_s    = vektor koeficijenata rješenja (na stranici)
  // r_s    = rezidual (na stranici)
  template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_boundary (const IG& ig, const LFSU& lfsu_s, const X& x_s,
                       const LFSV& lfsv_s, R& r_s) const
  {
    // tipovi
    typedef typename LFSU::Traits::FiniteElementType::Traits::LocalBasisType::Traits::DomainFieldType DF;
    typedef typename LFSU::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType RF;
    typedef typename LFSU::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeType Range;
    typedef typename LFSU::Traits::SizeType size_type;

    // dimenzije
    const int dim = IG::dimension;

    // integracijska formula na stranici //vrijednost rijesenja u integracijskoj tocki 
    auto gtface = ig.geometryInInside().type();
    const auto& rule = Dune::QuadratureRules<DF,dim-1>::rule(gtface,intorder);

    // petlja po svim integracijskim točkama
    for (auto it=rule.begin(); it!=rule.end(); ++it)
      {
        // Ako smo na Dirichletovoj granici preskačemo petlju
        if ( bctype.isDirichlet( ig, it->position() ) )
          continue;

        // pozicija int. točke u lokalnim koordinatam elementa
        auto local = ig.geometryInInside().global(it->position());

        // izračunaj bazne funkcije u integracijskoj tički
        std::vector<Range> phi(lfsu_s.size());
        lfsu_s.finiteElement().localBasis().evaluateFunction(local,phi);

        // rješenje u integracijskoj točki
        RF u=0.0;
        for (size_type i=0; i<lfsu_s.size(); ++i)
          u += x_s(lfsu_s,i)*phi[i];
		
		/* 
        //dodano na zadatak, računanje Neumannovog rubnog uvjeta
        Dune::FieldVector<RF,dim> globalpos = ig.geometry().global(it->position());
        RF j = 120; 

        // integracija
        RF factor = it->weight()*ig.geometry().integrationElement(it->position());

        for (size_type i=0; i<lfsu_s.size(); ++i)
          r_s.accumulate(lfsu_s,i, j*phi[i]*factor); */
      }
  }

private:
  const BCType& bctype;
  unsigned int intorder;
};
#endif
