#include <feel/feel.hpp>
int main( int argc, char** argv ) {
  using namespace Feel;
  po::options_description options ( "Membrane options ") ;
  options.add_options()
  ( "T", po::value<double>()->default_value(1.0), "Tension" )
  ( "x0", po::value<double>()->default_value(0.0), "x-coordinate" )
  ( "y0", po::value<double>()->default_value(0.0), "y-coordinate" )
  ( "sigma", po::value<double>()->default_value(0.02), "sigma" )
  ( "R", po::value<double>()->default_value(0.3), "Radius" )
  ( "theta", po::value<double>()->default_value(0), "Angle" )
  ( "A", po::value<double>()->default_value(1.0), "Amplitude" );
  Environment env( _argc=argc, _argv=argv,
                   _desc=options,
                   _about=about(_name="membrane",
                                _author="Feel++ Consortium",
                                _email="feelpp-devel@feelpp.org"));
  auto mesh = unitCircle();
  auto A=Environment::vm(_name="A").template as<double>();
  auto R=Environment::vm(_name="R").template as<double>();
  auto T=Environment::vm(_name="T").template as<double>();
  auto x0=Environment::vm(_name="x0").template as<double>();
  auto y0=Environment::vm(_name="y0").template as<double>();
  // $\mathbb{P}_{2}$ finite element space on
  // triangular elements
  auto Xh = Pch<2>(mesh);
  auto u = Xh->element();
  auto v = Xh->element();
  auto ue = expr("1+x*x+2*y*y:x:y");
  auto f = expr("-6:x:y");
  // $\int_\Omega f v$
  auto l = form1(_test=Xh);
  l= integrate( elements(mesh), f*id(v) );
  // $\int_\Omega \nabla u \cdot \nabla v$
  auto a = form2(_test=Xh,_trial=Xh);
  a = integrate( _range=elements(mesh),
                 _expr=gradt(u)*trans(grad(v)) );
  a+=on( boundaryfaces(mesh), u, l, ue );
  // solve a( u, v ) = l( v )
  a.solve( _solution=u, _rhs=l );

  std::cout << "|u-ue|= "
            << normL2(_range=elements(mesh),
                      _expr=(idv(u)-ue) )
            << "\n";

  auto e = exporter( _mesh=mesh ); 
  e->add( "u", u );
  // v interpolate ue
  v.on(_range=elements(mesh),_expr=ue );
  e->add( "ue", v );
  e->save();
}
