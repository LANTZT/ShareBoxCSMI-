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
  auto A=option(_name="A").as<double>();
  auto R=option(_name="R").as<double>();
  auto T=option(_name="T").as<double>();
  auto x0=option(_name="x0").as<double>();
  auto y0=option(_name="y0").as<double>();
  auto theta=option(_name="theta").as<double>();
  auto sigma=option(_name="sigma").as<double>();
    // $\mathbb{P}_{2}$ finite element space on
  // triangular elements
  auto Xh = Pch<2>(mesh);
  auto u = Xh->element();
  auto v = Xh->element();
  //auto ue = expr("1+x*x+2*y*y:x:y");
  //auto f = expr("-6:x:y");
    auto f=expr("4*exp(-0.5*(pow((R*x-x0)/sigma,2)) - 0.5*(pow((R*y-y0)/sigma,2))):x:y:R:T:sigma:x0:y0");
  // $\int_\Omega f v$
  auto l = form1(_test=Xh);
  l= integrate( _range=elements(mesh),_expr= f*id(v) );
    // $\int_\Omega \nabla u \cdot \nabla v$
  auto a = form2(_test=Xh,_trial=Xh);
  a = integrate( _range=elements(mesh),
                 _expr=gradt(u)*trans(grad(v)) );
  a+=on( _range=boundaryfaces(mesh),_element= u, _rhs=l, _expr=cst(0.) );
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
