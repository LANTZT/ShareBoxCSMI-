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
                                _author="CSMI"
                               ));
  auto mesh = unitCircle();
  // $\mathbb{P}_{2}$ finite element space on
  // triangular elements
  auto Xh = Pch<2>(mesh);
  auto u = Xh->element();
  auto v = Xh->element();
//  auto ue = expr("1+x*x+2*y*y:x:y");
//  auto f = expr("-6:x:y");
  auto f=expr("4*exp(-0.5*(pow((R*x-x0)/sigma,2)) - 0.5*(pow((R*y-y0)/sigma,2))):x:y:R:T:sigma:x0:y0"); 

  double T = option(_name="T").as<double>();
  double x0 = option(_name="x0").as<double>(); 
  double y0  = option(_name="y0").as<double>();
  double sigma = option(_name="sigma").as<double>();
  double R = option(_name="R").as<double>();
  double theta = option(_name="theta").as<double>();
  double A = option(_name="A").as<double>();

 auto F=vf::project( _space=Xh, _range=elements(mesh), _expr=4*exp(-0.5*(pow((R*Px()-x0)/sigma,2)) - 0.5*(pow((R*Py()-y0)/sigma,2))));


  // $\int_\Omega \nabla u \cdot \nabla v$
  auto a = form2(_test=Xh,_trial=Xh);
  a = integrate( _range=elements(mesh),
                 _expr=gradt(u)*trans(grad(v)) );
  
  // $\int_\Omega f v$
  auto l = form1(_test=Xh);
  l= integrate(_range=elements(mesh),
               _expr= f*id(v) );

//  a+=on( boundaryfaces(mesh), u, l, ue );
  
  a+=on(_range=boundaryfaces(mesh), _rhs=l, _element=u, _expr=cst(0.));

  // solve a( u, v ) = l( v )
  a.solve( _solution=u, _rhs=l );

//  std::cout << "|u-ue|= "
//            << normL2(_range=elements(mesh),
//                      _expr=(idv(u)-ue) )
//            << "\n";
  

  auto e = exporter( _mesh=mesh ); 
  e->add( "u", u );
  e->add("f", F);
  // v interpolate ue
 // v.on(_range=elements(mesh),_expr=ue );
 // e->add( "ue", v );
  e->save();
}
