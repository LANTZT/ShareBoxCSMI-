program chaleur_transfert
	implicit none
	integer, parameter :: nxmax=50,ntmax=50, cond=1
	real*8 :: bas(2:nxmax),haut(1:nxmax-1),diag(1:nxmax)
	integer :: i,j,n=1

	real*8 :: uold(nxmax) ! stockage de u_i^n 
	real*8 :: unew(nxmax) !stockage de u_i^(n+1)

	real*8 :: pold(nxmax) ! stockage de p_i^n 
	real*8 :: pnew(nxmax) !stockage de p_i^(n+1)

	real*8 :: xmil(nxmax)
	!gradJ : gradient de J, px : les valeurs de droite 
	real*8 :: gradJ,g(1:ntmax),px(1:ntmax), u_initiale=0.d0, u_desiree=10

	real*8 :: a,b,dx,dt,t,tmax,alpha,beta,eps=0.005d0, rho=1.0d0,tt
    
t=0.d0
tmax=1.0d0!2.d0!0.0001d0!
a=0
b=1
dx=(b-a)/nxmax
dt=tmax/ntmax

!ntmax=tmax/dt

alpha=0
beta=0

if(cond.eq.0) then
    g=0.0d0
    else
! récupération de la dernière donnée du flux
open (103, file='g')
do i=1,ntmax
	read(103,*) tt,g(i)
    end do
close(103)
end if

	do i=1,nxmax
	xmil(i)=a+i*dx-dx/2 
	call uini(xmil(i),u_initiale ,uold(i))
end do

diag=2
diag(1)=1
diag(nxmax)=1
diag=(dt/dx**2)*diag+1
bas=-dt/dx**2
haut=-dt/dx**2

!grande boucle : (en théorie) critere d'arret lorsque g(n) est stable
do j=1,1000

!1ere resolution en temps 
uold=0
do n=1,ntmax
	t=t+dt
	uold(1)=uold(1)+(dt/dx)*g(n)
	call factolu(nxmax,bas,diag,haut,uold,unew)
	uold=unew
!	write(*,*) 't= ', t, ' u',uold
	px(n)=uold(nxmax)
end do

pold=0

!2eme boucle
do n=ntmax,1,-1
	pold(nxmax)=pold(nxmax)+ (dx/dt)*(px(n)-u_desiree)
	call factolu(nxmax,bas,diag,haut,pold,pnew)
	pold=pnew
	g(n)=g(n) - rho*(pold(1) + eps*g(n))
end do
    end do

!traces
 open(101,file='num_transfert')
do i=1,nxmax
	write(101,*) xmil(i),unew(i)
end do
close(101)
open (102, file='u_droit')
    write(102,*) 0,0
do i=1,ntmax
	write(102,*) i*dt,px(i)
end do
close(102)
open (103, file='g')
do i=1,ntmax
	write(103,*) i*dt,g(i)
end do
close(103)	

call system('gnuplot plotcom')

end program chaleur_transfert

! -------------------------------------------------------------INI
subroutine uini(x,u_initiale,u)
	implicit none
    real*8 :: x, u, u_initiale
	u=u_initiale
end subroutine uini

! --------------------------------------------------------------LU
subroutine factolu(n,L,d,u,f,x)
  implicit none
  real*8 :: L(2:n),d(1:n),u(1:n-1),f(1:n),x(1:n)
  integer :: i,n
  logical,save :: facto = .true.


  if (facto) then
     facto = .false.
     write(*,*) 'factorisation'
     ! factorisation
     do i=2,n
        L(i)=L(i)/d(i-1)
        d(i)=d(i)-L(i)*u(i-1)
     end do
  end if

  ! descente
  x(1)=f(1)
  do i=2,n
     x(i)=f(i)-L(i)*x(i-1)
  end do

  ! remontÃ©e
  x(n)=x(n)/d(n)
  do i=n-1,1,-1
     x(i)=(x(i)-u(i)*x(i+1))/d(i)
  end do

end subroutine factolu

