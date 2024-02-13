program aa
implicit none
integer i,imax,k,n,kmax,nmax
parameter(imax=5000,kmax=100,nmax=530)
real fhb,d,del,dt,ksp,kgs,lam,lama,gam,x,xa,y,po,ax,axa,ay,m,s,t,sumi,sumr,om,amp
real xt(imax),fmax,ampom(kmax),fmaxmax(nmax)
parameter(gam=1.3,dt=1.e-2,lam=2.8,lama=10.,d=60.9,s=0.65,m=1.)
parameter(ksp=0.65,kgs=0.75,del=4.53)

open(1,file='aa') 
open(2,file='bb') 

do n=1,nmax
x=-120.;xa=-121.;y=1.
fmax=n*0.1

do i=1,imax
po=1./(1.+exp(41.4)*exp(-(x-xa)/del))            !open probability
!fhb=-lam*y-kgs*(x-xa-d*po)-ksp*x!applied force
!ay=-gam*y+fhb/m
ax=-(kgs/lam)*(x-xa-d*po)-(ksp/lam)*x
axa=(kgs/lama)*(x-xa-d*po)-(fmax/lama)*(1.-s*po)

y=y+ay*dt                                        !velocity of x
x=x+ax*dt                                        !position of hair bundle
xa=xa+axa*dt                                     !position of motor protein
t=i*dt
!write(1,*) t,x
xt(i)=x
enddo

!DFT
do k=1,kmax
om=6.28*k/(imax*dt)

t=0.
sumr=0.
sumi=0.

do i=1,imax
t=i*dt
sumr=sumr+xt(i)*cos(om*t)*dt
sumi=sumi-xt(i)*sin(om*t)*dt
enddo
amp=sqrt(sumi**2+sumr**2)
ampom(k)=amp
write(2,*) om,amp
enddo
fmaxmax=maxval(ampom)
write(1,*) fmax,fmaxmax
write(2,*) ''
enddo
end
