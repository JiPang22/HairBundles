program aa
implicit none
integer i,imax,k,n
parameter(imax=5000)
real fhb,fmax,d,del,dt,ksp,kgs,lam,lama,gam,x,xa,y,po,ax,axa,ay,m,s,t,sumi,sumr,om,amp,u1,u2,z1,z2
real xt(imax)
parameter(gam=1.3,dt=1.e-2,lam=2.8,lama=10.,d=60.9,s=0.65,m=1.,ksp=0.65,kgs=0.75,del=4.44)

open(1,file='aa') 
open(2,file='bb') 

!노이즈생성
call random_seed()  
do i = 1, imax / 2
call random_number(u1)
call random_number(u2)
z1 = sqrt(-2. * log(u1)) * cos(2.*3.14 * u2)
z2 = sqrt(-2. * log(u1)) * sin(2.*3.14 * u2)
noise(2 * i - 1) = z1
noise(2 * i) = z2
end do

!do n=1,1000
x=-101.;xa=-100.;y=1.
!fmax=(n-1.)*0.2
!fmax=200
!fmax=n
fmax=50.24

do i=1,imax

po=1./(1.+exp(10.)*exp(-(x-xa)/del))!open probability
fhb=-lam*y-kgs*(x-xa-d*po)-ksp*x!applied force
!!!!!!!!!!i!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ay=-gam*y+fhb/m+noise(i)
ax=y

axa=(kgs/lama)*(x-xa-d*po)-(fmax/lama)*(1.-s*po)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
y=y+ay*dt!velocity of x
x=x+ax*dt!position of hair bundle
xa=xa+axa*dt!position of motor protein
t=i*dt
!write(1,*) i*dt,x
write(1,*) t,x
xt(i)=x
enddo

!DFT
do k=1,40
om=6.28*k/(imax*dt)
!reset
t=0.
sumr=0.
sumi=0.
!sum
do i=1,imax
t=i*dt
sumr=sumr+xt(i)*cos(om*t)*dt
sumi=sumi-xt(i)*sin(om*t)*dt
enddo
amp=sqrt(sumi**2+sumr**2)
write(2,*) om,amp
enddo
!write(1,*) ''
!write(2,*) ''
!enddo
end
