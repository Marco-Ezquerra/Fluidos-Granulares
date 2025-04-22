 program simulacion
 implicit none
 integer elem, i, j, div, ncol, bin, dime
 real*8 pi, a, wijmax, accept, mu2, mu4, mu2cum, mu4cum
 parameter (pi=4.d0*datan(1.d0), elem=131.072, a=0.8d0,  ncol=100000, bin=10, dime=5*bin) !!elem 2**n

 real*8 vx(1:elem), vy(1:elem), vz(1:elem)
 real*8 hist, hx, hy, hz, v2, v4, v2cum, v4cum, a23, a25, a240
 real*8 mu20, mu21, mu40, mu41, a2262, a2264
 dimension hist(0:dime)
 dimension hx(-dime:dime), hy(-dime:dime), hz(-dime:dime)

 open(13, file='velocidades.txt')


 !Llamamos a la subrutina del número aleatorio
      CALL init_random_seed ()
 
 do i=1, elem
   call diraleat (vx(i),vy(i),vz(i))
   end do
      call reajuste(vx,vy,vz,elem)

! Empezamos a meter colisiones         c

 wijmax=1.d-10
 accept=0.d0
 v2cum=0.d0
 v4cum=0.d0
 mu2cum=0.d0
 mu4cum=0.d0 
  
 do i=1, int(ncol/2)
        call colisiones(vx,vy,vz,elem,a,accept,wijmax)
        call reajuste(vx,vy,vz,elem)
   end do

      call inithistog (hist,hx,hy,hz,dime)

 do i=1, int(ncol/2)
        call colisiones(vx,vy,vz,elem,a,accept,wijmax)
        call reajuste(vx,vy,vz,elem)
        call histog (vx,vy,vz,elem,hist,hx,hy,hz,bin,v2cum,v4cum, dime)

        call momentos(vx,vy,vz,elem,a,mu2,mu4)
   mu2cum=mu2cum+mu2
   mu4cum=mu4cum+mu4
   end do

 v2cum=v2cum/(ncol/2*elem)
 v4cum=v4cum/(ncol/2*elem)
 mu2cum=mu2cum/(ncol/2)
 mu4cum=mu4cum/(ncol/2)
 a23=4.d0/15.d0*v4cum-1.d0
 a25=16.d0*(1.d0-a)*(1.d0-2.d0*a**2)/ (81.d0-17.d0*a-30.d0*(1.d0-a)*a**2)

 a240=16.d0*(1.d0-a)*(1.d0-2.d0*a**2)/ (97.d0-33.d0*a-2.d0*(1.d0-a)*a**2)


 mu20=dsqrt(pi*2.d0)*(1.d0-a**2)
 mu21=3.d0/16.d0*mu20
 mu40=(9.d0/2.d0+a**2)*mu20
 mu41=(3.d0/32.d0*(69.d0+10.d0*a**2)+2.d0/(1.d0-a))*mu20
 a2262=(mu2cum-mu20)/mu21
 a2264=(mu4cum-mu40)/mu41 

 do i=0, dime
   write(13,*)(dble(i)+0.5d0)/bin,hist(i)*bin/(int(ncol/2)*elem),hx(i)*bin/(int(ncol/2)*elem),&
   hy(i)*bin/(int(ncol/2)*elem),hz(i)*bin/(int(ncol/2)*elem)
   end do

 print*," v2 ",v2cum
 print*," v4 ",v4cum
 print*," a2 (Eq.3) ",a23
 print*," a2 (Eq.5) ",a25
 print*," a2 (Eq.26.2) ",a2262
 print*," a2 (Eq.26.4) ",a2264
 print*," a2 (Eq.40) ",a240
 print*, " mu2 ", mu2cum
 print*, " mu4 ", mu4cum
   
 stop
 end program
!------------------------------------------subrutina que genera direcciones aleatorias

subroutine diraleat(vx,vy,vz) 
implicit none
integer i
real*8 vx, vy, vz, rnd, r, modus
parameter (r=1.d0)
   
modus=2.d0*abs(r)
do while (modus.gt.r)
       CALL RANDOM_NUMBER (rnd) 
  vx=rnd*2.d0-1.d0
       CALL RANDOM_NUMBER (rnd) 
  vy=rnd*2.d0-1.d0
       CALL RANDOM_NUMBER (rnd)
  vz=rnd*2.d0-1.d0
    
  modus=dsqrt(vx**2+vy**2+vz**2)
  vx=vx/modus
  vy=vy/modus
  vz=vz/modus
   
  end do

return
end subroutine

!------------------------------------------subrutina que ajusta la velocidad del centro de masas y hace el reescalado de las velocidades

subroutine reajuste(vx, vy, vz, elem)
implicit none
integer i, elem
real*8 vxm, vym, vzm, v2
real*8 vx(1:elem), vy(1:elem), vz(1:elem), modu(1:elem)
   

vxm=0.d0
vym=0.d0
vzm=0.d0
v2=0.d0

do i=1, elem
  vxm=vxm+vx(i)
  vym=vym+vy(i)         
  vzm=vzm+vz(i)
  end do
   
vxm=vxm/elem
vym=vym/elem
vzm=vzm/elem

do i=1, elem
  vx(i)=vx(i)-vxm
  vy(i)=vy(i)-vym
  vz(i)=vz(i)-vzm
  modu(i)=dsqrt(vx(i)**2+vy(i)**2+vz(i)**2)
  v2=v2+modu(i)**2
  end do
  
  v2=v2/elem
      
do i=1, elem
  vx(i)=dsqrt(1.5d0/v2)*vx(i)
  vy(i)=dsqrt(1.5d0/v2)*vy(i)
  vz(i)=dsqrt(1.5d0/v2)*vz(i)
  modu(i)=dsqrt(vx(i)**2+vy(i)**2+vz(i)**2)
  end do

return
end subroutine
 
!------------------------------------------subrutina que realiza las colisiones

subroutine colisiones(vx,vy,vz,elem,alpha,accept,wijmax) 
implicit none
integer i, j, a, b, elem, k
parameter (k=1024) ! intentos de colisión
real*8 vx, vy, vz, alpha ! velocidades de las partículas antes de las colisiones y coeficiente de restitución, alpha
real*8 sxij, syij, szij, modsij, wijmax, wij, accept, rnd, r ! variables auxiliares para calcular la probabilidad de colisión
parameter (r=1.d0) 
real*8 vxij, vyij, vzij ! velocidades relativas
  
dimension vx(1:elem),vy(1:elem),vz(1:elem)
  
do i=1, k

       call diraleat (sxij, syij, szij)

       call random_number(rnd)
  a=int(rnd*elem)+1
  b=a
  do while (a.eq.b)
         call random_number(rnd)
       b=int(rnd*elem)+1
    end do
          
  vxij=vx(b)-vx(a)
  vyij=vy(b)-vy(a)
  vzij=vz(b)-vz(a)
        
   wij=vxij*sxij+vyij*syij+vzij*szij
      
   if (abs(wij).gt.wijmax) then       !
    wijmax=abs(wij)             
    end if
       
!! Ahora se elige un número aleatorio entre 0 y 1
     call random_number(rnd)

!! Y se acepta la colisión con probabilidad wij/wijmax
  if (rnd.le.(wij/wijmax))then  !esto excluye probabilidades negativas
!! Calcula las nuevas coordenadas vx(p1), vx(p2) … :
    vx(b)=vx(b)-0.5d0*(1.d0+alpha)*wij*sxij    
    vy(b)=vy(b)-0.5d0*(1.d0+alpha)*wij*syij
    vz(b)=vz(b)-0.5d0*(1.d0+alpha)*wij*szij

    vx(a)=vx(a)+0.5d0*(1.d0+alpha)*wij*sxij
    vy(a)=vy(a)+0.5d0*(1.d0+alpha)*wij*syij
    vz(a)=vz(a)+0.5d0*(1.d0+alpha)*wij*szij
  
    accept=accept+1  !contador de colisiones aceptadas
  
    end if   
  end do

return
end subroutine
   
!-------------------------------------------------------  
subroutine inithistog (hist, hx, hy, hz, dime)
implicit none
integer dime, i
real*8 hist, hx, hy, hz, v2cum, v4cum
dimension hist(0:dime)
dimension hx(-dime:dime), hy(-dime:dime), hz(-dime:dime)



!debug
v2cum=0.d0
v4cum=0.d0

do i= 0, dime
  hist(i)=0.d0 
  end do

do i= -dime, dime
  hx(i)=0.d0
  hy(i)=0.d0
  hz(i)=0.d0
  end do

return
end subroutine


! ---------------------------------------------------------------------
!                            HISTOG
!      Acumula histogramas    
! ---------------------------------------------------------------------

  subroutine histog (vx,vy,vz,npar,hist,hx,hy,hz,bin,v2cum,v4cum, dime)        

  implicit none
  integer npar, dime, bin
  real*8 vx, vy, vz, rdime

  dimension vx(npar),vy(npar),vz(npar)
  real*8 hist,hx,hy,hz,v2cum,v4cum
  dimension hist(0:dime)
  dimension hx(-dime:dime), hy(-dime:dime), hz(-dime:dime)

  integer i,k
  real*8 velx,vely,velz,vel,v,v2,v4

  rdime=dble(dime)
  do i=1,npar
    velx=vx(i)
    vely=vy(i)
    velz=vz(i)

    v2= velx*velx+vely*vely+velz*velz
    v2cum= v2cum+v2

    v4= v2*v2
    v4cum= v4cum+v4

    v= dsqrt(v2)

    k= int(v*bin)
    if (k.gt.rdime) then
      print*," fuera de rango ",k
      stop
    else
      hist(k)=hist(k)+1.d0
      endif

    k= nint(velx*bin)
    if (abs(k).gt.rdime) then
      print*," fuera de rango ",k
      stop
    else
      hx(k)=hx(k)+1.d0
      endif

    k= nint(vely*bin)
    if (abs(k).gt.rdime) then
      print*," fuera de rango ",k
      stop
    else
      hy(k)=hy(k)+1.d0
      endif


    k= nint(velz*bin)
    if (abs(k).gt.rdime) then
      print*," fuera de rango ",k
      stop
    else
      hz(k)=hz(k)+1.d0
      endif

    end do

  return
  end subroutine


!-------------------------Subrutina para calcular momentos

 subroutine momentos(vx,vy,vz,n,alpha,mu2,mu4) 
 implicit none
 integer i,j,n,a,b
 real*8 vx,vy,vz,alpha,mu2,mu4,v0,v2
 real*8 vrx,vry,vrz,vsx,vsy,vsz,pe
 real*8 v12,u12,pi,rnd
 dimension vx(1:n),vy(1:n),vz(1:n)
 parameter (pi=4.d0*datan(1.d0))
   
 v2=0.d0
   
 do i=1, n
   v2=v2+vx(i)**2+vy(i)**2+vz(i)**2
   end do
  
 v2=v2/dble(n)
  
 v0=dsqrt(2.d0*v2/3.d0)
  
 do i=1, n
       call random_number(rnd)
   a=int(rnd*n)+1
   b=a
   do while (a.eq.b)
            call random_number(rnd)
      b=int(rnd*n)+1
     end do
   
   vrx=vx(a)-vx(b)
   vry=vy(a)-vy(b)
   vrz=vz(a)-vz(b)
   vsx=0.5d0*(vx(a)+vx(b)) 
   vsy=0.5d0*(vy(a)+vy(b))
   vsz=0.5d0*(vz(a)+vz(b))
   
   v12=dsqrt(vrx**2+vry**2+vrz**2)
   u12=dsqrt(vsx**2+vsy**2+vsz**2)
   pe=vsx*vrx+vsy*vry+vsz*vrz !producto escalar
   
   mu2=mu2+pi/8.d0*(1.d0-alpha**2)*v12**3
   mu4=mu4+pi/4.d0*v12*(5.d0/3.d0*(1.d0-alpha**2)*v12**2 *u12**2+(2.d0+alpha**2)*(1-alpha**2)/12.d0*v12**4+&
   (3.d0-alpha)*(1.d0+alpha)*(pe**2-1.d0/3.d0*v12**2*u12**2)) 
   end do
   
 mu2=mu2/(n*v0**3)
 mu4=mu4/(n*v0**5)

 return
 end subroutine
  
  
  !!------------RANDOM NUMBER-----------------------------------------------!!

subroutine init_random_seed()
    implicit none
    integer, allocatable :: seed(:)
    integer :: i, n, istat, dt(8), pid, t(2), s
    integer(8) :: count, tms

    ! Inicializar el tamaño de la semilla
    call random_seed(size = n)
    allocate(seed(n))

    ! Usar el reloj del sistema si no hay acceso a /dev/urandom
    call system_clock(count)
    if (count /= 0) then
        t = transfer(count, t)
    else
        call date_and_time(values=dt)
        tms = (dt(1) - 1970) * 365_8 * 24 * 60 * 60 * 1000 &
            + dt(2) * 31_8 * 24 * 60 * 60 * 1000 &
            + dt(3) * 24 * 60 * 60 * 60 * 1000 &
            + dt(5) * 60 * 60 * 1000 &
            + dt(6) * 60 * 1000 + dt(7) * 1000 &
            + dt(8)
        t = transfer(tms, t)
    end if

    ! Calcular una semilla usando el tiempo y el PID
    s = ieor(t(1), t(2))
    pid = 1099279  ! Puedes usar una constante si no puedes obtener el PID
    s = ieor(s, pid)

    ! Asignar valores a la semilla
    if (n >= 3) then
        seed(1) = t(1) + 36269
        seed(2) = t(2) + 72551
        seed(3) = pid
        if (n > 3) then
            seed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
        end if
    else
        seed = s + 37 * (/ (i, i = 0, n - 1 ) /)
    end if

    ! Inicializar el generador de números aleatorios con la semilla generada
    call random_seed(put=seed)

end subroutine init_random_seed

  
  
  
  
