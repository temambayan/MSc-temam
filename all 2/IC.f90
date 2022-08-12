
real*8 function integra(m,n,x,y)
  integer :: m,n 
  real*8  :: x(n),y(n)
  integra = 0.0 
  do i = 1,m
     integra = integra + 0.5 * (x(i+1)-x(i)) * (y(i)+y(i+1))
  enddo

end function integra


real*8 function dnfot(eCMB)

  real*8              ::  eCMB , const
  real*8              ::  dpi    = 3.1415926535897932
  real                ::  c     =    3e10	          !cm s^-1
  real                ::  h     =    4.14e-24             !GeV*s
  real                ::  T0    =    2.728                !K
  real                ::  k     =    8.62e-14             !GeV K^-1

  const  =    ( exp(eCMB/(k*T0)) - 1.0 )**(-1.0)
  dnfot   =    8.0*dpi*const*(eCMB**2)/(h*c)**3
 ! print *, dnfot, 'CMB ----------------------'

end function dnfot


real*8 function sigmakn(ex,eCMB,gel)
  real              ::  gammae,q,g ,ex,eCMB,gel
  real                 ::  mel     =   5.11e-4                     !GeV
  real                 ::  sigmat  =   6.652e-25                   !cm^2
  gammae  =   4.0*eCMB*gel/mel
  q       =   ex/(gammae*(gel*mel-ex))

  if  ( (q .lt. 1.0/(4.0*gel**2) ) .or. (q .gt. 1.0 ) ) then 
     sigmakn =  0.0d0
  else  
     g =  2.0*q*log(q)+(1.0+2.0*q)*(1.0-q) + (gammae*q)**2*(1.0-q)/(2.0*(1.+gammae*q))
     sigmakn =  3.0d0*sigmat*g/(4.0*eCMB*gel**2)
  endif
  if ( sigmakn .lt. 0.0) sigmakn = 0.0d0 
  !print *,  gammae, ex, gel, q, eCMB, sigmakn
  !print *, sigmakn, "sigma klein-Nishina cross-section"

end function sigmakn


real*8 function ICpower(ex,gel)
  real               ::  ex, gel
  real*8             ::  eCMB
  integer            ::  l
  integer            ::  e0m     =   100
  real               ::  c     =    3e10	          !cm s^-1
  real*8,allocatable :: assee0(:)
  real*8,allocatable :: ICp(:)
  real*8,allocatable :: sum1(:)
  allocate(assee0(0:e0m))
  allocate(ICp(0:e0m))
  allocate(sum1(0:e0m))

  do l = 0, e0m
     assee0(l)   = 10.0**(-17.0 + l*10.0/real(e0m))
     eCMB        = assee0(l)
     ICp(l)      = c*ex* dnfot(assee0(l)) * sigmakn(ex, assee0(l) ,gel)
 
     !print* , dnfot(eCMB), eCMB!; '    CMB ----------------------'
     !print* , sigmakn(ex,eCMB,gel)
     !print*, ICp
  enddo
  ! print *, eCMB , 'IC power integrand'

  !;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  !;;;;; NOW WE INTEGRATE OVER CMB energies to obtain the IC power    ;;;;;;;;
  !;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  !print *, assee0, 'the CMB energy range used'
  sum1 = 0.0
  do j = 1, e0m
     sum1(j) =  sum1(j) + 0.5* (assee0(j) -assee0(j-1))*(ICp(j)+ICp(j-1))
     !  print*, j, assee0(j), ICp(j), sum1(j)
  enddo
  ICpower   =   sum(sum1)
  !print *, ICpower, "IC power"
  deallocate(assee0, ICp, sum1)

end function ICpower




program IC
  !------------------------------------------------------------------------------
  ! gfortran -O3  -fdefault-real-8 ICsubahlos.f90 -o ICsubahlos


  implicit none 

  integer            :: i,j,k,l,m,n,ok,ll, ll2, ll3
  real*8             :: dpi = 3.1415926535897932

  ! constants used for CALCULATING electron density (el_den)
  real               ::  Msun     =  1.99e33         ! in g
  real               ::  Kpc      =  3.08568025e21   ! in cm
  real               ::  YHe      =  0.25
  real               ::  mpr      =  1.67262178e-24  !proton mass in g
  real               ::  factB 



  ! constants used for CALCULATING THE EQULIBRIUM SPECTRA---------------
  real               :: SigmaV   =  3.0e-26            !cm^3s^-1  
  real               :: z        =  0.028              !redshift
  real               :: r0       =  2.82e-13           !cm
  real               :: c        =  3.0e10             !cm/s
  real               :: me       =  0.511e-3           !GeV
  real               :: dl       =  96000.     !luminosity distance in Kpc(Coma)
  real               :: ome


  real               :: rad
  real               :: dx,dx3,odx,odx3
  integer            :: nn 
  integer            :: Ngrid, nstart
  integer            :: ncl,nrock !added the halo id as in rock star

  real               ::  dia
  real               ::  fcan   ! = Msun * Kpc**(-3)  /  ompr

  real,allocatable   :: dark(:,:,:),gas(:,:,:),gne(:,:,:),gz(:,:,:)
  real,allocatable   :: bmu(:,:,:),el_den(:,:,:),nu0(:,:,:),nup(:,:,:)

  real               :: A2(3)!, E_gamma


  ! SUSY -----------------------------------------------------------------

  integer              :: lsusy,lsum,lsup, len1
  real                 :: Mchi,tmp, tmp2
  real                 :: fact = 4.68150e+21

  real*8,allocatable :: dE(:),dNdE(:)
  real*8,allocatable :: E_gamma(:),dNdE_gamma(:)
  real*8,allocatable :: ngel(:)
  real*8,allocatable :: gel(:)
  real*8,allocatable :: gel2(:)
  real*8,allocatable :: fortrpz(:)
  real*8,allocatable :: S_IC(:,:,:),F_IC(:,:)
  real*8             :: inter 


  real*8,allocatable :: integel(:) 

  real*8             :: sfac
  real*8             :: d8,bic,bsyn,bcoul,bbrem,obtot
  real*8             :: o75
  real*8             :: bfac


  real*8,external    :: integra,dnfot,sigmakn,ICpower

  real*4 :: testvar
  real*4, allocatable :: test_array(:,:,:)
  !------------------------------------------------------------------------------

  character(2)       :: mode 
  character(3)       :: cg,cn,cnn,cMchi,cll2,cll3
  character(4)       :: cgrid
  character(5)       :: cncl,crad,crock
  character(8)       :: cE_gamma
  character(256)     :: path,name,path2

  !------------------------------------------------------------------------------

  read (*,'(a)'), path
  read (*,'(a)'), mode
  read (*,*), nn
  read (*,*), rad
  read (*,*), ncl
  read (*,*), nrock
  read (*,*), Ngrid

  
  print *,('-',i = 1,70)
  print *, 'path                          :', path 
  print *, 'mode                          :', mode 
  print *, 'cluster                       :', ncl 
  print *, 'halo ID                       :', nrock
  print *, 'grid must be 21, 101 or 201   :', Ngrid


  print *,('-',i = 1,70)

  write (cg ,'(I3.3)') Ngrid
  write (cll2 ,'(I3.3)') ll2
  write (cll3 ,'(I3.3)') ll2


!!!!!!!!!!!!!!!!!!

  Mchi = 35.0
  write(cMchi,'(I3.3)') nint(Mchi) 
  name = './gamma_bb_'//cMchi//'GeV.dat'

  print *, trim(name) 


  fact    = fact / Mchi**2.0     ! fact = 0.5*SigmaV /  (1.79 * 1.0e-24 )**2.0  
  Mchi    = Mchi*1.79*1.0e-24    ! in g

  open(                           & 
       unit   = 11,               & 
       file   = trim(name),       & 
       status = 'old'             & 
       ) 

  len1 = 0 
  do 

     read(11,*,iostat = ok) tmp2
     if ( ok .lt. 0 ) exit 
     len1 = len1 + 1

  enddo


  print *,('-',i = 1,70)
  print*, 'length of the susy data for pion decay for Mx = 35 GeV giving E_gamma vs dN/dE_gamma :', len1
  print *,('-',i = 1,70)

  rewind(11) 

  allocate (E_gamma(len1)) 
 allocate (dNdE_gamma(len1)) 

   do ll = 1,len1
      read(11,*)E_gamma(ll),dNdE_gamma(ll) 
      !  print *, ll, E_gamma(ll)
   enddo
   ! print *, ' gamma ray energy for 100 ',  E_gamma(100)
   print *,' min/max gamma ray energy' , minval(E_gamma), maxval(E_gamma)
   print *,' min/max number of gamma-rays per energy range' , minval(dNdE_gamma),maxval(dNdE_gamma)
   print * 

  !### Here we read the gamma-ray energy array


!!!!!!!!!!!!!!!!!!!
  ome      =  1.0/me
  fcan     =  4.0495e-8

  ! constants used for obtaining the magnetic field (bmu) based on Bonafede2010 et al. 
  factB    =  4.7/(3.44e-3)**0.5  !! factor for the magnetic field in micro gauss based on Bonafede2010 et al. 

  write(cnn  ,'(I3.3)'),nn
  write(crad ,'(I5.5)'),nint(rad)
  write(cgrid,'(I4.4)'),Ngrid
  write(cncl ,'(I5.5)'),ncl
  write(crock ,'(I5.5)'),nrock
  name  = mode//'.NN'//cnn//'.Rad'//crad//'.Grid'//cgrid//'.d'//cncl//'.'//crock


  open (                                &  
       unit = 11,                      & 
       file = trim(path)//trim(name),  &
       form = 'unformatted',           & 
       status = 'old'                  & 
       ) 

  read(11) j, testvar
  dia = testvar
  write(*,*) j, dia
  if ( j .ne. Ngrid ) stop 
  allocate (dark(Ngrid,Ngrid,Ngrid))
  allocate (test_array(Ngrid,Ngrid,Ngrid))
  read(11) test_array
  dark = test_array


  deallocate (test_array)
  read(11) j, testvar
  dia = testvar
  write(*,*) j, dia
  if ( j .ne. Ngrid ) stop 
  allocate (gas(Ngrid,Ngrid,Ngrid))
  allocate (test_array(Ngrid,Ngrid,Ngrid))
  read(11)test_array
  gas = test_array


  deallocate (test_array)
  read(11) j, testvar
  dia = testvar
  write(*,*) j, dia
  if ( j .ne. Ngrid ) stop 
  allocate (gne(Ngrid,Ngrid,Ngrid))
  allocate (test_array(Ngrid,Ngrid,Ngrid))
  read(11) test_array
  gne = test_array


  deallocate (test_array)
  read(11) j, testvar
  dia = testvar
  write(*,*) j, dia
  if ( j .ne. Ngrid ) stop 
  allocate (gz(Ngrid,Ngrid,Ngrid))
  allocate (test_array(Ngrid,Ngrid,Ngrid))
  read(11) test_array
  gz = test_array


  dx       = dia/real(Ngrid)
  dx3      = dx**3.0
  odx      = 1.0/dx     
  odx3     = odx**3.0

  close(11) 

  print *, trim(path)//trim(name) 

  print *,('-',i = 1,70)
  print *,'dia                           :', dia            ! diameter of the cluster in kpc
  print *,'dx                            :', dx               ! size of the grid in kpc
  print *,'dark min/max                  :',minval(dark),maxval(dark), '[Msun Kpc-3]'
  print *,'gas  min/max                  :',minval(gas),maxval(gas), '[Msun Kpc-3]'
  print *,'gne  min/max                  :',minval(gne),maxval(gne)
  print *,'gz   min/max                  :',minval(gz),maxval(gz)
  print *,('-',i = 1,70)

  !--------------Magnetic field in cubes using analytical formulation
  !with B0 = 5micro Gauss --------------------

  allocate(bmu(Ngrid,Ngrid,Ngrid))
  A2 = [ -0.229052,      3.98069  ,   -13.3836]    !for B_0 =  5 microgauss 

  bmu =  ( A2(1) * (log10(dark) - 1.0) )  + A2(2)  +  A2(3) / log10(dark)
  bmu = 10**bmu

  print *,'B min/max from dm for B_0 5micro gauss  :',minval(bmu), maxval(bmu), ' [ micro gausses ]'


  !------------- using gas density in each cube to calculate the electron density
  ! (el_den)  based on Sembolini et al.-------

  allocate(el_den(Ngrid,Ngrid,Ngrid))
  el_den   =  fcan * gne * gas * (1.0 - gz - YHe )

  print *, 


  print *, 'min/max of el_den', minval(el_den), maxval(el_den), '[cm-3]'

  !------------- using gas density to calculate the Bmu(magnetic field)
  !              in each cube based on Bonafede2010 et al. -------
  ! allocate(bmu(Ngrid,Ngrid,Ngrid))
  ! bmu      =  factB * (el_den )**0.5  
  ! print *, 'magnetic field from gas min/max:',minval(bmu), maxval(bmu),' [ micro gauss ]'

  print * 


  !-----------------SUSY: CALCULATING THE EQULIBRIUM SPECTRA---------------------

  Mchi = 35.0
  write(cMchi,'(I3.3)') nint(Mchi) 
  name = './pos_bb_'//cMchi//'GeV.dat'

  print *, trim(name) 

  open(                           & 
       unit   = 11,               & 
       file   = trim(name),       & 
       status = 'old'             & 
       ) 


  lsusy = 0 
  do 

     read(11,*,iostat = ok) tmp
     if ( ok .lt. 0 ) exit 
     lsusy = lsusy + 1

  enddo

  lsum = lsusy - 1
  lsup = lsusy + 1

  print *,('-',i = 1,70)
  print*, 'length of the susy data file, lsusy     :', lsusy
  print *,('-',i = 1,70)

  rewind(11) 

  allocate (dE(lsusy)) 
  allocate (dNdE(lsusy)) 

  do i = 1,lsusy
     read(11,*)dE(i),dNdE(i) 
  enddo

  print *,' min/max electron energy' , minval(dE), maxval(dE)
  print *,' min/max number of electron per energy range' , minval(dNdE), maxval(dNdE)
  print * 


  allocate(ngel(lsusy))
  allocate(gel(lsusy))
  allocate(gel2(lsusy))
  allocate(fortrpz(lsusy)) 
  allocate(S_IC(Ngrid, Ngrid, Ngrid)) 
  allocate(F_IC(Ngrid, Ngrid)) 
  allocate(integel(lsusy)) 



  sfac     = dx3    * Kpc / (4.0* dpi * dl**2 )
  o75      = 1.0 / 75.0 
  bfac     = 1.37e-20 * (1.+z)**4

  print *,

  dark     =  fcan * mpr  * dark ! converting the dark matter density from ' [ Msun Kpc^-3 ]' to '[g cm-3]'
  print *,'dark min/max  :',minval(dark),maxval(dark)  , '[g cm-3]'

  ngel     =  0.0d0
  inter    =  0.0d0

  do ll = 201,201

     print *, 'Energy of Gamma-ray emission', E_gamma(ll) , '[GeV]'
     print *, ll, nint(1e6*E_gamma(ll)), 'indice   & Gamma-ray energy [keV]' !! CHECK THIS ALLWAYS
     write(cE_gamma,'(I8.8)')   nint(1e6*E_gamma(ll))
     path2 = 'ICsubhalos/'//'/HaloID'//crock//'Egamma'//cE_gamma//'keV'
  
     open (                               &  
       unit = 12,                      & 
       file = trim(path2),             &
       form = 'formatted',             & 
       status = 'unknown'              & 
       ) 

     do l = 1, Ngrid
        do m = 1, Ngrid
           F_IC(l,m)  = 0.0d0   !! this is to initializ
           do n = 1, Ngrid
              do k = 1, lsum
                 gel(k)  = ome  * dE(k)
                 gel2(k) = (gel(k))**2.0  
                 bic     = bfac     * gel2(k)
                 bsyn    = 1.30e-21 * ( gel(k) * bmu(l,m,n) )**2.0 
                 bcoul   = 1.2e-12  * el_den(l,m,n) * (1.0 +  ( log( gel(k) / el_den(l,m,n) ) ) * o75 )
                 bbrem   = 1.51e-16 * el_den(l,m,n) * gel(k) * (log(gel(k)) + 0.36)

                !!! bsyn    = 1.30e-21 * ( gel(k) * maxval(bmu) )**2.0 
                !!! bcoul   = 1.2e-12  * maxval(el_den) * (1.0 +  ( log( gel(k) / maxval(el_den) ) ) * o75 )
                !!! bbrem   = 1.51e-16 * maxval(el_den) * gel(k) * (log(gel(k)) + 0.36)
                 d8      = dark(l,m,n)
                 !!!d8      = maxval(dark)
                 obtot   = d8**2.0 / (me * (bic+bsyn+bcoul+bbrem) ) ! dark is the dark matter density [g cm-3]
                 fortrpz(k) = ( 0.5 * ( dE(k+1) - dE(k) ) * ( dNdE(k) +dNdE(k+1) ) ) * fact * obtot  
              enddo

              do k = 1, lsum 
                 ngel(k)= sum(fortrpz(k:lsum)) ! THE EQULIBRIUM SPECTRA
              enddo

              !!------------- calculate the integerand for local emissivity----
              do k   = 1, lsum
                 integel(k) = 2.0*ngel(k)*ICpower(E_gamma(ll),gel(k))
               ! print *,"IC power", ICpower(E_gamma(ll),gel(k)),E_gamma,gel(k)
              enddo


              inter          = integra(lsum, lsusy, gel,integel) ! local emissivity
              !!------------ calculate the emissivity-------------------------
              S_IC    = sfac*inter

              F_IC(l,m) =   F_IC(l,m) + S_IC(l,m,n)

           enddo
           write(12, *) l, m, F_IC(l,m) 
        enddo
        !print *, ICpower( E_gamma, gel(k) ), E_gamma, gel(k)
     enddo


     print *, 'Neutralino mass used in g   :', Mchi
     print *, 'DM annihilation cross-section times relative velocity of DM particles in cm3 s-1 :', SigmaV
 !    print *, 'ngel', ngel

     !  print *, 'THE electron EQULIBRIUM SPECTRA',ngel(:), '[# cm-3 Gev-1]'
   !  print *, 'Gamma-ray individual[# cm-2 s-1]', S_IC

    print *, 'Gamma-ray l.o.s [# cm-2 s-1]', minval(F_IC),maxval(F_IC) 
     print *, 'Energy of Gamma-ray emission',E_gamma , '[GeV]'
     print *, 'indices', ll

   !  print *,"IC power", ICpower 
  enddo
  close(12)

end program IC

