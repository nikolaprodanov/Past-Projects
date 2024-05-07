module sub

contains
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!$
!!$  subroutine xmask(xarray,fmodulation)
!!$
!!$    ! *********************************
!!$    !
!!$    ! Genera una maschera per modulare il dato iniziale
!!$    !
!!$    ! *********************************
!!$
!!$    use nrtype
!!$    implicit none
!!$
!!$
!!$
!!$
!!$  end subroutine xmask


  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine initialXphon(fileId,xphon)

    ! *********************************
    !
    ! Carica il dato iniziale da file
    !
    ! *********************************

    use nrtype
    implicit none

    integer(i4b), intent(in) :: fileId
    real(dp), dimension(:),intent(out) :: xphon

    integer(i4b) :: ix,iselfcon, N


    N=size(xphon)


    ! Read old values of phonon to set the initial conditions for the self-consistent cicle
    ! debug

    do ix=1,N
       read(fileId,'(I4,E15.7)') iselfcon, xphon(ix)
       write(876,'(I4,E15.7)') iselfcon,xphon(ix)
    end do


  end subroutine initialXphon

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  subroutine updateXphon(damp,derr,diffsqr,xphonOld,xphonNew,xphon)

    use nrtype

    implicit none

    real(dp),intent(in) :: derr,damp,diffsqr
    real(dp),dimension(:),intent(in) :: xphonOld,xphonNew
    real(dp),dimension(:),intent(out) :: xphon

    real(dp) :: dampself

    dampself = damp*exp(-derr/diffsqr)

    xphon= dampself*xphonNew + (1.0d0-dampself)*xphonOld

  end subroutine updateXphon

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  subroutine print_section(iselfcon,pos,array)

    use nrtype
    implicit none

    integer(I4B),intent(IN) :: iselfcon,pos
    real(dp),dimension(:),intent(IN) :: array

    write(129,*) iselfcon,array(pos)
    flush(129)



  end subroutine print_section


  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine dgap(eig,dg,deig)

    use nrtype
    implicit none

    real(dp),dimension(:),intent(IN) :: eig
    real(dp),intent(OUT) :: dg,deig

    integer(I4B) :: icol,ncol,i,N


    i=1
    do
       if(eig(i) .gt. 0.0d0) EXIT
       i=i+1
    end do
    ncol = i-1

    dg=eig(ncol+1)-eig(ncol)
    deig=eig(ncol)-eig(ncol-1)


  end subroutine dgap


  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  subroutine check_solution(pos,xphonOld,xphonNew)

    use nrtype
    implicit none

    integer(I4B),intent(IN) :: pos
    real(dp),dimension(:),intent(IN) :: xphonOld,xphonNew

    if (xphonNew(pos)*xphonOld(pos) .gt. 0) return

    write(17,*) 'Cambio segno'



  end subroutine check_solution


  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine test_func(N,g_ep,K_el,t)

    use nrtype

    implicit none

    integer(I4B),intent(IN) :: N
    REAL(DP),intent(IN) :: g_ep,K_el,t

    REAL(DP) :: sum
    integer(I4B) :: i,numK,ix,nx
    REAL(DP) :: k,dx,x

    numK = int(dble(N)/2.0d0)

    nx = 100
    dx = 0.1d-9

    do ix=1,nx

       sum = 0.0d0

       x=dble(ix)*dx
       do i=0,numK-1

          k = 4.0d0*PI_D*dble(i)/dble(N)
          sum = sum + 1.d0/sqrt(g_ep**2*x**2 + 4.0d0*t**2 * cos(k/2.0d0)**2)
          !       write(17,*) k,func,cos(dble(k))
       end do
       sum = sum/dble(N) - K_el/g_ep**2
       write(286,*) x,sum
    end do





  end subroutine test_func

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!!$  subroutine polaron_local(eig,WN,xphon,localpol)
!!$
!!$    use nrtype
!!$    implicit none
!!$
!!$    real(dp),dimension(:),intent(IN) :: eig
!!$    real(dp),dimension(:,:),intent(IN) :: WN
!!$    real(dp),dimension(:),intent(IN) :: xphon
!!$    real(dp),dimension(:),intent(OUT) :: localpol
!!$
!!$    integer(I4B) :: icol,ncol,i,j,N
!!$    real(dp) :: avrgN
!!$
!!$    N=size(xphon)
!!$
!!$    i=1
!!$    do       
!!$       if(eig(i) .gt. 0.0d0) EXIT
!!$       i=i+1
!!$    end do
!!$    ncol = i-1
!!$
!!$    !Test
!!$    write(17,*) 'ncol_polaron_local=',ncol
!!$
!!$    charge = 0.0d0
!!$    do i=1,N
!!$       do icol=1,ncol
!!$          charge(i) = charge(i)+0.5d0*(WN(2*i-1,icol)**2 - WN(2*i,icol)**2)
!!$       end do
!!$       charge(i) = charge(i) + 0.5d0
!!$    end do
!!$
!!$    localpol = 0.0d0
!!$    localpol = charge*xphon
!!$
!!$  end subroutine polaron_local

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function bulkZeroEnergy(N,g_ep,K_el,t,delta)

    use nrtype
    implicit none

    integer(I4B),intent(IN) :: N
    REAL(DP), INTENT(IN) :: g_ep,K_el,t,delta
    REAL(DP) :: bulkZeroEnergy

    integer(I4B) :: i,numK
    REAL(DP) :: k,rad_k,delta_k,csi_k


    numK = int(dble(N)/2.0d0)
    bulkZeroEnergy = 0.0d0

    do i=0,numK-1

       k = 4.0d0*PI_D*dble(i)/dble(N)
       delta_k = 2.0d0*delta*sin(k/2.0d0)
       csi_k = 2.0d0*t*cos(k/2.0d0)
       rad_k=sqrt(delta_k**2 + csi_k**2)

       bulkZeroEnergy = bulkZeroEnergy - rad_k
    end do

    bulkZeroEnergy = bulkZeroEnergy/dble(N)



  end function bulkZeroEnergy
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function bulkEnergy(N,g_ep,K_el,t,dimer,delta)

    use nrtype
    implicit none

    integer(I4B),intent(IN) :: N
    REAL(DP), INTENT(IN) :: g_ep,K_el,t,dimer,delta
    REAL(DP) :: bulkEnergy

    integer(I4B) :: i,numK
    REAL(DP) :: k


    numK = int(dble(N)/2.0d0)
    bulkEnergy = 0.0d0

    do i=0,numK-1

       k = 4.0d0*PI_D*dble(i)/dble(N)
       bulkEnergy = bulkEnergy - 0.5d0*(sqrt((g_ep*dimer + 2.0d0*delta*sin(k/2.0d0))**2 + 4.0d0*t**2 * cos(k/2.0d0)**2) + & 
            & sqrt((g_ep*dimer - 2.0d0*delta*sin(k/2.0d0))**2 + 4.0d0*t**2 * cos(k/2.0d0)**2))
    end do

    bulkEnergy = bulkEnergy/dble(N) + 0.5d0*K_el*dimer**2

!!$    do k=1,numK
!!$       energy = energy - sqrt(g_ep**2*dimer**2 + 4.0d0*t**2 * cos(dble(k)/2.0d0)**2)
!!$    end do
!!$
!!$    energy = energy/dble(N)

  end function bulkEnergy



  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine bulkDimer(rtbis2,N,g_ep,K_el,t,delta,E0,dimer,metadimer,dflag)

    use nrtype
    use param
    implicit none

    integer(I4B),intent(IN) :: N
    REAL(DP), INTENT(IN) :: g_ep,K_el,t,delta,E0
    REAL(DP), INTENT(OUT) :: dimer,metadimer
    logical, intent(OUT) :: dflag

    REAL(DP) :: xacc,x1,x2,x3,x,dimer1,dimer2,E1,E2,Ed
    integer(I4B) :: i,nstep


!!$       FUNCTION func(x,N,g_ep,K_el,t)
!!$         USE nrtype
!!$         IMPLICIT NONE
!!$         REAL(DP), INTENT(IN) :: x
!!$         integer(I4B),intent(IN) :: N
!!$         REAL(DP),intent(IN) :: g_ep,K_el,t
!!$         REAL(DP) :: func
!!$       END FUNCTION func
    INTERFACE
       FUNCTION rtbis2(x1,x2,N,g_ep,K_el,t,delta,xacc)
         USE nrtype; USE nrutil, ONLY : nrerror
         USE nrtype
         IMPLICIT NONE
         REAL(DP), INTENT(IN) :: x1,x2,xacc
         integer(I4B),intent(IN) :: N
         REAL(DP), INTENT(IN) :: g_ep,K_el,t,delta
         REAL(DP) :: rtbis2

       END FUNCTION rtbis2
    END INTERFACE




    xacc=1.d-3
    !    x1= 1.0d-20
    x1= 1.0d-20
    x3= 2.5d0*g_ep/K_el
    nstep = int((x3-x1)/STEP_DX)

    ! Trovo i punti in cui la f cambia di segno (ho due soluzioni da isolare)
    x = x1
    x2 = x3
    do i=1,nstep
       x = x + STEP_DX
       if (funcdim2(x1,N,g_ep,K_el,t,delta)*funcdim2(x,N,g_ep,K_el,t,delta) .lt. ZERO) then
          x2=x
          write(17,*) 'x1,x2,x3',x1,x2,x3
          write(17,*) 'f1,f2,f3',funcdim2(x1,N,g_ep,K_el,t,delta),funcdim2(x2,N,g_ep,K_el,t,delta),funcdim2(x3,N,g_ep,K_el,t,delta)
          exit
       end if

    end do

    if(x2 .eq. x3) then ! Condizione per cui c'è solo uno zero
       dimer=0.0d0
       metadimer = 0.0d0
       dflag=.true.
       return
    end if

    dimer1= 0.0d0
    dimer1=rtbis2(x1,x2,N,g_ep,K_el,t,delta,xacc)
    E1 = bulkEnergy(N,g_ep,K_el,t,dimer1,delta)
    write(17,*) 'dimer1',dimer1,'E1',E1 !debug
    Ed=E1
    dimer=dimer1

    !    dimer2= 0.0d0
    dimer2= dimer1 !debug
    if(funcdim2(x2,N,g_ep,K_el,t,delta)*funcdim2(x3,N,g_ep,K_el,t,delta) .lt. ZERO) then
       dimer2=rtbis2(x2,x3,N,g_ep,K_el,t,delta,xacc)
       E2 = bulkEnergy(N,g_ep,K_el,t,dimer2,delta)
       write(17,*) 'dimer2',dimer2,'E2',E2 !debug

       ! Escludo il massimo dell'energia ???????
       if (E1 .lt. E2) then
          dimer = dimer1
          Ed = E1
       else
          dimer = dimer2
          Ed = E2
       end if
    end if

    write(17,*) 'E0=',E0,'Ed=',Ed !debug
    ! Separo soluzione stabile da soluzione metastabile
    if (E0 .lt. Ed) then
       metadimer = dimer
       dimer = 0.0d0
    else
       metadimer = 0.0d0
    end if
    !    dimer=rtbis(func,x1,x2,N,g_ep,K_el,t,xacc) ! Debug


  end subroutine bulkDimer


  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  FUNCTION rtbis(x1,x2,N,g_ep,K_el,t,xacc)
    USE nrtype; USE nrutil, ONLY : nrerror
    !    USE nrtype
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: x1,x2,xacc
    integer(I4B),intent(IN) :: N
    REAL(DP), INTENT(IN) :: g_ep,K_el,t
    REAL(DP) :: rtbis

    REAL(DP) :: func
    integer(I4B) :: i,numK
    REAL(DP) :: k




    INTEGER(I4B), PARAMETER :: MAXIT=100
    !Using bisection, find the root of a function func known to lie between x1 and x2. The
    !root, returned as rtbis, will be refined until its accuracy is ±xacc.
    !Parameter: MAXIT is the maximum allowed number of bisections.
    INTEGER(I4B) :: j
    REAL(DP) :: dx,f,fmid,xmid

    numK = int(dble(N)/2.0d0)
    !write(17,*) 'test_bis',x1,f,x2,fmid ! Debug
    fmid=0.0d0

    do i=0,numK-1

       k = 4.0d0*PI_D*dble(i)/dble(N)
       fmid = fmid + 1.d0/sqrt(g_ep**2*x2**2 + 4.0d0*t**2 * cos(k/2.0d0)**2)
       !       write(17,*) k,func,cos(dble(k))
    end do

    fmid = fmid/dble(N) - K_el/g_ep**2

    f = 0.0d0
    do i=0,numK-1

       k = 4.0d0*PI_D*dble(i)/dble(N)
       f = f + 1.d0/sqrt(g_ep**2*x1**2 + 4.0d0*t**2 * cos(k/2.0d0)**2)
       !       write(17,*) k,func,cos(dble(k))
    end do

    f = f/dble(N) - K_el/g_ep**2

    if (f*fmid > 0.0) call nrerror('rtbis: root must be bracketed')
    if (f < 0.0) then !Orient the search so that f>0 lies at x+dx.
       rtbis=x1
       dx=x2-x1
    else
       rtbis=x2
       dx=x1-x2
    end if
    do j=1,MAXIT !Bisection loop.
       dx=dx*0.5_dp
       xmid=rtbis+dx

       fmid=0.0d0

       do i=0,numK-1

          k = 4.0d0*PI_D*dble(i)/dble(N)
          fmid = fmid + 1.d0/sqrt(g_ep**2*xmid**2 + 4.0d0*t**2 * cos(k/2.0d0)**2)
       end do

       fmid = fmid/dble(N) - K_el/g_ep**2

       if (fmid <= 0.0) rtbis=xmid
       !if (abs(dx) < xacc .or. fmid == 0.0) RETURN
       !if (fmid == 0.0) RETURN
    end do
    !    call nrerror('rtbis: too many bisections')
    return

  END FUNCTION rtbis


  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine freezeXphon(NboundaryL,NboundaryR,dimer,xphon)

    use nrtype
    implicit none

    integer(I4B),intent(IN) :: NboundaryL,NboundaryR
    real(dp),intent(IN) :: dimer
    real(dp),dimension(:),intent(OUT) :: xphon

    integer(I4B) :: i,N

    N=size(xphon)


    do i=(NboundaryL+1),(N-NboundaryR)

       xphon(i)= (-1.0d0)**i * dimer

    end do

  end subroutine freezeXphon

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  subroutine print_energy(K_el,g_ep,eig,xphon,Etot)

    use nrtype
    implicit none

    real(dp),intent(IN) :: K_el,g_ep
    real(dp),dimension(:),intent(IN) :: eig,xphon
    real(dp),intent(OUT) :: Etot

    integer(I4B) :: icol,ncol,i,N

    N=size(xphon)


    i=1
    do
       if(eig(i) .gt. 0.0d0) EXIT
       i=i+1
    end do
    ncol = i-1

    !Test
    write(17,*) 'ncol_energy=',ncol

    Etot = 0.0d0
    do icol=1,ncol
       Etot = Etot + eig(icol)
    end do

    !    write(567,*) xphon**2,sum(xphon**2),Etot
    Etot = (Etot + 0.5d0*K_el*sum(xphon**2))/dble(N) !Controllare 
    !    Etot = (Etot)/N ! Debug
    !    Etot = (Etot + 0.5d0*g_ep*sum(xphon))/dble(N) !Debug


  end subroutine print_energy


  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine print_charge(mu,eig,WN)

    use nrtype
    use param
    implicit none

    real(dp),intent(IN) :: mu
    real(dp),dimension(:),intent(IN) :: eig
    real(dp),dimension(:,:),intent(IN) :: WN

    integer(I4B) :: icol,ncol,i,j,N

    real(dp),dimension(:),allocatable :: charge,Wvec
    real(dp),dimension(2,2) :: prova

    N = int(dble(size(eig))/dble(SMLLDIM))
    write(17,*) 'NprintCharge=',N !Test
    allocate(charge(N))

    i=1
    do       
       if(eig(i) .gt. 0.0d0) EXIT
       i=i+1
    end do
    ncol = i-1

    !Test
    write(17,*) 'ncol_charge=',ncol
!!$    ncol=count(eig .le. mu)
!!$    write(17,*) 'ncolTest=',ncol




!!$       Wvec= WN(:,icol)**2
!!$
!!$       do j=2,n4,2
!!$          Wvec(j) = - Wvec(j)
!!$       end do
!!$
!!$       prova = reshape(Wvec,(/N,2/))
!!$
!!$       do i=1,2
!!$          do j=1,2
!!$             write(17,*) i,j, prova(i,j)
!!$          end do
!!$       end do
!!$!       write(17,*) 'matrice', reshape(Wvec,(/N,2/))
!!$       charge = 0.5d0*sum(reshape(Wvec,(/N,2/)),dim=2) + 0.5d0
    charge = 0.0d0
    do i=1,N
       do icol=1,ncol
          charge(i) = charge(i)+0.5d0*(WN(2*i-1,icol)**2 - WN(2*i,icol)**2)
          !             write(17,*) 'autoval=',icol,'sito=',i,'u2=',WN(2*i-1,icol)**2,'v2=',WN(2*i,icol)**2
       end do
       charge(i) = charge(i) + 0.5d0
    end do


    do i=1,N
       write(19,'(I4,E15.7)') i,charge(i)
    end do

  end subroutine print_charge


  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine print_phonon(xphon)

    use nrtype

    implicit none

    real(dp),dimension(:),intent(IN) :: xphon

    integer(I4B) :: i,N

    N=size(xphon)

    do i=1,N
       write(18,'(I4,E15.7)') i,xphon(i)
    end do

  end subroutine print_phonon

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine sqrdiff(xphonOld,xphonNew,derr,Nl,Nr,diffsqr)

    use nrtype

    implicit none

    integer(I4B), intent(IN) :: Nl,Nr
    real(dp),intent(IN) :: derr
    real(dp),dimension(:),intent(IN) :: xphonOld,xphonNew
    real(dp),intent(OUT) :: diffsqr

    integer(I4B) :: N
    real(dp) :: diffsqr1,diffsqr2,eps

    N = size(xphonNew)
    eps = tiny(eps) + derr

    !    diffsqr = sum((xphonNew - xphonOld)**2/xphonOld**2)/dble(N)
    diffsqr1 = sum((xphonNew - xphonOld)**2)/(sum(xphonOld**2)+eps) / dble(N) !Prendo tutto il sistema
    diffsqr2 = (sum((xphonNew(1:Nl) - xphonOld(1:Nl))**2)/(sum(xphonOld(1:Nl)**2)+eps) + &
         & sum((xphonNew(N-Nr:N) - xphonOld(N-Nr:N))**2)/(sum(xphonOld(N-Nr:N)**2)+eps)) / dble(Nl+Nr) !Escludo il bulk
    if (diffsqr1 .gt. diffsqr2) then
       diffsqr =diffsqr1
    else
       diffsqr =diffsqr2
    end if
!    write(17,*) 'Nl,Nr',Nl,Nr
    write(17,*) '|xphonOld|**2',sum(xphonOld**2),'diff1',diffsqr1,'diff2',diffsqr2 !Debug
    !    diffsqr = sum(((xphonNew - xphonOld)/xphonOld)**2)/dble(N) !Debug

  end subroutine sqrdiff

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine setXphon(g_ep,K_el,mu,eig,WN,xphon)

    use nrtype

    implicit none

    real(dp),intent(IN) :: mu,g_ep,K_el
    real(dp),dimension(:),intent(IN) :: eig
    real(dp),dimension(:,:),intent(IN) :: WN
    real(dp),dimension(:),intent(OUT) :: xphon

    integer(I4B) :: icol,ncol,i,j,N
    real(dp) :: avrgN

    !  real(dp),dimension(:),allocatable :: Wvec


    !    nn4 = size(eig)
    N=size(xphon)
    !   allocate(Wvec(n4))


    i=1
    do       
       if(eig(i) .gt. 0.0d0) EXIT
       i=i+1
    end do
    ncol = i-1

    !Test
    write(17,*) 'ncol_phon=',ncol
!!$    ncol=count(eig .le. mu)
!!$    write(17,*) 'ncolTest=',ncol


!!$    do icol=1,ncol
!!$
!!$       Wvec= WN(:,icol)**2
!!$
!!$       do j=2,n4,2
!!$          Wvec(j) = - Wvec(j)
!!$       end do
!!$
!!$       xphon = g_ep/(2.0d0*K_el)*sum(reshape(Wvec,(/N,2/)),dim=2) + g_ep/(2.0d0*K_el)*WN(:,icol)**2
!!$
!!$    end do

    avrgN = 0.0d0
    do i=1,N
       do icol=1,ncol
          avrgN = avrgN + 0.50d0*(WN(2*i-1,icol)**2 - WN(2*i,icol)**2)
          !          write(17,*) 'sito=',i,'u2=',WN(2*i-1,icol)**2,'v2=',WN(2*i,icol)**2,'autoval=',icol
       end do
    end do

    avrgN = avrgN/dble(N) + 0.5d0
    write(17,*) 'avrgN',avrgN

    xphon = 0.0d0
    do i=1,N
       do icol=1,ncol
          xphon(i) = xphon(i)+(WN(2*i-1,icol)**2 - WN(2*i,icol)**2)*0.5d0
          !             write(17,*) 'autoval=',icol,'sito=',i,'u2=',WN(2*i-1,icol)**2,'v2=',WN(2*i,icol)**2
       end do
    end do
!    xphon = (xphon + 0.5d0 - avrgN)*(-2.0d0*g_ep)/K_el    !!!!### controllare il fattore 2
    xphon = (xphon + 0.5d0 - avrgN)*(-g_ep)/K_el    !!!!### controllare il fattore 2


    !    deallocate(Wvec)
  end subroutine setXphon

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine setXphon0(dimer,xphon,init_mode)

    use nrtype

    implicit none

    real(dp),intent(IN) :: dimer
    real(dp),dimension(:),intent(OUT) :: xphon
    real :: rnd
    integer(I4B) :: N,i
    character(12) init_mode

!    CALL RANDOM_INIT(.false., .false.)
    
    N = size(xphon)
!
! initialize random number generator
!    
    if (trim(init_mode).eq.'random') then 
     call RANDOM_SEED()
 !    CALL RANDOM_NUMBER(rnd)
 !    write(*,*) rnd
    end if
!    if (trim(init_mode).eq.'random') call init_random_seed()
    do i=1,N

	if (trim(init_mode).eq.'trimer') then
       xphon(i)= sin(2.*i*pi_D/3.d0) * dimer !Questo rompe la simmetria
      end if
      
	if (trim(init_mode).eq.'dimer') then
       xphon(i)= (-1.0d0)**i * dimer !Questo rompe la simmetria
      end if
	if (trim(init_mode).eq.'random') then
	 CALL RANDOM_NUMBER(rnd)
	 xphon(i) = dimer*(rnd-0.5)
	end if
	if (trim(init_mode).eq.'uniform') then
	 xphon(i) = dimer
	end if	
	write(98,*) xphon(i)
    end do

    !    xphon = 0.0d0

  end subroutine setXphon0

SUBROUTINE init_random_seed()
            INTEGER :: i, n, clock
            INTEGER, DIMENSION(:), ALLOCATABLE :: seed
          
            CALL RANDOM_SEED(size = n)
            ALLOCATE(seed(n))
          
            CALL SYSTEM_CLOCK(COUNT=clock)
          
            seed = clock + 37 * (/ (i - 1, i = 1, n) /)
            CALL RANDOM_SEED(PUT = seed)
          
            DEALLOCATE(seed)
END SUBROUTINE

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine setHe(n4,g_ep,xphon,Ham)

    use nrtype
    use param
    implicit none

    integer,intent(IN) :: n4
    real(dp),intent(IN) :: g_ep
    real(dp),dimension(:),intent(IN) :: xphon
    real(dp),dimension(:,:),intent(OUT) :: Ham

    real(dp),dimension(LRGDIM,LRGDIM) :: M
    real(dp),dimension(SMLLDIM,SMLLDIM) :: H1 !SMLLDIM=2
    integer :: isite,iHam,i,j,nHam

    nHam = 2*size(xphon) - SMLLDIM

    isite = 0
    do iHam=0,nHam,SMLLDIM

       H1=0.0d0
       M = 0.0d0

       isite = isite + 1
!
!	e-ph ??
!
       H1(1,1)=  g_ep*xphon(isite)/2.d0
       H1(2,2)= -g_ep*xphon(isite)/2.d0

       !       write(17,*) 'iHam=',iHam,'isite=',isite,'xphon(isite)=',xphon(isite)

       !       M(1:SMLLDIM,1:SMLLDIM) = H1
       !       M(LRGDIM-1:LRGDIM,LRGDIM-1:LRGDIM) = H1

       Ham(iHam+1:iHam+SMLLDIM,iHam+1:iHam+SMLLDIM)= Ham(iHam+1:iHam+SMLLDIM,iHam+1:iHam+SMLLDIM) + H1

       !       Ham(iHam+1:iHam+LRGDIM,iHam+1:iHam+LRGDIM)= Ham(iHam+1:iHam+LRGDIM,iHam+1:iHam+LRGDIM) + M

    end do

!!$    do i=1,n4+LRGDIM
!!$       do j=1,n4+LRGDIM
!!$          write(204,*) i,j,Ham(i,j)
!!$!          write(17,*) 'epporcaputtana'
!!$       end do
!!$    end do


  end subroutine setHe

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!!$  subroutine averageDos()
!!$
!!$    implicit none
!!$
!!$     en_grid = en_grid/dble(ndis)
!!$     en=en_start
!!$     do i=1,nGrid+1
!!$        write(16,*) en, en_grid(i)
!!$        en = en + delta_en
!!$     end do
!!$
!!$  end subroutine averageDos

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!!$  subroutine localDoS()
!!$
!!$    implicit none
!!$
!!$    if (idis .eq. 1) en_start=EIG(1)
!!$
!!$    en=en_start
!!$    delta_en = -2.0d0*en_start/dble(nGrid)
!!$    !     write(17,*) 'ciao',en_start,delta_en
!!$
!!$    do i=1,nGrid+1
!!$       !     write(18,*) en, dos(en,width,nn4,eig)
!!$       en_grid(i) = en_grid(i) + dos(en,width,nn4,eig)
!!$       en = en + delta_en
!!$    end do
!!$
!!$  end subroutine localDoS

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine printMajorana(Ham)

    use nrtype
    use param
    implicit none

    real(dp),dimension(:,:),intent(IN) :: Ham
    integer(I4B) :: jaut,i,j,n4
!    REAL(dp), ALLOCATABLE, DIMENSION (:):: MMSR, MMSI,MMS1, MMS2
    REAL(dp), ALLOCATABLE, DIMENSION (:) :: MMS1, MMS2

    n4 = size(Ham,dim=1)
    jaut = int(dble(n4)/2.0d0)

    allocate(MMS1(n4),MMS2(n4))

    !write(17,*) 'pollo',jaut,n4


101 FORMAT(I5,4F15.8)

    !MMs

    MMS1 = 0.0d0
    MMS2 = 0.0d0


    MMS1 = 0.5d0*(Ham(:,jaut) + Ham(:,jaut+1))**2
    MMS2 = 0.5d0*(Ham(:,jaut) - Ham(:,jaut+1))**2

!    MMSR= MMS1+cshift(MMS1,1)
!    MMSI=MMS2+cshift(MMS2,1)

    j=0
    do i=1,jaut
!       write(15,101) i, MMSR(i+j),MMSR(i+j+2),MMSI(i+j),MMSI(i+j+2)
       write(15,101) i, MMS1(i+j),MMS1(i+j+1),MMS2(i+j),MMS2(i+j+1)
       !write(15,*) 'i+j=',i+j,'i+j+2=',i+j+2
       j=j+1

    end do

    deallocate(MMS1,MMS2)
  end subroutine printMajorana
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine diag(nn4,Ham,eig)

    use nrtype

    implicit none

    integer,intent(IN) :: nn4
    real(dp),dimension(:,:),intent(INOUT) :: Ham
    real(dp),dimension(:),intent(OUT) :: eig

    integer(I4B) :: lwork,info,errore1
    REAL(dp), ALLOCATABLE, DIMENSION (:) :: WORK

    lwork=3*nn4-1

    allocate(WORK(lwork),STAT=errore1)

    !Utilizzo routine per calcolo autoval. e autovett. matrice real-simm
    call DSYEV('V','U',nn4,Ham,nn4,eig,WORK,lwork,INFO)
    if(info.ne.0)then
       write(6,*)'info=',info
       stop
    endif

    deallocate(WORK)
  end subroutine diag

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine setMajoranaBasis(UN,Ham)

    use nrtype

    implicit none

    real(dp),dimension(:,:),intent(IN) :: UN
    real(dp),dimension(:,:),intent(OUT) :: Ham

 !   real(dp),dimension(:,:),allocatable :: MN
!    integer(I4B) :: dim

!    dim = size(UN,dim=1)

!    allocate(MN(dim,dim))

!    MN = transpose(UN)
!    Ham = matmul(Ham,MN)
    Ham=matmul(UN,Ham) ! Passo alla base dei Majorana (forse)

!    deallocate(MN)


  end subroutine setMajoranaBasis

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!!$  subroutine setdis()
!!$
!!$    implicit none
!!$
!!$        !genero 4N numeri random e li memorizzo nell'array Rd
!!$
!!$        do j=1,nn4
!!$           CALL RANDOM_NUMBER(num)
!!$           Rd(j)=(num-0.50d0)*d
!!$        end do
!!$
!!$
!!$        forall(i=0:nlast:4)
!!$           Ham(i+1,i+1) = Ham(i+1,i+1) + Rd(i+1)
!!$           Ham(i+2,i+2) = Ham(i+2,i+2) + Rd(i+1)
!!$
!!$           Ham(i+3,i+3) = Ham(i+3,i+3) - Rd(i+1)
!!$           Ham(i+4,i+4) = Ham(i+4,i+4) - Rd(i+1)
!!$        end forall
!!$
!!$  end subroutine setdis

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine buildHam(n4,tt,mu,delta,Ham)

    use nrtype
    use param

    implicit none

    integer,intent(IN) :: n4
    real(dp),intent(IN) :: tt,mu,delta
    real(dp),dimension(:,:),intent(OUT) :: Ham

    real(dp),dimension(SMLLDIM,SMLLDIM) :: H1,H2 !SMLLDIM=2
    real(dp),dimension(LRGDIM,LRGDIM) :: M
    real(dp) :: t
    integer :: ilink,iHam,i,j


    Ham = 0.0d0
    ilink = 0


    do iHam=0,n4,SMLLDIM

       H1=0.0d0
       H2=0.0d0
       M = 0.0d0

       ilink = ilink + 1

       !       t = tt + dimer*(-1)**ilink
       t = tt
       !     write(17,*) 'ilink,t',ilink,t

       H1(1,1)=-mu/2.d0

       H2(1,1)=-t/2.d0
       H2(1,2)=delta/2.d0
                                !EDIT PRODANOV
       H1(2,2)=mu/2.d0
       H2(2,1)=-delta/2.d0
       H2(2,2)=t/2.d0


       forall(j=1:SMLLDIM,i=1:SMLLDIM,i.lt.j)
          H1(j,i)=H1(i,j)
       end forall

       M(1:SMLLDIM,1:SMLLDIM) = H1
       M(1:SMLLDIM,LRGDIM-1:LRGDIM) = H2
       M(LRGDIM-1:LRGDIM,1:SMLLDIM) = transpose(H2)
       M(LRGDIM-1:LRGDIM,LRGDIM-1:LRGDIM) = H1

       Ham(iHam+1:iHam+LRGDIM,iHam+1:iHam+LRGDIM)=M
       !     write(17,*) iHam
    end do

!    do i=1,8
!       do j=1,8
!          write(177,*) i,j,Ham(i,j)
!       end do
!    end do

  end subroutine buildHam
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine buildHamPBC(n4,tt,mu,delta,Ham)

    use nrtype
    use param

    implicit none

    integer,intent(IN) :: n4
    real(dp),intent(IN) :: tt,mu,delta
    real(dp),dimension(:,:),intent(OUT) :: Ham

    real(dp),dimension(SMLLDIM,SMLLDIM) :: H1,H2 !SMLLDIM=2
    real(dp),dimension(LRGDIM,LRGDIM) :: M
    real(dp) :: t
    integer :: iHam,i,j


    Ham = 0.0d0


    H1=0.0d0
    H2=0.0d0
    M = 0.0d0

    t = tt

    H1(1,1)=-mu/2.d0

    H2(1,1)=-t/2.d0
    H2(1,2)=delta/2.d0
                            !EDIT PRODANOV
    H1(2,2)=mu/2.d0
    H2(2,1)=-delta/2.d0
    H2(2,2)=t/2.d0


    forall(j=1:SMLLDIM,i=1:SMLLDIM,i.lt.j)
       H1(j,i)=H1(i,j)
    end forall

    M(1:SMLLDIM,1:SMLLDIM) = H1
    M(1:SMLLDIM,LRGDIM-1:LRGDIM) = H2
    M(LRGDIM-1:LRGDIM,1:SMLLDIM) = transpose(H2)
    M(LRGDIM-1:LRGDIM,LRGDIM-1:LRGDIM) = H1

    do iHam=0,n4,SMLLDIM
       Ham(iHam+1:iHam+LRGDIM,iHam+1:iHam+LRGDIM)=M
       !     write(17,*) iHam
    end do
    !
    !	PBC
    !
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    Ham(1:SMLLDIM,n4+SMLLDIM+1:n4+LRGDIM)=H2
!    Ham(n4+SMLLDIM+1:n4+LRGDIM,1:SMLLDIM)=transpose(H2)

    Ham(1:SMLLDIM,n4+SMLLDIM+1:n4+LRGDIM)=transpose(H2)
    Ham(n4+SMLLDIM+1:n4+LRGDIM,1:SMLLDIM)=H2

    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!$    do i=1,4
!!$       do j=1,4
!!$          write(17,*) i,j,Ham(i,j)
!!$       end do
!!$    end do

  end subroutine buildHamPBC
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  subroutine setarrayval(arr,val)

    use nrtype

    implicit none

    real(dp),dimension(:),intent(OUT) :: arr
    real(dp),intent(IN) :: val

    arr = val

  end subroutine setarrayval


  subroutine setId(M)

    use nrtype

    implicit none

    real(dp),dimension(:,:),intent(OUT) :: M

    integer :: dim,i

    dim = size(M,dim=1)

    M=0.0d0
    forall (i=1:dim)
       M(i,i) = 1.0d0
    end forall


  end subroutine setId

  subroutine setUN(U4,UN,nlast)

    use nrtype
    use param

    implicit none

    real(dp),dimension(:,:),intent(OUT) :: UN
    real(dp),dimension(:,:),intent(IN) :: U4
    integer,intent(IN) :: nlast

    integer :: i

    do i=0,nlast,SMLLDIM
       UN(i+1:i+SMLLDIM,i+1:i+SMLLDIM)=U4
    end do


  end subroutine setUN


  subroutine setU4(U4)

    use nrtype

    implicit none

    real(dp),dimension(:,:),intent(OUT) :: U4

    U4(1,1)=1.0d0
    U4(1,2)=1.0d0
    U4(2,1)=1.0d0
    U4(2,2)=-1.0d0

    U4=U4/sqrt(2.0d0)

  end subroutine setU4


  subroutine setmatrixval(Matrix,val)

    use nrtype

    implicit none

    real(dp),dimension(:,:),intent(OUT) :: Matrix
    real(dp),intent(IN) :: val

    Matrix = val

  end subroutine setmatrixval


  subroutine print_spectrum(eig,Ham,EFLAG,AFLAG,jaut,nn4)
    use nrtype
    implicit none

    real(dp),dimension(nn4),intent(IN) :: eig
    real(dp),dimension(nn4,nn4),intent(IN) :: Ham
    integer,intent(IN) :: EFLAG,AFLAG,jaut,nn4

    integer*4 :: i,j

    if (EFLAG .eq. 0) then
       !Scrittura su file autovalori
       do i=1,nn4
          write(13,*)i,eig(i)
       enddo
    end if

    if (AFLAG .eq. 0) then
       !Scrittura su file autovettori

!!$     do i=1,nn4
!!$        write(14,*) i,Ham(i,jaut),Ham(i,jaut+1) 
!!$     enddo

       do j=1,nn4
          do i=1,nn4
             write(14,'(1X,2I4,E15.7)') i,j,Ham(i,j)
          end do
       end do


    end if
    return
  end subroutine print_spectrum

  function dos(en,width,nn4,eig)
    use nrtype
    IMPLICIT NONE

    real(dp),intent(IN) :: en,width
    integer,intent(IN) :: nn4
    real(dp),dimension(nn4),intent(IN) :: eig

    real(dp) :: dos

    integer :: i
    real(dp) :: sum

    sum=0.0d0

    do i=1,nn4

       sum = sum + width/((en-eig(i))**2+width**2)

    end do

    dos = sum/(PI_D*dble(nn4))

  end function dos

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!!$  subroutine print_log()
!!$
!!$    implicit none
!!$
!!$
!!$  end subroutine print_log

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine open_file(fileId,inputFilename,filetype,date,time,outputFileName)

    use nrtype

    implicit none

    integer(I4B),intent(IN) :: fileId
    character(80),intent(IN) :: inputFileName,filetype
    character(8),intent(IN) :: date
    character(10),intent(IN) :: time
    character(80),intent(OUT) :: outputFileName


    write(outputFileName,'(a,a,a)') inputfilename,date,filetype
    open(unit=fileId,file=outputFileName,form='formatted',status='unknown')

    write(17,*) outputFileName

  end subroutine open_file

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function funcdim(x,N,g_ep,K_el,t,delta)
    USE nrtype
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: x
    integer(I4B),intent(IN) :: N
    REAL(DP),intent(IN) :: g_ep,K_el,t,delta
    integer(I4B) :: i,numK
    REAL(DP) :: k,delta_k,csi_k,radp_k,radm_k
    real(dp) :: funcdim

    numK = int(dble(N)/2.0d0)

    funcdim=0.0d0

    do i=0,numK-1

       k = 4.0d0*PI_D*dble(i)/dble(N)
       delta_k = 2.0d0*delta*sin(k/2.0d0)
       csi_k = 2.0d0*t*cos(k/2.0d0)
       radp_k=sqrt((g_ep*x + delta_k)**2 + csi_k**2)
       radm_k=sqrt((g_ep*x - delta_k)**2 + csi_k**2)

       funcdim = funcdim + (g_ep*x + delta_k)/radp_k + (g_ep*x - delta_k)/radm_k
       !       write(17,*) k,func,cos(dble(k))
    end do

    funcdim = (g_ep/(2.0d0*K_el*dble(N)))*funcdim - x

  end function funcdim

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function funcdim2(x,N,g_ep,K_el,t,delta)
    USE nrtype
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: x
    integer(I4B),intent(IN) :: N
    REAL(DP),intent(IN) :: g_ep,K_el,t,delta

    integer(I4B) :: i,numK
    REAL(DP) :: k,delta_k,csi_k,radp_k,radm_k
    real(dp) :: funcdim2

    numK = int(dble(N)/2.0d0)


    funcdim2=0.0d0

    do i=0,numK-1

       k = 4.0d0*PI_D*dble(i)/dble(N)
       delta_k = 2.0d0*delta*sin(k/2.0d0)
       csi_k = 2.0d0*t*cos(k/2.0d0)
       radp_k=sqrt((g_ep*x + delta_k)**2 + csi_k**2)
       radm_k=sqrt((g_ep*x - delta_k)**2 + csi_k**2)

!       funcdim2 = funcdim2 + (1.0d0/radp_k + 1.0d0/radm_k)*(1.0d0 - 2.0d0*delta_k**2/(radp_k*radm_k))
       funcdim2 = funcdim2 + 1.0d0/radp_k + 1.0d0/radm_k  - 4.0d0*delta_k**2/((radp_k*radm_k)*(radm_k + radp_k))
    end do

    funcdim2 = (g_ep**2/(dble(N)*2.0d0*K_el))*funcdim2 - 1.0d0

!!$    (g_ep*x+delta_k)/sqrt((g_ep*x + delta_k)**2 + csi_k**2) + & 
!!$            & (g_ep*x-delta_k)/sqrt((g_ep*x - delta_k)**2 + csi_k**2)

  end function funcdim2


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  FUNCTION rtbis2(x1,x2,N,g_ep,K_el,t,delta,xacc)
    USE nrtype; USE nrutil, ONLY : nrerror
!    USE nrtype
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: x1,x2,xacc
    integer(I4B),intent(IN) :: N
    REAL(DP), INTENT(IN) :: g_ep,K_el,t,delta
    REAL(DP) :: rtbis2

    INTEGER(I4B), PARAMETER :: MAXIT=100
    !Using bisection, find the root of a function func known to lie between x1 and x2. The
    !root, returned as rtbis, will be refined until its accuracy is ±xacc.
    !Parameter: MAXIT is the maximum allowed number of bisections.
    INTEGER(I4B) :: j
    REAL(DP) :: dx,f,fmid,xmid
    fmid=funcdim2(x2,N,g_ep,K_el,t,delta)
    f=funcdim2(x1,N,g_ep,K_el,t,delta)
    write(17,*) 'x1,f,x2,fmid',x1,f,x2,fmid ! Debug
    !if (f*fmid > 0.0) call nrerror('rtbis2: root must be bracketed')
    if (f*fmid > 0.0) return
!!$    if (f*fmid > 0.0) then
!!$       rtbis2 = 0.0d0 ! Unica soluzione in x=0 
!!$       return
!!$    end if
    !  if (f*fmid > 0.0) write(17,*) 'cazzo'
    if (f < 0.0) then !Orient the search so that f>0 lies at x+dx.
       rtbis2=x1
       dx=x2-x1
    else
       rtbis2=x2
       dx=x1-x2
    end if
    do j=1,MAXIT !Bisection loop.
       dx=dx*0.5_dp
       xmid=rtbis2+dx
       fmid=funcdim2(xmid,N,g_ep,K_el,t,delta)
       if (fmid <= 0.0) rtbis2=xmid
    end do
    !    call nrerror('rtbis: too many bisections')
    return

  END FUNCTION rtbis2
!
! strip spaces
!
    subroutine StripSpaces(string)
    character(len=*) :: string
    integer :: stringLen 
    integer :: last, actual

    stringLen = len (string)
    last = 1
    actual = 1

    do while (actual < stringLen)
        if (string(last:last) == ' ') then
            actual = actual + 1
            string(last:last) = string(actual:actual)
            string(actual:actual) = ' '
        else
            last = last + 1
            if (actual < last) &
                actual = last
        endif
    end do

    end subroutine
    

end module sub
