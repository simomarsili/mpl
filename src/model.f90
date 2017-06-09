! Copyright (C) 2015, 2016, Simone Marsili 
! All rights reserved.
! License: BSD 3 clause

module model
  use kinds
  use data, only: nv,ns,nd
  implicit none
  private

  public :: prm,grd

  public :: model_initialize
  public :: model_set_myv
  public :: model_put_myv
  public :: model_collect_prm
  public :: compute_pseudo_likelihood

  public :: mypl
  public :: etot
  public :: fields
  public :: couplings

  ! working variable for the processor
  integer(kint) :: myv 
  real(kflt), allocatable :: fields(:,:) ! NS x NV
  real(kflt), allocatable :: couplings(:,:,:) ! NS x NS x NV(NV-1)/2 
  real(kflt), allocatable :: p1(:,:) ! NS x NV
  real(kflt), allocatable :: p2(:,:,:) ! NS x NS x NV(NV-1)/2 

  ! arrays for fixed residue/sequence
  real(kflt), allocatable, save :: my_fields(:) ! NS
  real(kflt), allocatable, save :: my_couplings(:,:,:) ! NS x NV x NS
  real(kflt), allocatable, save :: my_p1(:) ! NS 
  real(kflt), allocatable, save :: my_p2(:,:,:) ! NS x NV x NS
  real(kflt), allocatable, save :: my_f1(:) ! NS 
  real(kflt), allocatable, save :: my_f2(:,:,:) ! NS x NV x NS

  ! regularization parameters 
  real(kflt) :: regularization_strength=0.01_kflt ! regularization strength; the default is l2 with regularization_strength=0.01

  ! "cost" functions
  real(kflt) :: mypl,ereg,etot

  ! gradients arrays
  real(kflt), allocatable, save :: grd(:) ! NS + NS x NV x NS
  real(kflt), allocatable, save :: prm(:) ! NS + NS x NV x NS

contains

  subroutine model_initialize(lambda)
    integer(kint) :: err
    real(kflt), intent(in) :: lambda

    regularization_strength = lambda

    allocate(p1(ns,nv),stat=err)
    allocate(p2(ns,ns,nv*(nv-1)/2),stat=err)
    allocate(fields(ns,nv),stat=err)
    allocate(couplings(ns,ns,nv*(nv-1)/2),stat=err)
    allocate(my_f1(ns),stat=err)
    allocate(my_f2(ns,nv,ns),stat=err)
    allocate(my_fields(ns),stat=err)
    allocate(my_couplings(ns,nv,ns),stat=err)
    allocate(my_p1(ns),stat=err)
    allocate(my_p2(ns,nv,ns),stat=err)
    allocate(grd(ns + ns*ns*nv),stat=err)
    allocate(prm(ns + ns*ns*nv),stat=err)


    p1 = 0.0_kflt
    p2 = 0.0_kflt
    fields = 0.0_kflt
    couplings = 0.0_kflt       

  end subroutine model_initialize

  subroutine model_zero()
    real(kflt), parameter :: small_number=1.e-2_kflt
    my_p1 = 0.0_kflt
    my_p2 = 0.0_kflt
    mypl = 0.0_kflt
    etot = 0.0_kflt
    ereg = - regularization_strength * (sum(prm(:ns)**2) + 0.5_kflt * sum(prm(ns+1:)**2))
  end subroutine model_zero

  subroutine model_set_myv(iv,err) ! my_couplings
    use data, only: nd,data_samples,w
    integer(kint), intent(in) :: iv
    ! make my_couplings given myv 
    ! must be called before looping on data
    integer(kint) :: id,jv,mys
    integer(kint) :: err
    integer, allocatable :: list(:)

    myv = iv
    my_p1 = 0.0_kflt
    my_p2 = 0.0_kflt
    my_f1 = 0.0_kflt
    my_f2 = 0.0_kflt
    my_fields = 0.0_kflt
    my_couplings = 0.0_kflt       
    mypl = 0.0_kflt
    etot = 0.0_kflt
    ereg = 0.0_kflt
    prm = 0.0_kflt
    grd = 0.0_kflt

    ! compute variable-specific arrays of frequencies
    allocate(list(nv),stat=err)
    my_f1 = 0.0_kflt
    my_f2 = 0.0_kflt
    do id = 1,nd
       list = data_samples(:,id)
       mys = list(myv)
       my_f1(mys) = my_f1(mys) + w(id)
       do jv = 1,nv
          if(jv /= myv) then 
             my_f2(list(jv),jv,mys) = my_f2(list(jv),jv,mys) + w(id)
          end if
       end do
    end do
    deallocate(list)

  end subroutine model_set_myv

  subroutine model_put_myv 
    ! store my_fields in fields 
    ! store my_couplings in couplings
    ! must be called before looping on data
    integer(kint) :: jv,k

    fields(:,myv) = my_fields
    
    !      --myv--
    !      x x x x
    !   |  1 x x x
    !  jv  2 4 x x
    !   |  3 5 6 x
    !
    ! lower triangle packing: jv > myv

    ! remove gauge before adding to couplings
    call model_gauge()

    fields(:,myv) = my_fields
    do jv = myv+1,nv ! jv > myv
       k = (myv - 1) * nv - myv * (myv + 1) / 2 + jv 
       couplings(:,:,k) = couplings(:,:,k) + my_couplings(:,jv,:)
    end do

    do jv = 1,myv-1 ! myv > jv: submatrices must be transposed 
       k = (jv - 1) * nv - jv * (jv + 1) / 2 + myv
       couplings(:,:,k) =  couplings(:,:,k) + transpose(my_couplings(:,jv,:))
    end do

  end subroutine model_put_myv

  subroutine model_gauge
    integer(kint) :: jv,is,js
    real(kflt) :: mat(ns,ns),arr(ns),marr
    real(kflt) :: rsum(ns),csum(ns),totsum

    arr = my_fields
    marr = sum(arr) / real(ns)
    arr = arr - marr
    do jv = 1,nv
       mat = my_couplings(:,jv,:)
       totsum = sum(mat)
       totsum = totsum / real(ns*ns)
       do is = 1,ns
          rsum(is) = sum(mat(is,:))
          csum(is) = sum(mat(:,is))
       end do
       rsum = rsum / real(ns)
       csum = csum / real(ns)
       if(jv /= myv) arr = arr + csum - totsum
       do js = 1,ns
          do is = 1,ns
             mat(is,js) = mat(is,js) - rsum(is) - csum(js) + totsum
          end do
       end do
       my_couplings(:,jv,:) = mat
    end do
    my_fields = arr
    
  end subroutine model_gauge

  subroutine compute_pseudo_conp()
    use data, only: data_samples,w,nd
    integer(kint) :: list(nv)
    real(kflt) :: conp(ns)
    real(kflt) :: r,rsum
    real(kflt) :: pp
    integer(kint):: mys
    integer :: id
    real(kflt) :: ww
    integer(kint) :: is, jv

    ! loop over data
    do id = 1,nd
       list = data_samples(:,id)
       ww = w(id)
       mys = list(myv)
       
       ! loop over the states of myv 
       do is = 1,ns
          r = my_fields(is) 
          do jv = 1,nv
             if(myv /= jv) then 
                r = r + my_couplings(list(jv),jv,is) 
             end if
          end do
          conp(is) = exp(r)
       end do
       
       rsum = sum(conp)
       conp = conp / rsum
       
       mypl = mypl + ww * log(conp(mys))
       
       ! update histograms 
       ! loop over the states of myv 
       do is = 1,ns
          pp = conp(is) * ww
          my_p1(is) = my_p1(is) + pp 
          do jv = 1,nv
             if(myv /= jv) then 
                my_p2(list(jv),jv,is) = my_p2(list(jv),jv,is) + pp
             end if
          end do
       end do
       
    end do
    
  end subroutine compute_pseudo_conp

  subroutine compute_pseudo_likelihood(it)
    integer(kint), intent(in) :: it
    real(kflt) :: etot0,de

    call model_parameters_unpack()

    etot0 = etot
    call model_zero()
    call compute_pseudo_conp()
    etot = mypl + ereg
    de = etot-etot0

    call model_parameters_pack()
    
  end subroutine compute_pseudo_likelihood

  subroutine model_parameters_pack
    integer(kint) :: dim
    real(kflt), parameter :: small_number=1.e-2_kflt

    dim = ns*ns*nv
    ! compute the gradient
    my_p1 = my_p1 - my_f1 
    my_p2 = my_p2 - my_f2 

    ! l2 regularization
    my_p1 = my_p1 + 2.0_kflt * regularization_strength * my_fields
    my_p2 = my_p2 + 2.0_kflt * 0.5_kflt * regularization_strength * my_couplings

    prm(1:ns) = my_fields
    prm(ns+1:) = reshape(my_couplings,(/dim/))

    grd(1:ns) = my_p1
    grd(ns+1:) = reshape(my_p2,(/dim/))

  end subroutine model_parameters_pack

  subroutine model_parameters_unpack

    my_fields = prm(1:ns) 
    my_couplings = reshape(prm(ns+1:),(/ns,nv,ns/))
    
  end subroutine model_parameters_unpack

  subroutine model_collect_prm

    couplings = 0.5_kflt * couplings

  end subroutine model_collect_prm
  

end module model
  
