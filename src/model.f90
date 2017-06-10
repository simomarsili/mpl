! Copyright (C) 2015, 2016, Simone Marsili 
! All rights reserved.
! License: BSD 3 clause

module model
  use kinds
  use data, only: nv,ns,nd
  implicit none
  private

  public :: prm,grd
  public :: initialize_model
  public :: model_set_myv
  public :: model_put_myv
  public :: model_collect_prm
  public :: update_gradient
  public :: cond_likelihood
  public :: etot
  public :: fields
  public :: couplings

  ! index of outcome variable (all the others nv - 1 variables are included in the set of explanatory variables)
  integer :: out_var
  
  real(kflt), allocatable :: fields(:,:)      ! fields (ns x nv)
  real(kflt), allocatable :: couplings(:,:,:) ! couplings (ns x ns x nv(nv-1)/2)

  ! arrays for fixed residue/sequence
  real(kflt), allocatable, save :: vfields(:)        ! out_var fields array (ns)
  real(kflt), allocatable, save :: vcouplings(:,:,:) ! out_var couplings (ns x nv x ns)
  real(kflt), allocatable, save :: model_f1(:)       ! model single-variable frequencies (ns)
  real(kflt), allocatable, save :: model_f2(:,:,:)   ! model frequencies for pairs of variables (ns x nv x ns)
  real(kflt), allocatable, save :: data_f1(:)        ! data single-variable frequencies (ns)
  real(kflt), allocatable, save :: data_f2(:,:,:)    ! data frequencies for pairs of variables (ns x nv x ns)

  ! regularization parameters 
  real(kflt) :: regularization_strength=0.01_kflt ! regularization strength; the default is l2 with regularization_strength=0.01

  ! "cost" function-related variables
  real(kflt) :: cond_likelihood,ereg,etot

  ! parameters and gradients 1D arrays for the optimization routines
  real(kflt), allocatable, save :: prm(:) ! 1D array of parameters (ns + ns x nv x ns)
  real(kflt), allocatable, save :: grd(:) ! 1D gradient array (ns + ns x nv x ns)

contains

  subroutine initialize_model(lambda)
    real(kflt), intent(in) :: lambda
    integer :: err

    regularization_strength = lambda

    allocate(fields(ns,nv),stat=err)
    allocate(couplings(ns,ns,nv*(nv-1)/2),stat=err)
    allocate(data_f1(ns),stat=err)
    allocate(data_f2(ns,nv,ns),stat=err)
    allocate(vfields(ns),stat=err)
    allocate(vcouplings(ns,nv,ns),stat=err)
    allocate(model_f1(ns),stat=err)
    allocate(model_f2(ns,nv,ns),stat=err)
    allocate(grd(ns + ns*ns*nv),stat=err)
    allocate(prm(ns + ns*ns*nv),stat=err)


    fields = 0.0_kflt
    couplings = 0.0_kflt       

  end subroutine initialize_model

  subroutine model_set_myv(iv,err) ! vcouplings
    use data, only: nd,data_samples,w
    integer, intent(in) :: iv
    ! make vcouplings given out_var 
    ! must be called before looping on data
    integer :: id,jv,mys
    integer :: err
    integer, allocatable :: list(:)

    out_var = iv
    model_f1 = 0.0_kflt
    model_f2 = 0.0_kflt
    data_f1 = 0.0_kflt
    data_f2 = 0.0_kflt
    vfields = 0.0_kflt
    vcouplings = 0.0_kflt       
    cond_likelihood = 0.0_kflt
    etot = 0.0_kflt
    ereg = 0.0_kflt
    prm = 0.0_kflt
    grd = 0.0_kflt

    ! compute variable-specific arrays of frequencies
    allocate(list(nv),stat=err)
    data_f1 = 0.0_kflt
    data_f2 = 0.0_kflt
    do id = 1,nd
       list = data_samples(:,id)
       mys = list(out_var)
       data_f1(mys) = data_f1(mys) + w(id)
       do jv = 1,nv
          if(jv /= out_var) then 
             data_f2(list(jv),jv,mys) = data_f2(list(jv),jv,mys) + w(id)
          end if
       end do
    end do
    deallocate(list)

  end subroutine model_set_myv

  subroutine model_put_myv 
    ! store vfields in fields 
    ! store vcouplings in couplings
    ! must be called before looping on data
    integer :: jv,k

    fields(:,out_var) = vfields
    
    !      --out_var--
    !      x x x x
    !   |  1 x x x
    !  jv  2 4 x x
    !   |  3 5 6 x
    !
    ! lower triangle packing: jv > out_var

    ! remove gauge before adding to couplings
    call fix_gauge()

    fields(:,out_var) = vfields
    do jv = out_var+1,nv ! jv > out_var
       k = (out_var - 1) * nv - out_var * (out_var + 1) / 2 + jv 
       couplings(:,:,k) = couplings(:,:,k) + vcouplings(:,jv,:)
    end do

    do jv = 1,out_var-1 ! out_var > jv: submatrices must be transposed 
       k = (jv - 1) * nv - jv * (jv + 1) / 2 + out_var
       couplings(:,:,k) =  couplings(:,:,k) + transpose(vcouplings(:,jv,:))
    end do

  end subroutine model_put_myv

  subroutine fix_gauge
    integer :: jv,is,js
    real(kflt) :: mat(ns,ns),arr(ns),marr
    real(kflt) :: rsum(ns),csum(ns),totsum

    arr = vfields
    marr = sum(arr) / real(ns)
    arr = arr - marr
    do jv = 1,nv
       mat = vcouplings(:,jv,:)
       totsum = sum(mat)
       totsum = totsum / real(ns*ns)
       do is = 1,ns
          rsum(is) = sum(mat(is,:))
          csum(is) = sum(mat(:,is))
       end do
       rsum = rsum / real(ns)
       csum = csum / real(ns)
       if(jv /= out_var) arr = arr + csum - totsum
       do js = 1,ns
          do is = 1,ns
             mat(is,js) = mat(is,js) - rsum(is) - csum(js) + totsum
          end do
       end do
       vcouplings(:,jv,:) = mat
    end do
    vfields = arr
    
  end subroutine fix_gauge

  subroutine update_model_averages()
    use data, only: data_samples,w,nd
    integer :: list(nv)
    real(kflt) :: conp(ns)
    real(kflt) :: r,rsum
    real(kflt) :: pp
    integer :: mys
    integer :: id
    real(kflt) :: ww
    integer :: is, jv

    ! loop over data
    do id = 1,nd
       list = data_samples(:,id)
       ww = w(id)
       mys = list(out_var)
       
       ! loop over the states of out_var 
       do is = 1,ns
          r = vfields(is) 
          do jv = 1,nv
             if(out_var /= jv) then 
                r = r + vcouplings(list(jv),jv,is) 
             end if
          end do
          conp(is) = exp(r)
       end do
       
       rsum = sum(conp)
       conp = conp / rsum
       
       cond_likelihood = cond_likelihood + ww * log(conp(mys))
       
       ! update histograms 
       ! loop over the states of out_var 
       do is = 1,ns
          pp = conp(is) * ww
          model_f1(is) = model_f1(is) + pp 
          do jv = 1,nv
             if(out_var /= jv) then 
                model_f2(list(jv),jv,is) = model_f2(list(jv),jv,is) + pp
             end if
          end do
       end do
       
    end do
    
  end subroutine update_model_averages

  subroutine update_gradient(it)
    ! update cost-related variables: etot, cond_likelihood, ereg and gradient grd
    integer, intent(in) :: it
    integer :: dim
    real(kflt) :: etot0,de

    call unpack_parameters()

    ! save old cost function value
    etot0 = etot
    
    ! reset averages and cost before looping over data samples
    model_f1 = 0.0_kflt
    model_f2 = 0.0_kflt
    cond_likelihood = 0.0_kflt
    etot = 0.0_kflt
    ereg = - regularization_strength * (sum(prm(:ns)**2) + 0.5_kflt * sum(prm(ns+1:)**2))

    ! take averages over model distribution
    call update_model_averages()

    ! add regularization term to conditional likelihood
    etot = cond_likelihood + ereg

    ! delta cost 
    de = etot - etot0

    ! update gradient 
    model_f1 = model_f1 - data_f1 + 2.0_kflt * regularization_strength * vfields
    model_f2 = model_f2 - data_f2 + 2.0_kflt * 0.5_kflt * regularization_strength * vcouplings

    dim = ns*ns*nv
    grd(1:ns) = model_f1
    grd(ns+1:) = reshape(model_f2,(/dim/))
    
  end subroutine update_gradient

  subroutine unpack_parameters()
    ! unpack field and couplings after updating parameters

    vfields = prm(1:ns) 
    vcouplings = reshape(prm(ns+1:),(/ns,nv,ns/))
    
  end subroutine unpack_parameters

  subroutine model_collect_prm()

    couplings = 0.5_kflt * couplings

  end subroutine model_collect_prm
  

end module model
  
