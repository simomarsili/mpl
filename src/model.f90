! Copyright (C) 2015-2017, Simone Marsili 
! All rights reserved.
! License: BSD 3 clause
! License: BSD 3 clause

module model
  use kinds
  use data, only: nv,ns,nd
  implicit none
  private

  public :: initialize_model
  public :: model_set_myv
  public :: fix_gauge
  public :: update_gradient
  public :: cond_likelihood
  public :: etot

  ! index of outcome variable (all the others nv - 1 variables are included in the set of explanatory variables)
  integer :: out_var
  
  ! arrays for fixed residue/sequence
  real(kflt), allocatable, save :: model_f1(:)       ! model single-variable frequencies (ns)
  real(kflt), allocatable, save :: model_f2(:,:,:)   ! model frequencies for pairs of variables (ns x nv x ns)
  real(kflt), allocatable, save :: data_f1(:)        ! data single-variable frequencies (ns)
  real(kflt), allocatable, save :: data_f2(:,:,:)    ! data frequencies for pairs of variables (ns x nv x ns)

  ! regularization parameters 
  real(kflt) :: regularization_strength=0.01_kflt ! regularization strength; the default is l2 with regularization_strength=0.01

  ! "cost" function-related variables
  real(kflt) :: cond_likelihood,ereg,etot

contains

  subroutine initialize_model(lambda)
    real(kflt), intent(in) :: lambda
    integer :: err

    regularization_strength = lambda

    allocate(data_f1(ns),stat=err)
    allocate(data_f2(ns,ns,nv),stat=err)
    allocate(model_f1(ns),stat=err)
    allocate(model_f2(ns,ns,nv),stat=err)

  end subroutine initialize_model

  subroutine model_set_myv(iv,vprm,grd,err) ! couplings
    use data, only: nd,data_samples,w
    integer, intent(in) :: iv
    real(kflt), intent(out) :: vprm(:)
    real(kflt), intent(out) :: grd(:)
    ! make couplings given out_var 
    ! must be called before looping on data
    integer :: id,jv,mys
    integer :: err
    integer, allocatable :: list(:)

    out_var = iv
    model_f1 = 0.0_kflt
    model_f2 = 0.0_kflt
    data_f1 = 0.0_kflt
    data_f2 = 0.0_kflt
    cond_likelihood = 0.0_kflt
    etot = 0.0_kflt
    ereg = 0.0_kflt
    vprm = 0.0_kflt
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
             data_f2(list(jv),mys,jv) = data_f2(list(jv),mys,jv) + w(id)
          end if
       end do
    end do
    deallocate(list)

  end subroutine model_set_myv

  subroutine fix_gauge(nv1,ns1,fields,couplings)
    integer, intent(in) :: nv1,ns1
    real(kflt), intent(inout) :: fields(ns1)
    real(kflt), intent(inout) :: couplings(ns1,ns1,nv1)
    integer :: jv,is,js
    real(kflt) :: mat(ns,ns),arr(ns),marr
    real(kflt) :: rsum(ns),csum(ns),totsum

    arr = fields
    marr = sum(arr) / real(ns)
    arr = arr - marr
    do jv = 1,nv
       mat = couplings(:,:,jv)
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
       couplings(:,:,jv) = mat
    end do
    fields = arr
    
  end subroutine fix_gauge

  subroutine update_model_averages(nv1,ns1,fields,couplings)
    use data, only: data_samples,w,nd
    integer, intent(in) :: nv1,ns1
    real(kflt), intent(in) :: fields(ns1)
    real(kflt), intent(in) :: couplings(ns1,ns1,nv1)
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
       conp = fields
       do jv = 1,nv
          if(out_var /= jv) then 
             do is = 1,ns
                conp(is) = conp(is) + couplings(list(jv),is,jv) 
             end do
          end if
       end do
       conp = exp(conp)
       
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
                model_f2(list(jv),is,jv) = model_f2(list(jv),is,jv) + pp
             end if
          end do
       end do
       
    end do
    
  end subroutine update_model_averages

  subroutine update_gradient(nv1,ns1,fields,couplings,grd1,grd2)
    ! update cost-related variables: etot, cond_likelihood, ereg and gradient grd
    integer, intent(in) :: nv1,ns1
    real(kflt), intent(in) :: fields(ns1)
    real(kflt), intent(in) :: couplings(ns1,ns1,nv1)
    real(kflt), intent(out) :: grd1(ns1)
    real(kflt), intent(out) :: grd2(ns1,ns1,nv1)
    real(kflt) :: etot0,de

    ! save old cost function value
    etot0 = etot
    
    ! reset averages and cost before looping over data samples
    model_f1 = 0.0_kflt
    model_f2 = 0.0_kflt
    cond_likelihood = 0.0_kflt
    etot = 0.0_kflt
    ereg = - regularization_strength * (sum(fields**2) + 0.5_kflt * sum(couplings**2))

    ! take averages over model distribution
    call update_model_averages(nv1,ns1,fields,couplings)

    ! add regularization term to conditional likelihood
    etot = cond_likelihood + ereg

    ! delta cost 
    de = etot - etot0

    ! update gradient 
    grd1 = model_f1 - data_f1 + 2.0_kflt * regularization_strength * fields
    grd2 = model_f2 - data_f2 + 2.0_kflt * 0.5_kflt * regularization_strength * couplings
    
  end subroutine update_gradient

end module model
  
