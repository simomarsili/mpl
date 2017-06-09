! Copyright (C) 2015, 2016, Simone Marsili 
! All rights reserved.
! License: BSD 3 clause

module dvmlm_wrapper
  use nrtype
  use model, only: prm, grd, etot
  ! wrapper to dvmlm subroutine 
  implicit none
  private 

  public :: dvmlm_minimize

contains

  subroutine dvmlm_min(accuracy,ndim,mstep,task,wa)
    integer(I4B), intent(in) :: accuracy,ndim,mstep
    character(60), intent(inout) :: task
    real(DP), intent(in) :: wa(:)
    integer(I4B) :: isave(5)
    real(DP) :: dsave(24)
    real(DP) :: f
    real(DP) :: frtol
    real(DP) :: fatol
    real(DP) :: fmin
    external dvmlm

    ! set prms for minimization
    ! frtol: desired relative error
    ! fatol: desired absolute error
    select case(accuracy)
    case(0)
       ! low accuracy
       frtol = 1.0e-2_DP
       fatol = 1.0e-3_DP
    case(1)
       ! moderate accuracy
       frtol = 1.0e-4_DP
       fatol = 1.0e-6_DP
    case(2)
       ! high accuracy
       frtol = 1.0e-8_DP
       fatol = 1.0e-12_DP
    end select

    fmin = -1.9e30_DP
    
    f = - etot
    call dvmlm(ndim,prm,f,grd,frtol,fatol,fmin,task,mstep,&
         wa(1),wa(ndim*mstep+1),wa(2*ndim*mstep+1),&
         isave,dsave,wa(2*ndim*mstep+mstep+1),wa(2*ndim*mstep+mstep+ndim+1))

  end subroutine dvmlm_min

  subroutine dvmlm_minimize(accuracy,iter,totiter)
    use model, only: compute_pseudo_likelihood, model_put_myv
    integer(I4B),intent(in) :: accuracy
    integer(I4B),intent(out) :: iter,totiter
    integer(I4B) :: err
    integer(I4B) :: ndim,mstep
    character(60) :: task
    real(DP), allocatable :: wa(:)
    integer :: lwa
    
    iter = 0
    totiter = 0
    task = 'START'

    ndim = size(prm)
    ! this is the number of steps for hessian approximation
    mstep = 100
    lwa = 2*ndim*mstep + 2*ndim + mstep
    allocate(wa(lwa),stat=err)
    

    call compute_pseudo_likelihood(iter)
    do 
       if(totiter > 100) then 
          write(0,*) 'warning: totiter > 100'
          flush(0)
       end if
       call dvmlm_min(accuracy,ndim,mstep,task,wa)
       if(task(1:2) == 'FG') then 
          ! update etot and gradient for line search
          totiter = totiter + 1
          call compute_pseudo_likelihood(iter)
       elseif(task(1:4) == 'NEWX') then
          ! start new line search
          iter = iter + 1
       elseif(task(1:4) == 'WARN') then 
          write(0,*) 'warning ', iter
          flush(0)
       elseif(task(1:4) == 'CONV') then 
          ! compute final values for likelihood 
          call compute_pseudo_likelihood(iter)
          ! put my prms back in fields and couplings arrays 
          call model_put_myv()
          exit
       end if
    end do

    deallocate(wa)

  end subroutine dvmlm_minimize

end module dvmlm_wrapper
