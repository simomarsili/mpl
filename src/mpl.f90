! Copyright (C) 2015-2017, Simone Marsili 
! All rights reserved.
! License: BSD 3 clause
! License: BSD 3 clause

program mpl
  use kinds
  use constants,     only: long_string
  use units
  use command_line,  only: read_args
  use data,          only: nd,nv,ns,data_read
  use model,         only: initialize_model, model_set_myv,fix_gauge
  use scrs,          only: compute_scores, print_scores
  use dvmlm_wrapper, only: dvmlm_minimize
  implicit none
  character(long_string) :: data_file
  character(long_string) :: scores_file
  real(kflt)             :: w_id
  real(kflt)             :: lambda
  logical                :: skip_gaps
  integer                :: accuracy
  integer                :: udata,uscrs
  integer                :: err,iv
  integer                :: niter,neval
  real(kflt_single)      :: finish,start,start_min,end_min
  character(long_string) :: syntax = 'syntax: mpl -i <data_file> -l <regularization_strength> [-w <weigths_file>] [-g]'
  real(kflt), allocatable :: prm(:,:) ! 1D array of parameters (nv, ns + ns x nv x ns, nv)
  real(kflt), allocatable :: grd(:) ! 1D gradient array (ns + ns x nv x ns)
  logical, parameter :: symmetrize=.true.
  real(kflt), allocatable :: fields(:,:) ! nv x ns
  real(kflt), allocatable :: couplings(:,:,:,:) ! nv x nv x ns x ns

  call units_initialize()

  ! get command line 
  call read_args(data_file,w_id,lambda,skip_gaps,accuracy,err)
  if(err /= 0) then 
     write(0,*) 'error: check syntax'
     write(0,*) trim(syntax)
     stop
  end if

  ! open data file
  call units_open(data_file,udata,'O',err)
  if(err /= 0) then 
     write(0,*) 'error: cant find data file'
     stop
  end if
  
  ! open scores file
  scores_file = trim(data_file)//'.scores'
  call units_open(scores_file,uscrs,'U',err)

  write(0,*) 'reading data...'
  call data_read(udata,w_id)

  write(0,*) 'initialize...'
  ! allocate parameters and gradient
  allocate(prm(ns + ns*ns*nv,nv),stat=err)
  allocate(grd(ns + ns*ns*nv),stat=err)
  call initialize_model(lambda)

  write(0,*) 'minimize...'
  call cpu_time(start_min)
  ! loop over features
  do iv = 1,nv
     call model_set_myv(iv,prm(:,iv),grd,err)
     niter = 0
     call cpu_time(start)
     call dvmlm_minimize(prm(:,iv),grd,size(prm(:,iv)),accuracy,niter,neval)
     call cpu_time(finish)
     write(0,'(a,i5,a,2i5,a,f8.3,a)') ' variable ', iv, &
          '  converged (niter,neval) ', niter, neval, ' in ', finish-start, ' secs'
     call fix_gauge(nv,ns,prm(:ns,iv),prm(ns+1:,iv))
  end do
  flush(0)
  call cpu_time(end_min)
  
  write(0,*) 'minimization total (secs): ', end_min-start_min
  flush(0)

  !call compute_scores(skip_gaps)
  !call compute_scores(nv,ns,prm,skip_gaps,symmetrize)
  !call print_scores(uscrs)

  ! reorder prm array into fields and couplings 
  allocate(fields(nv,ns),couplings(nv,nv,ns,ns),stat=err)
  call reshape_prm(nv,ns,prm,fields,couplings)
  deallocate(prm,grd)

  ! compute scores and print
  call go_scores(nv,ns,couplings,uscrs)

contains

  subroutine reshape_prm(nv,ns,prm,fields,couplings)
    ! reorder prm array into fields and couplings 
    implicit none
    integer, intent(in) :: nv,ns
    real(kflt), intent(in) :: prm(:,:)
    real(kflt), intent(out) :: fields(nv,ns)
    real(kflt), intent(out) :: couplings(nv,nv,ns,ns)
    integer :: err
    integer :: iv,jv,is,js,k1,k2,k

    k = 0
    do iv = 1,nv
       fields(iv,:) = prm(:ns,iv)
       k = ns
       do is = 1,ns
          do jv = 1,nv
             do js = 1,ns
                k = k + 1
                couplings(iv,jv,is,js) = prm(k,iv)
             end do
          end do
       end do
    end do
    
  end subroutine reshape_prm

  subroutine go_scores(nv,ns,couplings,uscrs)
    implicit none
    integer(kint), intent(in) :: nv,ns
    real(kflt), intent(in) :: couplings(nv,nv,ns,ns)
    integer, intent(in) :: uscrs
    integer :: iv,jv
    real(kflt) :: sij
    real(kflt) :: scores(nv,nv)
    real(kflt), allocatable :: sums(:)
    real(kflt) :: totsum

    scores = 0.0_kflt
    do iv = 1,nv-1
       do jv = iv+1,nv
          scores(iv,jv) = sum((couplings(iv,jv,:,:) + transpose(couplings(jv,iv,:,:)))**2)
          scores(jv,iv) = scores(iv,jv)
       end do
    end do
    scores = 0.5_kflt * sqrt(scores)
    
    allocate(sums(nv),stat=err)
    do iv = 1,nv
       sums(iv) = sum(scores(:,iv))
    end do
    totsum = sum(sums)

    do jv = 1,nv
       do iv = jv,nv
          scores(iv,jv) = scores(iv,jv) - sums(iv)*sums(jv)/totsum
          scores(jv,iv) = scores(iv,jv)
       end do
    end do

    deallocate(sums)

    do iv = 1,nv-1
       do jv = iv+1,nv
          sij = scores(iv,jv)
          write(uscrs,'(i5,1x,i5,1x,f8.5)') iv,jv,sij
       end do
    end do
    
  end subroutine go_scores

end program mpl
