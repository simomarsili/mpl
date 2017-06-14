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
  use scrs,          only: print_mat, compute_scores
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
  real(kflt), allocatable :: fields(:,:) ! ns x nv
  real(kflt), allocatable :: couplings(:,:,:,:) ! ns x ns x nv x nv
  real(kflt), allocatable :: scores(:,:) ! nv x nv
  character(1) :: scores_format

  call units_initialize()

  ! get command line 
  call read_args(data_file,w_id,lambda,skip_gaps,accuracy,scores_format,err)
  if(err /= 0) then 
     write(0,100) 
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
  allocate(scores(nv,nv),stat=err)
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

  ! reorder prm array into fields and couplings 
  allocate(fields(ns,nv),couplings(ns,ns,nv,nv),stat=err)
  call reshape_prm(nv,ns,prm,fields,couplings)
  deallocate(prm,grd)

  ! compute scores and print
  call compute_scores(nv,ns,couplings,scores)
  call print_mat(scores,uscrs,scores_format,err)
  if (err /= 0) stop

100 format(/&
         'mpl                                                            '/&
         '                                                               '/&
         'Usage:                                                         '/&
         '    mpl [options] -i <file>                                    '/&
         '                                                               '/&
         'Options:                                                       '/&
         '-h, --help                                                     '/&
         '    Display this help and exit.                                '/&
         '                                                               '/&
         '-l, --lambda <regularization_parameter> float, optional        '/&
         '    Controls L_2 regularization strength.                      '/&
         '    [default: 0.01]                                            '/&
         '                                                               '/&
         '-w, --reweight <percentage_identity> - float, optional         '/&
         '    Reweight data using a percentage identity threshold.       '/&
         '    [default: no reweight.]                                    '/&
         '                                                               '/&
         '-g, --no_gap - optional                                        '/&
         '    Do not consider state "1" in the calculation of the scores.'/&
         '                                                               '/&
         '-a, --accuracy <accuracy_level> - integer, optional            '/&
         '    Larger values correspond to increased accuracy and slower  '/&
         '    covergence. Possible values are {0, 1, 2}.                 '/&
         '    [default: 1.]                                              '/&
         '                                                               '/&
         '--scores_format <scores_matrix_format> - string, optional      '/&
         '    Possible values are {"rcv", "array", "coordinate"}.        '/&
         '    "rcv": (row_index, column_index, value) format             '/&
         '    "coordinate", "array": Matrix Market (MM) formats. See:    '/&
         '    https://people.sc.fsu.edu/~jburkardt/data/mm/mm.html       '/&
         '    [default: "rcv"].                                          '/)
  
  
contains

  subroutine reshape_prm(nv,ns,prm,fields,couplings)
    ! reorder prm array into fields and couplings 
    implicit none
    integer, intent(in) :: nv,ns
    real(kflt), intent(in) :: prm(:,:)
    real(kflt), intent(out) :: fields(ns,nv)
    real(kflt), intent(out) :: couplings(ns,ns,nv,nv)
    integer :: err
    integer :: iv,jv,is,js,k1,k2,k
    
    k = 0
    do iv = 1,nv
       fields(:,iv) = prm(:ns,iv)
       k = ns
       do jv = 1,nv
          do js = 1,ns
             do is = 1,ns
                k = k + 1
                couplings(js,is,jv,iv) = prm(k,iv)
             end do
          end do
       end do
    end do
  end subroutine reshape_prm

end program mpl
