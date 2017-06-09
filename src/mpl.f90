! Copyright (C) 2015, 2016, Simone Marsili 
! All rights reserved.
! License: BSD 3 clause

program mpl
  use kinds
  use units
  use command_line,  only: read_args
  use data,          only: nv, data_read
  use model,         only: model_initialize, model_set_myv, model_collect_prm
  use scrs,          only: compute_scores, print_scores
  use dvmlm_wrapper, only: dvmlm_minimize
  implicit none
  character(long_string) :: data_file
  character(long_string) :: scores_file
  real(kflt)             :: w_id
  real(kflt)             :: lambda
  logical                :: skip_gaps
  integer(kint)          :: accuracy
  integer(kint)          :: udata,uscrs
  integer(kint)          :: err,iv
  integer(kint)          :: niter,neval
  real(kflt_single)      :: finish,start,start_min,end_min
  character(long_string) :: syntax = 'syntax: mpl -i <data_file> -l <regularization_strength> [-w <weigths_file>] [-g]'
  

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
  call model_initialize(lambda)

  write(0,*) 'minimize...'
  call cpu_time(start_min)
  ! loop over features
  do iv = 1,nv
     call model_set_myv(iv,err)
     niter = 0
     call cpu_time(start)
     call dvmlm_minimize(accuracy,niter,neval)
     call cpu_time(finish)
     write(0,'(a,i5,a,2i5,a,f8.3,a)') ' variable ', iv, &
          '  converged (niter,neval) ', niter, neval, ' in ', finish-start, ' secs'
  end do
  flush(0)
  call cpu_time(end_min)
  
  write(0,*) 'minimization total (secs): ', end_min-start_min
  flush(0)
  call model_collect_prm()
  call compute_scores(skip_gaps)
  call print_scores(uscrs)

end program mpl
