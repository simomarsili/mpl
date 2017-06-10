! Copyright (C) 2015, 2016, Simone Marsili 
! All rights reserved.
! License: BSD 3 clause

module command_line
  use constants, only: long_string
  implicit none
  private
  public :: read_args
  character(long_string) :: syntax = 'syntax: mpl -i <data_file> -l <regularization_strength> [-w <weigths_file>] [-g]'
  character(1), parameter :: nl=achar(10)
  character(long_string) :: usage = &
       'SYNTAX'//nl//&
       nl//&
       '    mpl -i <data_file> [-l <regularization_parameter>] [-w <percentage_identity>] [-g] [-a <accuracy_level>]'//nl//&
       '    mpl -h                                                                                                  '//nl//&
       nl//&
       'OPTIONS'//nl//&
       '-i, --input <data_file>                                                                                     '//nl//&
       '    Data file                                                                                               '//nl//&
       nl//&
       '-l, --lambda <regularization_parameter> - float, optional                                                   '//nl//&
       '    Controls L_2 regularization strength. Default: 0.01                                                     '//nl//&
       nl//&
       '-w, --reweight <percentage_identity> - float, optional                                                      '//nl//&
       '    Percentage identity threshold to be used when reweighting data. Default: no reweight.                   '//nl//&
       nl//&
       '-g, --no_gap - optional                                                                                     '//nl//&
       '    Do not take into account state 1 into the calculation of interaction scores.                            '//nl//&
       nl//&
       '-a, --accuracy <accuracy_level> - integer, optional                                                         '//nl//&
       '    Controls thresholds for convergence. Possible values are {0, 1, 2}. Default: 1.                         '//nl//&
       '    Larger values correspond to accurate solutions but slower covergence.                                   '//nl

contains

  subroutine read_args(data_file,wid,lambda,skip_gaps,accuracy,nerrs)
    use kinds
    character(long_string), intent(out) :: data_file
    real(kflt),             intent(out) :: wid
    real(kflt),             intent(out) :: lambda
    logical,                intent(out) :: skip_gaps
    integer(kint),          intent(out) :: accuracy
    integer,                intent(out) :: nerrs
    integer :: iarg,nargs
    integer :: err
    character(long_string) :: cmd,arg


    call get_command(cmd)
    nargs = command_argument_count()
    if (nargs == 0) then
       write(0,'(a)') trim(usage)
       stop
    end if

    iarg = 1
    nerrs = 0
    data_file = ''
    lambda = 0.01_kflt
    wid = 0.0_kflt
    skip_gaps = .false.
    accuracy = 1
    do while(iarg <= nargs)
       call get_command_argument(iarg,arg)
       select case(trim(arg))
       case('-h')
          write(0,'(a)') trim(usage)
          stop
       case('-i','--input')
          ! input file
          iarg = iarg + 1
          call get_command_argument(iarg,arg)
          data_file = arg
          if(len_trim(data_file) == 0) then
             nerrs = nerrs + 1
             write(0,*) 'error: check data file'
          end if
          if(data_file(1:1) == '-') then
             nerrs = nerrs + 1
             write(0,*) 'error: check data file'
          end if
       case('-w','--reweight')
          iarg = iarg + 1
          call get_command_argument(iarg,arg)
          read(arg,*,iostat=err) wid
          if ( err/=0 ) then
             write(0,*) "error ! check %id (reweight)"
             nerrs = nerrs + 1
          end if
       case('-l','--lambda')
          iarg = iarg + 1
          call get_command_argument(iarg,arg)
          read(arg,*,iostat=err) lambda
          if ( err/=0 ) then
             write(0,*) "error ! check lambda"
             nerrs = nerrs + 1
          end if
       case('-g','--no_gap')
          ! remove contrib. from state 1 to the final scores
          skip_gaps = .true.
       case('-a','--accuracy')
          iarg = iarg + 1
          call get_command_argument(iarg,arg)
          read(arg,*,iostat=err) accuracy
          if ( err/=0 ) then
             write(0,*) "error ! check accuracy"
             nerrs = nerrs + 1
          end if
       case default
          write(0,'(a,1x,a)') 'error: unknown option',trim(arg)
          nerrs = nerrs + 1             
       end select
       iarg = iarg + 1
    end do

    ! check command line
    if (data_file == '') then
       write(0,*) 'error ! missing input data file'
       nerrs = nerrs + 1
    end if
    if (lambda < 1.e-5_kflt) then
       write(0,*) 'error ! regularization parameter too small (< 1.e-5)'
       nerrs = nerrs + 1
    end if
    if (.not. any([0,1,2] == accuracy)) then
       write(0,*) 'error ! possible accuracy values are 0,1,2'
       nerrs = nerrs + 1
    end if

  end subroutine read_args

end module command_line
