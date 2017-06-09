! Copyright (C) 2015, 2016, Simone Marsili 
! All rights reserved.
! License: BSD 3 clause

module command_line
  implicit none
  private
  public :: command_line_read
  character(2) :: opts(5) = ['-i','-w','-l','-g','-a']

contains

  subroutine command_line_read(data_file,wid,lambda,skip_gaps,accuracy,nerrs)
    use units, only: long_string
    use nrtype
    character(long_string), intent(out) :: data_file
    real(DP),                     intent(out) :: wid
    real(DP),                     intent(out) :: lambda
    logical,                      intent(out) :: skip_gaps
    integer(I4B),                 intent(out) :: accuracy
    integer,                      intent(out) :: nerrs
    integer :: iarg,nargs
    integer :: err
    character(long_string) :: cmd,arg


    call get_command(cmd)
    nargs = command_argument_count()
    write(0,*) trim(cmd)

    iarg = 1
    nerrs = 0
    data_file = ''
    lambda = 0.01_DP
    wid = 0.0_DP
    skip_gaps = .false.
    accuracy = 1
    do while(iarg <= nargs)
       call get_command_argument(iarg,arg)
       if (any(opts == trim(arg))) then
          select case(trim(arg))
          case('-i','--input')
             ! input file
             iarg = iarg + 1
             call get_command_argument(iarg,arg)
             data_file = arg
             if(len_trim(data_file) == 0) then
                nerrs = nerrs + 1
                write(0,*) 'error: no data file name'
             end if
             if(data_file(1:1) == '-') then
                nerrs = nerrs + 1
                write(0,*) 'error: no data file name'
             end if
          case('-w','--weigths')
             iarg = iarg + 1
             call get_command_argument(iarg,arg)
             read(arg,*,iostat=err) wid
             if ( err/=0 ) then
                write(0,*) "error ! check %id (reweight)"
                nerrs = nerrs + 1
             end if
          case('-l')
             iarg = iarg + 1
             call get_command_argument(iarg,arg)
             read(arg,*,iostat=err) lambda
             if ( err/=0 ) then
                write(0,*) "error ! check lambda"
                nerrs = nerrs + 1
             end if
          case('-g')
             ! remove contrib. from state 1 to the final scores
             skip_gaps = .true.
          case('-a')
             iarg = iarg + 1
             call get_command_argument(iarg,arg)
             read(arg,*,iostat=err) accuracy
             if ( err/=0 ) then
                write(0,*) "error ! check accuracy"
                nerrs = nerrs + 1
             end if
          end select
       else
          write(0,'(a,1x,a)') 'error: unknown option',trim(arg)
          nerrs = nerrs + 1
       end if
       iarg = iarg + 1
    end do

    ! check command line
    if (data_file == '') then
       write(0,*) 'error ! missing input data file'
       nerrs = nerrs + 1
    end if
    if (lambda < 1.e-5_DP) then
       write(0,*) 'error ! regularization parameter too small (< 1.e-5)'
       nerrs = nerrs + 1
    end if
    if (.not. any([0,1,2] == accuracy)) then
       write(0,*) 'error ! possible accuracy values are 0,1,2'
       nerrs = nerrs + 1
    end if

  end subroutine command_line_read

end module command_line
