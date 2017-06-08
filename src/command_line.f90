! Copyright (C) 2015, 2016, Simone Marsili 
! All rights reserved.
! License: BSD 3 clause

module command_line
  implicit none
  private
  public :: command_line_read
  character(2) :: opts(4) = ['-i','-w','-l','-g']
  character(1000) :: syntax = 'syntax: mpl -i <data_file> -l <regularization_strength> [-w <weigths_file>] [-g]'

contains

  subroutine command_line_read(nerrs,data_file,wid,regu,lambda,skip_gaps)
    use units, only: max_string_length
    use nrtype
    integer, intent(out) :: nerrs
    character(max_string_length), intent(out) :: data_file
    real(DP), intent(out) :: wid
    integer(I4B), intent(out) :: regu
    real(DP), intent(out) :: lambda
    logical, intent(out) :: skip_gaps
    character(max_string_length) :: cmd
    integer :: nargs
    character(max_string_length) :: arg
    integer :: iarg

    call get_command(cmd)
    nargs = command_argument_count()
    write(0,*) trim(cmd)

    iarg = 1
    nerrs = 0
    data_file = ''
    lambda = 0.01_DP
    wid = 0.0_DP
    regu = 2
    skip_gaps = .false.
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
             ! weights file
             iarg = iarg + 1
             call get_command_argument(iarg,arg)
             read(arg,*) wid
          case('-l')
             iarg = iarg + 1
             call get_command_argument(iarg,arg)
             regu = 2
             read(arg,*) lambda
          case('-g')
             ! remove contrib. from state 1 to the final scores
             skip_gaps = .true.
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
       write(0,*) trim(syntax)
       stop
    end if
    if (lambda < 1.e-5_DP) then
       write(0,*) 'error ! regularization parameter too small (< 1.e-5)'
       write(0,*) trim(syntax)
       stop
    end if

  end subroutine command_line_read

end module command_line
