! Copyright (C) 2015, 2016, Simone Marsili 
! All rights reserved.
! License: BSD 3 clause

module data
  use kinds
  implicit none
  private

  public :: nd ! number of data samples
  public :: nv ! number of variables
  public :: ns ! number of states

  public :: data_samples
  public :: w
  public :: neff

  public :: data_read

  integer :: nd   ! number of data samples
  integer :: nv   ! number of variables
  integer :: ns ! number of states per variable

  integer, allocatable :: data_samples(:,:)
  real(kflt), allocatable :: w(:)
  real(kflt) :: neff

contains

  subroutine data_initialize(udata)
    use constants, only: long_string
    use parser, only: parser_nfields
    integer, intent(in) :: udata
    integer :: nfields
    integer :: err
    character(long_string) :: line,newline

    ! read the first string
    read(udata,'(a)',iostat=err) line
    call parser_nfields(line,newline,nfields)
    rewind(udata)

    ! set nv
    nv = nfields

    ! count lines for nd
    nd = 0
    do
       read(udata,'(a)',iostat=err) line
       if(err < 0) exit
       nd = nd + 1
    end do
    rewind(udata)

    ! allocate memory for the data

    allocate(data_samples(nv,nd),stat=err)
    allocate(w(nd),stat=err)

    neff = nd

  end subroutine data_initialize

  subroutine data_read(udata,w_id)
    use constants, only: long_string
    use parser, only: parser_nfields
    integer,    intent(in) :: udata
    real(kflt), intent(in) :: w_id
    integer :: i
    integer :: nfields
    integer :: err
    character(long_string) :: line,newline

    call data_initialize(udata)

    do i = 1,nd
       if(mod(i,100) == 0)write(0,*) 'data: ', i
       read(udata,'(a)',iostat=err) line
       call parser_nfields(line,newline,nfields)
       read(newline,*,iostat=err) data_samples(:,i)
       if(err > 0) write(0,*) 'error: reading data'
    end do

    if (any(data_samples==0)) data_samples = data_samples + 1

    ns = maxval(data_samples)

    write(0,*) 'nd: ', nd
    write(0,*) 'nv: ', nv
    write(0,*) 'ns: ', ns
    flush(0)

    if(w_id > 1.E-10_kflt) then
       write(0,*) 'computing weights...'
       call data_reweight(w_id)
    else
       w = 1.0_kflt / real(nd)
    end if

  end subroutine data_read

  subroutine data_reweight(w_id)
    real(kflt), intent(in) :: w_id
    integer :: id,jd
    integer :: err
    integer, allocatable :: x(:),y(:)

    allocate(x(nv),y(nv),stat=err)

    w = 1.0_kflt
    do id = 1,nd-1
       if(mod(id,1000) == 0) then
          write(0,'(f5.1,a)') 100.0 * real(id) / real(nd),' %'
       end if
       x = data_samples(:,id)
       do jd = id+1,nd
          y = data_samples(:,jd)
          if(count(x == y) >= nint(nv * w_id * 0.01_kflt)) then
             w(id) = w(id) + 1.0_kflt
             w(jd) = w(jd) + 1.0_kflt
          end if
       end do
    end do

    w = 1.0_kflt / w
    neff = sum(w)
    w = w / neff
   write(0,*) 'neff: ', neff

  end subroutine data_reweight

end module data
