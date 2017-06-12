! Copyright (C) 2015-2017, Simone Marsili
! All rights reserved.
! License: BSD 3 clause
! License: BSD 3 clause

module scrs
  use kinds
  implicit none
  private

  public :: compute_scores2
  public :: print_scores

  real(kflt), allocatable, save :: scores(:,:)

contains

  real(kflt) function frobenius(a)
    real(kflt), intent(in) :: a(:,:)

    frobenius = sqrt(sum(a**2))

  end function frobenius

  subroutine apc_correction(a)
    real(kflt), intent(inout) :: a(:,:)
    real(kflt), allocatable :: sums(:)
    real(kflt) :: totsum
    integer :: err
    integer :: i,j
    integer :: n

    n = size(a,1)
    allocate(sums(n),stat=err)

    do i = 1,n
       sums(i) = sum(a(:,i))
    end do
    totsum = sum(sums)

    do j = 1,n
       do i = j,n
          a(i,j) = a(i,j) - sums(i)*sums(j)/totsum
          a(j,i) = a(i,j)
       end do
    end do

    deallocate(sums)

  end subroutine apc_correction

  subroutine compute_scores2(nv,ns,coups,skip_gaps)
    integer, intent(in) :: nv,ns
    real(kflt), intent(in) :: coups(ns,nv,ns,nv)
    logical, intent(in) :: skip_gaps
    integer :: iv,jv,is,js,k
    integer :: err
    real(kflt) :: fuffa

    ! at the very end of the run
    allocate(scores(nv,nv),stat=err)
    scores = 0.0_kflt

    do jv = 1,nv-1
       do iv = jv+1,nv
          fuffa = 0.0_kflt
          if (skip_gaps) then
             do js = 2,ns
                do is = 2,ns
                   fuffa = fuffa + (coups(is,iv,js,jv) + coups(js,jv,is,iv))**2
                end do
             end do
          else
             do js = 1,ns
                do is = 1,ns
                   fuffa = fuffa + (coups(is,iv,js,jv) + coups(js,jv,is,iv))**2
                end do
             end do
          end if
          fuffa = 0.5_kflt * sqrt(fuffa)
          scores(iv,jv) = fuffa
          scores(jv,iv) = fuffa
       end do
    end do

    call apc_correction(scores)

  end subroutine compute_scores2

  subroutine print_scores(uscrs)
    use units
    use data, only: nv
    integer, intent(in) :: uscrs
    integer :: iv,jv
    real(8) :: sij

    do iv = 1,nv-1
       do jv = iv+1,nv
          sij = scores(iv,jv)
          write(uscrs,'(i5,1x,i5,1x,f8.5)') iv,jv,sij
       end do
    end do

  end subroutine print_scores

end module scrs
