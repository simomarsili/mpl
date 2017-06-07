! Copyright (C) 2015, 2016, Simone Marsili 
! All rights reserved.
! License: BSD 3 clause

module scrs
  use nrtype
  implicit none
  private 

  public :: compute_scores
  public :: print_scores

  real(DP), allocatable, save :: scores(:,:)

contains

  subroutine compute_scores(skip_gaps)
    use data, only: nv,ns,neff
    use model, only: fields,couplings
    logical, intent(in) :: skip_gaps
    integer :: iv,jv,k,is,js
    real(DP) :: sums(nv),totsum
    integer(I4B) :: err
    real(DP) :: rpseudo=0.5_DP,qq,ll,sij
    real(DP) :: nrm,fpi(ns-1),fpj(ns-1),fpij(ns-1,ns-1),enti,entj,entij,mij

    ! at the very end of the run
    allocate(scores(nv,nv),stat=err)
    scores = 0.0_DP

    ! compute scores
    k = 0
    do jv = 1,nv-1
       do iv = jv+1,nv
          k = k + 1
          if (skip_gaps) then 
             scores(iv,jv) = frobenius(couplings(2:,2:,k))
          else
             scores(iv,jv) = frobenius(couplings(:,:,k))
          end if
          scores(jv,iv) = scores(iv,jv)
       end do
    end do
    
    !    call apc_correction(scores)

  contains

    real(DP) function frobenius(a)
      real(DP), intent(in) :: a(:,:)

      frobenius = sqrt(sum(a**2))

    end function frobenius

    subroutine apc_correction(a)
      real(DP), intent(inout) :: a(:,:)
      real(DP), allocatable :: sums(:)
      real(DP) :: totsum
      integer(I4B) :: err
      integer(I4B) :: i,j
      integer(I4B) :: n

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

  end subroutine compute_scores

  subroutine print_scores(uscrs)
    use units
    use data, only: nv,ns
    use model, only: fields,couplings
    integer(I4B), intent(in) :: uscrs
    integer(I4B) :: usc,uprm
    integer(I4B) :: err
    integer(I4B) :: iv,jv,is,js,k
    character(80) :: filename
    real(8) :: sij

    do iv = 1,nv-1
       do jv = iv+1,nv
          sij = scores(iv,jv)
          write(uscrs,'(i5,1x,i5,1x,f8.5)') iv,jv,sij
       end do
    end do

  end subroutine print_scores

end module scrs
