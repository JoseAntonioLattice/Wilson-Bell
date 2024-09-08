program main
  use iso_fortran_env, only: dp => real64, i4 => int32 
  implicit none
 
  real(dp), parameter :: pi = acos(-1.0_dp)
  integer(i4) :: i, j, k

  do i = 0, 2
     do j = i, 2
        do k = j, 2
           print('("[",3i0,"] ", f18.15)'), i,j,k, fourier([i,j,k],10,6.0_dp)
        end do
     end do
  end do
contains

  function omega_star(n,L)
    real(dp) :: omega_star
    integer(i4), dimension(3), intent(in) :: n
    integer(i4), intent(in) :: L
    integer(i4) :: i,l1,l2,l3
    integer(i4), dimension(3) :: l_vec
    real(dp) :: s_prod,q(3)
    integer, parameter :: m = 100

    q = 2*pi*n/L
    omega_star = 0.0_dp
    do l1 = -m, m
       do l2 = -m, m
          do l3 = -m, m
             l_vec = [l1,l2,l3]
             s_prod = 1.0_dp
             do i = 1, 3
                if( n(i) == 0 .and. l_vec(i) == 0) cycle
                s_prod = s_prod*(sin(q(i)/2)/(q(i)/2 + pi*l_vec(i)))**2
             end do
             omega_star = omega_star + s_prod/(norm2(q + 2*pi*l_vec))**2
          end do
       end do
    end do
  end function omega_star

  function rho_star(n,L,a_star)
    real(dp) :: rho_star
    real(dp), intent(in)  :: a_star
    integer(i4), dimension(3), intent(in) :: n
    integer(i4), intent(in) :: L

    if( n(1) == 0 .and. n(2) == 0 .and. n(3) == 0 )then
       rho_star = 0.0_dp
    else
       rho_star = a_star/(1.0_dp + a_star * omega_star(n,L))
    end if
    
  end function rho_star

  function fourier(r,L,a_star)
    real(dp), intent(in)  :: a_star
    integer(i4), intent(in) :: L
    integer(i4) :: n(3),n1,n2,n3
    real(dp) :: fourier, q(3)
    integer(i4), dimension(3), intent(in) :: r

    fourier = 0.0_dp
    do n1 = 0, L - 1
       do n2 = 0, L - 1
          do n3 = 0, L - 1
             n = [n1,n2,n3]
             q = 2*pi/L * n
             fourier = fourier + cos(dot_product(q,r)) * rho_star(n,L,a_star) 
          end do
       end do
    end do
    fourier = fourier / L**3
  end function fourier

  
end program main
