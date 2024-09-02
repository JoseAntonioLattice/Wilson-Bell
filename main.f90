program main
  use iso_fortran_env, only: dp => real64, i4 => int32 
  implicit none
 
  real(dp), parameter :: pi = acos(-1.0_dp)
  integer(i4), parameter :: d = 3
  real(dp), dimension(d) :: q
  integer(i4), dimension(d) :: n
  integer(i4), parameter :: L = 10
  integer(i4), dimension(d) :: l_vec
  real(dp) :: s, s_prod,rho_star, omega_star,a, b, a_star, z
  integer(i4) :: i,l1,l2,l3
  

  n = [1,1,1]
  omega_star = 0.0_dp
  do l1 = -L, L
     do l2 = -L, L
        do l3 = -L, L
           l_vec = [l1,l2,l3]
           q = 2*pi/L * n
           do i = 1, 3
              s_prod = 1.0_dp
              s = (sin(q(i)/2)/(q(i)/2 + pi*l_vec(i)))**2
              s_prod = s_prod*s
           end do
           omega_star = omega_star + s_prod/(norm2(q+2*pi*l_vec))**2
        end do
     end do
  end do

  z = 8.0_dp
  a = 8.0_dp
  b = 2**(-2.5_dp)
  a_star = a*(1.0_dp - 2**d * b**2)
  rho_star = a_star/(1.0_dp - (norm2(q))**2 * Omega_Star)
  
  print*, omega_star
  print*, rho_star
end program main
