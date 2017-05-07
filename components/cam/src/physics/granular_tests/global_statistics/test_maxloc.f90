program main

  real :: array(3)
  integer :: icol(1)

  array = (/1.,2.,3./)
  icol = maxloc( abs(array) )

  print*,  icol
end
