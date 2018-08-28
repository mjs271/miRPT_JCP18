module ADRE_mod
implicit none

! these can be used to get the component list and selected output headers
! from the commented code blocks near line 160 in dolo_ADRE.f90
! type component_list
!     character(:), allocatable :: comp
! end type
! type selectout_list
!     character(:), allocatable :: head
! end type

! global constants
double precision, parameter :: pi = 4.0d0 * atan(1.0d0)
integer, parameter          :: sp = kind(1.0), dp = kind(1.d0)
integer, parameter          :: nspec = 10

contains

! subroutine to initialize the random number generator seed from clock time
subroutine init_random_seed()
    integer                            :: i, n, clock
    integer, dimension(:), allocatable :: seed

    call random_seed(size = n)
    allocate(seed(n))
    call system_clock(count = clock)
    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    call random_seed(put = seed)
    deallocate(seed)
end subroutine init_random_seed

! this subroutine writes the PHREEQC input file
! NOTE: much of this is hard-coded and not general
subroutine phreeqc_input(calcite_in, na_in, mg_in, ca_in, cl_in, co2_in)
    double precision, intent(in   ) :: calcite_in, na_in, mg_in, ca_in,&
                                       cl_in, co2_in
    character*1                     :: tab

    tab = char(9)

    open (unit=11, file='dolomite_chem.in', action='write')
    write (11, *) '# Brine-CO2-Calcite-Quartz system'
    write (11, *) 'SOLUTION 0 Brine'
    write (11, *) tab, 'pH 7.0'
    write (11, *) tab, 'units mol/L'
    write (11, *) tab, 'temp 60.000000'
    write (11, *) tab, 'pressure 98.6923 atm' ! = 100 bar
    write (11, '(A, A, f8.6)') tab, 'Na ', Na_in
    write (11, '(A, A, f8.6)') tab, 'Mg ', mg_in
    write (11, '(A, A, f8.6)') tab, 'Ca ', ca_in
    write (11, '(A, A, f8.6, A)') tab, 'Cl ', cl_in, ' charge'
    write (11, *) 'EQUILIBRIUM_PHASES 0'
    write (11, '(A, A, f9.6)') tab, 'CO2(g) 2.0 ', co2_in
    write (11, *) 'SAVE solution 0'
    write (11, *) 'END'
    write (11, *) 'SOLUTION 1 Domain'
    write (11, *) tab, 'pH 7.0 charge'
    write (11, *) tab, 'temp 60.000000'
    write (11, *) tab, 'pressure 98.6923 atm' ! = 100 bar
    write (11, *) tab, 'Cl 10.0'
    write (11, *) 'EQUILIBRIUM_PHASES 1'
    write (11, '(A, A, f9.6)') tab, 'Calcite 0.000000', calcite_in
    write (11, *) tab, 'Dolomite 0.000000 0.000000'
    write (11, *) tab, 'Quartz 0.000000 22.000000'
    write (11, *) 'SAVE solution 1'
    write (11, *) 'SELECTED_OUTPUT'
    write (11, *) tab, '-simulation false'
    write (11, *) tab, '-state false'
    write (11, *) tab, '-solution false'
    write (11, *) tab, '-distance false'
    write (11, *) tab, '-time false'
    write (11, *) tab, '-step false'
    write (11, *) tab, '-ph true' ! get pH output
    write (11, *) tab, '-pe false'
    write (11, *) tab, '-equilibrium_phases Calcite Dolomite'
        ! get calcite and dolomite concentrations
    write (11, *) 'END'
    close (unit=11, status='keep')
end subroutine phreeqc_input

! advection subroutine, using explicit upwinding
subroutine advect(mat, v, dx, dt, ncell)
    double precision, intent(inout) :: mat(:, :)
    double precision, intent(in   ) :: v(:), dx, dt
    integer,          intent(in   ) :: ncell
    double precision                :: vmat(ncell - 1, nspec)
    integer                         :: i

    do i = 1, nspec
        vmat(:, i) = v(1 : ncell - 1)
    enddo

    mat(2 : ncell, :) = mat(2 : ncell, :) - ((vmat * dt)/dx) *&
                        (mat(2 : ncell, :) - mat(1 : ncell - 1, :))

end subroutine advect

! diffusion subroutine using explicit Euler
subroutine diffuse(mat, D, dx, dt, ncell)
    double precision, intent(inout) :: mat(:, :)
    double precision, intent(in   ) :: D(:), dx, dt
    integer,          intent(in   ) :: ncell
    double precision                :: Dmat(ncell - 2, nspec), temp(nspec)
    integer                         :: i

    ! stop concentrations from diffusing in across boundary by eliminating
        ! concentration gradient between cells 1 and 2, then replace at end
    temp = mat(1, :)
    mat(1, :) = mat(2, :)
    do i = 1, nspec
        Dmat(:, i) = D(2 : ncell - 1)
    enddo

    mat(2 : ncell - 1, :) = mat(2 : ncell - 1, :) + ((Dmat * dt)/(dx**2)) *&
                            (mat(3 : ncell, :) - 2 * mat(2 : ncell - 1, :)&
                                + mat(1 : ncell - 2, :))
    mat(1, :) = temp
end subroutine diffuse

end module ADRE_mod
