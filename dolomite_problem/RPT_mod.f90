module RPT_mod
use kdtree2_module
implicit none

! these can be used to get the component list and selected output headers
! from the commented code blocks near line 203 in RPT_dolo.f90
! type component_list
!     character(:), allocatable :: comp
! end type
! type selectout_list
!     character(:), allocatable :: head
! end type

integer, parameter          :: sp = kind(1.0), dp = kind(1.0d0)
double precision, parameter :: pi = 4.0d0 * atan(1.0d0), kappaD = 0.5d0,&
                               kappa = (1.0d0 - kappaD), kappaM = 0.5d0,&
                               kappaI = 1.0d0 - kappaM
integer, parameter          :: nspec = 10, numsd = 5, num_alloc = 1000
    ! NOTE: nspec (the number of chemical species in the model) is hard-coded for
        ! this specific problem
    ! numsd is distances over which particle reactions are considered
    ! num_alloc is the maximum number of nearby particles expected to be found

! mobile particle type
! NOTE: these are only for a 1D problem
type mparticle
    double precision :: loc ! real-valued spatial location
    double precision :: concs(nspec) ! vector of chemical concentrations
    logical          :: active
        ! indicates whether particle is active and within the domain
    integer          :: bin
        ! used to indicate which spatial grid point a particle resides within
end type

! immobile particle type
type iparticle
    double precision :: loc
    double precision :: concs(nspec)
end type

! a couple of derived types for the kD tree search
! holds indices of nearby particles
type index_array
    integer, allocatable :: indices(:)
end type
! holds the distances to the corresponding particle held by index_array
type dist_array
    double precision, allocatable :: dists(:)
end type

contains

! subroutine to initialize the random number generator seed from clock time
subroutine init_random_seed()
    integer              :: i, n, clock
    integer, allocatable :: seed(:)

    call random_seed(size = n)
    allocate (seed(n))
    call system_clock(count = clock)
    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    call random_seed(put = seed)
    deallocate(seed)
end subroutine init_random_seed

! this generates the PHREEQC input file for this specific problem
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

! moves active particles via advection
subroutine advect(p, v, dt, alive)
    type(mparticle),  intent(inout) :: p(:) ! mobile particle array
    double precision, intent(in   ) :: v(:), dt ! velocity grid and time step
    integer,          intent(in   ) :: alive(:)
        ! array of indices of active particles

    ! use the velocity of the relevant cell, based on bin value
    p(alive)%loc = p(alive)%loc + v(p(alive)%bin) * dt
end subroutine advect

! moves active particles via diffusion
subroutine diffuse(p, np, D, dt, alive)
    type(mparticle),  intent(inout) :: p(:) ! mobile particle array
    double precision, intent(in   ) :: D(:), dt
        ! diffusion coefficient grid and time step
    integer,          intent(in   ) :: np, alive(:)
    ! number and array of indices of active particles
    double precision                :: normvec(np)
        ! vector which will hold Normal(0, 1) values

    ! call N(0, 1) generator
    call box_mullerp(np, normvec)

    ! use the diffusion coeff of the relevant cell, based on bin value
    p(alive)%loc = p(alive)%loc + sqrt(2.0d0 * kappaD * D(p(alive)%bin) * dt) * normvec
end subroutine diffuse

! reflective lower boundary
subroutine reflectlow(p, low, alive)
    type(mparticle),  intent(inout) :: p(:) ! mobile particle array
    double precision, intent(in   ) :: low ! lower spatial boundary
    integer,          intent(in   ) :: alive(:)
        ! array of indices of active particles

    ! if particle has exited lower boundary, flip the negative location
    ! to the same positive value
    where (p(alive)%loc < low) p(alive)%loc = -p(alive)%loc
end subroutine reflectlow

! absorbing upper boundary
subroutine absorbhigh(p, high, alive, na, n)
    type(mparticle),  intent(inout) :: p(:) ! mobile particle array
    double precision, intent(in   ) :: high ! upper spatial boundary
    integer,          intent(in   ) :: alive(:), na
    integer,          intent(  out) :: n
    logical                         :: gone(na)

    ! if particle has exited upper boundary, make the particle inactive
    ! also, make bin and loc -999 to catch any errors
    gone = p(alive)%loc > high
    n = count(gone)
    where (gone) p(alive)%active = .false.
    where (gone) p(alive)%bin = -999
    where (gone) p(alive)%loc = -999
end subroutine absorbhigh

! since PHREEQCRM can't accept the particle array as input, this subroutine
! assigns the values in the 2D concs array to its corresponding immobile
! particle
subroutine concs_to_iparts(c, p, np)
    double precision, intent(in   ) :: c(:, :)
        ! 2D concentration array used by PHREEQCRM
    type(iparticle),  intent(inout) :: p(:) ! immobile particle array
    integer,          intent(in   ) :: np ! number of immobile particles
    integer                         :: i ! iteration variable

    do i = 1, np
        p(i)%concs = c(i, :)
    enddo
end subroutine concs_to_iparts

! this does the opposite of above
subroutine iparts_to_concs(c, p, np)
    double precision, intent(inout) :: c(:, :)
        ! 2D concentration array used by PHREEQCRM
    type(iparticle),  intent(in   ) :: p(:) ! immobile particle array
    integer,          intent(in   ) :: np ! number of immobile particles
    integer                         :: i ! iteration variable

    do i = 1, np
        c(i, :) = p(i)%concs
    enddo
end subroutine iparts_to_concs

! this subroutine builds the distance matrix for all mobile-immobile pairwise
! distances less than the cutoff radius
subroutine build_Distmat(ip, mp, na, alive, ni, D, dt, Distmat)
    type(iparticle),  intent(in   ) :: ip(:) ! immobile particle array
    type(mparticle),  intent(in   ) :: mp(:) ! mobile particle array
    integer,          intent(in   ) :: na, alive(:), ni
        ! number and array of indices of active mobile particles and number
        ! of immobile particles
    double precision, intent(in   ) :: D(:), dt
    double precision, intent(  out) :: Distmat(ni, na)
        ! diffusion coefficient array, time step, and domain length
    type(kdtree2), pointer          :: tree ! this is the KD tree
    integer                         :: ntot, dim = 1, bindex, i, j
        ! total number of particles (active mobile + immobile)
        !****Note: hard coded one spatial dimension
    real(kdkind)                    :: locs(na + ni), r2
        ! array holding locations of immobile and active immobile particles
        ! and value of squared search radius for KD search
    type(index_array), allocatable  :: closeguys(:)
        ! this holds the indices of nearby particles
    type(dist_array), allocatable   :: close_dists(:)
        ! this holds the distances to the corresponding nearby particle

    Distmat = 0.0d0

    ! calculate total number of particles to be considered for mass balance
    ntot = na + ni
    ! build locs array--immobile particles will be at the beginning of the
    ! array, and mobile will be at the end
    locs(1 : ni) = real(ip%loc, kdkind)
    locs(ni + 1 : ntot) = real(mp(alive)%loc, kdkind)
    ! calculate interaction distance to be numsd standard deviations of the
    ! Brownian Motion process--r2 is this distance squared
    ! ****NOTE: numsd is a global variable that is hard-coded above
    r2 = (real(numsd, kdkind) * sqrt(4.0_kdkind * maxval(real(D, kdkind)) *&
                                     real(dt, kdkind)))**2

    ! build the KD tree and search it
    ! ****NOTE: num_alloc is a global variable that is hard-coded above
    call maketree(tree, dim, ntot, locs)

    allocate (closeguys(ni), close_dists(ni))
    ! this finds the closest mobile particles to each immobile particle
    call search(1, ni, tree, r2, num_alloc, closeguys, close_dists)
    ! NOTE: this search returns the SQUARED distance between two points
        ! also, the point itself is included in the closeguys list
    call kdtree2_destroy(tree)

    ! loop over immobile particles to build distance matrix
    do i = 1, ni ! immobile particle loop
        do j = 1, size(closeguys(i)%indices) ! mobile particle loop
            ! this is the mobile particle loop
            ! note that the closeguys array is indexed to the loc array,
            ! and thus its true index in the mp array is calculated below

            ! current mobile particle's index in locs array
            if (closeguys(i)%indices(j) <= ni) cycle
                ! if B is an immobile particle, skip this loop
                ! also prevents distance with self
            bindex = closeguys(i)%indices(j) - ni
            ! NOTE: Distmat is indexed to the mobile alive array
            Distmat(i, bindex) = close_dists(i)%dists(j)
        enddo
    enddo
    deallocate (closeguys, close_dists)
end subroutine build_Distmat

! algorithm set forth in mobile-immobile JCP paper using explicit matrix forward
! solve
subroutine immobile2mobile(ip, mp, na, alive, ni, D, dt, nc, Distmat)
    type(iparticle),  intent(inout) :: ip(:) ! immobile particle array
    type(mparticle),  intent(inout) :: mp(:) ! mobile particle array
    integer,          intent(in   ) :: na, alive(:), ni, nc
        ! number and array of indices of active mobile particles and number
        ! of immobile particles
    double precision, intent(in   ) :: D(:), dt, Distmat(ni, na)
        ! diffusion coefficient array, time step, and domain length
    integer                         :: i
    double precision                :: WImat(na, ni), denom(na, ni), DTmat(na, ni),&
                                       colsum(ni)
    ! denom is used in v(s) calculation but is pre-calculated for efficiency
    ! v_s is encounter density for a mobile/immobile particle pair

    DTmat = transpose(Distmat)

    do i = 1, ni
        denom(:, i) = -kappa * kappaI * 4.0d0 *  D(mp(alive)%bin) * dt
    enddo

    where (DTmat /= 0.0d0) WImat = exp(DTmat / denom)

    colsum = sum(WImat, 1)

    do i = 1, ni
        WImat(:, i) = WImat(:, i) / colsum(i)
    enddo

    do i = 1, nc
        mp(alive)%concs(i) = matmul(WImat, ip%concs(i))
    enddo
end subroutine immobile2mobile

! this does the opposite of the above
subroutine mobile2immobile(ip, mp, na, alive, ni, D, dt, nc, Distmat)
    type(iparticle),  intent(inout) :: ip(:) ! immobile particle array
    type(mparticle),  intent(inout) :: mp(:) ! mobile particle array
    integer,          intent(in   ) :: na, alive(:), ni, nc
        ! number and array of indices of active mobile particles and number
        ! of immobile particles
    double precision, intent(in   ) :: D(:), dt, Distmat(ni, na)
        ! diffusion coefficient array, time step, and domain length
    integer                         :: i
    double precision                :: WMmat(ni, na), denom(ni, na)
    ! denom is used in v(s) calculation but is pre-calculated for efficiency
    ! v_s is encounter density for a mobile/immobile particle pair

    do i = 1, ni
        denom(i, :) = -kappa * kappaM * 4.0d0 * D(mp(alive)%bin) * dt
    enddo

    where (Distmat /= 0.0d0) WMmat = exp(Distmat / denom)

    do i = 1, na
        WMmat(:, i) = WMmat (:, i) / sum(WMmat(:, i), 1)
    enddo

    do i = 1, nc
        ip%concs(i) = matmul(WMmat, mp(alive)%concs(i))
    enddo
end subroutine mobile2immobile

! this builds a KD tree
subroutine maketree(tree2, d, n, locs)
    type(kdtree2), pointer, intent(  out) :: tree2 ! this is the KD tree
    integer,                intent(in   ) :: d, n
        ! number of spatial dimensions, number of particles
    real(kdkind),           intent(in   ) :: locs(d, n)
        ! location array for particles, with dimension d x n (number of
        ! spatial dimensions x number of particles)

    ! build the tree
    tree2 => kdtree2_create(locs, dim=d, sort=.false., rearrange=.true.)
        ! currently don't see a need to sort, as false is quicker, while
        ! rearrange = true is quicker
end subroutine maketree

! this searches an already built KD tree
subroutine search(start, end, tree, r2, num_alloc, closeguys, close_dists)
    integer,                        intent(in   ) :: start, end, num_alloc
        ! number of particles and how large to to preallocate results array
        ! within KD tree module
    type(kdtree2), pointer,         intent(in   ) :: tree ! the KD tree
    real(kdkind),                   intent(in   ) :: r2 ! squared search radius
    type(index_array),  intent(  out) :: closeguys(:)
        ! this holds the indices of nearby particles
    type(dist_array),   intent(  out) :: close_dists(:)
        ! this holds the distances to the corresponding nearby particle
    integer                                       :: i, n, nf
        ! loop iterator and number of particles found by search
    type(kdtree2_result), allocatable             :: results(:)
        ! results array from KD tree module

    allocate (results(num_alloc))
    n = end - start + 1

    ! loop over all particles
    do i = 1, n
        ! the type of search used here finds all the particles within
        ! squared distance r2 from the i^th particle in the list
        ! the hard-coded 0 is the 'correlation time' of the search
        call kdtree2_r_nearest_around_point(tree, i + start - 1, 0, r2, nf, num_alloc, results)

        ! allocate these based on how many nearby particles were found
        allocate (closeguys(i)%indices(nf), close_dists(i)%dists(nf))

        closeguys(i)%indices = results(1 : nf)%idx
        close_dists(i)%dists = results(1 : nf)%dis
    end do

    deallocate (results)
end subroutine search

! these next two subroutines use the Box-Muller transform to generate N(0,1)
! random numbers from U(0,1)
! https://goo.gl/DQgmMu
! Note: this polar formulation seems to be consistently ~20% faster than the
! version below that uses trig functions
! reference for polar version (and standard version):
! https://www.taygeta.com/random/gaussian.html
subroutine box_mullerp(n, z)
    integer,          intent(in   ) :: n ! size of random vector to be generated
    double precision, intent(  out) :: z(n)
    integer                         :: j
    double precision                :: w, x1, x2

    call init_random_seed()

    do j = 1, n/2
        w = 1.0d0
        do while (w >= 1.0d0)
            x1 = 2.0d0 * rand() - 1.0d0
            x2 = 2.0d0 * rand() - 1.0d0
            w = x1**2 + x2**2
        enddo
        w = sqrt((-2.0d0 * log(w)) / w)
        z(2 * j - 1 : 2 * j) = (/x1 * w, x2 * w/)
    enddo

    if (mod(n, 2) /= 0) then
        w = 1.0d0
        do while (w >= 1.0d0)
            x1 = 2.0d0 * rand() - 1.0d0
            x2 = 2.0d0 * rand() - 1.0d0
            w = x1**2 + x2**2
        enddo
        w = sqrt((-2.0d0 * log(w)) / w)
        z(n) = x1 * w
    endif
end subroutine box_mullerp

! subroutine box_muller()
!     integer, parameter :: n = 1e8
!     integer            :: i, j
!     double precision   :: x1, x2, y1, y2, z(n)

!     call init_random_seed()
!     i = 1

!     do j = 1, n/2
!         x1 = rand()
!         x2 = rand()
!         y1 = sqrt(-2.0d0 * log(x1)) * cos(2.0d0 * pi * x2)
!         y2 = sqrt(-2.0d0 * log(x1)) * sin(2.0d0 * pi * x2)
!         z(i : i + 1) = (/y1, y2/)
!         i = i + 2
!     enddo
! end subroutine box_muller

end module RPT_mod
