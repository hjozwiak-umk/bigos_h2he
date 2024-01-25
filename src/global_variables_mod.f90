module global_variables_mod
   !! This module defines global variables used throughout the code.
   !---------------------------------------------------------------------------!
   use, intrinsic :: iso_fortran_env, only: int32, sp => real32, dp => real64
   !---------------------------------------------------------------------------!
   implicit none
   !---------------------------------------------------------------------------!
   ! mathematical constants, converters and tolerance control variables
   !---------------------------------------------------------------------------!
   real(dp), parameter :: pi                  = dacos(-1.0_dp)
      !! \\(\Pi\\) value
   real(dp), parameter :: amu_to_au           = 1822.8884862_dp
      !! converter from atomic mass units to atomic units of mass
   real(dp), parameter :: bohr_to_angstrom    = 0.5291772109_dp
      !! converter from atomic units of length to angstrom
   real(dp), parameter :: hartree_to_cm       = 219474.631363_dp
      !! converter from atomic units of energy to cm\\(^{-1}\\)
   real(dp), parameter :: unitary_tolerance   = 1e-4_dp
      !! tolerance on the unitarity check, see Eq. (13) in
      !! "Solution of the coupled equations" section
   !---------------------------------------------------------------------------!
   ! variables read on input (namelist "input")
   !---------------------------------------------------------------------------!
   character(len = 80) :: label = "Test scattering calculations"
      !! user-defined label
   character(len = 80) :: coupling_terms_file_name = "RadialTerms.dat"
      !! name of the file with tabulated radial coupling terms of the PES
   character(len = 8)  :: coupling_terms_r_unit = "--------"
      !! unit used for intermolecular distances in the radial
      !! coupling terms file; only "bohr" and "angstrom" are accepted
   character(len = 80) :: partial_xs_file_name = "partial_xs_file_name.dat"
      !! name of the file holding partial state-to-state cross-sections
   character(len = 80) :: s_matrix_file_name = "s_matrix_file_name.dat"
      !! name of the S-matrix file
   integer(int32) :: relative_energy_flag = 0
      !! if set to 0, "energy" is interpreted as the total energy,
      !! if set to 1, "energy" is interpreted as kinetic energy calculated with
      !! respect to selected "initial" level
   integer(int32) :: jtot_min = 0, jtot_max = -1, jtot_step = 1
      !! range of total angular momenta, for which coupled equations are solved
   integer(int32) :: steps = 10
      !! number of steps per half-de Broglie wavelength of the scattering
      !! system; if provided on input, the number of \\(R\\)-grid points is
      !! determined as
      !! \\(N = \frac{R\_{max}-R\_{min}}{\pi(k\_{max}+k\_{potential\_depth})}\mathrm{steps}\\),
      !! where \\(k\_{max}\\) is the largest wavevector in a given block,
      !! and \\(k\_{potential\_depth}\\) is the correction due to the
      !! depth of the potential, see "potential_depth"
   integer(int32) :: consecutive_blocks_threshold = 1
      !! number of consecutive blocks for which the threshold condition
      !! on cross-sections needs to be fulfilled to terminate calculations
      !! if "jtot_max" = -1
   integer(int32) :: initial_level = -1
      !! if "relative_energy_flag" = 1, "initial" points to a specific level
      !! in the basis indicating the initial molecular level
   integer(int32) :: number_of_basis_levels = -1
      !! number of rovibrational levels in the basis
   integer(int32) :: number_of_r_points = -1
      !! number of grid points for the radial coupling terms in
      !! the potential expansion
   integer(int32) :: number_of_legendre_indices = -1
      !! number of the radial coupling terms in the potential expansion
   integer(int32) :: total_number_of_coupling_terms = 1
      !! total number of coupling terms provided in "coupling_terms_file_name"
   integer(int32) :: n_skip_lines = 0
      !! number of lines that at the beginning of "coupling_terms_file_name"
      !! that are to be skipped
   integer(int32) :: print_level = 2
      !! print level control; should be >= 0
   logical :: print_partial_cross_sections = .false.
      !! if .true. partial cross-sections will be saved
      !! to "partial_xs_file_name"
   real(dp) :: reduced_mass = -1.0_dp
      !! reduced mass of the scattering system
   real(dp) :: energy = -1.0_dp
      !! if "relative_energy_flag" = 0, "energy" is the total energy
      !! if "relative_energy_flag" = 1, "energy" is the kinetic energy
   real(dp) :: r_min = -1.0_dp, r_max = -1.0_dp
      !! range of the propagation
   real(dp) :: r_step = -1.0_dp
      !! if >=0, r_step is used to determine number of steps on the \\(R\\) grid;
      !! otherwise, grid is determined from "steps"
   real(dp) :: potential_depth = 0.0_dp
      !! the absolute value of the depth of the potential, included
      !! in the determination of the step size of the propagator through
      !! \\(k\_{potential\_depth} = \sqrt{2\mu(\mathrm{potential\_depth})}\\)
   real(dp) :: elastic_xs_threshold = 0.1_dp
      !! threshold condition on elastic cross-sections used in "jtot_max"=-1
   real(dp) :: inelastic_xs_threshold = 0.1_dp
      !! threshold condition on elastic cross-sections used in "jtot_max"=-1
   !---------------------------------------------------------------------------!
   ! variables read on input (namelist "basis")
   !---------------------------------------------------------------------------!
   integer(int32), allocatable :: vib_levels(:)
      !! array holding vibrational quantum numbers,
      !! of number_of_basis_levels size
   integer(int32), allocatable :: rot_levels(:)
      !! array holding rotational quantum numbers,
      !! of number_of_basis_levels size,
   real(dp), allocatable :: internal_energies(:)
      !! array holding energies of rovibrational levels of the molecule,
      !! corresponding to vib_levels and rot_levels elements;
      !! of number_of_basis_levels size
   !---------------------------------------------------------------------------!
   ! variables read on input (namelist "potential")
   !---------------------------------------------------------------------------!
   integer(int32), allocatable :: legendre_indices(:)
      !! array holding legendre indices of the PES expansion, \\(\lambda\\),
      !! of number_of_legendre_indices size;
      !! see Eq. 2 in the "Coupling Matrix" section
   integer(int32), allocatable :: vib_couplings(:), rot_couplings(:),          &
      vib_prime_couplings(:), rot_prime_couplings(:)
      !! arrays holding quantum numbers of radial coupling terms:
      !! \\(\eta = v, j\\), and \\(\eta' = v', j'\\), of
      !! "total_number_of_coupling_terms" size;
      !! see Eq. 2 in the "Coupling Matrix" section
   !---------------------------------------------------------------------------!
   ! input/output units
   !---------------------------------------------------------------------------!
   integer(int32), parameter :: input_unit                = 5
      !! unit number for reading the input file
   integer(int32), parameter :: coupling_terms_file_unit  = 8
      !! unit number for reading the coupling terms file
   integer(int32), parameter :: s_matrix_unit             = 11
      !! unit number for writing S-matrix to an external file
   integer(int32), parameter :: partial_file_unit         = 12
      !! unit number for writing partial cross-sections to an external file
   !---------------------------------------------------------------------------!
   ! additional global variables
   !---------------------------------------------------------------------------!
   integer(int32) :: minimal_number_of_coupling_terms
      !! minimal number of coupling terms based on levels provided in the basis
   logical :: units_converted = .false.
      !! if .true. mass and energy units are converted
   real(dp) :: radial_term_distance_converter
      !! converter for the units of length used during coupling terms read
   real(dp) :: radial_term_energy_converter
      !! converter for the units of energy used during coupling terms read
   integer(int32), allocatable :: reduced_vib_couplings(:),                    &
      reduced_rot_couplings(:), reduced_vib_prime_couplings(:),                &
      reduced_rot_prime_couplings(:)
      !! arrays holding quantum numbers of the _necessary_ radial coupling terms
      !! (based on levels provided in the basis):
      !! \\(\eta = v, j\\), and \\(\eta' = v', j'\\), of
      !! "minimal_number_of_coupling_terms" size;
      !! see Eq. 2 in the "Coupling Matrix" section
   real(dp), allocatable :: r_grid(:)
      !! array holding \\(R\\) grid points, filled while reading radial coupling
      !! terms from the "coupling_terms_file_name"; of "number_of_r_points" size
   real(dp), allocatable :: tabulated_coupling_terms(:,:,:)
      !! array holding read radial coupling terms of the PES;
      !! the 3 dimensions correspond to \\(R, \lambda, \eta \eta'\\) grids
      !! and are of (number_of_r_points, number_of_legendre_indices,
      !! total_number_of_coupling_terms) size
   real(dp), allocatable :: coupling_terms(:,:,:)
      !! array holding _necessary_ radial coupling terms of the PES
      !! (based on levels provided in the basis);
      !! the 3 dimensions correspond to \\(R, \lambda, \eta \eta'\\) grids
      !! and are of (number_of_r_points, number_of_legendre_indices,
      !! minimal_number_of_coupling_terms) size
   real(dp), allocatable :: coupling_terms_b_coeffs(:,:,:),                    &
      coupling_terms_c_coeffs(:,:,:), coupling_terms_d_coeffs(:,:,:)
      !! arrays holding spline coefficients interpolating coupling_terms
   !---------------------------------------------------------------------------!
end module global_variables_mod
