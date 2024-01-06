module data_mod
   !! ...
   use, intrinsic :: iso_fortran_env, only: int32, sp => real32, dp => real64
!------------------------------------------------------------------------------!
   implicit none
!------------------------------------------------------------------------------!
   character(len = 80) :: label = "Test scattering calculations"
      !! user-defined label
   character(len = 80) :: potentialfile = "RadialTerms.dat"
      !! name of the file with tabulated radial coupling terms of the PES
   character(len = 80) :: smatrixfile = "SmatrixFile.dat"
      !! name of the S-matrix file
   character(len = 80) :: partialfile = "PartialFile.dat"
      !! name of the file holding partial state-to-state cross-sections
   integer(int32) :: relative_energy_flag = 0
      !! if set to 0, "energy" is interpreted as the total energy,
      !! if set to 1, "energy" is interpreted as kinetic energy calculated with
      !! respect to selected "initial" level
   integer(int32) :: jtotmin = 0, jtotmax = -1, jtotstep = 1
      !! range of total angular momenta, for which coupled equations are solved
   integer(int32) :: steps = 10
      !! number of steps per half-de Broglie wavelength of the scattering system;
      !! if provided on input, the number of \\(R\\)-grid points is determined
      !! as \\(N = \frac{R\_{max}-R\_{min}}{\pi(k\_{max}+k\_{vdepth})}\mathrm{steps}\\),
      !! where \\(k\_{max}\\) is the largest wavevector in a given block,
      !! and \\(k\_{vdepth}\\) is the correction due to the depth of the potential,
      !! see "vdepth"
   integer(int32) :: consecutive_blocks_threshold = 1
      !! number of consecutive blocks for which the threshold condition
      !! on cross-sections needs to be fulfilled to terminate calculations
      !! if "jtotmax" = -1
   integer(int32) :: initial = -1
      !! if "relative_energy_flag" = 1, "initial" points to a specific level
      !! in the basis indicating the initial molecular level
   integer(int32) :: nlevel = -1
      !! number of rovibrational levels in the basis
   integer(int32) :: nr = -1
      !! number of grid points for the radial coupling terms in the potential expansion
   integer(int32) :: nterms = -1
      !! number of the radial coupling terms in the potential expansion
   integer(int32) :: total_number_of_coupling_terms = 1
      !! total number of coupling terms provided in "potentialfile"
   integer(int32) :: n_skip_lines = 0
      !! number of lines that at the beginning of "potentialfile" that are to be skipped
   integer(int32) :: iunits = 0
      !! TO BE CORRECTED
   integer(int32) :: prntlvl = 2
      !! print level control; should be >= 0
   logical :: print_partial_cross_sections = .false.
      !! if .true. partial cross-sections will be saved to "partialfile"
   real(dp) :: reduced_mass = -1.0_dp
      !! reduced mass of the scattering system
   real(dp) :: energy = -1.0_dp
      !! if "relative_energy_flag" = 0, "energy" is the total energy (kinetic+internal)
      !! if "relative_energy_flag" = 1, "energy" is the kinetic energy
   real(dp) :: rmin = -1.0_dp, rmax = -1.0_dp
      !! range of the propagation
   real(dp) :: dr = -1.0_dp
      !! if >=0, dr is used to determine number of steps on the \\(R\\) grid;
      !! otherwise, grid is determined from "steps"
   real(dp) :: vdepth = 0.0_dp
      !! the absolute value of the depth of the potential, included
      !! in the determination of the step size of the propagator through
      !! \\(k\_{vdepth} = \sqrt{2\mu(\mathrm{vdepth})}\\)
   real(dp) :: elastic_xs_threshold = 0.1_dp
      !! threshold condition on elastic cross-sections used in "jtotmax"=-1
   real(dp) :: inelastic_xs_threshold = 0.1_dp
      !! threshold condition on elastic cross-sections used in "jtotmax"=-1
   !---------------------------------------------------------------------------!
   integer(int32) :: minimal_number_of_coupling_terms
      !! minimal number of coupling terms based on levels provided in the basis
   real(dp) :: radial_term_distance_converter, radial_term_energy_converter
   !---------------------------------------------------------------------------!
   integer(int32), parameter :: input_unit        = 5
   integer(int32), parameter :: pes_file_unit     = 8
   integer(int32), parameter :: s_matrix_unit     = 11
   integer(int32), parameter :: partial_file_unit = 12
   !---------------------------------------------------------------------------!
   real(dp), parameter :: amutoau           = 1822.8884862d0
   real(dp), parameter :: bohrtoangstrom    = 0.5291772109d0
   real(dp), parameter :: hartreetocm       = 219474.631363d0
   real(dp), parameter :: pi                = dacos(-1.d0)
   real(dp), parameter :: unitary_tolerance = 1e-6_dp
   !---------------------------------------------------------------------------!
   integer(int32), allocatable :: v1array(:)
      !! array holding vibrational quantum numbers, of nlevel size
   integer(int32), allocatable :: j1array(:)
      !! array holding rotational quantum numbers, of nlevel size,
   integer(int32), allocatable :: l1tab(:)
      !! array holding legendre indices of the PES expansion, \\(\lambda\)),
      !! of nterms size; see Eq. 2 in the "Coupling Matrix" section
   integer(int32), allocatable :: v1pes(:), j1pes(:), v1ppes(:), j1ppes(:)
      !! arrays holding quantum numbers of radial coupling terms:
      !! \\(\eta = v, j\\), and \\(\eta' = v', j'\\), of "total_number_of_coupling_terms"size;
      !! see Eq. 2 in the "Coupling Matrix" section
   integer(int32), allocatable :: reduced_v1pes(:), reduced_j1pes(:),          &
      reduced_v1ppes(:), reduced_j1ppes(:)
      !! arrays holding quantum numbers of the _necessary_ radial coupling terms
      !! (based on levels provided in the basis):
      !! \\(\eta = v, j\\), and \\(\eta' = v', j'\\), of "minimal_number_of_coupling_terms" size;
      !! see Eq. 2 in the "Coupling Matrix" section
   real(dp), allocatable :: elevel(:)
      !! array holding energies of rovibrational levels of the molecule, corresponding to
      !! v1array and j1array elements; of nlevel size
   !---------------------------------------------------------------------------!
   real(dp), allocatable :: rmat(:)
      !! array holding \\(R\\) grid points, filled while reading radial coupling
      !! terms from the "potentialfile"; of "nr" size
   real(dp), allocatable :: read_vmat3D(:,:,:)
      !! array holding read radial coupling terms of the PES;
      !! the 3 dimensions correspond to \\(R, \lambda, \eta \eta'\\) grids
      !! and are of (nr, nterms, total_number_of_coupling_terms) size
   real(dp), allocatable :: vmat3D(:,:,:)
      !! array holding _necessary_ radial coupling terms of the PES
      !! (based on levels provided in the basis);
      !! the 3 dimensions correspond to \\(R, \lambda, \eta \eta'\\) grids
      !! and are of (nr, nterms, minimal_number_of_coupling_terms) size
   real(dp), allocatable :: bmat3D(:,:,:), cmat3D(:,:,:), dmat3D(:,:,:)
      !! arrays holding spline coefficients interpolating vmat3D
   logical :: pes_file_exists = .false.
      !! a flag for checking existence of "potentialfile" -> MOVE
   logical :: units_converted = .false.
      !! if .true. mass and energy units are converted
   !---------------------------------------------------------------------------!
end module data_mod
