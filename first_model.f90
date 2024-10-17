
! This is the main code file for Hanif Kawousi's Msc-project.

!  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!    NOTE: This is the 2Gene Model WITH BACKGROUND MORTALITY AND RISK SIGNALLING
!  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

!------------------------------------
!--------PARAMETER MODULE------------
!------------------------------------
module params
    implicit none

! NON-OBJECT: change these to access or avoid data that will affect compilation.

    logical, public, parameter :: IS_DEBUG_SCREEN = .FALSE.
    logical, public, parameter :: IS_DEBUG_ARRAY_READ = .FALSE.
    logical, public, parameter :: IS_SAVE_POPULATION = .TRUE.
    
    logical, public, parameter :: IS_SAVE_POPULATION_GENES = .FALSE. !remember: IS_ZIP_DATA and IS_DEBUG_DATA must be .TRUE.

    !Both must be the same: .TRUE. or .FALSE., otherwise error-message
    logical, public, parameter :: IS_ZIP_DATA =  .FALSE. ! compress large data before it takes a lot of space
    logical, public, parameter :: IS_DEBUG_DATA= .FALSE. !Run all recording of a value at each timestep: weight, emotion_state, food_availability, food_eaten

    character(*), public, parameter :: ZIP_PROG = 'gzip' !program can be changed in the future.

    ! Numerical code for missing value
    real, parameter, public :: MISSING = -9999.0
    integer, parameter, public :: UNDEFINED = -9999

! ============= BIRD ====================
    integer, parameter, public :: HUNGER_MIN = 0, HUNGER_MAX = 10
    integer, parameter, public :: FEAR_MIN = 0, FEAR_MAX = 10

    real, parameter, public :: BIRD_INITIAL_WEIGHT = 10.0 !20 grams 
    real, parameter, public :: BIRD_MAXIMUM_WEIGHT_MULTIPLICATOR = 4 !The number we multiply initial weight with. This defines birds maximum weight.
    real, parameter, public :: BIRD_MINIMUM_WEIGHT_MULTIPLICATOR = 0.4 !The number we multiply initial weight with. This defines birds minimum weight. 
    real, parameter, public :: WEIGHT_REDUCTION_CONSTANT = 0.0015 !Only for when the bird encounters predator, but escapes. The added stress of escape will manifest in excess weight reduction.
    real, parameter, public :: METABOLIC_COST_CONSTANT = 0.0001 ! The cost of life to be multiplied with birds weight.

    real, parameter, public :: FEAR_DECAY_RATE = 0.02 ! Adjust this value as needed
    real, parameter, public :: FEAR_INCREMENT_AFTER_ATTACK = 0.3 !This is the increased fear of the bird in attack situations.

    real, parameter, public :: HUNGER_PLASTICITY = 1 ! Param for variation in hunger between birds. 
    real, parameter, public :: FEAR_PLASTICITY = 1

    real, parameter, public :: EXPLORATION_PROBABILITY = 0.001 !Should be between 0.0 and 1.0 (most likely below 0.1)
   
   
    !params for fraction in fuction bird_eat_fraction_when_fear
    real, parameter, public :: FRAC_S1_FEAR = 1 !max value 
    real, parameter, public :: FRAC_S2_FEAR = 0.3 !min val
    real, parameter, public :: SIGM_K_FEAR= 1.5

    !params for fraction in fuction bird_eat_fraction_when_hunger
    real, parameter, public :: FRAC_S1_HUNGER = .9 !max value 
    real, parameter, public :: FRAC_S2_HUNGER = 0 !min val
    real, parameter, public :: SIGM_K_HUNGER = 1.25

    !bird limit_food_intake params:
    real, parameter, public :: SIGMOID_STEEPNESS = 1
    real, parameter, public :: SIGMOID_MIDPOINT = 0.5
    real, parameter, public :: VIGILANCE_FACTOR = 0.01

! ========================================

! =========== PREDATOR ====================
   real, parameter, public :: FREQUENCY_OF_PREDATOR_MIN = 0.01
   real, parameter, public :: FREQUENCY_OF_PREDATOR_MAX = 0.15

   ! For signal generations only:
   real, parameter, public :: FEAR_SIGNAL_MULTIPLIER = 1.0 ! if 1 then we are in evolutionary gens. 

!parameter for predation_risk, risk of the bird to meet predator in a particular environment.
! if the bird sees the predator, the predator sees the bird.

  real, parameter, public :: ATTACK_RATE = 0.15!risk of bird being eaten by predator if attacked. 0.1 if Is_Evolutionary_Gens
  
!  ========================================

! =========== ENVIRONMENT ====================
    integer, parameter, public :: ENVIRONMENT_SIZE = 100 !size of envornment, nr of cells = 100

    real, parameter, public :: FOOD_AVAILABILITY_MEAN = 3, FOOD_AVAILABILITY_VARIANCE = 2 !parameter for food, measured in weight grams added to the birds mass
    real, parameter, public :: FOOD_AVAILABILITY_MIN = 1, &
        FOOD_AVAILABILITY_MAX = FOOD_AVAILABILITY_MEAN + FOOD_AVAILABILITY_VARIANCE 
! ===========================================


!! =========== OUTPUT ====================
    ! File name that keeps all output data from all generations of the model
    character(len=*), parameter, public :: MODEL_OUTPUT_FILE = "my_model_output_2G.csv"

    character(len=*), parameter, public :: FOOD_PRED_DISTRIBUTION_FILE = "food_pred_distribution_2G.csv"
! =========================================

!!=== CREATING DIRECTORIES FOR OUTPUT FILES =====

    character(*), parameter, public :: GENERATION_SUBFOLDER = "generations/"

    character(*), parameter, public :: FACTOR_EXPERIMENT_SUBFOLDER = "factor_experiment/"

    character(*), parameter, public :: WEIGHT_HISTORY_SUBFOLDER = "weight_history_2G/"

    character(*), parameter, public :: EMOTION_STATE_SUBFOLDER = "emotion_history_2G/"

    character(*), parameter, public :: EMOTION_GENE_SUBFOLDER = "gene_history_2G/"

!================================================

!======= READ "FACTOR" FILE FOR EXPERIMENT =======
    character(len=*), parameter, public :: FACTOR_FILE = "signal_values.csv"
!=================================================

!=================== GA ===================
    !Population size
    integer, parameter, public :: GENERATIONS = 300
   ! 30 for testing, 300 for sim. 
    integer, parameter, public :: SIGNAL_ONLY_GENERATIONS = 1 !10 for testing fear only, 
                                                               !5 for theoretical experiment with factor from file and no mutation. 

    !Generations for when attack_rate is turned off and only signal remains
    
    integer, parameter, public :: EASY_GENERATIONS = 20 !10 for testing, 20 for sim

    integer, parameter, public :: POP_SIZE = 5000



    real, parameter, public :: GA_REPRODUCE_THRESHOLD = 0.75 !minimum fitness for reproduction
    !Proportion of the best reproducing birds of selection
    real, parameter :: GA_REPRODUCE_PR = 0.7
    !exact number of birds reproducing
    integer, parameter, public :: GA_REPRODUCE_N = int(POP_SIZE * GA_REPRODUCE_PR)

    real, parameter, public :: GA_PROB_REPR_MAX = 0.85
    real, parameter, public :: GA_PROB_REPR_MIN = 0.0
    ! integer, parameter, public :: NUM_EGGS_RANGE_MAX = 9
    ! integer, parameter, public :: NUM_EGGS_RANGE_MIN = 1

    integer, parameter, public :: SELECTION_MULTIPLICATOR = 1

    integer, parameter, public :: TOTAL_TIME_STEP = 365

    real, parameter, public :: MUTATION_PROBABILITY = 0.01 !0.02 for evo, 0.0 when Is_Evolutionary_Generations is .not. and/or Is_Factor_Experiments = .TRUE.
    real, parameter, public :: BACKGROUND_MORTALITY = 0.001 !if Is_Evolutionary_Generations = True
    !Background mortality is a common word birds killed by events as disease, hostile weather-conditions, parasitism, etc...
!==========================================

!==================== Global Variables ===========================
    integer, public :: Current_Generation 
    logical, public :: Is_Evolutionary_Generations 
    logical, public :: Is_Factor_Experiments = .FALSE. ! TRUE IF EXPERIMENT WITH FACTOR FROM EXTERNAL CSV FILE
    
    character(len=:), allocatable, public :: Genome_File_Name

    real, public :: Full_Pop_Factor 
    real, public :: Last_Full_Pop_Factor
    ! the real value one must multiply
    ! with survivors in order to return to POP_SIZE
    ! or size of the population at it's full. 
! ================================================================


    contains

  !-----------------------------------------------------------------------------
  !> Force a value within the range set by the vmin and vmax dummy parameter
  !! values. If the value is within the range, it does not change, if it
  !! falls outside, the output force value is obtained as
  !! min( max( value, FORCE_MIN ), FORCE_MAX )
  !! @param[in] value_in Input value for forcing transformation.
  !! @param[in] vmin minimum value of the force-to range (lower limit), if
  !!            not present, a lower limit of 0.0 is used.
  !! @param[in] vmax maximum value of the force-to range (upper limit)
  !! @returns   an input value forced to the range.
  !! @note      Note that this is the **real** precision version of the
  !!            generic `within` function.
  elemental function within(value_in, vmin, vmax) result (value_out)
    real, intent(in) :: value_in
    real, optional, intent(in) :: vmin
    real, intent(in) :: vmax
    real :: value_out

    ! Local copies of optionals.
    real :: vmin_here

    ! Check optional minimum value, if absent, set to a default value 0.0.
    if (present(vmin)) then
      vmin_here =  vmin
    else
      vmin_here = 0.0
    end if

    value_out = min( max( value_in, vmin_here ), vmax )

  end function within

  !-----------------------------------------------------------------------------
  !> @brief    Rescale a real variable with the range A:B to have the new
  !!           range A1:B1.
  !! @details  Linear transformation of the input value `value_in` such
  !!           `k * value_in + beta`, where the `k` and `beta` coefficients
  !!           are found by solving a simple linear system:
  !!           @f$  \left\{\begin{matrix}
  !!                A_{1}= k \cdot A + \beta;   \\
  !!                B_{1}= k \cdot B + \beta
  !!                \end{matrix}\right. @f$. It has this solution:
  !!           @f$ k=\frac{A_{1}-B_{1}}{A-B},
  !!               \beta=-\frac{A_{1} \cdot B-A \cdot B_{1}}{A-B} @f$
  !! @warning  The function does not check if `value_in` lies within [A:B].
  !! @note     Code for wxMaxima equation solve:
  !! @code
  !!           solve(  [A1=A*k+b, B1=B*k+b] ,[k,b] );
  !! @endcode
  elemental function rescale(value_in, A, B, A1, B1) result(rescaled)
    real, intent(in) :: value_in
    real, intent(in) :: A, B, A1, B1
    real :: rescaled

    ! Local variables
    real :: ck, cb

    !> ### Implementation details ###

    !> First, find the linear coefficients `ck` and `cb`
    !! from the simple linear system.
    ck = (A1-B1) / (A-B)
    cb = -1.0 * ((A1*B - A*B1) / (A-B))

    !> Second, do the actual linear rescale of the input value.
    rescaled = value_in*ck + cb

  end function rescale


  !-----------------------------------------------------------------------------
  !> Calculate an average value of a real array, excluding MISSING values.
  !! @param vector_in The input data vector
  !! @param missing_code Optional parameter setting the missing data code,
  !!        to be excluded from the calculation of the mean.
  !! @param undef_ret_null Optional parameter, if TRUE, the function returns
  !!        zero rather than undefined if the sample size is zero.
  !! @returns The mean value of the vector.
  !! @note This is a real array version.
  pure function average (array_in, missing_code, undef_ret_null)            &
                                                              result (mean_val)

    ! @param vector_in The input data vector
    real, dimension(:), intent(in) :: array_in

    ! @param missing_code Optional parameter setting the missing data code,
    !        to be excluded from the calculation of the mean.
    real, optional, intent(in) :: missing_code

    ! @param undef_ret_null Optional parameter, if TRUE, the function returns
    ! zero rather than undefined if the sample size is zero.
    logical, optional, intent(in) :: undef_ret_null

    ! @returns The mean value of the vector.
    real :: mean_val

    ! Local missing code.
    real :: missing_code_here

    ! Local sample size, N of valid values.
    integer :: count_valid

    ! Define high precision kind for very big value
    integer, parameter :: HRP = selected_real_kind(33, 4931)

    ! Big arrays can result in huge sum values, to accommodate them,
    ! use commondata::hrp and commondata::long types
    real(HRP) :: bigsum, bigmean

    !> ### Implementation details ###

    !> Check if missing data code is provided from dummy input.
    !! If not, use global parameter.
    if (present(missing_code)) then
      missing_code_here = missing_code
    else
      missing_code_here = MISSING
    end if

    !> Fist, count how many valid values are there in the array.
    count_valid = count(array_in /= missing_code_here)

    !> If there are no valid values in the array, mean is undefined.
    if (count_valid==0) then
      if (present(undef_ret_null)) then
        if (undef_ret_null) then
          mean_val = 0.0    !> still return zero if undef_ret_null is TRUE.
        else
          mean_val = MISSING
        end if
      else
        mean_val = MISSING
      end if
      return
    end if

    bigsum = sum( real(array_in, HRP), array_in /= missing_code_here )
    bigmean = bigsum / count_valid

    mean_val =  real(bigmean)

  end function average


!   pure function update_pop_regain_size_factor(pop_survived, pop_full) result(Full_Pop_Factor)
!   integer, dimension(:), intent(in) :: pop_survived
!   real, dimension(:), intent(out) :: pop_full ! Assuming pop_full is meant to be returned as well
!   real :: Full_Pop_Factor


  !> Calculates the factor to regain the full population size from the survived population.
  !!
  !! This function takes the number of individuals that survived and the full population size,
  !! and calculates the factor that can be used to scale up the survived population to the
  !! full population size.
  !!
  !! @param pop_survived The number of individuals that survived.
  !! @param pop_full The full population size.
  !! @return The factor to scale up the survived population to the full population size.
  pure function update_pop_regain_size_factor(pop_survived, pop_full) result(factor)
  integer, intent(in) :: pop_survived
  integer, intent(in) :: pop_full ! Assuming pop_full is meant to be returned as well
  real :: factor

    factor = real(pop_full) / real(pop_survived)


  end function update_pop_regain_size_factor


  !For finding probability of reproduction
  !> Calculates the probability of reproduction for an individual based on its current 
  !! mass (m) and the mass thresholds for reproduction (m_0 and m_max).
  !!
  !! This function uses a linear interpolation between the minimum and maximum 
  !! probabilities of reproduction (GA_PROB_REPR_MIN and GA_PROB_REPR_MAX) to determine the probability of reproduction for the given mass.
  !!
  !! @param m The current mass of the individual.
  !! @param m_0 The minimum mass threshold for reproduction.
  !! @param m_max The maximum mass threshold for reproduction.
  !! @return The probability of reproduction for the individual.

  function prob_repr(m, m_0, m_max) result(prob)

  real, intent(in) :: m, m_0, m_max
  real  :: prob
  real ::  k, b
  real, parameter :: P  = GA_PROB_REPR_MAX
  real, parameter :: P_0 = GA_PROB_REPR_MIN

  call solve_linear( k, b, m_0, P_0, m_max, P )

  prob = k * m + b

  prob = within(prob, P_0, P)

  end function prob_repr

  
  ! General solver for linear equation
  !
  ! Given the x_min, x_max and y_min y_max,
  !   determine the coefficients for the linear
  !   equation k and b
  ! y_max +       *
  !       |     * .    y = k x + b
  !       |   *   .        ?     ?
  !       | *     .
  ! y_min +-------+
  !     x_min   x_max
  
  !> Solves for the coefficients k and b of a linear equation y = kx + b, 
  !! given the minimum and maximum values of x and y.
  !!
  !! This subroutine takes the minimum and maximum values of x and y, 
  !! and calculates the coefficients k and b for the linear equation y = kx + b 
  !! that passes through those points.
  !!
  !! @param k The slope of the linear equation.
  !! @param b The y-intercept of the linear equation.
  !! @param x_min The minimum value of x.
  !! @param y_min The minimum value of y.
  !! @param x_max The maximum value of x.
  !! @param y_max The maximum value of y.
  pure subroutine solve_linear(k, b, x_min, y_min, x_max, y_max)
    real, intent(out) :: k,b
    real, intent(in)  :: x_min, y_min, x_max, y_max

    k = (y_min - y_max) / (x_min-x_max)
    b = -1 * (x_max*y_min - x_min*y_max) / ( x_min - x_max )

  end subroutine solve_linear


  !function and subsequent subroutine for calculating median value. 
  ! https://rosettacode.org/wiki/Averages/Median#Fortran
  function median (array_in) result(median_value)
     implicit none
     real, dimension(:), intent(in) :: array_in
  
     real, dimension(:), allocatable :: sorted_array
     real :: median_value
     integer :: i, n, mid 

     integer :: n_valid, index_valid


     ! make a smaller array from the input array that excludes missing values: 
     ! Make size of sorted array to 
 
     !Copy array to avoid modifying original data
     n_valid = 0

     do i=1, size(array_in)
       if (array_in(i) /= MISSING) n_valid = n_valid + 1
     end do

     allocate(sorted_array(n_valid))
     

     
     index_valid = 0 
     do i=1, size(array_in)
        if (array_in(i) /= MISSING) then 
          index_valid = index_valid + 1 
          sorted_array(index_valid) = array_in(i)
        end if
     end do
     
     n = size(sorted_array)

     call sort(sorted_array)

     !Calculate the median 
     
     if (mod(n, 2) == 0) then 
       mid = n / 2 
       median_value = (sorted_array(mid) + sorted_array(mid + 1)) / 2.0
     else
       mid = (n + 1) / 2
       median_value = sorted_array(mid)
     end if 

  end function median
  
 
  subroutine sort(array)
    implicit none
    real, intent(inout) :: array(:)
    integer :: n
    integer :: i, j
    real :: temp

    n = size(array)

  ! Making a bubble sort algorithm: 
  ! Bubble sort algorithm, also known as sinking sort,
  !  is the simplest sorting algorithm that runs through the list repeatedly,
  !   compares adjacent elements, and swaps them if they are out of order.
  ! Code inspired by https://rosettacode.org/wiki/Sorting_algorithms/Bubble_sort#Fortran 
  !
   do i = 1, n - 1
     do j = i + 1, n
       if (array(i) > array(j)) then
         temp = array(i)
         array(i) = array(j)
         array(j) = temp
       end if
     end do
   end do

  end subroutine sort


  
  !-----------------------------------------------------------------------------
  !> Calculate standard deviation using trivial formula:
  !! @f[ \sigma=\sqrt{\frac{\sum (x-\overline{x})^{2}}{N-1}} . @f]
  !! @note This is a real array version.
  function std_dev(array_in, missing_code, undef_ret_null) result (stddev)
    !> @param vector_in The input data vector
    real, dimension(:), intent(in) :: array_in
    !> @param missing_code Optional parameter setting the missing data code,
    !!        to be excluded from the calculation of the mean.
    real, optional, intent(in) :: missing_code
    !> @param undef_ret_null Optional parameter, if TRUE, the function returns
    !! zero rather than undefined if the sample size is zero.
    logical, optional, intent(in) :: undef_ret_null
    !> @returns The standard deviation of the data vector.
    real :: stddev

    ! Local missing code.
    real :: missing_code_here

    ! Local sample size, N of valid values.
    integer :: count_valid

    ! Minimum sample size resulting in real calculation, everythin less
    ! than this returns invalid.
    integer, parameter :: MIN_N = 4

    ! An array of squared deviations
    real, dimension(size(array_in)) :: dev2

    ! Mean value
    real :: sample_mean

    !> Check if missing data code is provided from dummy input.
    !! If not, use global parameter.
    if (present(missing_code)) then
      missing_code_here = missing_code
    else
      missing_code_here = MISSING
    end if

    count_valid = count(array_in /= missing_code_here)

    !> If there are no valid values in the array, std. dev. is undefined.
    if (count_valid <  MIN_N) then
      if (present(undef_ret_null)) then
        if (undef_ret_null) then
          stddev = 0.0    !> still return zero if undef_ret_null is TRUE.
        else
          stddev = MISSING
        end if
      else
        stddev = MISSING
      end if
      return
    end if

    sample_mean = average ( array_in, missing_code_here )

    where ( array_in /= missing_code_here )
      dev2 = ( array_in - sample_mean ) ** 2
    elsewhere
      dev2 = missing_code_here
    end where

    stddev = sqrt( sum( dev2, dev2 /= missing_code_here ) / (count_valid - 1) )

  end function std_dev

  !-----------------------------------------------------------------------------
  !> Converts logical to standard (kind SRP) real, .FALSE. => 0, .TRUE. => 1
  !! @note Note that this function is required to place logical data
  !!       like the survival status (alive) into the reshape array that is
  !!       saved as the CSV data file.
  !!       **Example: **
  !!       @code
  !!          call CSV_MATRIX_WRITE ( reshape(                              &
  !!                     [ habitat_safe%food%food%x,                        &
  !!                       habitat_safe%food%food%y,                        &
  !!                       habitat_safe%food%food%depth,                    &
  !!                       conv_l2r(habitat_safe%food%food%eaten),          &
  !!                       habitat_safe%food%food%size],                    &
  !!                     [habitat_safe%food%number_food_items, 5]),         &
  !!                     "zzz_food_s" // MODEL_NAME // "_" // MMDD //       &
  !!                       "_gen_" // TOSTR(generat, GENERATIONS) // csv,   &
  !!                     ["X   ","Y   ", "D   ", "EATN", "SIZE"]            &
  !!                     ) @endcode
  elemental function l2r(flag, code_false, code_true) result (num_out)
    logical, intent(in) :: flag
    real, optional, intent(in) :: code_false
    real, optional, intent(in) :: code_true
    real :: num_out
    ! Local copies of optionals.
    real :: code_false_loc, code_true_loc
    ! Local real parameters for the default FALSE and TRUE.
    real, parameter :: FALSE_DEF=0.0, TRUE_DEF=1.0

    !> First, check optional parameters.
    if (present(code_false)) then
      code_false_loc = code_false
    else
      code_false_loc = FALSE_DEF
    end if

    if (present(code_true)) then
      code_true_loc = code_true
    else
      code_true_loc = TRUE_DEF
    end if

    ! Second, do the actual conversion.
    if (flag .eqv. .TRUE.) then
      num_out = code_true_loc
    else
      num_out = code_false_loc
    end if
  end function l2r



! Function for calculating linear increase in predation to ensure soft pressure from the environment in the initial 20 generations.

  function approx_easy(generation, min_y, max_y) result (y_approx)
    real :: y_approx
    integer, intent(in) :: generation
    real, intent(in) :: min_y, max_y
    !when M is at it's highest value, the attack_rate is also at it's highest,
    !when M is lower than it's highest value, the attack rate is somewhere along
    !the linear slope that goes from minimum attack rate (min_y) to maximum attack rate
    !(max_y).

    real :: k, b ! parameters of the linear equation

    integer, parameter :: M = EASY_GENERATIONS !Parameter for the generations where the environment is less harsh and allows for an individual's mistakes.
    

     k = -1 * (min_y - max_y) / (M - 1)
     b = (M*min_y - max_y) / (M - 1)
  
     if (generation < M) then !if the generation nr is less than M(=20)
       y_approx = k * real(generation) + b !y_approx, the environments predation is k * (real)generation.
     else                                  !"Real" due to temporal continuity within and between generation (slope of the linear eq)
       y_approx = max_y !if generation is 20 or above, continue with universal parameters.
     end if


  end function approx_easy



  subroutine get_global_runtime_params()
! Purpose
!   The subroutine get_global_runtime_params is designed to retrieve and 
! process command-line arguments passed to the program. It sets global 
! parameters based on the presence and values of these arguments.
! Description:
!   This subroutine fetches the number of command-line arguments 
! and stores them in an allocatable array. It then processes the 
! first argument, if available, to set the values of two global 
! parameters: Genome_File_Name and Is_Evolutionary_Generations.
! Usage:
!   This subroutine does not take any input arguments nor return any output. 
! It relies on command-line arguments passed when the Fortran program is executed.
! Example, in terminal, if we run gfortran -o example example.f90:
    !./example
    ! expected output should be:      Is_Evolutionary_Generations:  T
    !if we run ./example my_genome_file.txt
    ! expected output:                Is_Evolutionary_Generations:  F

  integer :: i, n_args
  character(100), dimension(:), allocatable :: cmd_args

  ! Get the number of command-line arguments
  n_args = command_argument_count()
  
  ! Allocate array to hold command-line arguments
  allocate(cmd_args(n_args))
  
  ! Retrieve each command-line argument and store it in cmd_args
  do i = 1, n_args
    call get_command_argument(number = i, value = cmd_args(i))
  end do

  ! Process the command-line arguments to set global parameters
  if (n_args > 0) then
    Genome_File_Name = trim(cmd_args(1))
    Is_Evolutionary_Generations = .FALSE.
  else
    Genome_File_Name = " "
    Is_Evolutionary_Generations = .TRUE.
  end if

end subroutine get_global_runtime_params

  

end module params

!------------------------------------
!------ENVIRONMENT MODULE------------
!------------------------------------

module environment

use params
implicit none

! Spatial location at one cell, in the initial model it is
! just a single cell, integer number
! Note: defining an elementary spatial unit makes it easy to
!       extend the model to 2d or 3d
type, public :: location
    integer :: x
    contains
    procedure, public :: put_random => location_place_random
    procedure, public :: place => location_place_object_location !for bird and predator
    procedure, public :: walk => location_random_walk !only for bird


end type location

! Each spatial cell has additional data characteristics:
! - food availability
! - predation risk
type, extends(location) :: env_cell
  real :: food_availability
  real :: predator_frequency
  contains
  procedure, public :: init => cell_init
end type env_cell

! The whole environment is a collection of cells
! Note: note that the size of the environment is an allocatable (dynamic)
!       array.
type whole_environ
  type(env_cell), allocatable, dimension(:) :: point
  contains
  procedure, public :: init => env_init
  procedure, public :: save => environemnt_save_csv
  procedure, public :: load_env_csv => load_environment_from_csv

end type whole_environ

contains

! This subroutine places the basic spatial object (location class) to a
! random place.
subroutine location_place_random(this, min_pos, max_pos)
    use BASE_RANDOM
    class(location), intent(inout) :: this
    ! Optional parameters defining a range of positions to place the spatial
    ! object, but by default the position is from 1 to ENVIRONMENT_SIZE
    integer, optional, intent(in) :: min_pos, max_pos

    integer :: min_loc, max_loc

    if (present(min_pos)) then
        min_loc = min_pos
    else
        min_loc = 1
    end if

    if (present(max_pos)) then
        max_loc = max_pos
    else
        max_loc = ENVIRONMENT_SIZE
    end if

    this%x = RAND(min_loc, max_loc)

end subroutine location_place_random

! Place a spatial object to a specific location (cell) within the environment
subroutine location_place_object_location(this, where)
    class(location), intent(inout) :: this
    integer, intent(in) :: where

    this%x = where

end subroutine location_place_object_location



subroutine location_random_walk(this, min_pos, max_pos)
  use BASE_RANDOM
  class(location), intent(inout) :: this

  integer, optional, intent(in) :: min_pos, max_pos

  !deside if left or right-movement
  logical :: is_going_left

  integer :: min_loc, max_loc

  ! Process optional parameters using the min_pos and max_pos from
  ! the subroutine call or default values
  if (present(min_pos)) then
        min_loc = min_pos
  else
        min_loc = 1
  end if

  if (present(max_pos)) then
        max_loc = max_pos
  else
        max_loc = ENVIRONMENT_SIZE
  end if

  if ( RAND() > 0.5 ) then
    is_going_left = .TRUE.
  else
    is_going_left = .FALSE.
  end if

  if ( is_going_left .and. this%x == min_loc ) then
    is_going_left = .FALSE.
  elseif ( .NOT. is_going_left .and. this%x == max_loc) then
    is_going_left = .TRUE.
  end if

  if (is_going_left) then
    this%x = this%x - 1
  else
    this%x = this%x + 1
  end if

end subroutine location_random_walk

! Initialize a single cell within the environment
subroutine cell_init(this, x, env)
  use BASE_RANDOM
  class(env_cell), intent(out) :: this
  integer, intent(in) :: x
  type(whole_environ), optional, intent(inout) :: env

  this%x = x ! set location
  this%food_availability = within(RNORM(FOOD_AVAILABILITY_MEAN, FOOD_AVAILABILITY_VARIANCE), &
                                  FOOD_AVAILABILITY_MIN, FOOD_AVAILABILITY_MAX)
  if (Is_Evolutionary_Generations) then
      this%predator_frequency = RAND(FREQUENCY_OF_PREDATOR_MIN, FREQUENCY_OF_PREDATOR_MAX)
  else
      ! if (present(env)) then
      !     call env%load_env_csv()
      !     this%predator_frequency = env%point(x)%predator_frequency
      !  else
    this%predator_frequency = RAND(FREQUENCY_OF_PREDATOR_MIN * FEAR_SIGNAL_MULTIPLIER,     &
        FREQUENCY_OF_PREDATOR_MAX * FEAR_SIGNAL_MULTIPLIER)
      ! end if
      ! this%predator_frequency = this%predator_frequency * FEAR_SIGNAL_MULTIPLIER
      !  end if
  end if
end subroutine cell_init




! Initialize the whole environment
subroutine env_init(this, max_size) !max size of the environ, but optional. Made for testing with shorter arrays.
  class(whole_environ), intent(inout) :: this
  integer, intent(in), optional :: max_size

  integer :: i, max_size_loc

  ! process optional parameter defining the maximum size of the environment
  if (present(max_size)) then
    max_size_loc = max_size
  else
    max_size_loc = ENVIRONMENT_SIZE
  end if

  ! First, allocate the dynamic environment array with its size
  ! Note that we test if the array has already been allocated to guard against
  !      possible error, e.g. if init is called more than once
  if (.not. allocated(this%point) ) then
    allocate (this%point(max_size_loc))
  else
    ! Repeated initialization of the environment
    ! would be strange so we report this
    write(*,*) "WARNING: repeated initialization of the environment detected"
    deallocate (this%point)
    allocate (this%point(max_size_loc))
  end if

  do i=1, size(this%point)
    call this%point(i)%init(i, this) ! we initialize the food and predation data
  end do

end subroutine env_init



!Subroutine for saving info on food-distribution across the one-dimensional 
!matrix that is our birds array of habitats.
!Subroutine logic:
!Using only CSV_ARRAY_WRITE from CSV_IO, an array is created. This array 
!makes an output array(output_array) which is basically a one dimensional 
!array of cells that correspond the habitats our bird walks in. This array
!consist of however many cells we want ("this"). When we loop over it from 
!i = 1 until total array size, for each loop the food that is within each 
!cell is recorded. This only needs to be done once, since the food in cells 
!are a non-changing value. This is then saved to output_array, that gets sent 
!to FOOD_PRED_DISTRIBUTION_FILE
subroutine environemnt_save_csv(this)
  use CSV_IO, only : CSV_MATRIX_WRITE
  class(whole_environ), intent(in) :: this
  integer  :: i
  real, allocatable, dimension(:,:) :: output_array


  allocate(output_array(size(this%point), 2))

  do i=1, size(this%point)
   output_array(i,1) = this%point(i)%food_availability
   output_array(i,2) = this%point(i)%predator_frequency
  end do

  call CSV_MATRIX_WRITE(output_array, FOOD_PRED_DISTRIBUTION_FILE, colnames = ["FOOD_AVAILABILITY", "RISK_PREDATION   "])

end subroutine environemnt_save_csv



!> Loads the environment data from a CSV file and initializes the environment cells.
!>
!> If the CSV file is successfully opened, the function reads the food availability 
!  and predator frequency for each cell in the environment and stores them in the 
!  corresponding `this%point` elements. If the file does not contain data for all cells, 
!  the remaining cells are initialized using the `cell_init` subroutine.
!>
!> If the CSV file cannot be opened, the function generates a new environment 
!  by initializing all cells using the `cell_init` subroutine.
!>
!> @param this The `whole_environ` object representing the environment.
subroutine load_environment_from_csv(this)
  class(whole_environ), intent(inout) :: this
  integer :: i, io_status
  character(len=100) :: line
  
  open(unit=10, file=FOOD_PRED_DISTRIBUTION_FILE, status='old', action='read', iostat=io_status)
  
  if (io_status == 0) then
      read(10, *, iostat=io_status) ! Skip header line
      
      do i = 1, size(this%point)
          read(10, *, iostat=io_status) this%point(i)%food_availability !, this%point(i)%predator_frequency
          if (io_status /= 0) exit
          this%point(i)%x = i
      end do
      
      close(10)
      
      if (i <= size(this%point)) then
          print*, "Warning: Not all environment data loaded. Generating remaining cells."
          do i = i, size(this%point)
              call cell_init(this%point(i), i)
          end do
      end if
  else
      print*, "Error opening CSV file. Generating new environment."
      do i = 1, size(this%point)
          call cell_init(this%point(i), i)
      end do
  end if
end subroutine load_environment_from_csv

end module environment



!--------------------------------------
!-------GENOME MODULE------------------
!--------------------------------------
module organism

use params
use environment
use BASE_RANDOM
implicit none

!describe how we define genes
!make an integer range that describes fear vs hunger

!--------------------------------
!-----------LOCATION/GENE--------
!--------------------------------
type, extends(location), public :: GENOME
    integer :: gene_fear !the higher the gene_fear the lower the gene, opposite proportional
    integer :: gene_hunger !DOCUMENT ME! 
    contains
    procedure, public :: init_genome => genome_init_all
    procedure, public :: mutate => gene_mutate


end type GENOME


!--------------------------------
!-------------BIRD/GENOME----------
!--------------------------------
type, extends(GENOME), public :: BIRD
    real :: weight  
    !In this model set by the params BIRD_INITAL_WEIGHT and BIRD_MAXIMUM_WEIGHT_MULTIPLICATOR/BIRD_MINIMUM_WEIGHT_MULTIPLICATOR
    real, dimension(TOTAL_TIME_STEP) :: weight_history 
    !recording weight in matrix by timestep and generation
    real, dimension(TOTAL_TIME_STEP) :: fear_history 
    ! recording emotion in matrix by timestep and generation - previously emotion_history
    real, dimension(TOTAL_TIME_STEP) :: hunger_history
    ! recording genes in matrix by timestep and generation 
    real, dimension(TOTAL_TIME_STEP) :: hunger_gene_history
    real, dimension(TOTAL_TIME_STEP) :: fear_gene_history
    real :: state_fear !values between 0 and 10
    real :: state_hunger !values between 0 and 10
    logical :: is_alive 
    logical :: is_killed_by_predator
    integer :: bird_meets_predator_counter
    


    !integer :: death_count !this relates to subroutine count_the_dead_birds.
    contains
    procedure, public :: init => bird_init  !updated
    procedure, public :: killed_from_starvation => bird_dead_from_starvation !updated
    procedure, public :: is_starved => bird_is_starved  !updated
    procedure, public :: fly => bird_do_fly ! updated
    procedure, public :: do_feed => bird_feeds !updated
    procedure, public :: do_explore => bird_do_explore

    procedure, public :: is_killed_by_background => bird_dead_from_background_mortality !NEW! SUBROUTINE
    procedure, public :: environment_kills => background_mortality_probability !NEW! FUNCTION
    procedure, public :: decay_fear => bird_decay_fear


    procedure, public :: genetic_lower_hunger => bird_hunger_min_genetic_limit !NEW!
    procedure, public :: genetic_upper_hunger => bird_hunger_max_genetic_limit !NEW!

    procedure, public :: genetic_lower_fear => bird_fear_min_genetic_limit   !NEW!
    procedure, public :: genetic_upper_fear => bird_fear_max_genetic_limit   !NEW!

    procedure, public :: hunger_genetic_limit => bird_force_hunger_within_genetic_limits !updated!
    procedure, public :: fear_genetic_limit => bird_force_fear_within_genetic_limits   !NEW!

    procedure, public :: add_to_history => bird_weight_add_to_history !UPDATED!
    procedure, public :: emo_state_to_history => bird_emotion_state_add_to_history !UPDATED!
    procedure, public :: pay_for_life => bird_subtract_metabolic_cost  !Keep as is!

    procedure, public :: is_within => bird_is_within_environment_check 
    !checking if bird is within the environment with correct number of cells 

end type BIRD

!--------------------------------
!-----------PREDATOR-------------
!--------------------------------
type, public :: PREDATOR
    real :: risk !If PREDATOR sees the bird according to the risk_of_the, the risk of attack is this


  contains

    procedure, public :: init => predator_init_new !keep as is!

    procedure, public :: predator_attack_bird

end type PREDATOR


!--------------------------------
!-----------POPULATION-----------
!--------------------------------

type, public :: POPULATION
!population size:
  type(BIRD), dimension(:), allocatable :: birds
  integer :: num_dead !global counter of dead birds
  logical :: is_alive

!any other properties of the population?
!-----
  contains
  !subroutines and functions that apply for the population
  procedure, public :: init_pop => population_init_new
  procedure, public :: dead_count => population_get_dead_count
  procedure, public :: alive_count => population_get_alive_count
  procedure, public :: time_steps => population_time_steps
  procedure, public :: sort => sort_by_fitness
  procedure, public :: save => population_save_csv
  procedure, public :: load => population_load_genome_csv
  procedure, public :: smallest_weight => population_minimum_weight_alive
  procedure, public :: biggest_weight => population_maximum_weight_alive





end type POPULATION



!--------------------------------
!---------SUBROUTINES-----------
!--------------------------------
contains
!subroutines for POPULATION


!----------------------------------------

! The "population_init_new" subroutine initializes a population of birds within a 
! POPULATION class object. It allows for the optional specification of an 
! alternative population size (alt_size) and an alternative gene initialization 
! (is_gene_init) for testing and debugging purposes.

! "this": A POPULATION class object passed by reference (intent(out)). This parameter 
! is the target for the population initialization. "alt_size": An optional integer 
! parameter (intent(in)) that specifies an alternative population size. If not provided, 
!the default population size (POP_SIZE) defined in the parameters module is used.

! "is_gene_init": An optional logical parameter (intent(in)) that indicates whether 
! to use an alternative gene initialization. If not provided, the default value is .FALSE..!


subroutine population_init_new(this, alt_size, is_gene_init)
  class(POPULATION), intent(out) :: this
  integer, optional, intent(in) :: alt_size! if we want other value than POP_SIZE
  logical, optional, intent(in) :: is_gene_init !or testing and debugging

  integer :: size
  integer :: i

  logical :: is_gene_init_here

  if (present(is_gene_init)) then
    is_gene_init_here = is_gene_init
  else
    is_gene_init_here = .FALSE.
  end if

  !process for optional parameter where we can define alternative pop_size:
  if (present(alt_size)) then
    size = alt_size
  else
    size = POP_SIZE !POP_SIZE defined in parameters mod.
  end if

  allocate (this%birds(size)) !allocate the array to a specific size

  do i = 1, size
    call this%birds(i)%init(is_gene_init_here) !initialize array of birds

  end do

end subroutine population_init_new

! Creating a csv-file with data under the column names seen below. 
subroutine population_save_csv(this, file_name) 
  use CSV_IO
  use params
  
  class(POPULATION), intent(in) :: this
  character(len=*), intent(in) :: file_name
  
  integer, parameter :: GENPOP_COL_NUMBER = 9
  character(len=*), dimension(GENPOP_COL_NUMBER), parameter :: GENPOP_COLUMNS = &
    ["BIRD                 ",                      &  ! 1
     "GENE_FEAR            ",                      &  ! 2
     "GENE_HUNGER          ",                      &  ! 3
     "IS_ALIVE             ",                      &  ! 4 !0 = FALSE, 1 = TRUE
     "WEIGHT               ",                      &  ! 5
     "BIRD_MEETS_PRED_COUNT",                      &  ! 6
     "STATE_FEAR           ",                      &  ! 7
     "STATE_HUNGER         ",                      &  ! 8
     "REPRODUCTION_PROB.   "]                        ! 9
    


  real, dimension(POP_SIZE, GENPOP_COL_NUMBER) :: out_data

  integer :: i

   do i = 1 ,POP_SIZE
     out_data(i,1) = i
     out_data(i,2) = real(this%birds(i)%gene_fear)
     out_data(i,3) = real(this%birds(i)%gene_hunger)
     out_data(i,4) = l2r(this%birds(i)%is_alive) !0 = FALSE, 1 = TRUE
     out_data(i,5) = this%birds(i)%weight
     out_data(i,6) = real(this%birds(i)%bird_meets_predator_counter)
     out_data(i,7) = real(this%birds(i)%state_fear)
     out_data(i,8) = real(this%birds(i)%state_hunger)
     out_data(i,9) = (prob_repr(this%birds(i)%weight, this%smallest_weight(), this%biggest_weight()))
    ! out_data(i,10)= prob_repr_to_egg(out_data(i,9))
   end do

  call CSV_MATRIX_WRITE(out_data, file_name, colnames=GENPOP_COLUMNS)

end subroutine population_save_csv




!This subroutine loads a population of birds from a CSV file. 
!It reads the gene values and alive/dead status for each bird, 
!and initializes the population accordingly. 
!If the 'IS_KILL_BIRDS_DEAD_IN_FILE' parameter is set to .TRUE., 
!any birds that were dead in the original simulation will be 
!marked as dead when loaded from the file.subroutine population_load_genome_csv(this, file_name)
subroutine population_load_genome_csv(this, file_name)    
use CSV_IO
    class(POPULATION), intent(inout) :: this
    character(len=*), intent(in) :: file_name

    logical, parameter :: IS_KILL_BIRDS_DEAD_IN_FILE = .FALSE.
    !Start with false, TRUE if time for further experimentation
    ! if true, all the birds from the data are alive.

    integer, parameter :: GENPOP_COL_NUMBER = 9
    integer, parameter :: FEAR_GENE_COL = 2, HUNGER_GENE_COL = 3, ALIVE_COL = 4
    integer, parameter :: WEIGHT_COL = 5

    logical :: is_success
    integer :: inputdata_rows, inputdata_cols, i
    real, allocatable, dimension(:,:) :: in_data

    in_data = CSV_MATRIX_READ(file_name, is_success, .FALSE.)
    !see documentation at AHAmodel https://ahamodel.uib.no/
    !errorflag is by default .TRUE. => everything works as it should.

    inputdata_rows = size(in_data, 1)
    inputdata_cols = size(in_data, 2)

    if (.NOT. is_success) then
        print*, "ERROR: file error for ", file_name
        stop
    end if

    if(size(in_data,1) /= POP_SIZE) then
        print*, "WARNING: input data wrong N of rows: n", size(in_data,1), &
                ";Population size = ", POP_SIZE
    end if

    if (inputdata_cols /= GENPOP_COL_NUMBER) then
        print*, "WARNING: input data wrong N of columns: n", inputdata_cols, &
                ", must be ", GENPOP_COL_NUMBER
        stop
    end if

    this%birds(1:inputdata_rows)%gene_fear = int(in_data(:,FEAR_GENE_COL))
    this%birds(1:inputdata_rows)%gene_hunger = int(in_data(:,HUNGER_GENE_COL))
    this%birds(1:inputdata_rows)%weight = in_data(:,WEIGHT_COL)

    if (IS_KILL_BIRDS_DEAD_IN_FILE) then
        do i=1, inputdata_rows
            if (in_data(i, ALIVE_COL) == 0.0) this%birds(i)%is_alive = .FALSE.
        end do
    end if

    this%birds(inputdata_rows+1 : POP_SIZE)%is_alive = .FALSE.
    !make code to exclude dead birds
    !make warning if incorrect data is running

end subroutine population_load_genome_csv



! Functions for counting dead birds, updating it. Iterating over the array.
! 
function population_get_dead_count(this) result(get_dead_count)
!this = population, num_dead is the number of dead birds.
!this subroutine initialize a population of birds within a POPULATION class object. It allows
!for the optional specification of an alternative population size (alt_size), and an alternative
!gene initialization (is_gene_init) for testing and debugging.

  class(POPULATION), intent(in) :: this
  integer :: get_dead_count

  get_dead_count = count(.not.this%birds%is_alive) !birdS = population of birds.

end function population_get_dead_count


function population_get_alive_count(this) result(get_alive_count)! this = population, 
!num_dead is the number of dead birds.
!this subroutine initialize a population of birds within a POPULATION class object. It allows
!for the optional specification of an alternative population size (alt_size), and an alternative
!gene initialization (is_gene_init) for testing and debugging.

  class(POPULATION), intent(in) :: this
  integer :: get_alive_count

  get_alive_count = count(this%birds%is_alive) 

end function population_get_alive_count


function population_minimum_weight_alive(this) result(min_size)
  class(POPULATION), intent(in) :: this
  real :: min_size
  integer :: i 

  i = this%alive_count()
  min_size = this%birds( i )%weight
  
end function population_minimum_weight_alive

function population_maximum_weight_alive(this) result(max_size)
  class(POPULATION), intent(in) :: this
  real :: max_size

  max_size = this%birds( 1 )%weight
  
end function population_maximum_weight_alive


! Subroutines for PREDATOR
subroutine predator_init_new(this, generation)

! This subroutine sets a condition for the attack rate: If the time past is less than 20
! generations, the risk is lower - somewhere on a linear slope - towards it's max which is in
! turn defined by ATTACK_RATE parameter. The conditions given are defined by the function aprox_y.

  class(PREDATOR), intent(inout) :: this
  integer, intent(in) :: generation
  real :: attack_rate_signal
    
    if (Is_Evolutionary_Generations) then 
      this%risk = approx_easy(generation, 0.0, ATTACK_RATE)
    else 
      attack_rate_signal = FEAR_SIGNAL_MULTIPLIER
      this%risk = attack_rate_signal
    end if
   
end subroutine predator_init_new



function bird_is_within_environment_check(this, environment_in) result(is_within)

  class(BIRD), intent(in) :: this
  class(whole_environ), intent(in) :: environment_in
  logical :: is_within
  integer :: env_size

  env_size = size(environment_in%point)

  if (this%location%x < env_size) then
    is_within = .TRUE. 
  else
   is_within = .FALSE. 
  end if

end function bird_is_within_environment_check


subroutine predator_attack_bird(this, bird_prey, environment_in, predator_is_present, prey_is_killed)
  class(PREDATOR), intent(in) :: this
  class(BIRD), intent(inout) :: bird_prey
  class(whole_environ), intent(in) :: environment_in
  logical, optional, intent(out) :: predator_is_present
  logical, optional, intent(out) :: prey_is_killed
  logical :: p_is_present, p_prey_dies
  real :: escape_chance, fear_gene_factor

  p_is_present = .FALSE.
  p_prey_dies = .FALSE.

  if (bird_prey%is_alive) then
    p_is_present = (RAND() < environment_in%point(bird_prey%location%x)%predator_frequency)
    
    if (p_is_present) then
      if (Is_Evolutionary_Generations) then
        fear_gene_factor = real(bird_prey%gene_fear) / real(FEAR_MAX)
        escape_chance = 0.2 * fear_gene_factor  ! Adjust the 0.2 factor as needed
        p_prey_dies = (RAND() < (this%risk * (1.0 - escape_chance)))
      end if

      if (p_prey_dies) then
        bird_prey%is_alive = .FALSE.
        bird_prey%bird_meets_predator_counter = 0
        bird_prey%is_killed_by_predator = .TRUE.
      else
        bird_prey%is_alive = .TRUE.
        bird_prey%bird_meets_predator_counter = bird_prey%bird_meets_predator_counter + 1
        bird_prey%weight = bird_prey%weight - bird_prey%weight * &
          (WEIGHT_REDUCTION_CONSTANT**(bird_prey%bird_meets_predator_counter))
        
        bird_prey%state_fear = within( &
          bird_prey%state_fear + bird_prey%state_fear * FEAR_INCREMENT_AFTER_ATTACK, &
          bird_prey%genetic_lower_fear(), bird_prey%genetic_upper_fear() )
        call bird_prey%fear_genetic_limit()      
      end if
    end if
  end if

  if (present(predator_is_present)) predator_is_present = p_is_present
  if (present(prey_is_killed)) prey_is_killed = p_prey_dies
end subroutine predator_attack_bird




!****************************
!subroutine for initial gene
!****************************

!subroutine for initializing gene.
!The initial gene for hunger is randomly generated from 
!the parameters HUNGER_MIN, HUNGER_MAX, set to 0 and 10.
!The initial gene for fear is randomly generated from the
!parameters FEAR_MIN and FEAR_MAX, set to 0 and 10.
subroutine genome_init_all(this, gene_num)
    use BASE_RANDOM !importing this.
    class(GENOME), intent(out) :: this
    integer, optional, intent(in) :: gene_num 
    !gene number made for gene seperately initialized in gene_mutate(), so that genes are
    ! mutated seperately, independently of each other. 
    integer :: gene_num_loc

    if (present(gene_num)) then
      gene_num_loc = gene_num
    else
      gene_num_loc =  0
    end if

    if (gene_num_loc == 1) then
      this%gene_fear = RAND(FEAR_MIN, FEAR_MAX +1)       ! 1
    else if (gene_num_loc == 2) then
      this%gene_hunger = RAND(HUNGER_MIN, HUNGER_MAX +1)     ! 2
    else
    !classname%attribute_name
    this%gene_fear = RAND(FEAR_MIN, FEAR_MAX +1)       ! 1
    this%gene_hunger = RAND(HUNGER_MIN, HUNGER_MAX +1)     ! 2
    end if
end subroutine genome_init_all



! The gene given to each individual bird can be mutated. 
! Since we are working with a single gene, the overall 
! chance of mutations are high (set by parameter "MUTATION_PROBABILITY").
! Altered randomly (by using function RAND()).
subroutine gene_mutate(this, is_mutate)
  class(GENOME), intent(inout) :: this
  logical, optional, intent(out) :: is_mutate 
  !this local parameter is useful for debugging and testing purposes. 
  logical :: is_mutate_loc

  is_mutate_loc = .FALSE. 
  ! The standard logical value for mutaten is "false", 
  !if mutation occurs by the logic explained below then it will be "true". 

  if (RAND() < MUTATION_PROBABILITY) then 
    ! if probability for mutation is higher than the 
    !randomly generated value (RAND()), then mutation 
    !of the gene occurs. 
    call this%init_genome(1)
    is_mutate_loc = .TRUE.
  end if

  if (RAND() < MUTATION_PROBABILITY) then 
    ! if probability for mutation is higher than the 
    !randomly generated value (RAND()), then mutation of the gene occurs. 
    call this%init_genome(2)
    is_mutate_loc = .TRUE.
  end if

  if (present(is_mutate)) is_mutate = is_mutate_loc

end subroutine gene_mutate



subroutine bird_init(this, is_gene_init)

!The bird_init subroutine initializes a BIRD class object with 
!default attributes and optionally initializes its genetic attributes.
!It sets various properties of the bird, such as its weight, 
!fear/hunger state, survival status, and other relevant attributes.

!this: A BIRD class object passed by reference (intent(out)). 
!This parameter represents the bird to be initialized.
!is_gene_init: An optional logical parameter (intent(in)) 
!that indicates whether to initialize the bird's genetic attributes. 
!If not provided, the default value is .FALSE..


  class(BIRD), intent(out) :: this 
  !by calling class(BIRD), we inherit all qualities of BIRD and 
  !supply it with new atributes through new subrout.
  logical, optional, intent(in) :: is_gene_init

  logical :: is_gene_init_here

  if (present(is_gene_init)) then
    is_gene_init_here = is_gene_init
  else
    is_gene_init_here = .FALSE.
  end if

  if (is_gene_init_here) then
    call this%init_genome()
   ! print*, "GENE INIT"
  end if



  call this%put_random()
  this%weight = BIRD_INITIAL_WEIGHT
  this%weight_history = MISSING 
  !When initial weight is = 20 grams, the whole history of the array = MISSING,
  ! until the value is filled by the time_step, 
  !revealing the development of the bird by time_step and generations.
  this%state_fear = real(this%gene_fear)
  this%state_hunger = real(this%gene_hunger)
  this%is_alive = .TRUE.
  this%is_killed_by_predator = .FALSE. 
  this%bird_meets_predator_counter = 0
  this%fear_history = MISSING
  this%hunger_history = MISSING
  this%fear_gene_history = MISSING
  this%hunger_gene_history = MISSING


end subroutine bird_init


function bird_is_starved(this, threshold) result (is_dead) 

!The bird_is_starved function checks if a BIRD class object 
!has died of starvation by comparing its weight to a specified threshold. 
!If the bird's weight drops below the threshold, the function returns .TRUE.,
! indicating that the bird has died of starvation.


!this: A BIRD class object passed by reference (intent(in)). 
!This parameter represents the bird whose starvation status is to be determined.
!threshold: An optional real parameter (intent(in)) 
!that specifies the weight threshold below which the bird is 
!considered to have died of starvation. 
!If not provided, the default threshold is calculated as 
!BIRD_INITIAL_WEIGHT * BIRD_MINIMUM_WEIGHT_MULTIPLICATOR.

  class(BIRD), intent(in) :: this
  real, optional, intent(in) :: threshold
  logical :: is_dead

  real :: local_threshold

  if (present(threshold)) then
    local_threshold = threshold
  else
    local_threshold = BIRD_INITIAL_WEIGHT * BIRD_MINIMUM_WEIGHT_MULTIPLICATOR
  end if

  if (this%weight < local_threshold) then
    is_dead = .TRUE.
  else
    is_dead = .FALSE.
  end if

end function bird_is_starved





function background_mortality_probability(this) result(is_killed_by_background_mortality)
 class(BIRD), intent(inout) :: this
 logical :: is_killed_by_background_mortality

 is_killed_by_background_mortality = .FALSE.

 if (RAND() < BACKGROUND_MORTALITY) then
   is_killed_by_background_mortality = .TRUE. 
 end if
 
end function background_mortality_probability




subroutine bird_feeds(this, in_environment)
  class(BIRD), intent(inout) :: this
  class(whole_environ), intent(in) :: in_environment
  real :: current_cell_food, intake_factor, fear_factor, hunger_factor

  if (this%is_alive) then 
    current_cell_food = in_environment%point(this%location%x)%food_availability
    intake_factor = limit_food_intake(this)
    fear_factor = bird_eat_fraction_when_fear(this)
    hunger_factor = bird_eat_fraction_when_hunger(this)
    this%weight = this%weight + current_cell_food * intake_factor * fear_factor * hunger_factor
    this%state_hunger = this%state_hunger - this%state_hunger *      &
      ((intake_factor * fear_factor * hunger_factor) + (1 - hunger_factor))

    call this%hunger_genetic_limit()
  end if
end subroutine bird_feeds



function bird_eat_fraction_when_fear(this) result(fract)
  class(BIRD), intent(in) :: this
  real :: fract
  real :: midpoint, emotion, fear_gene_factor
  real :: adjusted_sigm_k

  emotion = this%state_fear
  midpoint = (this%genetic_lower_fear() + this%genetic_upper_fear()) / 2.0

  fear_gene_factor = real(this%gene_fear) / real(FEAR_MAX)
  adjusted_sigm_k = SIGM_K_FEAR * (1.0 + 0.5 * fear_gene_factor)  ! Adjust the 0.5 factor as needed

  fract = 1.0 / (1.0 + exp(adjusted_sigm_k * (emotion - midpoint)))
  fract = FRAC_S2_FEAR + (FRAC_S1_FEAR - FRAC_S2_FEAR) * fract
end function bird_eat_fraction_when_fear


function bird_eat_fraction_when_hunger(this) result(fract)
  class(BIRD), intent(in) :: this
  real :: fract
  real :: midpoint, emotion


  emotion = this%state_hunger
!SIGM_K IS STEEPNESS OF CURVE
  midpoint = (this%genetic_lower_hunger() + this%genetic_upper_hunger()) / 2.0


  if (this%weight < BIRD_INITIAL_WEIGHT * (BIRD_MAXIMUM_WEIGHT_MULTIPLICATOR/2)) then
    fract = 1.0
  else
    fract = 1.0 / (1.0 + exp(SIGM_K_HUNGER * (emotion - midpoint)))
    fract = FRAC_S2_FEAR + (FRAC_S1_HUNGER - FRAC_S2_HUNGER) * fract

  end if

 ! fract = within(fract, 0.0, 1.0)
end function bird_eat_fraction_when_hunger




subroutine bird_dead_from_starvation(this)

!The bird_dead_from_starvation subroutine checks if a BIRD class object has died
!from starvation by calling the is_starved method on the BIRD object. 
!If the bird is starved, the subroutine sets the bird's is_alive attribute to .FALSE., 
!indicating that the bird has died.

  class(BIRD), intent(inout) :: this

  if (this%is_starved()) this%is_alive = .FALSE.

end subroutine bird_dead_from_starvation


subroutine bird_dead_from_background_mortality(this)
 
 class(BIRD), intent(inout) :: this
 integer :: environment_kill_counter
 
  environment_kill_counter = 0
  if (this%environment_kills()) then
    this%is_alive = .FALSE. 
    environment_kill_counter = environment_kill_counter + 1

  end if

end subroutine bird_dead_from_background_mortality



subroutine bird_subtract_metabolic_cost(this)

!The bird_subtract_metabolic_cost subroutine adjusts the weight of a BIRD class 
!object by subtracting a fraction of its current weight,which represents the 
!metabolic cost of the bird's activities in relation to it's current weight. 
!This subroutine simulates the energy 
!expenditure of the bird, which can affect its survival and overall health. 
!The cost is applied at every timestep.

  class(BIRD), intent(inout) :: this

  this%weight = this%weight - this%weight * METABOLIC_COST_CONSTANT

end subroutine bird_subtract_metabolic_cost



function limit_food_intake(this) result(intake_factor)
  class(BIRD), intent(inout) :: this
  real :: intake_factor, sigmoid_factor
  real :: weight_ratio, max_weight, min_weight
  real :: gene_influence
  real, parameter :: GENE_SCALE_FACTOR = 0.0002 ! Adjust as needed

  max_weight = BIRD_INITIAL_WEIGHT * BIRD_MAXIMUM_WEIGHT_MULTIPLICATOR
  min_weight = BIRD_INITIAL_WEIGHT * BIRD_MINIMUM_WEIGHT_MULTIPLICATOR

  weight_ratio = (this%weight - min_weight) / (max_weight - min_weight)

  ! Calculate gene influence on stomach size
  gene_influence = 1 + GENE_SCALE_FACTOR * (real(this%gene_hunger) / real(HUNGER_MAX))

  ! Sigmoid function
  sigmoid_factor = 1.0 / (1.0 + exp(SIGMOID_STEEPNESS * (weight_ratio - SIGMOID_MIDPOINT)))

  ! Scale the sigmoid output to the desired range and apply gene influence
  intake_factor = 0.01 * sigmoid_factor * gene_influence

  ! Ensure the intake factor is within the desired range
  intake_factor = max(0.0, min(0.02, intake_factor))

  ! Increase fear based on weight ratio
  this%state_fear = this%state_fear + (weight_ratio * VIGILANCE_FACTOR)

  ! Ensure state_fear stays within genetic limits
  call this%fear_genetic_limit()
end function limit_food_intake


! function limit_food_intake(this) result(intake_factor)
!   class(BIRD), intent(inout) :: this
!   real :: intake_factor
!   real :: weight_ratio, max_weight, min_weight
!   real, parameter :: MAX_INTAKE = 0.1
!   real, parameter :: MIN_INTAKE = 0.05
!   real, parameter :: VIGILANCE_FACTOR = 0.5 ! Adjust this value to control fear increase

!   max_weight = BIRD_INITIAL_WEIGHT * BIRD_MAXIMUM_WEIGHT_MULTIPLICATOR
!   min_weight = BIRD_INITIAL_WEIGHT * BIRD_MINIMUM_WEIGHT_MULTIPLICATOR

!   weight_ratio = (this%weight - min_weight) / (max_weight - min_weight)

!   ! Linear function
!   intake_factor = MAX_INTAKE * (1.0 - weight_ratio) + MIN_INTAKE

!   ! Ensure the intake factor is within the desired range
!   intake_factor = max(MIN_INTAKE, min(MAX_INTAKE + MIN_INTAKE, intake_factor))

!   ! Increase state_fear based on food intake
!   this%state_fear = this%state_fear + (intake_factor * VIGILANCE_FACTOR)

!   ! Ensure state_fear stays within genetic limits
!   call this%fear_genetic_limit()

! end function limit_food_intake



subroutine bird_do_fly(this, in_environment)
!This subroutine is designed to simulate the behavior of a bird (this) 
!in response to its environment (in_environment). The bird's behavior 
  !is influenced by its hunger state, specifically its fear and hunger levels.
!The subroutine first checks if the bird is alive by evaluating this%is_alive. 
!If the bird is not alive, the subroutine ends without performing any actions.
!If the bird is alive, it calculates min_loc and max_loc based on the size of 
! the point array in the in_environment object.

!The subroutine then checks the bird's fear and hunger state by evaluating 
!this%state_fear. If the birds state of fear has a value higher than 
!fear_fly_threshold (range 0-10), the bird will fly and not eat. 
! This is simulated by calling the walk method of the BIRD class
!with min_loc and max_loc as arguments.

!If the bird's fear and hunger state is not within the specified range, 
!the subroutine does not perform any action, implying that the bird will 
!eat instead of moving. 

!This subroutine is designed to simulate the behavior of a bird (this) 
!in response to its environment (in_environment). The bird's behavior 
  !is influenced by its hunger state, specifically its fear and hunger levels.
!The subroutine first checks if the bird is alive by evaluating this%is_alive. 
!If the bird is not alive, the subroutine ends without performing any actions.
!If the bird is alive, it calculates min_loc and max_loc based on the size of 
! the point array in the in_environment object.

!The subroutine then checks the bird's fear and hunger state by evaluating 
!this%state_fear. If the birds state of fear has a value higher than 
!fear_fly_threshold (range 0-10), the bird will fly and not eat. 
! This is simulated by calling the walk method of the BIRD class
!with min_loc and max_loc as arguments.

!If the bird's fear and hunger state is not within the specified range, 
!the subroutine does not perform any action, implying that the bird will 
!eat instead of moving. 

  class(BIRD), intent(inout) :: this
  class(whole_environ), intent(in) :: in_environment

  real :: fear_fly_threshold

  integer :: min_loc, max_loc

  fear_fly_threshold = this%genetic_upper_fear() * 0.9

  
  if (this%is_alive) then

    min_loc = 1 
    max_loc = size(in_environment%point) 

    if (this%state_fear > fear_fly_threshold ) then
      call this%walk(min_loc, max_loc)
    end if

  end if

end subroutine bird_do_fly

subroutine bird_do_explore (this, in_environment)
  class(BIRD), intent (inout) :: this
  class(whole_environ), intent(in) :: in_environment

  integer :: min_loc, max_loc

  min_loc = 1
  max_loc = size(in_environment%point)

  if (this%is_alive) then
    if(RAND() < EXPLORATION_PROBABILITY) then
      call this%walk(min_loc, max_loc)
    end if
  end if


end subroutine bird_do_explore



!The two functions below calculate the minimum and maximum genetic limits for a bird's hunger, 
!which is used to determine the bird's hunger state. 
!The function is designed to be called on a BIRD class object and returns a 
!real number representing the lower limit of the bird's genetic hunger range.

!this: A class object of type BIRD. This parameter is passed with the intent(in) attribute, 
!indicating that the function can only read from this object but cannot modify it.

!lower_limit: A real number representing the minimum genetic limit for the bird's hunger.

!range: A real number representing the difference between the maximum and minimum hunger values 
!(HUNGER_MAX and HUNGER_MIN).

!The function calculates the range by subtracting HUNGER_MIN from HUNGER_MAX.
!It then calculates the lower_limit by subtracting the product of range and HUNGER_PLASTICITY from the 
!bird's genetic fear of hunger (this%gene_fear).
!The lower_limit is then adjusted to ensure it falls within the defined range of hunger values 
!(HUNGER_MIN and HUNGER_MAX) using a hypothetical within function. This function is not defined 
!within the provided code snippet, so its specific implementation and behavior are not detailed here.
!Finally, the lower_limit is returned as the result of the function.

!The two functions below calculate the minimum and maximum genetic limits for a bird's hunger, 
!which is used to determine the bird's hunger state. 
!The function is designed to be called on a BIRD class object and returns a 
!real number representing the lower limit of the bird's genetic hunger range.

!this: A class object of type BIRD. This parameter is passed with the intent(in) attribute, 
!indicating that the function can only read from this object but cannot modify it.

!lower_limit: A real number representing the minimum genetic limit for the bird's hunger.

!range: A real number representing the difference between the maximum and minimum hunger values 
!(HUNGER_MAX and HUNGER_MIN).

!The function calculates the range by subtracting HUNGER_MIN from HUNGER_MAX.
!It then calculates the lower_limit by subtracting the product of range and HUNGER_PLASTICITY from the 
!bird's genetic fear of hunger (this%gene_fear).
!The lower_limit is then adjusted to ensure it falls within the defined range of hunger values 
!(HUNGER_MIN and HUNGER_MAX) using a hypothetical within function. This function is not defined 
!within the provided code snippet, so its specific implementation and behavior are not detailed here.
!Finally, the lower_limit is returned as the result of the function.
function bird_hunger_min_genetic_limit(this) result(lower_limit)
  class(BIRD), intent(in) :: this
  real :: lower_limit

  real :: range

  range = 1

  lower_limit = this%gene_hunger - range * HUNGER_PLASTICITY

  lower_limit = within(lower_limit, real(HUNGER_MIN), real(HUNGER_MAX))


end function bird_hunger_min_genetic_limit


function bird_hunger_max_genetic_limit(this) result(upper_limit)
  class(BIRD), intent(in) :: this
  real :: upper_limit

  real :: range

  range = 1

  upper_limit = this%gene_hunger + range * HUNGER_PLASTICITY

  upper_limit = within(upper_limit, real(HUNGER_MIN), real(HUNGER_MAX))

end function bird_hunger_max_genetic_limit


!! FEAR
!The function is designed to be called on a BIRD class object and returns a 
!real number representing the lower limit of the bird's genetic fear range.

!this: A class object of type BIRD. This parameter is passed with the intent(in) attribute, 
!indicating that the function can only read from this object but cannot modify it.

!lower_limit: A real number representing the minimum genetic limit for the bird's fear.

!range: A real number representing the difference between the maximum and minimum fear values 
!(FEAR_MAX and FEAR_MIN).

!The function calculates the range by subtracting FEAR_MIN from FEAR_MAX.
!It then calculates the lower_limit by subtracting the product of range and FEAR_PLASTICITY from the 
!bird's genetic fear (this%gene_fear).
!The lower_limit is then adjusted to ensure it falls within the defined range of fear values 
!(FEAR_MIN and FEAR_MAX) using a hypothetical within function. This function is not defined 
!within the provided code snippet, so its specific implementation and behavior are not detailed here.
!Finally, the lower_limit is returned as the result of the function.
function bird_fear_min_genetic_limit(this) result(lower_limit)
  class(BIRD), intent(in) :: this
  real :: lower_limit

  real :: range

  range = 1 

  lower_limit = this%gene_fear - range * FEAR_PLASTICITY

  lower_limit = within(lower_limit, real(FEAR_MIN), real(FEAR_MAX))


end function bird_fear_min_genetic_limit



!The function is designed to be called on a BIRD class object and returns a 
!real number representing the upper limit of the bird's genetic fear range.

!this: A class object of type BIRD. This parameter is passed with the intent(in) attribute, 
!indicating that the function can only read from this object but cannot modify it.

!upper_limit: A real number representing the maximum genetic limit for the bird's fear.

!range: A real number representing the difference between the maximum and minimum fear values 
!(FEAR_MAX and FEAR_MIN).

!The function calculates the range by subtracting FEAR_MIN from FEAR_MAX.
!It then calculates the upper_limit by adding the product of range and FEAR_PLASTICITY to the 
!bird's genetic fear (this%gene_fear).
!The upper_limit is then adjusted to ensure it falls within the defined range of fear values 
!(FEAR_MIN and FEAR_MAX) using a hypothetical within function. This function is not defined 
!within the provided code snippet, so its specific implementation and behavior are not detailed here.
!Finally, the upper_limit is returned as the result of the function.
function bird_fear_max_genetic_limit(this) result(upper_limit)
  class(BIRD), intent(in) :: this
  real :: upper_limit

  real :: range

  range = 1 

  upper_limit = this%gene_fear + range * FEAR_PLASTICITY

  upper_limit = within(upper_limit, real(FEAR_MIN), real(FEAR_MAX))


end function bird_fear_max_genetic_limit






subroutine bird_force_fear_within_genetic_limits(this)
! setting a limited variance for genetic dispositions in birds.
! This way there will be heterogeneity in the genetic compositions.
 class(BIRD), intent(inout) :: this

  this%state_fear = within(this%state_fear,                         &
                                  this%genetic_lower_fear(),                    &
                                  this%genetic_upper_fear() )

end subroutine bird_force_fear_within_genetic_limits



subroutine bird_force_hunger_within_genetic_limits(this)
! setting a limited variance for genetic dispositions in birds.
! This way there will be heterogeneity in the genetic compositions.
 class(BIRD), intent(inout) :: this

  this%state_hunger = within(this%state_hunger,                         &
                                  this%genetic_lower_hunger(),                    &
                                  this%genetic_upper_hunger() )

end subroutine bird_force_hunger_within_genetic_limits

subroutine bird_decay_fear(this)
  class(BIRD), intent(inout) :: this
  this%state_fear = this%state_fear - this%state_fear * FEAR_DECAY_RATE
end subroutine bird_decay_fear


!Subroutine for history array

subroutine bird_weight_add_to_history(this, time_step)
  class(BIRD), intent(inout) :: this
  integer, intent(in) :: time_step

  associate ( A => this%weight_history, W => this%weight, i => time_step)
   A(i) = W
  end associate

end subroutine bird_weight_add_to_history





!> Adds the current fear and hunger state of a BIRD object to its history arrays.
!>
!> @param this The BIRD object whose fear and hunger state history is being recorded.
!> @param time_step The current time step, used as the index for the history arrays.
subroutine bird_emotion_state_add_to_history(this, time_step)
  class(BIRD), intent(inout) :: this
  integer, intent(in) :: time_step

  this%fear_history(time_step) = this%state_fear
  this%hunger_history(time_step) = this%state_hunger

end subroutine bird_emotion_state_add_to_history

subroutine bird_emotion_gene_add_to_history(this, time_step)
  class(BIRD), intent(inout) :: this
  integer, intent(in) :: time_step
  this%fear_gene_history(time_step) = this%gene_fear
  this%hunger_gene_history(time_step) = this%gene_hunger

end subroutine bird_emotion_gene_add_to_history


!> This subroutine simulates the time steps for a population of birds in an environment, including their behavior and interactions with predators.
!>
!> @param this The POPULATION object containing the birds.
!> @param environment_in The WHOLE_ENVIRON object representing the environment the birds are in.
!> @param predator_in The PREDATOR object representing the predators in the environment.
!>
!> The subroutine iterates through each bird in the population and performs the following actions for each time step:
!> - Checks if the bird is alive
!> - Prints debug information if the IS_DEBUG_SCREEN flag is set
!> - Checks if the bird has starved and kills it if so
!> - Calls the bird's fly, feed, pay_for_life, do_explore, and is_killed_by_background subroutines
!> - Calls the predator's predator_attack_bird subroutine to simulate predator attacks
!> - Adds the bird's weight and emotional state to its history if the IS_DEBUG_DATA flag is set
!>
!> The subroutine also keeps track of the number of birds killed by background mortality.
subroutine population_time_steps(this, environment_in, predator_in) !bird is in this environment = environment_in
  class(POPULATION), intent(inout) :: this
  class(whole_environ), intent(in) :: environment_in
  class(PREDATOR), intent(inout) :: predator_in

  integer :: current_step
  integer :: bird_current !this current bird
  integer :: environment_kill_counter


  logical :: debug_check , debug_ckeck2
   
  environment_kill_counter = 0


  do bird_current = 1, size(this%birds)

    do current_step=1, TOTAL_TIME_STEP

      if (this%birds(bird_current)%is_alive) then


        DEBUG_GET_INFO: if (IS_DEBUG_SCREEN) then
                        print*, " "
                        print *, "Bird", bird_current, "at step", current_step
                        if (this%birds(bird_current)%is_starved()) then
                          print*, "STARVED"
                        end if
                        print *, "Location:", this%birds(bird_current)%x, ",   food",                   &
                                  environment_in%point( this%birds(bird_current)%x )%food_availability

                        print *, "Location:", this%birds(bird_current)%x, ",   risk",           &
                                  environment_in%point( this%birds(bird_current)%x )%predator_frequency


                        print *, "state fear   ", this%birds(bird_current)%state_fear, &
                                          "[", this%birds(bird_current)%genetic_lower_fear(),     &
                                              this%birds(bird_current)%genetic_upper_fear(), "]"


                        print *, "state hunger ", this%birds(bird_current)%state_hunger, &
                                          "[", this%birds(bird_current)%genetic_lower_hunger(),     &
                                              this%birds(bird_current)%genetic_upper_hunger(), "]"

                        print *, " mass", this%birds(bird_current)%weight
                        
                      
                        if (this%birds(bird_current)%environment_kills()) then
                          print*, "*************************************************************************************"
                          print*, "                  BIRD IS KILLED BY BACKGROUND MORTALITY                             "
                          print*, "*************************************************************************************"
                          environment_kill_counter = environment_kill_counter + 1
                        end if
                        print*, "TOTAL BACKGROUND MORTALITY KILL: ", environment_kill_counter

                        

                    end if DEBUG_GET_INFO

        call this%birds(bird_current)%killed_from_starvation() !is the bird dead or alive from starvation?
        !what birds do if is alive
       
        call this%birds(bird_current)%fly(environment_in) !bird is flying

        call predator_in%predator_attack_bird(this%birds(bird_current), environment_in,                   &
              predator_is_present=debug_ckeck2, prey_is_killed=debug_check)! predator attacks bird

        DEBUG_PRINT_KILL: if (IS_DEBUG_SCREEN) then
                                  print *, "Predator ", debug_ckeck2
                                  if (debug_check) print *, "KILLED by predator"
                                end if DEBUG_PRINT_KILL
        call this%birds(bird_current)%do_feed(environment_in)!birds is feeding, and weight is incrementing
        call this%birds(bird_current)%pay_for_life() !bird looses some weight per time step due to metabolism.
        call this%birds(bird_current)%do_explore(environment_in) !bird may explore by moving to another habitat. 
        call this%birds(bird_current)%is_killed_by_background()!background mortality    
        call this%birds(bird_current)%decay_fear()
        
        !Ensure genetic limits for fear and hunger
        call this%birds(bird_current)%hunger_genetic_limit()
        call this%birds(bird_current)%fear_genetic_limit()




        DEBUG_GET_HISTORY: if (IS_DEBUG_DATA) then
                                    call this%birds(bird_current)%add_to_history(current_step)
                                    call this%birds(bird_current)%emo_state_to_history(current_step) !both fear and hunger within emotion
                            end if DEBUG_GET_HISTORY

        else
                         if (IS_DEBUG_SCREEN) print *, bird_current , "DEAD"
      end if


    end do
  end do

end subroutine population_time_steps

subroutine sort_by_fitness(this)!sorting by bodymass
  class(POPULATION), intent(inout) :: this

  !-----------------------------------------------------------------------------
  ! Sort parents by fitness:

  ! This subroutine sorts the birds in a POPULATION object by their fitness, 
  ! which is represented by their body mass. The sorting is performed in descending
  ! order, so that the fittest birds (those with the highest body mass) are placed 
  ! first in the array.

   ! Initial Sorting: The subroutine begins by calling QsortC on the birds array of the POPULATION object. 
   ! This initiates the QuickSort algorithm, which sorts the birds based on their weight.

   ! Reversing the Order: After the initial sorting, the subroutine reverses the order of the birds array. 
   ! This is done to ensure that the fittest birds (those with the highest body mass) are placed first in the array.

   ! QuickSort Algorithm: The QsortC subroutine is a recursive implementation of the QuickSort algorithm, 
   ! adapted for sorting BIRD objects based on their weight. It uses a helper subroutine, Partition, to partition 
   ! the array around a pivot point, which is the weight of the first bird in the array.

   ! Partition Subroutine: The Partition subroutine partitions the array of BIRD objects around a pivot point, 
   ! which is the weight of the first bird in the array. It uses two indices, i and j, to track the boundaries 
   ! of the partitions. The subroutine iterates through the array, swapping elements as necessary to ensure that 
   ! all elements less than the pivot are to its left, and all elements greater than the pivot are to its right.
  !-----------------------------------------------------------------------------

  call QsortC(this%birds)

  !Reverse the sorting so that the fittest birds(highest mass) go first (10 -> 1).
  this%birds = this%birds(size(this%birds):1:-1)

  contains
    ! The two subroutines below are a variant of the recursive QuickSort
    ! algorithm adapted for BIRD fitness (real)

    recursive subroutine QsortC(A)

      type(BIRD), intent(in out), dimension(:) :: A
      integer :: iq

      if(size(A) > 1) then
        call Partition(A, iq)
        call QsortC(A(:iq-1))
        call QsortC(A(iq:))
      endif

    end subroutine QsortC

    subroutine Partition(A, marker)

      type(BIRD), intent(in out), dimension(:) :: A
      integer, intent(out) :: marker
      integer :: i, j
      type(BIRD) :: temp
      real :: x      ! pivot point

      x = A(1)%weight
      i= 0
      j= size(A) + 1

      do
        j = j-1
        do
            if (A(j)%weight <= x) exit
            j = j-1
        end do
        i = i+1
        do
            if (A(i)%weight >= x) exit
            i = i+1
        end do
        if (i < j) then
            ! exchange A(i) and A(j)
            temp = A(i)
            A(i) = A(j)
            A(j) = temp
        elseif (i == j) then
            marker = i+1
            return
        else
            marker = i
            return
        endif
      end do

    end subroutine Partition

end subroutine sort_by_fitness

end module organism




module factor_reader

  use params
  implicit none

  ! Define the Fecundity_factor_def type
  type :: Fecundity_factor_def
      real :: signal
      real, dimension(:), allocatable :: ffactor
  end type Fecundity_factor_def

  ! Global variable to store fecundity factors
  type(real), dimension(:), allocatable, public :: Glob_Fecundity_Signal
  type(real), dimension(:,:), allocatable, public :: Glob_Fecundity_Factor



contains

  ! Function to get fecundity factor based on signal and generation
  !> Retrieves the fecundity factor for the given signal and generation.
  !!
  !! This function searches the global `Glob_Fecundity_Signal` and `Glob_Fecundity_Factor` 
  ! arrays to find the fecundity factor that corresponds to the provided signal and generation. 
  ! It uses a tolerance of `TOLER_DIFF` to match the signal value.
  !!
  !! @param signal The signal value to match.
  !! @param generation The generation index to retrieve the fecundity factor for.
  !! @return The fecundity factor value, or -9999.9 if no match is found.

  function get_ffactor(signal, generation) result(ffact_out)
      real, intent(in) :: signal
      integer, intent(in) :: generation
      real :: ffact_out
      integer :: i, correct_index_col
      real, parameter :: TOLER_DIFF = 0.01

      ffact_out = -9999.9


      do i= 1, size (Glob_Fecundity_Signal)
        if (abs(signal - Glob_Fecundity_Signal(i)) < TOLER_DIFF ) then 
          correct_index_col = i 
          ffact_out = Glob_Fecundity_Factor(generation, correct_index_col)
          exit
        end if 
      end do 



   !   print*, "ZZZZ", "Index col: ", correct_index_col
  end function get_ffactor


  
  

 !Subroutine to read matrix data from CSV files
 !> Reads a matrix of fecundity factor data from a CSV file.
 !!
 !! This subroutine reads a matrix of fecundity factor data from a CSV file specified by
  ! the `filename_ffact` argument. The data is stored in the global `Glob_Fecundity_Signal` 
  !and `Glob_Fecundity_Factor` arrays.
 !!
 !! @param filename_ffact The filename of the CSV file containing the fecundity factor data.
 subroutine read_matrix(filename_ffact)
  use CSV_IO

    character(len=*), intent(in) ::  filename_ffact
    real, dimension(:,:), allocatable :: tmp_signal, tmp_ffact
    integer :: i
    logical :: errorflag
    integer, parameter :: SKIP_COLS = 1
	
      integer :: fact_rows, fact_cols

      if (IS_DEBUG_ARRAY_READ) then 
        print *, "DEBUG: Entering read_matrix subroutine"
        print *, "DEBUG: filename_ffact = ", filename_ffact
      end if

      
      ! Read fecundity factor data
      if (IS_DEBUG_ARRAY_READ) print *, "DEBUG: Reading fecundity factor file"
      tmp_ffact = CSV_MATRIX_READ(filename_ffact, errorflag)
      if (errorflag .neqv. .TRUE.) then
          print *, "ERROR: Error reading fecundity factor file"
          return
      end if


      fact_rows = size(tmp_ffact,1)-1
      fact_cols = size(tmp_ffact,2)-1

      print *, "READING ", fact_rows, fact_cols

      if (allocated(Glob_Fecundity_Signal)) deallocate(Glob_Fecundity_Signal)
      if (allocated(Glob_Fecundity_Factor)) deallocate(Glob_Fecundity_Factor)

      allocate(Glob_Fecundity_Signal(fact_cols))
      allocate(Glob_Fecundity_Factor(fact_rows,fact_cols))

      do i = 1, fact_cols
        Glob_Fecundity_Signal(i)= tmp_ffact(1,i+SKIP_COLS)
       if (IS_DEBUG_ARRAY_READ) print *,"XXX", Glob_Fecundity_Signal(i)
      end do

      if (IS_DEBUG_ARRAY_READ) then 
        print *, "DEBUG: Fecundity factor file read successfully"

        ! Populate glob_fecundity_factor with read data
        print *, "DEBUG: Populating glob_fecundity_factor, ", fact_rows, " rows"
      end if
      do i = 1, fact_rows
          Glob_Fecundity_Factor(i,1:fact_cols) = tmp_ffact(i+1,1+SKIP_COLS:fact_cols+SKIP_COLS)
          ! accessing the i-th row and columns 1 to 5
          if (IS_DEBUG_ARRAY_READ) print*, tmp_ffact(i+1,1:fact_cols)
          if (IS_DEBUG_ARRAY_READ) print *, ">>>", Glob_Fecundity_Factor(i,:)
      end do


    ! Deallocate arrays
     ! deallocate(Glob_Fecundity_Signal)
  end subroutine read_matrix



end module factor_reader



module GA !genetic algorythm
use params
use environment
use BASE_RANDOM
use CSV_IO
use organism
use factor_reader
use, intrinsic :: ISO_FORTRAN_ENV, only : OUTPUT_UNIT, ERROR_UNIT
implicit none


type(whole_environ) :: habitat

type(POPULATION),target :: generation_one
type(POPULATION),target :: generation_two

type(POPULATION), pointer :: parent_generation
type(POPULATION), pointer :: offspring_generation

type(PREDATOR) :: predator_in_habitat


!---------------------------------------------
! Making the output data:
!---------------------------------------------
! collumns
! 1: generation nr (int)
! 2: Num alive (int)
! 3: average mass all (real)
! 4: average mass alive (real)
! 5: best mass(real)
! 6: average mass top 25% (real)
! 7: Average gene alive (real)
! 8: Average emotion alive (real)
! 8: standard deviation of mass for alive (real))
! 9: Frequency of predator (real)



integer, parameter :: CSV_COL_NUMBER = 25
character(len=*), dimension(CSV_COL_NUMBER), parameter :: CSV_COLUMNS =    &
  ["GENERATION                   ",                 &  !1
   "NUM-ALIVE                    ",                 &  !2
   "AVERAGE-MASS-ALL             ",                 &  !3
   "AVERAGE-MASS-ALIVE           ",                 &  !4
   "BEST-MASS                    ",                 &  !5
   "AVERAGE-MASS-DOMINANTS       ",                 &  !6
   "MEDIAN MASS ALIVE            ",                 &  !7
   "AVERAGE-GENE-FEAR-ALIVE      ",                 &  !8
   "AVERAGE-GENE-HUNGER-ALIVE    ",                 &  !9
   "AVERAGE-FEAR-STATE-ALIVE     ",                 &  !10
   "AVERAGE-HUNGER-STATE-ALIVE   ",                 &  !11
   "STD-DEV-MASS-ALIVE           ",                 &  !12
   "BIRD-MEETS-PREDATOR-AVERAGE  ",                 &  !13
   "NUM-KILLED-BY-PREDATOR       ",                 &  !14
   "NUM-KILLED-BY-PREDATOR-AVG   ",                 &  !15
   "MUTATION PROBABILITY         ",                 &  !16
   "AVERAGE-GENE-FEAR-DEAD       ",                 &  !17
   "AVERAGE-GENE-HUNGER-DEAD     ",                 &  !18
   "AVERAGE-STATE-FEAR-DEAD      ",                 &  !19
   "AVERAGE-STATE-HUNGER-DEAD    ",                 &  !20
   "STD-DEV-GENE_FEAR-ALIVE      ",                 &  !21
   "STD-DEV-GENE_HUNGER-ALIVE    ",                 &  !22
   "AVEREAGE-GENE-FEAR-ALL       ",                 &  !23
   "AVERAGE-PROB-REPRODUCTION    ",                 &  !24
   "Regain Pop when multiply with"]



real, dimension(GENERATIONS, CSV_COL_NUMBER) :: csv_generation_output

private :: CSV_COL_NUMBER, CSV_COLUMNS, csv_generation_output




contains


subroutine save_gene_history(population_birds, gen_n)
  use params
  use BASE_UTILS, only : TOSTR
  use CSV_IO

  class(POPULATION), intent(in) :: population_birds
  integer, intent(in) :: gen_n

  character(len=:), allocatable :: file_name_fear_gene
  character(len=:), allocatable :: file_name_hunger_gene

  real, dimension(POP_SIZE) :: fear_gene_population
  real, dimension(POP_SIZE) :: hunger_gene_population

  integer :: i

  file_name_fear_gene = EMOTION_GENE_SUBFOLDER // "fear_gene_gen_" // TOSTR(gen_n) // ".csv"
  file_name_hunger_gene = EMOTION_GENE_SUBFOLDER // "hunger_gene_gen_" // TOSTR(gen_n) // ".csv"

  do i = 1, POP_SIZE
    fear_gene_population(i) = population_birds%birds(i)%gene_fear
    hunger_gene_population(i) = population_birds%birds(i)%gene_hunger
  end do

  call CSV_MATRIX_WRITE(matrix = reshape(fear_gene_population, [POP_SIZE, 1]), &
                        csv_file_name = file_name_fear_gene, &
                        colnames = ["FEAR_GENE"])

  call CSV_MATRIX_WRITE(matrix = reshape(hunger_gene_population, [POP_SIZE, 1]), &
                        csv_file_name = file_name_hunger_gene, &
                        colnames = ["HUNGER_GENE"])

end subroutine save_gene_history



subroutine genetic_algorithm()
!This subroutine simulates the evolution of a population of birds over multiple 
!generations using a genetic algorithm. The algorithm involves initializing the 
!population, running through time steps, sorting by fitness, selecting the best 
!individuals for reproduction, and swapping generations.

!Initialization: The subroutine begins by initializing the random seed, the habitat, 
!and the predator presence within the habitat. It also initializes two generations 
!of birds, parent_generation and offspring_generation, and creates directories for 
!saving weight history and emotion state data.


!Main Loop: The subroutine enters a loop that runs for a specified number of generations. 
  !In each generation, it performs the following steps:
  ! Initialization: If it's not the first generation, the parent generation is re-initialized.
  ! Time Steps: The birds in the parent generation go through time steps, interacting with the habitat
  !  and predators.
  ! Predator Initialization: The predator presence is updated for the current generation.
  ! Sorting by Fitness: The birds are sorted by their fitness, which is determined by their body mass.
  ! Data Saving: If specified, the current generation's data is saved to a CSV file.
  ! Output Data Building: Characteristics of the parent generation are recorded for output.
  ! Check for Surviving Birds: The subroutine checks if there are enough surviving birds. If not, 
  ! it writes a warning to the error unit and stops the simulation.
  ! Selection and Reproduction: The best individuals are selected for reproduction, 
  ! and the next generation is produced.
  ! Generation Swap: The roles of the parent and offspring generations are swapped for the next generation.

!Finalization: After all generations have been processed, the subroutine writes the final 
  !generation's data to a CSV file 
!   optionally saves the habitat's nutritional value data if debugging is enabled.


  use params

  use BASE_UTILS, only : TOSTR

  real, dimension(POP_SIZE) :: fitness

  integer :: n_alive_get, n_of_generations
  
  csv_generation_output = MISSING

  call RANDOM_SEED_INIT()

  call habitat%init()

  Full_Pop_Factor = MISSING
  Last_Full_Pop_Factor = MISSING

  !!!!!call predator_in_habitat%init(Current_Generation)
  
  Is_Evolutionary_Generations = .TRUE. 
  !Comment in if doing experiments. 

  call get_global_runtime_params()

  if (Is_Evolutionary_Generations) then
    n_of_generations = GENERATIONS
  else
    n_of_generations = SIGNAL_ONLY_GENERATIONS
  end if

  parent_generation => generation_one
  offspring_generation => generation_two


  call parent_generation%init_pop(is_gene_init = .TRUE.)
  call offspring_generation%init_pop(is_gene_init = .TRUE.)


!First generation of experimental generations
  if (.NOT. Is_Evolutionary_Generations) then
    call parent_generation%load(Genome_File_Name)
    !call offspring_generation%load(Genome_File_Name
  end if


  call FS_MKDIR(WEIGHT_HISTORY_SUBFOLDER)
  call FS_MKDIR(EMOTION_STATE_SUBFOLDER)
  call FS_MKDIR(EMOTION_GENE_SUBFOLDER)

! Save gene history for the first generation
  call save_gene_history(parent_generation, Current_Generation)
  do Current_Generation = 1, n_of_generations
    
   ! if (.NOT. Is_Evolutionary_Generations) call parent_generation%load(Genome_File_Name)

    write(OUTPUT_UNIT, *) "Generation: ", Current_Generation

    call predator_in_habitat%init(Current_Generation)


    if (Current_Generation < 1 .AND. Is_Evolutionary_Generations) then     
      call parent_generation%init_pop(is_gene_init = .TRUE.)
      ! Saving gene history for the generations > 1
      call save_gene_history(parent_generation, Current_Generation)
    end if


    call parent_generation%time_steps(habitat, predator_in_habitat)
      
    !sorting by fitness
    call parent_generation%sort()
    
    SAVE: block

    character(len=500) :: generation_data_file
    character(*), parameter :: GENERATION_SUBFOLDER = "generations/"
    call FS_MKDIR(GENERATION_SUBFOLDER)
    generation_data_file = GENERATION_SUBFOLDER // "gen_" // TOSTR(Current_Generation) // ".csv"
    if (IS_SAVE_POPULATION) call parent_generation%save(generation_data_file)

    end block SAVE


    !record parent characteristics to output file
    call build_output_data(Current_Generation)

    n_alive_get = parent_generation%alive_count()

    write(OUTPUT_UNIT, *)                                                                   &
        "|Num alive: ", n_alive_get,                                                        &
        "|  avg. mass alive:         ", csv_generation_output(Current_Generation, 4),       &
        "|  avg. gene fear alive:    ", csv_generation_output(Current_Generation, 8),       &
        "|  avg. gene hunger alive:  ", csv_generation_output(Current_Generation, 9),      &                                                     
        "|  avg. fear-state alive:   ", csv_generation_output(Current_Generation, 10),       &
        "|  avg. hunger-state alive: ", csv_generation_output(Current_Generation, 11)!,    &
    !    "     Mutation probability: ", MUTATION_PROBABILITY,                                &
    !    "     attack rate: ", predator_in_habitat%risk
 

      call save_history_array(parent_generation, Current_Generation)


    if (n_alive_get < 10) then
      write(ERROR_UNIT, *) "WARNING: Too few birds remain in the generation"

      call CSV_MATRIX_WRITE(csv_generation_output,                      &
                            MODEL_OUTPUT_FILE,colnames=CSV_COLUMNS)
      stop
    end if

    call select_and_reproduce()


    if (Is_Evolutionary_Generations) then
      call generation_swap()
    else
      ! If experimental generations 
      parent_generation = reset_generation()
      call minimum_weight_exp_gens
    end if

  end do

call CSV_MATRIX_WRITE(csv_generation_output,                      &
                      MODEL_OUTPUT_FILE,colnames=CSV_COLUMNS)


! calling function to save nutritional value of each habitat upon environment initiation.
                      ! currently restarting with both predation and nutritional cloning, 
                      ! effectively canceling the signal-increase in experimental generations.
                      ! Therefore habitat cloning is not used during experimental generations 
                      ! per august 23rd 2024.
if (IS_DEBUG_DATA) then
  call habitat%save()
end if


print*, "SIMULATION COMPLETE"

end subroutine genetic_algorithm



! Renamed from save_history_array_weight, since the former is less general.
! We did this since we wanted a general subroutine for saving "history" of
! birds walk and all its updated parameters for each step.
subroutine save_history_array(population_birds, gen_n)
!This subroutine is designet to save the weight and emotion numerical data depicting 
!the birds mass and emotional state at each time step. It takes two arguments: 
!  population_birds, which is of the POPULATION class, and gen_n which is an integer 
!  value representing the generation. 

! 1. The subroutine begins by using the BASE_UTILS module for string conversion (TOSTR) 
! and declaring several variables and arrays. These include: 
!   var_name: An array of character strings to hold the names of the time steps
!   file_name_weight_history and file_name_emo_history: Allocatable character strings 
!   for the file names where the weight and emotion history will be saved. 
!   weigh_history_population and emotion_history_population: 2D arrays to store the 
!   weight and emotion history of each bird in the population. 
!
! 2. Debug Data Check: The subroutine check if IS_DEBUG_DATA is true. If it is, the 
!  subroutine proceeds to generate file names for the weight and emotion history files, 
!  appending the generation number to the file names. 
!
! 3. Generating Variable Arrays: For each time step, a variable name is generated in the 
!  format STEP_ followed by the time step number. This is done using the TOSTR function, 
!  which is presumably defined in the BASE_UTILS module.
!
! 4. Population History Arrays: The subroutine then iterates over each bird in the population 
!   and each time step, populating the weight_history_population and emotion_history_population 
!   arrays with the corresponding data from the population_birds object. 
!
! 5. Writing to CSV Files: the CSV_MATRIX_WRITE subroutine (not shown in this subroutine) 
!    is called to write the weight and emotion history arrays to CSV files. The file name 
!    and column names (generated variable names) are passed as arguments. 
!
! 6. File Compression: If IS_ZIP_DATA is true, the subroutine compresses the generated 
!    CSV-files using a command line tool specified by ZIP_PROG. WAIT = .FALSE. simply 
!    means that the subroutine is allowed to continue regardless of the compression done. 

  use BASE_UTILS, only : TOSTR
  use params

  class(POPULATION), intent(in) :: population_birds
  integer, intent(in) :: gen_n

  !                 STEP_100000
  !                 12345678901
  character(len=11), dimension(TOTAL_TIME_STEP) :: var_name
  character(len=:), allocatable :: file_name_weight_history
  character(len=:), allocatable :: file_name_fear_history
  character(len=:), allocatable :: file_name_hunger_history
  character(len=:), allocatable :: file_name_fear_gene_history
  character(len=:), allocatable :: file_name_hunger_gene_history

  real, dimension(POP_SIZE,TOTAL_TIME_STEP) :: weight_history_population
  real, dimension(POP_SIZE, TOTAL_TIME_STEP) :: fear_history_population, hunger_history_population
  real, dimension(POP_SIZE, TOTAL_TIME_STEP) :: fear_gene_history_population, hunger_gene_history_population
  
  integer :: i, j
 if (IS_DEBUG_DATA) then
     file_name_weight_history = WEIGHT_HISTORY_SUBFOLDER// "weight_history_gen_" // TOSTR(gen_n) // ".csv"
     file_name_fear_history = EMOTION_STATE_SUBFOLDER // "fear_history_gen_" // TOSTR(gen_n) // ".csv"
     file_name_hunger_history = EMOTION_STATE_SUBFOLDER // "hunger_history_gen_" // TOSTR(gen_n) // ".csv"
    !  file_name_fear_gene_history = EMOTION_GENE_SUBFOLDER // "fear_gene_history_gen_" // TOSTR(gen_n) // ".csv"
    !  file_name_hunger_gene_history = EMOTION_GENE_SUBFOLDER // "hunger_gene_history_gen_" // TOSTR(gen_n) // ".csv"


    do i = 1, TOTAL_TIME_STEP
       var_name(i) = "STEP_" // TOSTR(i, TOTAL_TIME_STEP) ! add second parameter that defines the number format
     end do

     do i = 1, POP_SIZE
       do j = 1, TOTAL_TIME_STEP
         weight_history_population(i,j) = population_birds%birds(i)%weight_history(j)
         fear_history_population(i,j) = population_birds%birds(i)%fear_history(j)
         hunger_history_population(i,j) = population_birds%birds(i)%hunger_history(j)
         if (IS_SAVE_POPULATION_GENES .eqv. .TRUE.) then
           fear_gene_history_population(i,j) = population_birds%birds(i)%fear_gene_history(j)
           hunger_gene_history_population(i,j) = population_birds%birds(i)%hunger_gene_history(j)
         end if
       end do
     end do

     call CSV_MATRIX_WRITE(matrix = weight_history_population,                     &
                           csv_file_name = file_name_weight_history,               &
                           colnames = var_name)

     if (IS_ZIP_DATA) call execute_command_line(ZIP_PROG//' '//file_name_weight_history, WAIT= .FALSE.)
     ! compresses everything if IS_ZIP_DATA = .TRUE.

     call CSV_MATRIX_WRITE(matrix = fear_history_population,                        &
                           csv_file_name = file_name_fear_history,                  &
                           colnames = var_name)

     call CSV_MATRIX_WRITE(matrix = hunger_history_population,                        &
                           csv_file_name = file_name_hunger_history,                  &
                           colnames = var_name)

     call CSV_MATRIX_WRITE(matrix = hunger_gene_history_population,                        &
                      csv_file_name = file_name_hunger_gene_history,                  &
                      colnames = var_name)

     call CSV_MATRIX_WRITE(matrix = fear_gene_history_population,                        &
                      csv_file_name = file_name_fear_gene_history,                  &
                      colnames = var_name)



     if (IS_ZIP_DATA) call execute_command_line(ZIP_PROG//' '//file_name_fear_history, WAIT= .FALSE.)
     ! compresses everything if IS_ZIP_DATA = .TRUE.
     if (IS_ZIP_DATA) call execute_command_line(ZIP_PROG//' '//file_name_hunger_history, WAIT= .FALSE.)
     ! compresses everything if IS_ZIP_DATA = .TRUE.
     if (IS_ZIP_DATA) call execute_command_line(ZIP_PROG//' '//file_name_hunger_gene_history, WAIT= .FALSE.)
     ! compresses everything if IS_ZIP_DATA = .TRUE.
     if (IS_ZIP_DATA) call execute_command_line(ZIP_PROG//' '//file_name_fear_gene_history, WAIT= .FALSE.)
     ! compresses everything if IS_ZIP_DATA = .TRUE.

  end if

end subroutine save_history_array

! collumns
! 1: generation nr (int)
! 2: Num alive (int)
! 3: average mass all (real)
! 4: average mass alive (real)
! 5: best mass(real)
! 6: average mass top 25% (real)
subroutine build_output_data(row)
 ! This subroutine compiles and organizes various statistics and metrics from a simulation
 ! of a population of birds for a given generation.

  character(len = 15) :: value_str  ! Temporary string for formatting real values
  real :: value_real                ! Temporary real variable for storing formatted values

! these values represent the values for each single generation.
  integer, intent(in) :: row        ! 1  
  integer :: n_alive                ! 2  
  real :: average_mass_all          ! 3  
  real :: average_mass_alive        ! 4  
  real :: best_mass                 ! 5  
  real :: average_mass_dominants    ! 6  
  real :: average_gene_fear_alive   ! 7     
  real :: average_gene_hunger_alive ! 8
  real :: average_hunger_alive      ! 9
  real :: average_fear_alive        ! 10
  real :: std_dev_mass_alive        ! 11 
!  real :: freq_predation            ! 12 
  real :: meet_predator             ! 12
  integer :: n_is_killed_by_pred    ! 13 
  real :: n_is_killed_average       ! 14 
  real :: MUTATION_PROBABILITY      ! 15 
  real :: average_gene_fear_dead    ! 16
  real :: average_gene_hunger_dead  ! 17    
  real :: average_fear_dead         ! 18
  real :: average_hunger_dead       ! 19
  real :: std_dev_gene_fear_alive   ! 20
  real :: std_dev_gene_hunger_alive ! 21
  real :: average_gene_fear_all     ! 22
  real :: median_mass_alive         ! 23
  real :: avg_repr_prob             ! 24
  !real :: Full_Pop_Factor           ! 25





  ! temporary array keeping mass of birds that are alive, mass of dead birds
  ! is equal to MISSING
  real, dimension(POP_SIZE) :: mass_alive

  ! temporary array to keep population genes
  integer, dimension(POP_SIZE) :: gene_fear_alive, gene_hunger_alive, gene_fear_all

  !temporary array to keep the fear state (state_fear) of alive birds.
  real, dimension(POP_SIZE) :: fear_alive

  !temporary array to keep the hunger state (state_hunger) of alive birds.
  real, dimension(POP_SIZE) :: hunger_alive


  !for figuring out average gene of dead birds.
  real, dimension(POP_SIZE) :: gene_fear_dead, gene_hunger_dead

  !for figuring out average emotion of dead birds.
  real, dimension(POP_SIZE) :: fear_state_dead, hunger_state_dead

  !For figuring out how many birds have met predator

  integer, dimension(POP_SIZE) :: meet_predator_count

 ! Initialize temporary arrays with default values

 real, dimension(POP_SIZE) :: probability_repro

  mass_alive = MISSING

  gene_fear_dead = MISSING
  gene_hunger_dead = MISSING
  

  gene_fear_alive = -9999
  gene_hunger_alive = -9999
  

  fear_alive = -9999
  fear_state_dead = -9999

  hunger_alive = -9999

  meet_predator_count = UNDEFINED




! Calculate various stats and metrics. 
!NUMBER OF LIVING BIRDS AND THEIR AVERAGE MASS
  n_alive = count(parent_generation%birds%is_alive)
  average_mass_all = average(parent_generation%birds%weight)


  where(parent_generation%birds%is_alive)
    mass_alive = parent_generation%birds%weight
  end where


! AVERAGE MASS ALIVE
  average_mass_alive = average(mass_alive)
! BEST MASS
  best_mass = parent_generation%birds(1)%weight
!AVERAGE MASS DOMINANTS
  average_mass_dominants = average(parent_generation%birds(1:GA_REPRODUCE_N)%weight)

!MEDIAN OF MASS ALIVE BIRDS: 
! Calculate median of mass of alive birds

  median_mass_alive = median(mass_alive)

!AVG REPRODUCTION PROBABILITY (FITNESS)
  probability_repro = prob_repr(average_mass_alive,         &
  parent_generation%birds(n_alive)%weight, mass_alive(1))

  avg_repr_prob = average(probability_repro)



!AVERAGE GENE FEAR/HUNGER OF ALIVE BIRDS
  where(parent_generation%birds%is_alive)
    gene_fear_alive = parent_generation%birds%gene_fear
    gene_hunger_alive = parent_generation%birds%gene_hunger
  end where

  average_gene_fear_alive = average( real(gene_fear_alive) )
  average_gene_hunger_alive = average( real(gene_hunger_alive) )




!AVERAGE GENE FEAR OF DEAD BIRDS

  where(.not.parent_generation%birds%is_alive)
    gene_fear_dead = parent_generation%birds%gene_fear
    gene_hunger_dead = parent_generation%birds%gene_hunger
  end where

  average_gene_fear_dead = average( real(gene_fear_dead) )
  average_gene_hunger_dead = average( real(gene_hunger_dead) )

! AVERAGE GENE FEAR ALL
  gene_fear_all = real(parent_generation%birds%gene_fear)
  average_gene_fear_all = average( real(gene_fear_all) )


!AVERAGE EMOTION OF DEAD BIRDS

  where(.not.parent_generation%birds%is_alive)
    fear_state_dead = parent_generation%birds%state_fear
    hunger_state_dead = parent_generation%birds%state_hunger
    
  end where
  average_fear_dead = average( real( fear_state_dead ) )
  average_hunger_dead = average( real( hunger_state_dead ) )

!AVERAGE STATE FEAR /HUNGER OF ALIVE BIRDS
  where(parent_generation%birds%is_alive)
      fear_alive = parent_generation%birds%state_fear
      hunger_alive = parent_generation%birds%state_hunger
  end where

  average_fear_alive = average( real(fear_alive) )
  average_hunger_alive = average( real(hunger_alive) )


! STANDARD DEVIATION OF MASS
  std_dev_mass_alive = std_dev(mass_alive)


! STANDARD DEVIATION OF GENE

  std_dev_gene_fear_alive = std_dev(real(gene_fear_alive))
  std_dev_gene_hunger_alive = std_dev(real(gene_hunger_alive))

! MEET THE PREDATOR AVERAGE
! Calculate the average number of alive birds' number of predator encounters
   where(parent_generation%birds%is_alive)
     meet_predator_count = ( parent_generation%birds%bird_meets_predator_counter )

   elsewhere
     meet_predator_count = UNDEFINED !for all birds that are dead: assign UNDEFINED value

   end where

   meet_predator = average( real(meet_predator_count) ) 


! Update and display global variable
   Full_Pop_Factor = update_pop_regain_size_factor(n_alive, POP_SIZE)

  ! Update Last_Full_Pop_Factor if this is the last generation
  if (.not. Is_Evolutionary_Generations .and.                           &
  row == SIGNAL_ONLY_GENERATIONS - (SIGNAL_ONLY_GENERATIONS + 1)) then
    Last_Full_Pop_Factor = update_pop_regain_size_factor(n_alive,       &
       size(parent_generation%birds))
    Full_Pop_Factor = Last_Full_Pop_Factor
  end if

    
!PREDATOR KILLS BIRD

 n_is_killed_by_pred = count(parent_generation%birds%is_killed_by_predator) 
 !count works with arrays, and therefore we do not call birds(1)
 n_is_killed_average = (n_is_killed_by_pred / meet_predator) 

!
write(value_str, '(f6.3)') MUTATION_PROBABILITY
write(value_str, '(f6.3)') n_is_killed_average


 ! Add calculated values to the output array

  csv_generation_output(row, 1)  = row                       
  csv_generation_output(row, 2)  = n_alive                   
  csv_generation_output(row, 3)  = average_mass_all          
  csv_generation_output(row, 4)  = average_mass_alive        
  csv_generation_output(row, 5)  = best_mass                 
  csv_generation_output(row, 6)  = average_mass_dominants    
  csv_generation_output(row, 7)  = median_mass_alive
  csv_generation_output(row, 8)  = average_gene_fear_alive   
  csv_generation_output(row, 9)  = average_gene_hunger_alive 
  csv_generation_output(row, 10) = average_fear_alive        
  csv_generation_output(row, 11) = average_hunger_alive     
  csv_generation_output(row, 12) = std_dev_mass_alive               
  csv_generation_output(row, 13) = meet_predator            
  csv_generation_output(row, 14) = n_is_killed_by_pred      
  csv_generation_output(row, 15) = n_is_killed_average      
  csv_generation_output(row, 16) = MUTATION_PROBABILITY     
  csv_generation_output(row, 17) = average_gene_fear_dead   
  csv_generation_output(row, 18) = average_gene_hunger_dead 
  csv_generation_output(row, 19) = average_fear_dead        
  csv_generation_output(row, 20) = average_hunger_dead      
  csv_generation_output(row, 21) = std_dev_gene_fear_alive  
  csv_generation_output(row, 22) = std_dev_gene_hunger_alive
  csv_generation_output(row, 23) = average_gene_fear_all
  csv_generation_output(row, 24) = avg_repr_prob
  !csv_generation_output(row, 25) = Full_Pop_Factor 


end subroutine build_output_data

subroutine select_and_reproduce()
  use iso_fortran_env, only: real64
  implicit none

  integer :: bird_index, offspring_index, selected_bird_index, alive_counter, reproducing_counter
  integer :: remaining_slots, i
  type(BIRD), allocatable :: reproducing_birds(:)

  ! Initialize new offspring
  offspring_index = 1

  ! Calculate the number of alive birds
  alive_counter = parent_generation%alive_count()

  ! Allocate array to store reproducing birds
  allocate(reproducing_birds(POP_SIZE))
  reproducing_counter = 0

  ! Selection process
  do bird_index = 1, parent_generation%alive_count()
    associate(m => parent_generation%birds(bird_index)%weight, &
              m_0 => parent_generation%smallest_weight(), &
              m_max => parent_generation%biggest_weight())

      if (Is_Evolutionary_Generations) then
        if (GA_REPRODUCE_THRESHOLD < prob_repr(m, m_0, m_max)) then
          reproducing_counter = reproducing_counter + 1
          reproducing_birds(reproducing_counter) = parent_generation%birds(bird_index)

          ! Inline the selection process
          do selected_bird_index = 1, SELECTION_MULTIPLICATOR
            offspring_generation%birds(offspring_index) = parent_generation%birds(bird_index)
            offspring_index = offspring_index + 1
            if (offspring_index > POP_SIZE) exit
          end do
        end if
      else
        if (prob_repr(m, m_0, m_max) > GA_PROB_REPR_MIN) then
          reproducing_counter = reproducing_counter + 1
          reproducing_birds(reproducing_counter) = parent_generation%birds(bird_index)

          ! Inline the selection process
          do selected_bird_index = 1, SELECTION_MULTIPLICATOR
            offspring_generation%birds(offspring_index) = parent_generation%birds(bird_index)
            offspring_index = offspring_index + 1
            if (offspring_index > POP_SIZE) exit
          end do
        end if
      end if
    end associate
  end do

 ! Fill the remaining slots in the offspring array with randomly chosen reproducing birds
  remaining_slots = POP_SIZE - (offspring_index - 1)
  do bird_index = 1, remaining_slots
    selected_bird_index = RAND_I(1, reproducing_counter)
    offspring_generation%birds(offspring_index) = reproducing_birds(selected_bird_index)
    offspring_index = offspring_index + 1
  end do

  ! Deallocate the reproducing_birds array
  deallocate(reproducing_birds)

  ! Perform gene recombination for birds from index 101 to POP_SIZE
  call recombine_genes()
  ! Add this after the offspring generation is created
  do i = 1, POP_SIZE
    offspring_generation%birds(i)%weight = BIRD_INITIAL_WEIGHT
  end do

  ! Save gene history for the offspring generation
  call save_gene_history(offspring_generation, Current_Generation + 1)

!> Performs gene recombination for the offspring generation by exchanging the gene_fear and gene_hunger values between pairs of birds.
!> The recombination is performed for birds with indices from 101 to the end of the population size (POP_SIZE).
!> The recombination is done by associating each bird with the bird that is 100 indices behind it, and then exchanging the gene_fear and gene_hunger values between the pair.
contains

subroutine recombine_genes()
  integer :: current_bird_index, partner_bird_index_fear, partner_bird_index_hunger
  integer :: gene_swap_distance_fear, gene_swap_distance_hunger, total_alive_birds
  integer, parameter :: MIN_SWAP_DISTANCE = 50, MAX_SWAP_DISTANCE = 150
  type(BIRD) :: temp_bird_storage

  total_alive_birds = count(offspring_generation%birds%is_alive)

  if (total_alive_birds < MIN_SWAP_DISTANCE + 1) return

  do current_bird_index = 1, total_alive_birds - MIN_SWAP_DISTANCE
    ! Generate different swap distances for fear and hunger genes
    gene_swap_distance_fear = INT(RAND() * (MAX_SWAP_DISTANCE - MIN_SWAP_DISTANCE + 1)) + MIN_SWAP_DISTANCE
    gene_swap_distance_hunger = INT(RAND() * (MAX_SWAP_DISTANCE - MIN_SWAP_DISTANCE + 1)) + MIN_SWAP_DISTANCE

    ! Ensure different partner birds for fear and hunger genes
    partner_bird_index_fear = min(current_bird_index + gene_swap_distance_fear, total_alive_birds)
    partner_bird_index_hunger = min(current_bird_index + gene_swap_distance_hunger, total_alive_birds)

    ! Recombine fear genes
    associate(current_bird => offspring_generation%birds(current_bird_index), &
              partner_bird_fear => offspring_generation%birds(partner_bird_index_fear))
      temp_bird_storage%gene_fear = current_bird%gene_fear
      current_bird%gene_fear = partner_bird_fear%gene_fear
      partner_bird_fear%gene_fear = temp_bird_storage%gene_fear
    end associate

    ! Recombine hunger genes
    associate(current_bird => offspring_generation%birds(current_bird_index), &
              partner_bird_hunger => offspring_generation%birds(partner_bird_index_hunger))
      temp_bird_storage%gene_hunger = current_bird%gene_hunger
      current_bird%gene_hunger = partner_bird_hunger%gene_hunger
      partner_bird_hunger%gene_hunger = temp_bird_storage%gene_hunger
    end associate
  end do
end subroutine recombine_genes




end subroutine select_and_reproduce



subroutine save_new_popsize_fact_info(new_popsize_fact)
  use params
  use BASE_UTILS, only : TOSTR
  implicit none

  integer, intent(in) :: new_popsize_fact
  character(len=500) :: new_popsize_fact_file

  call FS_MKDIR(FACTOR_EXPERIMENT_SUBFOLDER)
  new_popsize_fact_file = FACTOR_EXPERIMENT_SUBFOLDER // "new_popsize_fact_gen_" //       &
    TOSTR(Current_Generation) // ".csv"
  
  call save_subset(offspring_generation, new_popsize_fact_file, 1, new_popsize_fact)

end subroutine save_new_popsize_fact_info



subroutine save_subset(self, csv_file_name, start_index, end_index)
  use CSV_IO
  implicit none

  class(POPULATION), intent(in) :: self
  character(len=*), intent(in) :: csv_file_name
  integer, intent(in) :: start_index, end_index

  integer :: i
  character(len=500) :: record

  open(unit=10, file=csv_file_name, status='replace')

  ! Write header
  write(10, '(a)') "Bird,Weight,Gene_Fear,Gene_Hunger,State_Fear,State_Hunger,Is_Alive"

  ! Write data for each bird in the subset
  do i = start_index, end_index
    write(record, '(i0,",",f10.6,",",i0,",",i0,",",f10.6,",",f10.6,",",l1)')     &
      i, self%birds(i)%weight, self%birds(i)%gene_fear, self%birds(i)%gene_hunger, &
      self%birds(i)%state_fear, self%birds(i)%state_hunger, self%birds(i)%is_alive
    write(10, '(a)') trim(record)
  end do

  close(10)

end subroutine save_subset

subroutine minimum_weight_exp_gens()

  integer :: i
  real :: m, m_0, m_max

  associate(  m =>  parent_generation%birds(i)%weight,          &
              m_0 => parent_generation%smallest_weight(),       &
              m_max => parent_generation%biggest_weight())


  ! need minimum threshold weight that is dynamically set by smallest surviving mother
  ! that produces eggs from generation "Genome_File_Name"
   if (.NOT.((parent_generation%birds(i)%is_alive) .and.   &
      (prob_repr(m, m_0, m_max) > GA_PROB_REPR_MIN ))) then  !make a parameter 

 
   end if

  end associate
end subroutine minimum_weight_exp_gens



subroutine generation_swap()

  if(associated(parent_generation, target=generation_one)) then
    parent_generation => generation_two
    offspring_generation => generation_one
  else
    parent_generation => generation_one
    offspring_generation => generation_two
  end if

end subroutine generation_swap



function reset_generation() result(popout)

 type(POPULATION) :: popout

 !This subroutine resets the parent generation to its initial state
 call popout%init_pop(is_gene_init = .TRUE.)
 call popout%load(Genome_File_Name) 

end function reset_generation

end module