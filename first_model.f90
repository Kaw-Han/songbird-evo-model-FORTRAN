
!
!------------------------------------
!--------PARAMETER MODULE------------
!------------------------------------
module params
    implicit none

    logical, public, parameter :: IS_DEBUG = .FALSE.

    logical, public, parameter :: IS_SAVE_POPULATION = .TRUE.

    ! Numerical code for missing value
    real, parameter, public :: MISSING = -9999.0

    integer, parameter, public :: HUNGER_MIN = -10, HUNGER_MAX = 10 !HUNGER_MIN = maximum fear, HUNGER_MAX = maximum boldness
    integer, parameter, public :: ENVIRONMENT_SIZE = 100 !size of envornment, nr of cells = 100

    real, parameter, public :: FOOD_AVAILABILITY_MEAN = 0.4, FOOD_AVAILABILITY_VARIANCE = 0.39 !parameter for food, measured in weight grams added to the birds mass
    real, parameter, public :: FOOD_AVAILABILITY_MIN = 0.01, &
        FOOD_AVAILABILITY_MAX = FOOD_AVAILABILITY_MEAN + FOOD_AVAILABILITY_VARIANCE *1.5

    real, parameter, public :: FREQUENCY_OF_PREDATOR_MIN = 0.0, FREQUENCY_OF_PREDATOR_MAX = 0.005 !parameter for predation_risk, risk of the bird to meet predator in a particular environment.
! if the bird sees the predator, the predator sees the bird.
    real, parameter, public :: BIRD_INITIAL_WEIGHT = 20.0 !20 grams
    real, parameter, public :: WEIGHT_DEATH_THRESHOLD = 15.0! if birds weight falls to or below 75% of BIRD_INITIAL_WEIGHT, the bird is dead.
    real, parameter, public :: BIRD_WEIGHT_COST_OF_ATTACK = 0.1 !The birds weight loss after survival of an attack.
    real, parameter, public :: WEIGHT_REDUCTION_CONSTANT = 0.999
    real, parameter, public :: ATTACK_RATE = 0.01 !risk of bird being eaten by predator if attacked.

    real, parameter, public :: HUNGER_DECREMENT_MULTIPLICATOR = 0.3 !This relates to the increased hunger of the bird in feeding situations.
    real, parameter, public :: FEAR_INCREMENT_AFTER_ATTACK = 0.3 !This relates to the increased fear of the bird in attack situations.
    real, parameter, public :: HUNGER_STEP_FRACTION = 0.1! This speaks to function for emotion_step

    real, parameter, public :: HUNGER_PLASTICITY = 0.5 ! Param for variation in hunger between birds.


    integer, parameter, public :: TOTAL_TIME_STEP = 200

    real, parameter, public :: MUTATION_PROBABILITY = 0.0019 !https://en.wikipedia.org/wiki/Genomic_evolution_of_birds


    ! File name that keeps all output data from all generations of the model
    character(len=*), parameter, public :: MODEL_OUTPUT_FILE = "my_model_output.csv"


    !Population size
    integer, parameter, public :: GENERATIONS = 100
    integer, parameter, public :: POP_SIZE = 100000

    !Proportion of the best reproducing birds of selection
    real, parameter :: GA_REPRODUCE_PR = 0.25
    !exact number of birds reproducing
    integer, parameter, public :: GA_REPRODUCE_N = int(POP_SIZE * GA_REPRODUCE_PR)


    !Variables for fraction in fuction bird_eat_fraction_when_fear
    real, parameter, public :: FRAC_S1 = 0.01
    real, parameter, public :: FRAC_S2 = 1.0



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
  real :: risk_of_predation
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
subroutine cell_init(this, x) !this = foraging area, x = x-coordinate
    use BASE_RANDOM
    class(env_cell), intent(out) :: this
    integer, intent(in) :: x

    this%x = x ! set location
    this%food_availability =                                              &
        within(RNORM(FOOD_AVAILABILITY_MEAN, FOOD_AVAILABILITY_VARIANCE), &
                     FOOD_AVAILABILITY_MIN, FOOD_AVAILABILITY_MAX ) !food variance, normal
    this%risk_of_predation = RAND(FREQUENCY_OF_PREDATOR_MIN, FREQUENCY_OF_PREDATOR_MAX) !risk predation, uniform

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
    call this%point(i)%init(i) ! we initialize the food and predation data
  end do

end subroutine env_init



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
type, extends(location), public :: GENE
    integer :: gene_fear_hunger !the higher the gene_fear_hunger the lower the gene, opposite proportional
    contains
    procedure, public :: init_gene => gene_init
    procedure, public :: mutate => gene_mutate
end type GENE


!--------------------------------
!-------------BIRD/GENE----------
!--------------------------------
type, extends(GENE), public :: BIRD
    real :: weight
    real :: state_fear_hunger
    logical :: is_alive
    integer :: bird_meets_predator_counter
    !integer :: death_count !this relates to subroutine count_the_dead_birds.
    contains
    procedure, public :: init => bird_init
    procedure, public :: killed_from_starvation => bird_dead_from_starvation !now we can use this in the future in the folloing manner: this%killed_from_starvation
    procedure, public :: is_starved => bird_is_starved !for logical boolean functions "is" is good to use to refer to the true/false_statement
    procedure, public :: fly => bird_do_fly
    procedure, public :: do_feed => bird_feeds
    procedure, public :: hunger_update => bird_update_emotion
    procedure, public :: genetic_lower_hunger => bird_hunger_min_genetic_limit
    procedure, public :: genetic_upper_hunger => bird_hunger_max_genetic_limit
    procedure, public :: hunger_genetic_limit => bird_force_hunger_within_genetic_limits


end type BIRD

!--------------------------------
!-----------PREDATOR-------------
!--------------------------------
type, public :: PREDATOR
    real :: risk !If PREDATOR sees the bird according to the risk_of_the, the risk of attack is this

  contains

    procedure, public :: init => predator_init_new
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
  procedure, public :: time_steps => population_time_steps
  procedure, public :: sort => sort_by_fitness
  procedure, public :: save => population_save_csv


end type POPULATION



!--------------------------------
!---------SUBROUTINES-----------
!--------------------------------
contains
!subroutines for POPULATION
subroutine population_init_new(this, alt_size)
  class(POPULATION), intent(out) :: this! intent(out) because the initialization will be only an output.
  integer, optional, intent(in) :: alt_size! if we want some other value than POP_SIZE

  integer :: size
  integer :: i

  !process for optional parameter where we can define alternative pop_size:
  if (present(alt_size)) then
    size = alt_size
  else
    size = POP_SIZE !POP_SIZE defined in parameters mod.
  end if

  allocate (this%birds(size)) !allocate the array to a specific size

  do i = 1, size
    call this%birds(i)%init() !initialize array of birds

  end do


end subroutine population_init_new


subroutine population_save_csv(this, file_name)
  use CSV_IO
  class(POPULATION), intent(in) :: this
  character(len=*), intent(in) :: file_name

  integer, parameter :: GENPOP_COL_NUMBER = 6
  character(len=*), dimension(GENPOP_COL_NUMBER), parameter :: GENPOP_COLUMNS = &
    ["BIRD                 ",                      &  ! 1
     "GENE                 ",                      &  ! 2
     "IS_ALIVE             ",                      &  ! 3
     "WEIGHT               ",                      &  ! 4
     "BIRD_PRED_COUNT      ",                      &  ! 5
     "STATE_FEAR_HUNGER    " ]                        ! 6

  real, dimension(POP_SIZE, GENPOP_COL_NUMBER) :: out_data

  integer :: i

  do i = 1, POP_SIZE
    out_data(i,1) = i
    out_data(i,2) = real(this%birds(i)%gene_fear_hunger)
    out_data(i,3) = l2r(this%birds(i)%is_alive)
    out_data(i,4) = this%birds(i)%weight
    out_data(i,5) = real(this%birds(i)%bird_meets_predator_counter)
    out_data(i,6) = this%birds(i)%state_fear_hunger
  end do

  call CSV_MATRIX_WRITE(out_data, file_name, colnames=GENPOP_COLUMNS)

end subroutine population_save_csv



! Functions for counting dead birds, updating it. Iterating over the array.
!
function population_get_dead_count(this) result(get_dead_count)! this = population, num_dead is the number of dead birds.
  class(POPULATION), intent(in) :: this
  integer :: get_dead_count

  get_dead_count = count(.not.this%birds%is_alive) !birdS = population of birds.

end function population_get_dead_count




! Subroutines for PREDATOR
subroutine predator_init_new(this)
  class(PREDATOR), intent(inout) :: this

  this%risk = ATTACK_RATE !chance of prey being successfully attacked and eaten by pred.

end subroutine predator_init_new

subroutine predator_attack_bird(this, bird_prey, environment_in,              &
                                predator_is_present, prey_is_killed)!defines how predator attacked bird
  class(PREDATOR), intent(in) :: this
  class(BIRD), intent (inout) :: bird_prey !the bird enters, and predator attacks the bird. The bird changes as a result of this. How the bird changes depends as the result of this
  class(whole_environ), intent(in) :: environment_in

  ! The two optional parameters are useful for testing and sebuygging:
  ! optional output indicator if the predator was present in the bird's cell
  logical, optional, intent(out) :: predator_is_present
  ! optional indicator that predator actually did attack the prey
  logical, optional, intent(out) :: prey_is_killed

  ! local copies of optional parameters
  logical :: p_is_present, p_prey_dies

  p_is_present = .FALSE.
  p_prey_dies = .FALSE.
  !***********************************************************
  !if statement that engages predator if bird is present here.
  !***********************************************************

!Function for counting the number of times the bird has met the Predator and survived. If Bird dies, result should be "is_alive = .FALSE.".
!If Bird survives meeting with predator, and in the next cell does not meet predator, the counter should return to 0.
!The counter will therefore be equivelent to the exponent in the number that defines birds cost in weight for meeting a predator and fleeing,
!and therefore not getting to eat itself. The weight reduction can be 0.9^x of the birds weight before interaction with predator.
!Ex: if a bird meets predator three times in a row, the weight reduction of the bird will be "birds_weight * 0.9**3", where 3 = x.
!if a birds weight falls to or below 75% of it's initital weight, the bird is dead. 75% of inital weight corresponds to WEIGHT_DEATH_THRESHOLD (params mod).

  if (bird_prey%is_alive) then !execute if-statement ONLY if bird is alive.
    !we check if the predator is present in the cell, given the risk of predation.
    p_is_present = (RAND() < environment_in%point( bird_prey%location%x )%risk_of_predation)
    if ( p_is_present ) then
      !We check if the predator will attack given it is in the same cell as the prey.
      p_prey_dies = ( RAND() < this%risk )
      if (p_prey_dies) then
        bird_prey%is_alive = .FALSE. !if random is larger then this%risk, the bird is killed by predator
        bird_prey%bird_meets_predator_counter = 0  !if bird is killed, count for exponent in weight equation is 0.
      else
        bird_prey%is_alive = .TRUE. !if bird survives attack and flees then...
        bird_prey%bird_meets_predator_counter = bird_prey%bird_meets_predator_counter + 1!counts how many times in a row bird meets predator and survives to escape
        bird_prey%weight = bird_prey%weight * (WEIGHT_REDUCTION_CONSTANT**(bird_prey%bird_meets_predator_counter))!the exponent = how many times in a row a bird survives, but does not eat.
        bird_prey%state_fear_hunger = within( bird_prey%state_fear_hunger -                               &
                                  bird_prey%state_fear_hunger * FEAR_INCREMENT_AFTER_ATTACK,              &
                                  bird_prey%genetic_lower_hunger(), bird_prey%genetic_upper_hunger() )

      end if
    else
     ! bird_prey%bird_meets_predator_counter = 0
    end if
  end if

  if (present(predator_is_present)) predator_is_present = p_is_present
  if (present(prey_is_killed)) prey_is_killed = p_prey_dies

end subroutine predator_attack_bird


!****************************
!subroutine for initial gene
!****************************
subroutine gene_init(this)
    use BASE_RANDOM !importing this.
    class(GENE), intent(out) :: this !in = input, out = output, inout = both

    !classname%attribute_name
    this%gene_fear_hunger = RAND(HUNGER_MIN, HUNGER_MAX)


end subroutine gene_init


subroutine gene_mutate(this, is_mutate)
  class(GENE), intent(inout) :: this
  logical, optional, intent(out) :: is_mutate

  logical :: is_mutate_loc

  is_mutate_loc = .FALSE.

  if (RAND() < MUTATION_PROBABILITY) then
    call this%init_gene()
    is_mutate_loc = .TRUE.
  end if

  if (present(is_mutate)) is_mutate = is_mutate_loc


end subroutine gene_mutate


!Subroutine to initalize the bird

subroutine bird_init(this)

  class(BIRD), intent(out) :: this !by calling class(BIRD), we inherit all qualities of BIRD and supply it with new atributes through new subrout.

  call this%init_gene()
  call this%put_random()
  this%weight = BIRD_INITIAL_WEIGHT
  this%state_fear_hunger = real(this%gene_fear_hunger)
  this%is_alive = .TRUE.
  this%bird_meets_predator_counter = 0

end subroutine bird_init




!determine if bird is dead or not. We use a function for doing this, so that we can declare an output with result-syntax.
function bird_is_starved(this, threshold) result (is_dead) !if the bird is starved and weight drops below threshold, then the bird dies of starvation.
  class(BIRD), intent(in) :: this
  real, optional, intent(in) :: threshold
  logical :: is_dead

  real :: local_threshold

  if (present(threshold)) then
    local_threshold = threshold
  else
    local_threshold = WEIGHT_DEATH_THRESHOLD
  end if

  if (this%weight < local_threshold) then
    is_dead = .TRUE.
  else
    is_dead = .FALSE.
  end if

end function bird_is_starved



!Subroutine for weight gain according to cells nutritional value.
subroutine bird_feeds(this, in_environment)
  class(BIRD), intent(inout) :: this
  class(whole_environ), intent(in) :: in_environment
  real :: current_cell_food

if (this%is_alive) then

    current_cell_food = in_environment%point(this%location%x)%food_availability !this exact location where the bird is, has a certain food_availability specific to the cell.


  if (this%state_fear_hunger < 0.0) then
    this%state_fear_hunger = this%state_fear_hunger + HUNGER_STEP_FRACTION*this%state_fear_hunger
    this%weight = this%weight + current_cell_food * &
                        bird_eat_fraction_when_fear(this%state_fear_hunger, &
                        this%genetic_lower_hunger(), this%genetic_upper_hunger() )
  else
  ! Increment the bird weight
    this%weight = this%weight + current_cell_food * &                           !because array, the parenthesis is used
                        bird_eat_fraction_when_fear(this%state_fear_hunger, &
                        this%genetic_lower_hunger(), this%genetic_upper_hunger() )
    ! Hunger is decremented == Increment fear
    this%state_fear_hunger = within(this%state_fear_hunger -                                            &
                              this%state_fear_hunger * HUNGER_DECREMENT_MULTIPLICATOR,                  &
                                                      this%genetic_lower_hunger(),                      &
                                                      this%genetic_upper_hunger() )
  end if
  call this%hunger_genetic_limit()
end if

end subroutine bird_feeds

! Calculate the proportion of food availability eaten at specific state of
! fear (state_hunger_hunger<0)
! The proportion of food eaten is 0.0 at state_fear_hunger = HUNGER_MIN
!                             and 1.0 at state_fear_hunger = HUNHER_MAX
!
!  S = k * emotion + b
!
!   | S1 = k * H_MIN + b
!   | S2 = k * M_MAX + b
!
!   solution:   k = (S2 - S1) / (H_MIN - H_MAX)
!               b = (S2 * H_MIN - S1 * H_MAX) / ( H_MIN-H_MAX )
!
function bird_eat_fraction_when_fear(emotion, min, max) result(fract)
  real, intent(in) :: emotion
  real, intent(in), optional :: min, max

  real :: fract
  real :: k, b
  real, parameter :: S1 = FRAC_S1, S2 = FRAC_S2

  real :: min_loc, max_loc

  if (present(min)) then
    min_loc = min
  else
    min_loc = HUNGER_MIN
  end if

  if (present(max)) then
    max_loc = max
  else
    max_loc = HUNGER_MAX
  end if

  k = (S2 - S1) / (max_loc - min_loc)
  b = (S2 * min_loc - S1 * max_loc) / ( min_loc-max_loc )

  fract = within(k * emotion + b, 0.0, 1.0)

end function bird_eat_fraction_when_fear


subroutine bird_dead_from_starvation(this)
  class(BIRD), intent(inout) :: this

  if (this%is_starved()) this%is_alive = .FALSE.

end subroutine bird_dead_from_starvation



subroutine bird_do_fly(this, in_environment) !dependent on the condition based on hunger-state of the bird.
  class(BIRD), intent(inout) :: this
  class(whole_environ), intent(in) :: in_environment

    ! Optional parameters defining a range of positions to place the spatial
    ! object, but by default the position is from 1 to ENVIRONMENT_SIZE

  integer :: min_loc, max_loc

  if (this%is_alive) then



    min_loc = 1 !min location = 1
    max_loc = size(in_environment%point) !max location is the size of environment. This way we can make the value dynamic


    if (this%state_fear_hunger > -5.0 .and.               &
          this%state_fear_hunger < 0.0) then !if bird is in the state of fear then it will fly, if not it will eat
      call this%walk(min_loc, max_loc)
    end if

  end if

end subroutine bird_do_fly

subroutine bird_update_emotion(this)
  class(BIRD), intent(inout) :: this


  this%state_fear_hunger = this%state_fear_hunger + emotion_step(this%weight, &
                  min=this%genetic_lower_hunger(), max=this%genetic_upper_hunger() )

  call this%hunger_genetic_limit()



end subroutine bird_update_emotion

function bird_hunger_min_genetic_limit(this) result(lower_limit)
  class(BIRD), intent(in) :: this
  real :: lower_limit

  real :: range

  range = HUNGER_MAX - HUNGER_MIN

  lower_limit = this%gene_fear_hunger - range * HUNGER_PLASTICITY

  lower_limit = within(lower_limit, real(HUNGER_MIN), real(HUNGER_MAX))


end function bird_hunger_min_genetic_limit


function bird_hunger_max_genetic_limit(this) result(upper_limit)
  class(BIRD), intent(in) :: this
  real :: upper_limit

  real :: range

  range = HUNGER_MAX - HUNGER_MIN

  upper_limit = this%gene_fear_hunger + range * HUNGER_PLASTICITY

  upper_limit = within(upper_limit, real(HUNGER_MIN), real(HUNGER_MAX))

end function bird_hunger_max_genetic_limit


subroutine bird_force_hunger_within_genetic_limits(this)
! setting a limited variance for genetic dispositions in birds.
! This way there will be heterogeneity in the genetic compositions.
 class(BIRD), intent(inout) :: this

  this%state_fear_hunger = within(this%state_fear_hunger,                         &
                                  this%genetic_lower_hunger(),                    &
                                  this%genetic_upper_hunger() )

end subroutine bird_force_hunger_within_genetic_limits

! This function calculates the emotion state increment dependent on the
! bird body mass.
!  1. If body mass is high (up to 25), then hunger state is reduced
!     up to zero, while fear increases up to a maximum value
!  2. If body mass is small, then hunger increases while fear decreases,
!
! Here
! - hunger decrement reaches minimum at mass 25.0 (fear is maximum at
!   this mass)
! - hunger reaches maximum at the starvation threshold mass 15, while fear
!   reaches minimum
! - when ther bird mass is equal to the starting mass 20g, no change
!   in change imotion state occurs
elemental function emotion_step(mass, frac, min, max) result (step)
  real, intent(in) :: mass
  real, optional, intent(in) :: frac !***** what is frac? *****
  !Fraction of the increment the birds emotional value based on its level of fear (y-axis).
  !A fraction of decrement or increment responds to the amount of food the bird eats.
  real, optional, intent(in) :: min, max
  real :: step

  real :: frac_loc, min_loc, max_loc
  real :: k, b
  real :: h_min, h_max
  real, parameter :: T_0 = WEIGHT_DEATH_THRESHOLD, W_M =  BIRD_INITIAL_WEIGHT * 1.25

  if (present(frac)) then
    frac_loc = frac
  else
    frac_loc = HUNGER_STEP_FRACTION
  end if

  if (present(min)) then
    min_loc = min
  else
    min_loc = HUNGER_MIN
  end if

  if (present(max)) then
    max_loc = max
  else
    max_loc = HUNGER_MAX
  end if

  h_min = min_loc+1.0
  h_max = max_loc

  ! We need to solve the simple linear system
  ! ----
  ! \left\{\begin{matrix}
  ! H_{max} = k T_0 + b \\
  ! H_{min} = k W_m + b
  ! \end{matrix}\right.
  ! ----
  ! Solution of this system is:
  ! k = \frac { H_{max} - H_{min} } { T_0 - W_m }
  ! b = \frac { H_{min} T_0 - H_{max} W_m } { T_0 - W_m }
  k = (h_max - h_min) / (T_0 - W_M)
  b = (h_min*T_0 - h_max * W_M) / (T_0 - W_M)

  step = k * mass + b
  step = step * frac_loc

end function emotion_step



subroutine population_time_steps(this, environment_in, predator_in) !bird is in this environment = environment_in
  class(POPULATION), intent(inout) :: this
  class(whole_environ), intent(in) :: environment_in
  class(PREDATOR), intent(inout) :: predator_in

  integer :: current_step
  integer :: bird_current !this current bird


  logical :: debug_check , debug_ckeck2



  do bird_current = 1, size(this%birds)

    do current_step=1, TOTAL_TIME_STEP

      if (this%birds(bird_current)%is_alive) then


        DEBUG_GET_INFO: if (IS_DEBUG) then
                        print *, "bird", bird_current, "at step", current_step
                        if (this%birds(bird_current)%is_starved()) then
                          print*, "STARVED"
                        end if
                        print *, "Location:", this%birds(bird_current)%x, ",   food",                   &
                                  environment_in%point( this%birds(bird_current)%x )%food_availability
                        print *, "fear versus hunger ", this%birds(bird_current)%state_fear_hunger, &
                                          "[", this%birds(bird_current)%genetic_lower_hunger(),     &
                                              this%birds(bird_current)%genetic_upper_hunger(), "]"
                        print *, " mass", this%birds(bird_current)%weight
                    end if DEBUG_GET_INFO

        call this%birds(bird_current)%killed_from_starvation() !is the bird dead or alive from starvation?
        !what birds do if is alive
        call this%birds(bird_current)%fly(environment_in) !bird is flying
        call predator_in%predator_attack_bird(this%birds(bird_current), environment_in, &
              predator_is_present=debug_ckeck2, prey_is_killed=debug_check)! predator attacks bird

        DEBUG_PRINT_KILL: if (IS_DEBUG) then
                                  print *, "Predator ", debug_ckeck2
                                  if (debug_check) print *, "KILLED by predator"
                                end if DEBUG_PRINT_KILL
          call this%birds(bird_current)%do_feed(environment_in)!birds is feeding, and weight is incrementing
          call this%birds(bird_current)%hunger_update() !updates the hunger of the bird, next is behaviour
        else
                         if (IS_DEBUG) print *, bird_current , "DEAD"
      end if
    end do
  end do

end subroutine population_time_steps

subroutine sort_by_fitness(this)!sorting by bodymass
  class(POPULATION), intent(inout) :: this

  !-----------------------------------------------------------------------------
  ! Sort parents by fitness
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


!Module for managing simulations
module GA !genetic algorythm
use params
use environment
use BASE_RANDOM
use CSV_IO
use organism
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
! 7: standard deviation of mass for alive (real))
integer, parameter :: CSV_COL_NUMBER = 8
character(len=*), dimension(CSV_COL_NUMBER), parameter :: CSV_COLUMNS = &
  ["GENERATION                 ",                      &  !1
   "NUM-ALIVE                  ",                      &  !2
   "AVERAGE-MASS-ALL           ",                      &  !3
   "AVERAGE-MASS-ALIVE         ",                      &  !4
   "BEST-MASS                  ",                      &  !5
   "AVERAGE-MASS-DOMINANTS     ",                      &  !6
   "AVERAGE-GENE-ALIVE         ",                      &  !7
   "STD-DEV-MASS-ALIVE         "]                         !8





real, dimension(GENERATIONS, CSV_COL_NUMBER) :: csv_generation_output

private :: CSV_COL_NUMBER, CSV_COLUMNS, csv_generation_output


contains

subroutine genetic_algorith()

  use BASE_UTILS, only : TOSTR

  integer :: current_generation

  real, dimension(POP_SIZE) :: fitness

  integer :: n_alive_get

  csv_generation_output = MISSING

  call RANDOM_SEED_INIT()

  call habitat%init()

  call predator_in_habitat%init()



  parent_generation => generation_one
  offspring_generation => generation_two

  call parent_generation%init_pop()
  call offspring_generation%init_pop()

  do current_generation = 1, GENERATIONS

    write(OUTPUT_UNIT, *) "Generation: ", current_generation


    call parent_generation%time_steps(habitat, predator_in_habitat)

    !sorting by fitness
    call parent_generation%sort()

    SAVE: block

    character(len=500) :: generation_data_file
    character(*), parameter :: GENERATION_SUBFOLDER = "generations/"
    call FS_MKDIR(GENERATION_SUBFOLDER)
    generation_data_file = GENERATION_SUBFOLDER // "gen_" // TOSTR(current_generation) // ".csv"
    if (IS_SAVE_POPULATION) call parent_generation%save(generation_data_file)

    end block SAVE


    !record parent characteristics to output file
    call build_output_data(current_generation)

    n_alive_get = csv_generation_output(current_generation, 2)

    write(OUTPUT_UNIT, *) "    Num alive:n ", n_alive_get, &
        "        avg.mass alive ", csv_generation_output(current_generation, 4)


    if (n_alive_get < 10) then
      write(ERROR_UNIT, *) "WARNING: Too few birds remain in the generation"

      call CSV_MATRIX_WRITE(csv_generation_output,                      &
                            MODEL_OUTPUT_FILE,colnames=CSV_COLUMNS)
      stop
    end if

    call selection_of_dominant() !selecting best fit - guaranteed mating.

    call select_and_reproduce()

    call generation_swap()


  end do

call CSV_MATRIX_WRITE(csv_generation_output,                      &
                      MODEL_OUTPUT_FILE,colnames=CSV_COLUMNS)

end subroutine genetic_algorith





! collumns
! 1: generation nr (int)
! 2: Num alive (int)
! 3: average mass all (real)
! 4: average mass alive (real)
! 5: best mass(real)
! 6: average mass top 25% (real)
subroutine build_output_data(row)
! these values represent the values for each single generation.
  integer, intent(in) :: row     ! 1
  integer :: n_alive             ! 2
  real :: average_mass_all       ! 3
  real :: average_mass_alive     ! 4
  real :: best_mass              ! 5
  real :: average_mass_dominants ! 6
  real :: average_gene_alive     ! 7
  real :: std_dev_mass_alive     ! 8


  ! temporary array keeping mass of birds that are alive, mass of dead birds
  ! is equal to MISSING
  real, dimension(POP_SIZE) :: mass_alive

  ! temporary array to keep population genes
  integer, dimension(POP_SIZE) :: gene_alive

  mass_alive = MISSING

  gene_alive = -9999

  n_alive = count(parent_generation%birds%is_alive)
  average_mass_all = average(parent_generation%birds%weight)


  where(parent_generation%birds%is_alive)
    mass_alive = parent_generation%birds%weight
  end where

  average_mass_alive = average(mass_alive)

  best_mass = parent_generation%birds(1)%weight

  average_mass_dominants = average(parent_generation%birds(1:GA_REPRODUCE_N)%weight)

  where(parent_generation%birds%is_alive)
    gene_alive = parent_generation%birds%gene_fear_hunger
  end where

  average_gene_alive = average( real(gene_alive) )

  std_dev_mass_alive = std_dev(mass_alive)




  csv_generation_output(row, 1) = row                   ! 1
  csv_generation_output(row, 2) = n_alive               ! 2
  csv_generation_output(row, 3) = average_mass_all      ! 3
  csv_generation_output(row, 4) = average_mass_alive    ! 4
  csv_generation_output(row, 5) = best_mass             ! 5
  csv_generation_output(row, 6) = average_mass_dominants! 6
  csv_generation_output(row, 7) = average_gene_alive    ! 7
  csv_generation_output(row, 8) = std_dev_mass_alive    ! 8





end subroutine build_output_data



subroutine selection_of_dominant() !selection of genes from parents to offspring.

  offspring_generation%birds(:GA_REPRODUCE_N) = parent_generation%birds(:GA_REPRODUCE_N)

end subroutine selection_of_dominant

subroutine select_and_reproduce() !25% of best fitness, 25% of remaining less fit, both base for new offspring

integer :: i, i1, i2, spos

  do i = GA_REPRODUCE_N + 1, POP_SIZE
    i1 = RAND_I(1, POP_SIZE / 2 )

    offspring_generation%birds(i) = parent_generation%birds(i1)

    !mutate randomly...

    call offspring_generation%birds(i)%mutate()


  end do
end subroutine select_and_reproduce

subroutine generation_swap()

  if(associated(parent_generation, target=generation_one)) then
    parent_generation => generation_two
    offspring_generation => generation_one
  else
    parent_generation => generation_one
    offspring_generation => generation_two
  end if

end subroutine generation_swap


end module GA
