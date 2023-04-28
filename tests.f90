!simulation is done in this code
! targeting for testing will also work, creating a test-model.
!
program first_model

use params
use organism

implicit none !always after use-statements

character(*), parameter :: DASHES=repeat("-", 72)

integer :: i, j

type(whole_environ) :: our_test_habitat

type(BIRD) :: our_bird

type(PREDATOR) :: evil_predator

type(POPULATION) :: bird_pop, bird_pop2

real :: weight_loss, hunger_change




! Logical flags to test predator attack code
logical :: p_is_here, prey_is_killed

real :: p_escape


! test hunger emotion-related feeding fraction
print *, DASHES
print *, "Testing emotion-related feeding fraction"
print *, "at -10", bird_eat_fraction_when_fear(-10.0)
print *, "at   0", bird_eat_fraction_when_fear(  0.0)
print *, "at  10", bird_eat_fraction_when_fear( 10.0)


! test emotion increment-decrement function
! Bird dies when mass < 15 grams, and is saturated when mass >= 25 grams. 
print *, DASHES
print *, "Test emotion increment-decrement function"
print *, "Minimum body mass threshold (15grams), raw", emotion_step(15.0, 1.0)
print *, "Maximum body mass threshold (25grams), raw", emotion_step(25.0, 1.0)
print *, ""
print *, "Test hunger change term for different body mass:"
print *, " - test mass= 5.0,  0.8", emotion_step( 5.0, 0.80)
print *, " - test mass=15.0*, 0.8", emotion_step(15.0, 0.80)
print *, " - test mass=20.0*, 0.8", emotion_step(20.0, 0.80)
print *, " - test mass=25.0*, 0.8", emotion_step(25.0, 0.80)
print *, " - test mass=45.0,  0.8", emotion_step(45.0, 0.80)

print*, DASHES
print*, ""



! test the within function
print *, DASHES

!Test the equation function of emotion increment-decrement
print*, "Testing emotion_step(mass, frac)"
print*, emotion_step(20.0, 5.0)
print*, emotion_step(22.0, 5.0)

print *, "Test the within function"
print *, "Must be 1.0", within(0.5, 1.0, 2.0)
print *, "Must be 2.0", within(2.5, 1.0, 2.0)
print *, "Must be 1.7", within(1.7, 1.0, 2.0)

!Test if the enviornment can be created
print *, DASHES
print *, "Test if the enviornment can be created"

call our_test_habitat%init(10)
print*, "Coordinates of habitats", our_test_habitat%point%x
print*, "Food availability", our_test_habitat%point%food_availability
print*, "Risk of predation", our_test_habitat%point%risk_of_predation
print*, "First cell location", our_test_habitat%point(10)%x

! ---------------------------------------------------------------------------
print*, DASHES
print*, "Test if our population is created"

call bird_pop%init_pop()
print*, "The number of birds are: ", size(bird_pop%birds) !size gives us the number of birds.
print*, "The number of dead birds are: ", bird_pop%dead_count()
print*, "Testing the killing of one bird:"
bird_pop%birds(3)%is_alive = .FALSE.
print*, bird_pop%dead_count()
!stop
! ! ---------------------------------------------------------------------------

! testing if bird is created correctly
print*, DASHES
print *, "TEST: testing if bird is created correctly"

call our_bird%init()
print*, "The bird attributes are: ", our_bird
!print*, "The bird is dead?", our_bird%is_alive() !True is false
print*, "Location: ", our_bird%x
print*, "Gene: ", our_bird%gene_fear_hunger
print*, "Weight: ", our_bird%weight
print*, "State fear vs hunger: ", our_bird%state_fear_hunger

if (our_bird%state_fear_hunger /= our_bird%gene_fear_hunger) then
    print*, "ERROR: state fear vs hunger is wrongly initialized"
print*, DASHES
end if

! ! ---------------------------------------------------------------------------
! ! Test how a bird can do random walks
! print *, DASHES
! print *, "Test how a bird can do random walks"
! call our_bird%place(5)

! do i = 1, 50

!     call our_bird%walk(1, 10)
!     print*, "Walk ", i, our_bird%x
!     print*, "Weight: ", our_bird%weight



! end do

! stop

! ---------------------------------------------------------------------------
!Test how a bird can do random flight
print *, DASHES
print *, "TEST: Test how a bird can do random flight"

our_bird%state_fear_hunger = -0.5
do i = 1, 50

    call our_bird%fly(our_test_habitat)
    print*, "Fly ", i, our_bird%x
end do


! ---------------------------------------------------------------------------
! Test if predator is working
print *, DASHES
print *, "TEST: Test if predator is working"

call evil_predator%init()     ! create a predator


print *, "Predator risk level: ", evil_predator%risk

print *, "Bird is alive: ", our_bird%is_alive
call our_bird%place(10)       ! place bird to the same cell 10
print *, "Bird is in cell ", our_bird%x, ", predation risk here is ",         &
        our_test_habitat%point(our_bird%x)%risk_of_predation

print *, "Do attacks ..."
j=1 ! counter of successful escapes

do i=1, 10
  print *, "TEST attack ", i
  weight_loss = our_bird%weight
  hunger_change = our_bird%state_fear_hunger
  call evil_predator%predator_attack_bird(our_bird, our_test_habitat,         &
                     p_is_here, prey_is_killed ) ! predator attacks the bird
  if (our_bird%is_alive) j=j+1  ! count the number of successful escapes = alive

  ! print *, "Bird is alive after attack ", i, our_bird%is_alive

  ! force the bird alive again because we need to calculate success rate
  ! fith new birds for testing
  print*, " MAX HUNGER = MIN FEAR = 10"
  print*, " MIN HUNGER = MAX FEAR = -10"
  print*, " Predator is in the bird's cell (PRESENT):", p_is_here
  print*, " Predator did kill the bird (KILLED = T, still alive = F):", prey_is_killed
  print*, " Weight before attack: ", weight_loss
  print*, " Weight after attack: ", our_bird%weight
  print*, DASHES
  print*, " Hunger before attack: ", hunger_change, "(genetic value: ", our_bird%gene_fear_hunger, ")"
  print*, " Hunger after attack: ", our_bird%state_fear_hunger
  print*, " Bird is alive: ", our_bird%is_alive

  if (our_bird%state_fear_hunger < our_bird%genetic_lower_hunger() ) print*, "WARNING: below hunger limit!"
  if (our_bird%state_fear_hunger > our_bird%genetic_upper_hunger() ) print*, "WARNING: over hunger limit!"

  our_bird%is_alive = .TRUE.
end do

print *, "Escape rate is ", real(j)/real(i-1), "(in ", i-1, "runs)."
associate ( p1 => our_test_habitat%point( our_bird%x )%risk_of_predation,     &
            p2 => evil_predator%risk )
    ! Calculate the probability of prey durvival given
    !          p1: Probability of predator to occur in the cell
    !          p2: Probability of successful attack when the prey dies
    p_escape = 1-p1*p2
end associate
print *, "  Note: It mus agree with theoretical survival", p_escape




! ! ---------------------------------------------------------------------------
print*, DASHES
print *, "testing  large number of birds attacked"
call our_test_habitat%init()
call bird_pop2%init_pop(100)

print*, DASHES

print*, "testing large number of birds attacked"


print*, "Testing how timesteps are running"
call our_test_habitat%init()
call bird_pop2%init_pop()

print *, "Do attacks ..."
j=1 ! counter of successful escapes

do i=1, 10

    call evil_predator%predator_attack_bird(bird_pop2%birds(i), our_test_habitat )
    !predator attacks the bird

    !call bird_pop2%birds(i)%hunger_genetic_limit()

    if (bird_pop2%birds(i)%state_fear_hunger < bird_pop2%birds(i)%genetic_lower_hunger() ) print*,                &
                                    "WARNING: below hunger limit!", bird_pop2%birds(i)%state_fear_hunger
    if (bird_pop2%birds(i)%state_fear_hunger > bird_pop2%birds(i)%genetic_upper_hunger() ) print*,                &
                                    "WARNING: over hunger limit!", bird_pop2%birds(i)%state_fear_hunger


end do

print*, DASHES
! ! ---------------------------------------------------------------------------

print*, "Testing how timesteps are running"
call our_test_habitat%init(1000)

print *, "Food availability (minval), (max):", minval(our_test_habitat%point%food_availability), &
         maxval(our_test_habitat%point%food_availability)

call bird_pop%init_pop(1000)
print*, "Number of birds alive before: ", count(bird_pop%birds%is_alive)
call bird_pop%time_steps(our_test_habitat, evil_predator)
print*, "Number of birds alive: ", count(bird_pop%birds%is_alive)
print*, "Number of birds dead: ", count(bird_pop%birds%is_alive .eqv. .FALSE.)

! ! ---------------------------------------------------------------------------
print*, DASHES


block !blocks confines the code within the block, so that values does not effect other non-related code.
    integer, parameter :: NUMBER_REPS = 10000
    logical, dimension (NUMBER_REPS) :: mutation
    print*, "Test mutation: "



    do i = 1, NUMBER_REPS
        call our_bird%mutate(is_mutate=mutation(i))
    end do
    print*, "Number of mutations: ", count(mutation), ", rate = ",      &
                                    real(count(mutation)) /NUMBER_REPS
                                  !count() counts the amount of logical happenings
end block

! ! ---------------------------------------------------------------------------
stop
print*, DASHES
print *, "Test sorting birds by mass"

print *, bird_pop%birds%weight
print *, "Min: ", minval(bird_pop%birds%weight)
print *, "Max: ", maxval(bird_pop%birds%weight)

call bird_pop%sort()

print *, bird_pop%birds%weight







end program first_model

