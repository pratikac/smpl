General RRT* planner inside the ros framework
=============================================

dynamical_system.h:

  state_c:
    Contains the state of the system. This can be a logical or physical state.
    Copy constructors and operator need to implemented as shown in the example
    class. Also implement a method that calculates |x - y| between two states
    x and y called "dist".

  control_c:
    Inherits from state_c. Used to store the control values for the system.
    It reimplements constructors because c++ does not inherit non-default
    constructors from the base class.

  trajectory_c:
    Contains a vector of states and controls. Also contains the length of the
    trajectory in the variable called total_variation. The necessary methods
    for this class are:
      clear:
        clears the trajectory and/or its memory
      append:
        appends two trajectories
      reverse:
        reverses the vector of states and controls (used in RRT*)

  optimization_data_c:
    This is an abstract class for any optimization data that one might want to
    store to quickly calculate the steering function between two states. It can
    contain, e.g., turning_radius (as shown in the example) or polynomial coefficients
    for spline interpolation or any other thing.

  dynamical_system_c:
    This is an abstract class for a dynamical system. The user will typically inherit
    from this and implement the "extend_to" and "evaluate_extend_cost" functions.

dubins.h
  
  This implements a class dubins_optimization_data_c which is inherited from the
  general optimization_data_c class. It also implements dubins_c which is derived
  from dynamical_system_c.

map.h:
  This is an abstract class for maps. The user will typicall derive from this and
  implement the three methods:
    sample_free_space:
      Generates a random sample of dimension N from free space
    is_in_collision
      Returns true if the sample is in collision, false otherwise
    get_state_cost:
      Returns a float which is the cost for a particular state, e.g., higher cost
      for going into the left lane while driving.
