TODO: 

* maybe we need to change the policy base class for explicit methods. 
Right now each method has its own base policy, but it seems to me 
that all expicit methods have the same policy base class. 
So it makes no sense to have separate ones. 
NOTE THAT THIS DOES NOT HOLD TRUE FOR IMPLICIT METHODS WHERE THINGS 
CHANGE FOR EVERY METHOD.

* fix how time is passed around, because now it is defaulted to 
ode::details::time_type, which is a double.

* fix initialization of the states and containers because this can create 
issues. Find a way to have a consistent initializataion




June 2nd: 

Need to think carefully about how to solve the problem with Euler implicit 
when we need to modify the jacobian and state vector using the left-sing vectors from the svd. 
Because i believe theere has to be a way to do this without duplicating the method 
all over whenever we need to do a new projection. 
Ideas: decorator. policy. Or maybe just to SFINAE.

I think one solution is the have a standard euler (or whatever) method in the 
ode package, and a modified version inside the ROM package to fit need of reduced case. 
IT IS KEY that we separate interface from implementation. Because the actual 
signature implementatation of the residual and jacobian for the objective function is the same, 
whether we have a full state or a reduced one. I mean, the math is the same, the only thing 
that changes is how we prepare the data to supply to these functions. 
The standard class has regular stuff for doing regular integration. 
The modified class, which exists inside the rom, is an actual stepper class but it also 
contains methods to preprocess the states and res and jacob, i mean projecting them 
using the basis. This class does not know anything about HOW this is done, it simply 
calls the target functor provided which is an object of one of the ROM methods that knows 
how to project a state vector or residual or jacobian because it owns the basis.

if we do this, we need to have only one class defining a target ROM algorithm, which can 
then be used in association with any other time stepper as long as there is a time stepper 
class for doing time integration for the ``reduced'' case using a particular time marching scheme.

Why is this better than the other scenario where instead we define a class for a target 
ROM method for any possible time marching scheme? 
say that we have LSPG. In this case, we can have different ways to rescale 
the residual because we can do collocation, identity or gnat, etc. 
So if we were to couple one of these to any time marching scheme, we would 
end up writing a million combinations of (a) time scheme and (b) particular ROM algorithm. 
Instead, we can generalize the time marching scheme by making it general enough so that 
it provided the signature of the steps one should follow to do that particular scheme. 
It is then up to the concrete ROM class to provide different implementations. 

Looks like policy-based design is promising: basically, we refactor the ODE stepper 
classes (maybe just the implicit ones) so that they take 3 policies: 
(a) policy for how to compute the residual, (b) how to compute jacobian


