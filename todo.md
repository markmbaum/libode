## TODO
#### Improvements to make, if it's ever worth the time:

1. Good and robust adaptive step sizing for implicit methods!
2. Handling non-convergence of Newton's method(s) with implicit methods. The Newton solver can return a code indicating success/failure. Code was modified to do that in other repositories. If failure, do what though? Take a smaller step?
3. Initial guess of k values for nonlinear solves in implicit methods.
4. L-stable Rosenbrock method.
5. Look into extra steps taken when snapping in an adaptive solves.
5. Tridiagonal Jacobian support? The matrix routine is already in `linalg.cc` but there isn't a way to access it with implicit methods.
