# Internal Details

This is a page for detailing some of the inner workings to help future
contributors to the library.

## Observables and Variable Elimination

In the variable "elimination" algorithms, what is actually done is that variables
are removed from being states and equations are moved into the `observed` category
of the system. The `observed` equations are explicit algebraic equations which
are then substituted out to completely eliminate these variables from the other
equations, allowing the system to act as though these variables no longer exist.

However, as a user may have wanted to interact with such variables, for example,
plotting their output, these relationships are stored and are then used to
generate the `observed` equation found in the `SciMLFunction` interface, so that
`sol[x]` lazily reconstructs the observed variable when necessary. In this sense,
there is an equivalence between observables and the variable elimination system.
