# Deciding the existence of solutions

**HasRealSolutions**(eqs, positive, nonzero, opts)
where 
- eqs, positive and nonzero are lists of polynomials with rational
  coefficients
    - eqs stands for _equations_
    - positive stands for _positivity contraints_
    - nonzer stands for _non-vanishing constraints_
- opts is a set indicating some options (this is an optional argument)

**HasRealSolutions** decides if the semi-algebraic set defined by 
- the vanishing of all polynomials in eqs,
- the positivity of all polynomials in positive,
- the non-vanishing of all polynomials in nonzero
is empty or not. 

In case of emptiness, it returns an empty list,
otherwise it returns a list of witness points of non-emptiness.

Each such point is given by an isolating box.

# Computing at least one point per connected component

**PointsPerComponents**(eqs, positive, nonzero, opts)
where 
- eqs, positive and nonzero are lists of polynomials with rational
  coefficients
    - eqs stands for _equations_
    - positive stands for _positivity contraints_
    - nonzer stands for _non-vanishing constraints_
- opts is a set indicating some options (this is an optional argument)

**PointsPerComponents** computes at least one point in each connected
component of the semi-algebraic set defined by 
- the vanishing of all polynomials in eqs,
- the positivity of all polynomials in positive,
- the non-vanishing of all polynomials in nonzero.

It returns a list of points meeting all connected components of the
semi-algebraic set under study.  
Each such point is given by an isolating box.

# Dependency

SemiAlgebraicSolve requires the [msolve](https://msolve.lip6.fr)
library to be installed with its interface to Maple. See
[https://msolve.lip6.fr](https://msolve.lip6.fr) for documentation on
[msolve](https://msolve.lip6.fr) and
[https://github.com/algebraic-solving/msolve](https://github.com/algebraic-solving/msolve)
for its source code and installation instructions.

# Options

Possible options are: 

- "verb"=INT (default value is 0)

  This option controls the verbosity. 

  Setting INT to 1 will display some global information on the
  computation, hiding details of the computations performed by
  msolve. 

  Setting INT to 2 will display both global informations on the
  computation and will activate the verbosity of msolve, hence
  displaying some data on the computations it performs.

- "nthreads"=INT (default value is 1)

  This option tells msolve how many threads it can use. At the
  moment, RAGlib does not support multi-threading; this is delegated
  to msolve computations.

- "isbounded"=INT (default value is 0)

  When it is known in advance that the semi-algebraic under study is
  bounded (for instance because some inequalities define a box or a
  ball), setting INT is set to 1, will activate a dedicated
  algorithm that take advantage of this information (and is usually
  faster than more general purpose algorithms). 

  This option must be used with care. It should not be used if there
  is a doubt on the boundedness of the semi-algebraic set under
  study. 

