# Installation script for RAGlib procedures.
# Run from the repository root with:
# maple -q install.mm

homedir:=kernelopts(homedir):
savelibname:=cat(homedir, "/libs/"):
mladirname:=cat(savelibname, "rag.mla"):
ssystem(cat("rm -f ", mladirname)):
libname:=savelibname,libname:

try
  with(MSolve):
catch:
  printf("Error: failed to load MSolve. Please install the file interface first.\n"):
  quit:
end try:

raglib_root:=currentdir():
raglib_src_dir:=cat(raglib_root, "/src"):
kernelopts(includepath=raglib_src_dir):

raglib_src:=cat(raglib_src_dir, "/rag.mm"):
try
  read raglib_src:
catch:
  printf("Error: failed to load RAGlib sources. Please run the installation script from the repository root.\n"):
  quit:
end try:

rag_procs:=[
  IsMinimal, MinimalGeneratorsDichotomy, MinimalGenerators, WitnessLift,
  MonomialSupport, TableCoeffsSinglePoly, OldTableCoeffs, TableCoeffs,
  LiftPolynomials, NewValuesWitness, GeneratorsLift, ElimModSatIntersect,
  ElimModSatIntersectLM, ModSatIntersect, ModSatIntersectLM,
  ElimSaturateIntersect, SaturateIntersect, ComputeMaximalMinors, IsRegular,
  HaveFiniteIntersections, GoodFiberValue_svars, GoodFiberValue,
  DegreeTruncate, Increment, NextForms, TestGenericLineDegreeSingular,
  TestGenericLineDegreeRegular, HasFiniteCriticalLocus, FindGenericLineRegular,
  FindGenericLineSingular, CoeffDeform_eps, CoeffDeform, FindGenericLine,
  SmallMidRational, ConstructFibers, HasOverLapCoupleOfIntervals, HasOverLap,
  ComputeBoundsRegular, ManageOverLapComputeBoundsSingular,
  ElimComputeBoundsSingular, ComputeBoundsSingular, ComputeBounds,
  ModularLimits, LimitsDeformedCriticalPoints, CriticalPointsSingular,
  SplitSystem_cstr, SplitSystem, CriticalPoints, ExactSolSelection,
  GenerateCriticalPointsFamilies, FamCriticalPoints, UnivariateSolveFamily,
  AdmissibleSolutions, GenerateDeformedFamilies_eps,
  OLDGenerateDeformedFamilies_eps, UnboundedComponents,
  DegenerateDeformedSystem, InfiniteBranches, ZeroDimBoundaries, SolveFamily,
  SemiAlgebraicSolveIterateOnFamilies, PointsPerComponentsAlgebraic,
  SemiAlgebraicSolve, HasRealSolutions, PointsPerComponents
]:

march(`create`, mladirname):

for p in rag_procs do
  savelib(p):
end do:

printf("Installed RAGlib procedures (total %d) into %s\n", nops(rag_procs), mladirname):
