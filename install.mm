# Installation script for the RAG package.
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

if not(evalb(whattype(eval(RAG)) = package)) then
  printf("Error: failed to load RAG package.\n"):
  quit:
end if:

march(`create`, mladirname):
savelib(`RAG`):
printf("Installed RAG package into %s\n", mladirname):
