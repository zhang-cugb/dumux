# additional macros
include(AddGstatFileLinks)
include(AddInputFileLinks)
include(DumuxDoxygen)
include(DumuxTestMacros)

find_package(GLPK)
find_package(Gnuplot)
set(HAVE_GNUPLOT ${GNUPLOT_FOUND})
find_package(Gstat)
find_package(Gmsh)
find_package(NLOPT)
find_package(PTScotch)
find_package(PVPython)
find_package(Valgrind)
