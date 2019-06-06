dune_add_test(NAME test_md_mortar_sharp_darcy1p_stokes1p
              LABELS multidomain multidomain_boundary mortar stokesdarcy 1p stokes
              SOURCES main.cc
              CMAKE_GUARD "( dune-foamgrid_FOUND AND HAVE_UMFPACK )"
              COMMAND ./test_md_mortar_sharp_darcy1p_stokes1p
              CMD_ARGS ../params.imput)
