dune_symlink_to_source_files(FILES "params.input")

dumux_add_test(NAME test_fem_algebraic
              SOURCES test_fem_algebraic.cc
              LABELS unit discretization fem)
set(CMAKE_BUILD_TYPE Debug)
