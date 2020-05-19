dune_symlink_to_source_files(FILES "params.input")

dumux_add_test(SOURCES test_timestepping.cc
               LABELS unit timestepping)
set(CMAKE_BUILD_TYPE Release)
