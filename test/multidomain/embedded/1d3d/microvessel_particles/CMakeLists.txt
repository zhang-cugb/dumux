dune_add_test(NAME test_md_embedded_1d3d_microvessel_particles
              SOURCES main.cc)

dune_symlink_to_source_files(FILES "params.input")
