dumux_add_test(SOURCES test_particles_rotating_disc.cc
               LABELS particles
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --script fuzzy
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_particles_rotating_disc"
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_particles_rotating_disc-reference.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/test_particles_rotating_disc-00050.vtp)
