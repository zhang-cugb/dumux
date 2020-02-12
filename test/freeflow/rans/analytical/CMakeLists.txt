add_input_file_links()

dumux_add_test(NAME test_ff_rans_analytical_komega
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=AnalyticalKOmega
              LABELS freeflow rans
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/a.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/a.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_rans_analytical_komega params.input
                        -Problem.Name test_ff_rans_analytical_komega")
