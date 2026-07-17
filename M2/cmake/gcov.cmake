# gcc/gcov code coverage for the engine (configure with -DGCOV=ON).
# Instrumentation flags are added in configure.cmake; this module only provides
# convenience targets to clear and report on the accumulated .gcda data.  The data
# itself is generated automatically whenever a binary built in coverage mode exits
# (the unit tests, ctest, or the M2 binary run on your own code).
#
# Usage:
#   cmake -S . -B BUILD/cov -GNinja -DGCOV=ON -DCMAKE_BUILD_TYPE=Debug -DBUILD_TESTING=ON
#   cmake --build BUILD/cov --target M2-engine M2-unit-tests
#   cmake --build BUILD/cov --target coverage-reset   # start from a clean slate
#   ctest --test-dir BUILD/cov -R unit-tests          # (or run M2 however you like)
#   cmake --build BUILD/cov --target coverage-report  # writes coverage/index.html

if(GCOV)
  find_program(GCOVR NAMES gcovr)
  set(_engine_objdir ${CMAKE_BINARY_DIR}/Macaulay2/e/CMakeFiles/M2-engine.dir)
  set(_coverage_dir  ${CMAKE_BINARY_DIR}/coverage)

  # coverage-reset clears all accumulated data; the M2 binary may have written
  # .gcda anywhere in the tree, so this is deliberately not engine-scoped.
  add_custom_target(coverage-reset
    COMMENT "Deleting accumulated .gcda coverage counters"
    COMMAND find ${CMAKE_BINARY_DIR} -name "*.gcda" -delete)

  if(NOT GCOVR)
    message(WARNING "GCOV is ON but gcovr was not found; the 'coverage-report' target will not be created. Install gcovr (e.g. 'pip install gcovr').")
  else()
    # Report is scoped to the engine by default; run gcovr by hand for other scopes.
    add_custom_target(coverage-report
      COMMENT "Generating gcov/gcovr coverage report for the engine"
      COMMAND ${CMAKE_COMMAND} -E make_directory ${_coverage_dir}
      COMMAND ${GCOVR}
        --root ${CMAKE_SOURCE_DIR}/Macaulay2/e
        --object-directory ${_engine_objdir}
        --exclude ".*/unit-tests/.*"
        --html-details ${_coverage_dir}/index.html
        --print-summary
      COMMAND ${CMAKE_COMMAND} -E echo "Coverage report: ${_coverage_dir}/index.html"
      USES_TERMINAL)
  endif()
endif()
