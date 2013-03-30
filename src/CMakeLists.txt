# Simple CMake file

# esac algorithms

include_directories("${PROJECT_SOURCE_DIR}/stxxl/include" "${PROJECT_BINARY_DIR}/stxxl/include")

if(BUILD_ALL)

  add_executable(esactest-skew esactest.cc)
  target_link_libraries(esactest-skew stxxl ${ZLIB_LIBRARIES})
  set_property(TARGET esactest-skew PROPERTY COMPILE_DEFINITIONS "COMPILE_ALGORITHM=skew" "LCP_CALC=0")

  add_executable(esactest-skewlcp esactest.cc)
  target_link_libraries(esactest-skewlcp stxxl ${ZLIB_LIBRARIES})
  set_property(TARGET esactest-skewlcp PROPERTY COMPILE_DEFINITIONS "COMPILE_ALGORITHM=skew" "LCP_CALC=1")

  add_executable(esactest-eSAIS esactest.cc)
  target_link_libraries(esactest-eSAIS stxxl ${ZLIB_LIBRARIES})
  set_property(TARGET esactest-eSAIS PROPERTY COMPILE_DEFINITIONS "COMPILE_ALGORITHM=esais" "ESAIS_DISCARD_UNIQUES=0" "LCP_CALC=0")

  add_executable(esactest-eSAISdisc esactest.cc)
  target_link_libraries(esactest-eSAISdisc stxxl ${ZLIB_LIBRARIES})
  set_property(TARGET esactest-eSAISdisc PROPERTY COMPILE_DEFINITIONS "COMPILE_ALGORITHM=esais" "ESAIS_DISCARD_UNIQUES=1" "LCP_CALC=0")

  add_executable(esactest-eSAISlcpint esactest.cc)
  target_link_libraries(esactest-eSAISlcpint stxxl ${ZLIB_LIBRARIES})
  set_property(TARGET esactest-eSAISlcpint PROPERTY COMPILE_DEFINITIONS "COMPILE_ALGORITHM=esais" "ESAIS_DISCARD_UNIQUES=0" "LCP_CALC=1" "ESAIS_LCP_CALC_INT=1")

  add_executable(esactest-eSAISlcpext esactest.cc)
  target_link_libraries(esactest-eSAISlcpext stxxl ${ZLIB_LIBRARIES})
  set_property(TARGET esactest-eSAISlcpext PROPERTY COMPILE_DEFINITIONS "COMPILE_ALGORITHM=esais" "ESAIS_DISCARD_UNIQUES=0" "LCP_CALC=1" "ESAIS_LCP_CALC_EXT=1")

  if(WITH_MALLOC_COUNT)

    add_definitions(-DMALLOC_COUNT)

    add_library(malloc_count SHARED tools/malloc_count.c)
    set_target_properties(malloc_count PROPERTIES COMPILE_FLAGS "-Wstrict-prototypes -fPIC -D_GNU_SOURCE")
    target_link_libraries(malloc_count dl)

    target_link_libraries(esactest-skew malloc_count)
    target_link_libraries(esactest-skewlcp malloc_count)
    target_link_libraries(esactest-eSAIS malloc_count)
    target_link_libraries(esactest-eSAISdisc malloc_count)
    target_link_libraries(esactest-eSAISlcpint malloc_count)
    target_link_libraries(esactest-eSAISlcpext malloc_count)

  endif(WITH_MALLOC_COUNT)

else()

  add_executable(esactest esactest.cc)
  target_link_libraries(esactest stxxl ${ZLIB_LIBRARIES})

  if(WITH_MALLOC_COUNT)

    add_definitions(-DMALLOC_COUNT)

    add_library(malloc_count SHARED tools/malloc_count.c)
    set_target_properties(malloc_count PROPERTIES COMPILE_FLAGS "-Wstrict-prototypes -fPIC -D_GNU_SOURCE")
    target_link_libraries(malloc_count dl)

    target_link_libraries(esactest malloc_count)

  endif(WITH_MALLOC_COUNT)

endif()

#set_property(TARGET esac PROPERTY COMPILE_DEFINITIONS "STXXL_VERBOSE_LEVEL=0")
