project(md_unittest VERSION 0.1 LANGUAGES C)

add_executable(md_unittest unittest.c)

target_compile_definitions(md_unittest PRIVATE MD_UNITTEST_DATA_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/../test_data\" MD_UNITTEST_SOURCE_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}\" UTEST_USE_OLD_QPC ${MD_DEFINES})
target_compile_options(md_unittest PRIVATE ${MD_FLAGS_EXT} ${MD_FLAGS_AVX} $<$<CONFIG:Debug>:${MD_FLAGS_DEB}> $<$<CONFIG:Release>:${MD_FLAGS_REL}>)
target_link_options(md_unittest PRIVATE ${MD_LINK_FLAGS} $<$<CONFIG:Debug>:${MD_LINK_FLAGS_DEB}> $<$<CONFIG:Release>:${MD_LINK_FLAGS_REL}>)
target_link_libraries(md_unittest PRIVATE mdlib ${MD_LIBS})

set (SRC_FILES
    test_allocator.c
    test_array.c
    test_atomic.c
    test_bitop.c
    test_grid.c
    test_bitset.c
    test_bonds.c
    test_chem.c
    test_contact.c
    test_hbond.c
    test_water.c
    test_nonbonded.c
    test_infer_atomic_number.c
    test_str.c
    test_hash.c
    test_handle.c
    test_edr.c
    test_cube.c
    test_dcd.c
    test_fifo.c
    test_intrinsics.c
    test_simd.c
    test_simd_ops.c
    test_vec_math.c
    test_os.c
    test_pdb.c
    test_attributes.c
    test_gro.c
    test_itp.c
    test_tpr.c
    test_trr.c
    test_mmcif.c
    test_match.c
    test_ring.c
    test_lipid.c
    test_monomer.c
    test_util.c
    test_secondary_structure.c
    test_unit.c
    test_script.c
    test_filter.c
    test_vial_spec.c
    test_vial_parse.c
    test_vial_fuzz.c
    test_spatial_hash.c
    test_writer.c
    test_csv.c
    test_base64.c
    test_xtc.c
    test_xyz.c
    test_ase_traj.c
    test_pot.c
    test_parse.c
    test_json.c
    test_vector_graphics.c
    test_flow.c
    test_gisaxs.c
    test_xvg.c
    rmsd.c
    test_lammps.c
    test_system_invariants.c
    test_run_cache.c
    test_molden.c
    test_gto_int.c
)

if (MD_ENABLE_HDF5)
    set (SRC_FILES ${SRC_FILES} test_vlx.c test_gto.c test_gto_iso.c test_gto_dmat.c test_trexio.c test_h5md.c)
endif()

if (MD_ENABLE_GPU AND (MD_GPU_BACKEND STREQUAL "VULKAN" OR MD_GPU_BACKEND STREQUAL "METAL"))
  set (SRC_FILES ${SRC_FILES} test_gpu.c test_gpu_tlsf.c test_gpu_render.c)
  # The kernels and shaders the GPU tests launch. They belong to the test
  # executable: nothing in mdlib uses them.
  compile_gpu_shaders(gpu_test_shaders.inl
    TARGET    md_unittest
    NAMESPACE md_shader
    SOURCE    ${CMAKE_CURRENT_SOURCE_DIR}/shaders/gpu_test.slang
    ENTRIES   fill scale_add sum_reduce tex_write tex_read bump layout_probe tex_probe sample_read spin
  )
  compile_gpu_shaders(gpu_render_test_shaders.inl
    TARGET    md_unittest
    NAMESPACE md_shader
    SOURCE    ${CMAKE_CURRENT_SOURCE_DIR}/shaders/gpu_render_test.slang
    ENTRIES   read_tex gen_pos
    VERTEX    basic_vs indirect_vs point_vs
    FRAGMENT  color_fs mrt_fs id_fs depth_fs
  )
  if (MD_ENABLE_HDF5)
    # test_gto_iso.c: gridless GTO isosurface extraction (prototype)
    compile_gpu_shaders(gto_iso_shaders.inl
      TARGET    md_unittest
      NAMESPACE md_shader
      SOURCE    ${CMAKE_CURRENT_SOURCE_DIR}/shaders/gto_iso.slang
      ENTRIES   iso_grid iso_octree_level iso_brick_hash iso_corners iso_cells iso_quads iso_project iso_relax iso_triangulate
      DEPENDS   ${CMAKE_CURRENT_SOURCE_DIR}/../src/shaders/md_gpu.slang
    )
  endif()
  compile_gpu_shaders(gpu_imgui_test_shaders.inl
    TARGET    md_unittest
    NAMESPACE md_shader
    SOURCE    ${CMAKE_CURRENT_SOURCE_DIR}/shaders/gpu_imgui_test.slang
    VERTEX    imgui_vs
    FRAGMENT  imgui_fs
  )
endif()

target_sources(md_unittest PRIVATE ${SRC_FILES})

target_include_directories(md_unittest
  PRIVATE
  ../ext/xtc/
  ${mdlib_BINARY_DIR}/gen   # generated shader headers of mdlib (test_gto_iso.c includes md_gto_iso.c); not
                            # CMAKE_BINARY_DIR, which is the parent's when mdlib is built as a subproject
)

# adapted from https://github.com/Kitware/CMake/blob/master/Modules/GoogleTest.cmake
# and https://github.com/Kitware/CMake/blob/master/Modules/GoogleTestAddTests.cmake
function(discover_tests TARGET)
  set(options)
  set(oneValueArgs
    TARGET
    WORKING_DIRECTORY
    DISCOVERY_TIMEOUT
  )
  set(multiValueArgs
    PROPERTIES
  )
  cmake_parse_arguments(
    ""
    "${options}"
    "${oneValueArgs}"
    "${multiValueArgs}"
    ${ARGN}
  )

  if(NOT _WORKING_DIRECTORY)
    set(_WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
  endif()
  set(_TEST_LIST ${TARGET}_TESTS)
  if(NOT _DISCOVERY_TIMEOUT)
    set(_DISCOVERY_TIMEOUT 5)
  endif()

  set(_DISCOVERY_MODE POST_BUILD)

  get_property(
    has_counter
    TARGET ${TARGET}
    PROPERTY CTEST_DISCOVERED_TEST_COUNTER
    SET
  )
  if(has_counter)
    get_property(
      counter
      TARGET ${TARGET}
      PROPERTY CTEST_DISCOVERED_TEST_COUNTER
    )
    math(EXPR counter "${counter} + 1")
  else()
    set(counter 1)
  endif()
  set_property(
    TARGET ${TARGET}
    PROPERTY CTEST_DISCOVERED_TEST_COUNTER
    ${counter}
  )

  # Define rule to generate test list for aforementioned test executable
  set(ctest_file_base "${CMAKE_CURRENT_BINARY_DIR}/${TARGET}[${counter}]")
  set(ctest_include_file "${ctest_file_base}_include.cmake")
  set(ctest_tests_file "${ctest_file_base}_tests.cmake")

  add_custom_command(
    TARGET ${TARGET} POST_BUILD
    BYPRODUCTS "${ctest_tests_file}"
    COMMAND "${CMAKE_COMMAND}"
            -D "TEST_TARGET=${TARGET}"  # maybe not needed?
            -D "TEST_EXECUTABLE=$<TARGET_FILE:${TARGET}>"
            -D "TEST_WORKING_DIR=${_WORKING_DIRECTORY}"
            -D "TEST_PROPERTIES=${_PROPERTIES}"
            -D "CTEST_FILE=${ctest_tests_file}"
            -D "TEST_DISCOVERY_TIMEOUT=${_DISCOVERY_TIMEOUT}"
            -P "${CMAKE_CURRENT_LIST_DIR}/add_tests.cmake"
    VERBATIM
  )

  file(WRITE "${ctest_include_file}"
    "if(EXISTS \"${ctest_tests_file}\")\n"
    "  include(\"${ctest_tests_file}\")\n"
    "else()\n"
    "  add_test(${TARGET}_NOT_BUILT ${TARGET}_NOT_BUILT)\n"
    "endif()\n"
  )

  # Add discovered tests to directory TEST_INCLUDE_FILES
  set_property(DIRECTORY
    APPEND PROPERTY TEST_INCLUDE_FILES "${ctest_include_file}"
  )
endfunction()

discover_tests(md_unittest
  PROPERTIES
    LABELS "unit"
  )
