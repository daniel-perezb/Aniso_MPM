#Generate the local Makefile

if(CREATE_LOCAL_MAKEFILES)
  file(RELATIVE_PATH ZIRAN_RELATIVE_PATH ${ZIRAN_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
  configure_file(${ZIRAN_SOURCE_DIR}/Scripts/Makefile.in ${CMAKE_CURRENT_SOURCE_DIR}/Makefile @ONLY)
endif()

add_ziran_executable(anisofracture anisofracture.cpp)
target_compile_options(anisofracture PRIVATE ${ZIRAN_CXXFLAGS_LTO})
target_link_libraries(anisofracture ziran mpmdev tbbmalloc tbbmalloc_proxy ${ZIRAN_CXXFLAGS_LTO})

