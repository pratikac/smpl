# Create a shared library lib${POD_NAME}.so with all source files
file(GLOB cc_files *.cc) 
file(GLOB cc_files *.c) 
file(GLOB cpp_files *.cpp)
set(all_sources ${cpp_files} ${cc_files})
add_library(${POD_NAME} SHARED ${all_sources})

add_executable(test_main test_main.cpp kdtree.c)

# c++11 experimental flags
set(CMAKE_CXX_FLAGS "-std=c++0x -o3")

# uncomment these lines to link against another library via pkg-config
#set(REQUIRED_PACKAGES pkg_a pkg_b)
#pods_use_pkg_config_packages(${POD_NAME} ${REQUIRED_PACKAGES})

# uncomment these lines to link against another library
#set(REQUIRED_LIBS lib_a)
#target_link_libraries (${POD_NAME} ${REQUIRED_LIBS})


# install all headers to include/${POD_NAME}
# make an aggregate header
file(GLOB h_files *.h)
file(GLOB hpp_files *.hpp)
file(GLOB hxx_files *.hxx)
set(all_headers ${h_files} ${hpp_files} ${hxx_files})

set(aggregate_header_name ${CMAKE_CURRENT_SOURCE_DIR}/${POD_NAME}.h)
string(TOUPPER "__${aggregate_header_name}__" aggregate_header_guard) 
file(WRITE ${aggregate_header_name}
   "/** THIS IS AN AUTOMATICALLY GENERATED FILE.  DO NOT MODIFY BY HAND! */\n" 
   "#ifndef ${aggregate_header_guard}\n"
   "#define ${aggregate_header_guard}\n\n")

foreach(header ${all_headers})
    file(RELATIVE_PATH rel_header ${CMAKE_CURRENT_SOURCE_DIR} ${header})
    if (NOT ${header} STREQUAL ${aggregate_header_name})
        file(APPEND ${aggregate_header_name} "#include \"${rel_header}\"\n")
    endif()
endforeach()      
file(APPEND ${aggregate_header_name} "\n#endif //${aggregate_header_guard} \n")


# install all the headers
pods_install_headers(${aggregate_header_name} ${all_headers} 
                        DESTINATION ${POD_NAME})

# make the library public
pods_install_libraries(${POD_NAME})


# create a pkg-config file for the library
pods_install_pkg_config_file(${POD_NAME}
    CFLAGS
    LIBS -l${POD_NAME} ${REQUIRED_LIBS}
    REQUIRES ${REQUIRED_PACKAGES}
    VERSION 0.0.1)
