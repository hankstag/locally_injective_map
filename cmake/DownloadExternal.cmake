################################################################################
include(DownloadProject)

# With CMake 3.8 and above, we can hide warnings about git being in a
# detached head by passing an extra GIT_CONFIG option
if(NOT (${CMAKE_VERSION} VERSION_LESS "3.8.0"))
    set(FLATTEN_EXTRA_OPTIONS "GIT_CONFIG advice.detachedHead=false")
else()
    set(FLATTEN_EXTRA_OPTIONS "")
endif()

# Shortcut function
function(flatten_download_project name)
    download_project(
        PROJ         ${name}
        SOURCE_DIR   ${FLATTEN_EXTERNAL}/${name}
        DOWNLOAD_DIR ${FLATTEN_EXTERNAL}/.cache/${name}
        QUIET
        ${FLATTEN_EXTRA_OPTIONS}
        ${ARGN}
    )
endfunction()

################################################################################

# ## libigl LGPL - triangle - tetgen
# function(flat_download_libigl)
#     polyfem_download_project(libigl
#         GIT_REPOSITORY https://github.com/libigl/libigl.git
#         GIT_TAG        c7c06e3735cdf6188bd17507403362065c4ae9dc
#     )
# endfunction()

# ## Json MIT
# function(flat_download_json)
#     flat_download_project(json
#         GIT_REPOSITORY https://github.com/jdumas/json
#         GIT_TAG        0901d33bf6e7dfe6f70fd9d142c8f5c6695c6c5b
#     )
# endfunction()

## CppNumericalSolvers MIT
#  function(polyfem_download_CppNumericalSolvers)
#     polyfem_download_project(CppNumericalSolvers
#         GIT_REPOSITORY https://github.com/PatWie/CppNumericalSolvers.git
#         GIT_TAG        7eddf28fa5a8872a956d3c8666055cac2f5a535d
#     )
# endfunction()



endfunction()