if (MKL_cmake_included)
    return()
endif ()
set(MKL_cmake_included true)

function(detect_mkl LIBNAME)
    if (HAVE_MKL)
        return()
    endif ()

    message(STATUS "Detecting Intel(R) MKL: trying ${LIBNAME}")

    find_path(MKLINC mkl_cblas.h
            HINTS ${MKLROOT}/include $ENV{MKLROOT}/include)
    if (NOT MKLINC)
        file(GLOB_RECURSE MKLINC
                ${CMAKE_CURRENT_SOURCE_DIR}/external/*/mkl_cblas.h)
        if (MKLINC)
            # if user has multiple version under external/ then guess last
            # one alphabetically is "latest" and warn
            list(LENGTH MKLINC MKLINCLEN)
            if (MKLINCLEN GREATER 1)
                list(SORT MKLINC)
                list(REVERSE MKLINC)
                # message(STATUS "MKLINC found ${MKLINCLEN} files:")
                # foreach(LOCN IN LISTS MKLINC)
                #     message(STATUS "       ${LOCN}")
                # endforeach()
                list(GET MKLINC 0 MKLINCLST)
                set(MKLINC "${MKLINCLST}")
                # message(WARNING "MKLINC guessing... ${MKLINC}.  "
                #     "Please check that above dir has the desired mkl_cblas.h")
            endif ()
            get_filename_component(MKLINC ${MKLINC} PATH)
        endif ()
    endif ()
    if (NOT MKLINC)
        return()
    endif ()

    get_filename_component(__mklinc_root "${MKLINC}" PATH)
    find_library(MKLLIB NAMES ${LIBNAME}
            HINTS ${MKLROOT}/lib ${MKLROOT}/lib/intel64
            $ENV{MKLROOT}/lib $ENV{MKLROOT}/lib/intel64
            ${__mklinc_root}/lib ${__mklinc_root}/lib/intel64)
    if (NOT MKLLIB)
        return()
    endif ()

    if (WIN32)
        set(MKLREDIST ${MKLINC}/../../redist/)
        find_file(MKLDLL NAMES ${LIBNAME}.dll
                HINTS
                ${MKLREDIST}/mkl
                ${MKLREDIST}/intel64/mkl
                ${__mklinc_root}/lib)
        if (NOT MKLDLL)
            return()
        endif ()
    endif ()

    if (NOT CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
        get_filename_component(MKLLIBPATH ${MKLLIB} PATH)
        find_library(MKLIOMP5LIB
                NAMES "iomp5" "iomp5md" "libiomp5" "libiomp5md"
                HINTS ${MKLLIBPATH}
                ${MKLLIBPATH}/../../lib
                ${MKLLIBPATH}/../../../lib/intel64
                ${MKLLIBPATH}/../../compiler/lib
                ${MKLLIBPATH}/../../../compiler/lib/intel64)
        if (NOT MKLIOMP5LIB)
            return()
        endif ()
        if (WIN32)
            find_file(MKLIOMP5DLL
                    NAMES "libiomp5.dll" "libiomp5md.dll"
                    HINTS ${MKLREDIST}/../compiler ${__mklinc_root}/lib)
            if (NOT MKLIOMP5DLL)
                return()
            endif ()
        endif ()
    else ()
        set(MKLIOMP5LIB)
        set(MKLIOMP5DLL)
    endif ()

    get_filename_component(MKLLIBPATH "${MKLLIB}" PATH)
    string(FIND "${MKLLIBPATH}" ${CMAKE_CURRENT_SOURCE_DIR}/external __idx)
    if (${__idx} EQUAL 0)
        if (WIN32)
            install(PROGRAMS ${MKLDLL} DESTINATION lib)
        else ()
            install(PROGRAMS ${MKLLIB} DESTINATION lib)
        endif ()
        if (MKLIOMP5LIB)
            if (WIN32)
                install(PROGRAMS ${MKLIOMP5DLL} DESTINATION lib)
            else ()
                install(PROGRAMS ${MKLIOMP5LIB} DESTINATION lib)
            endif ()
        endif ()
    endif ()

    if (WIN32)
        # Add paths to DLL to %PATH% on Windows
        get_filename_component(MKLDLLPATH "${MKLDLL}" PATH)
        set(CTESTCONFIG_PATH "${CTESTCONFIG_PATH}\;${MKLDLLPATH}")
        set(CTESTCONFIG_PATH "${CTESTCONFIG_PATH}" PARENT_SCOPE)
    endif ()

    # TODO: cache the value
    set(HAVE_MKL TRUE PARENT_SCOPE)
    set(MKLINC ${MKLINC} PARENT_SCOPE)
    set(MKLLIB "${MKLLIB}" PARENT_SCOPE)

    if (WIN32)
        set(MKLDLL "${MKLDLL}" PARENT_SCOPE)
    endif ()
    if (MKLIOMP5LIB)
        set(MKLIOMP5LIB "${MKLIOMP5LIB}" PARENT_SCOPE)
    endif ()
    if (WIN32 AND MKLIOMP5DLL)
        set(MKLIOMP5DLL "${MKLIOMP5DLL}" PARENT_SCOPE)
    endif ()
endfunction()

detect_mkl("mklml_intel")
detect_mkl("mklml")
detect_mkl("mkl_rt")

if(HAVE_MKL)
    add_definitions(-DUSE_MKL -DUSE_CBLAS)
    include_directories(AFTER ${MKLINC})
    list(APPEND MKL_LINKER_LIBS ${MKLLIB})

    set(MSG "Intel(R) MKL:")
    message(STATUS "${MSG} include ${MKLINC}")
    message(STATUS "${MSG} lib ${MKLLIB}")
    if (MKLIOMP5LIB)
        message(STATUS "${MSG} OpenMP lib ${MKLIOMP5LIB}")
    else ()
        message(STATUS "${MSG} OpenMP lib provided by compiler")
    endif ()
    if (WIN32)
        message(STATUS "${MSG} dll ${MKLDLL}")
        if (MKLIOMP5DLL)
            message(STATUS "${MSG} OpenMP dll ${MKLIOMP5DLL}")
        else ()
            message(STATUS "${MSG} OpenMP dll provided by compiler")
        endif ()
    endif ()
else ()
    if (DEFINED ENV{FAIL_WITHOUT_MKL} OR DEFINED FAIL_WITHOUT_MKL)
        set(SEVERITY "FATAL_ERROR")
    else ()
        set(SEVERITY "WARNING")
    endif ()
    message(${SEVERITY}
            "Intel(R) MKL not found. Please get a full version from "
            "https://software.intel.com/en-us/intel-mkl")
endif ()

