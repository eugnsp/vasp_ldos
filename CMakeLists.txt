cmake_minimum_required(VERSION 3.13)
project(vasp_ldos VERSION 1.0)

set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/external/findFFTW")

add_executable(vasp_ldos src/vasp_ldos.cpp)

target_compile_features(vasp_ldos PUBLIC cxx_std_17)
target_compile_options(vasp_ldos PRIVATE -Wall -Wpedantic -Wextra -Werror=return-type -march=native $<$<CONFIG:DEBUG>:-g>)

find_package(FFTW COMPONENTS FLOAT_LIB DOUBLE_LIB)

if(FFTW_FOUND)
	message("Using FFTW")
	target_include_directories(vasp_ldos PUBLIC FFTW_INCLUDE_DIRS)
	target_link_directories(vasp_ldos PUBLIC FFTW_LIBRARIES)
	target_link_libraries(vasp_ldos PUBLIC m fftw3 fftw3f)
else()
	message("Using Intel MKL")
	target_include_directories(vasp_ldos PUBLIC "$ENV{MKLROOT}/include")
	target_compile_definitions(vasp_ldos PUBLIC MKL_ILP64 USE_MKL_FFT)
	target_link_directories(vasp_ldos PUBLIC "$ENV{MKLROOT}/lib/intel64")
	target_link_libraries(vasp_ldos PUBLIC mkl_intel_ilp64 mkl_sequential mkl_core m dl)
endif()
