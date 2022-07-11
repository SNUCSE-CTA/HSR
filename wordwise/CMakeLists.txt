cmake_minimum_required(VERSION 2.8)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Header path of installed HEaan.STAT lib
include_directories(/opt/cryptolab/heaan/include/)
# Path of installed HEaan.STAT lib
link_directories(/opt/cryptolab/heaan/lib/)

set (CMAKE_CXX_STANDARD 14)

add_compile_options(-Ofast)

add_compile_options(-pthread)

ADD_EXECUTABLE( main
                main.cpp enc.cpp eval.cpp)

TARGET_LINK_LIBRARIES(main
                       HEAAN HEAAN_UTILS ntl gmp pthread crypto stdc++fs)