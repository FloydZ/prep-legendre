cmake_minimum_required(VERSION 3.20)
project(LegendrePrecompute)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_CXX_STANDARD 20)

# possible sanitzers = -fsanitize=[address,leak,thread,memory,undefined]
set(SANITIZE_FLAG "")
#set(SANITIZE_FLAG "-fsanitize=undefined")

set(DEBUG_FLAGS "-g -DDEBUG -DUSE_AVX2 -march=native ${SANITIZE_FLAG}")
set(OPTIMIZE_FLAGS "-O3 -march=native -flto -DUSE_AVX2 ${SANITIZE_FLAG}")
set(SANITIZER_FLAGS -fsanitize=address -fno-omit-frame-pointer)
set(GCC_HELPER_FLAGS -fvar-tracking-assignments)
set(LINK_FLAGS -pthread gmp gmpxx ${SANITIZE_FLAG})

set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} ${DEBUG_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${DEBUG_FLAGS}")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} ${OPTIMIZE_FLAGS}")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${OPTIMIZE_FLAGS}")

add_executable(code main.cpp)
target_link_libraries(code ${LINK_FLAGS})

add_executable(code_inner main_inner.cpp)
target_link_libraries(code_inner ${LINK_FLAGS})

add_executable(code_noprep main_noprep.cpp)
target_link_libraries(code_noprep ${LINK_FLAGS})
