cmake_minimum_required(VERSION 2.8.5)
cmake_policy(VERSION 2.8.5)

set(LOCAL_PROJECT_NAME InteractiveRegistration)
set (Tutorial_VERSION_MAJOR 0)
set (Tutorial_VERSION_MINOR 1)

project(${LOCAL_PROJECT_NAME})

# Find ITK.
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

add_subdirectory(PolyAffineTransform)

