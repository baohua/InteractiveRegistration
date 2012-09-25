cmake_minimum_required(VERSION 2.8.5)
cmake_policy(VERSION 2.8.5)

project(itkPolyAffineTransformTest)
include_directories ("${CMAKE_SOURCE_DIR}/ThirdParty")

# add the executable
add_executable(itkPolyAffineTransformTest itkPolyAffineTransformTest.cxx
               itkPolyAffineTransform.h itkPolyAffineTransform.hxx
               itkDebugHelper.h itkDebugHelper.hxx)
target_link_libraries(itkPolyAffineTransformTest ${ITK_LIBRARIES})

# add the executable
add_executable(itkLocalAffineTransformTest itkLocalAffineTransformTest.cxx
               itkLocalAffineTransform.h itkLocalAffineTransform.hxx)
target_link_libraries(itkLocalAffineTransformTest ${ITK_LIBRARIES})

# add package
set(CPACK_PACKAGE_NAME "InteractiveRegistration")
set(CPACK_PACKAGE_VENDOR "CMake.org")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "InteractiveRegistration")
set(CPACK_PACKAGE_VERSION "0.1.x")
set(CPACK_PACKAGE_VERSION_MAJOR "0")
set(CPACK_PACKAGE_VERSION_MINOR "1")
set(CPACK_PACKAGE_VERSION_PATCH "x")
set(CPACK_PACKAGE_INSTALL_DIRECTORY "InteractiveRegistration")
include(CPack)

