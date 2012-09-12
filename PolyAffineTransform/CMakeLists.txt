cmake_minimum_required(VERSION 2.8.5)
cmake_policy(VERSION 2.8.5)

project(antsPolyAffineTransformTest)
 
# add the executable
add_executable(antsPolyAffineTransformTest antsPolyAffineTransformTest.cxx)

target_link_libraries(antsPolyAffineTransformTest ${ITK_LIBRARIES})

