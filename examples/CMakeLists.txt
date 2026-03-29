add_executable(example example.cpp)
target_include_directories(example PRIVATE ${HOMLIB_INCLUDE_DIRS})
target_link_libraries(example LINK_PUBLIC HomLib Eigen3::Eigen PoseLib::PoseLib)

add_executable(example_ransac example_ransac.cpp)
target_include_directories(example_ransac PRIVATE ${HOMLIB_INCLUDE_DIRS})
target_include_directories(example_ransac PRIVATE ../thirdparty/RansacLib)
target_link_libraries(example_ransac LINK_PUBLIC HomLib Eigen3::Eigen PoseLib::PoseLib)