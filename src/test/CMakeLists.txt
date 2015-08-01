add_executable(test_main test_main.cpp ../kdtree.c)
pods_use_pkg_config_packages(test_main ${POD_NAME})

add_executable(test_dubins test_dubins.cpp)
pods_use_pkg_config_packages(test_dubins ${POD_NAME})

