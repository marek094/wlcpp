# List all your source files
set(SOURCE_FILES main.cpp genevalformulas.cpp hom2logic.cpp freetypes.cpp justpaths.cpp tables.cpp complex_approx.cpp)

# Find the required Boost components
find_package(Boost REQUIRED COMPONENTS serialization)

# Loop over the source files
foreach(source_file ${SOURCE_FILES})
    # Get the filename without the extension
    get_filename_component(executable_name ${source_file} NAME_WE)

    # Add an executable for this source file
    add_executable(${executable_name} ${source_file})

    # Link the Boost libraries to the executable
    target_link_libraries(${executable_name} Boost::boost Boost::serialization)  # Replace Boost::boost with the required Boost components
endforeach()