#This shit could totally be wrong lol, just a single C and header file dude.
cmake_minimum_required(VERSION 3.4)
project(DFC)

set(CMAKE_C_STANDARD 17)

option(BUILD_SHARED "ON: DFC is built as a static library | OFF: DFC is built as a shared library" OFF)
if(BUILD_SHARED)
	add_definitions(-D DFC_BUILD_STATIC)
else()
	add_definitions(-D DFC_BUILD_DLL)
endif()

# Add codec source directory
include_directories(${DFC_SOURCE_DIR}/Codec)
# Add all codec source files
file(GLOB DFC_SRC
	"${DFC_SOURCE_DIR}/Codec/*.h"
	"${DFC_SOURCE_DIR}/Codec/*.c"
)

if(BUILD_SHARED)
	add_library(DFC SHARED ${DFC_SRC})
else()
	add_library(DFC STATIC ${DFC_SRC})
endif()