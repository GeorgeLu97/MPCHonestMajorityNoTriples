cmake_minimum_required(VERSION 2.8)
project(MPCHonestMajorityNoTriples)


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O3 -std=c++11 -maes -mpclmul -mbmi2 -pthread")

INCLUDE_DIRECTORIES($ENV{HOME} $ENV{HOME}/libscapi/install/include
        $ENV{HOME}/libscapi/install/include/OTExtensionBristol)

set(SOURCE_FILES main.cpp
        ProtocolParty.cpp ProtocolParty.h ProtocolTimer.cpp ProtocolTimer.h
        BAParty.cpp BAParty.h ECC.cpp ECC.h LinearParty.hpp LinearParty.cpp
        ZpKaratsubaElement.cpp ZpKaratsubaElement.h
        HashEncrypt.cpp HashEncrypt.h HashAbstract.h TemplateField.cpp)


link_directories(/usr/ssl/lib/ $ENV{HOME}/libscapi/install/lib )


add_executable(MPCHonestMajorityNoTriples ${SOURCE_FILES})

# /home/hanjun/libscapi/build/NTL/src
# TARGET_LINK_LIBRARIES(MPCHonestMajorityNoTriples $ENV{HOME}/ntl-9.10.0/src/ntl.a
TARGET_LINK_LIBRARIES(MPCHonestMajorityNoTriples $ENV{HOME}/libscapi/install/lib/libntl.a
        $ENV{HOME}/libscapi/libscapi.a gmp gmpxx
        libboost_system.a libboost_thread.a pthread crypto dl ssl z)
