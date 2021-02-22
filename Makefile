HTSLIB_PATH := ../htslib/htslib

# C++ Compiler
CXX=g++
INCLUDE_DIRS=-I . -I $(HTSLIB_PATH)/htslib
CXXFLAGS=-O3 -g -Wall -std=c++17 $(INCLUDE_DIRS) $(CXXEXTRAFLAGS)
# Linker
LD=g++
LIBS=-lpthread -lhts
LDFLAGS=-O3 -L $(HTSLIB_PATH)
# Debugger
XDB=gdb

# Project specific :
TARGET := console_app
SOURCES := console_app.cpp
OBJS := $(SOURCES:.cpp=.o)
CPP_SOURCES := $(wildcard *.cpp)
DEPENDENCIES := $(CPP_SOURCES:.cpp=.d)

# Rules
all : $(TARGET) $(DEPENDENCIES)

# Make and open in debugger
# TODO : Handle the debug flags !
debug : $(TARGET)
	$(XDB) ./$(TARGET)

# Link the target
$(TARGET) : $(OBJS)
	$(LD) $(LDFLAGS) $^ $(LIBS) -o $@

# Do not include the depency rules for "clean"
ifneq ($(MAKECMDGOALS),clean)
-include $(DEPENDENCIES)
endif

# Compile
%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Rule to generate the dependency files
%.d : %.cpp
	$(CXX) $(INCLUDE_DIRS) -MG -MP -MM -MT '$(@:.d=.o)' $< -MF $@

# Remove artifacts
clean :
	rm -f $(OBJS) $(TARGET) $(DEPENDENCIES)

# Rules that don't generate artifacts
.PHONY :
	all clean debug