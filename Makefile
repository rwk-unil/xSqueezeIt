HTSLIB_PATH := ./htslib/

# C++ Compiler
CXX=g++
INCLUDE_DIRS=-I include -I $(HTSLIB_PATH)/htslib -I $(BENCHMARK_PATH)/include
#EXTRA_FLAGS=-fsanitize=address -fsanitize=undefined -fsanitize=pointer-subtract -fsanitize=pointer-compare -fno-omit-frame-pointer -fstack-protector-all -fcf-protection
CXXFLAGS=-O3 -g -Wall -std=c++17 $(INCLUDE_DIRS) $(CXXEXTRAFLAGS) $(EXTRA_FLAGS)
# Linker
LD=g++
LIBS=-lpthread -lhts -lzstd
LDFLAGS=-O3 $(EXTRA_FLAGS) -L $(HTSLIB_PATH)
# Debugger
XDB=gdb

# Project specific :
TARGET := console_app
SOURCE := console_app.cpp
OBJ := $(SOURCE:.cpp=.o)
CPP_SOURCES := $(wildcard *.cpp)
CPP_OBJS := $(CPP_SOURCES:.cpp=.o)
CPP_OBJS := $(CPP_OBJS:.c=.o)
OBJS := xcf.o bcf_traversal.o Accessor.o $(OBJ)
DEPENDENCIES := $(CPP_SOURCES:.cpp=.d)
DEPENDENCIES := $(DEPENDENCIES:.c=.d)

BENCHMARK_EXECUTABLE := main_bench
BENCHMARK_SOURCE := main_bench.cpp
BENCHMARK_OBJ := $(BENCHMARK_SOURCE:.cpp=.o)

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
%.o : %.c
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Rule to generate the dependency files
%.d : %.cpp
	$(CXX) $(INCLUDE_DIRS) -MG -MP -MM -MT '$(@:.d=.o)' $< -MF $@
%.d : %.c
	$(CXX) $(INCLUDE_DIRS) -MG -MP -MM -MT '$(@:.d=.o)' $< -MF $@

# Remove artifacts
clean :
	rm -f $(OBJS) $(TARGET) $(DEPENDENCIES)

# Rules that don't generate artifacts
.PHONY :
	all clean debug