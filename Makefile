HTSLIB_PATH := ./htslib
ZSTD_PATH := ./zstd/lib

# C++ Compiler
CXX=g++
INCLUDE_DIRS=-I include -I $(HTSLIB_PATH)/htslib -I $(ZSTD_PATH)
#EXTRA_FLAGS=-fsanitize=address -fsanitize=undefined -fsanitize=pointer-subtract -fsanitize=pointer-compare -fno-omit-frame-pointer -fstack-protector-all -fcf-protection
CXXFLAGS=-O3 -g -Wall -std=c++17 $(INCLUDE_DIRS) $(CXXEXTRAFLAGS) $(EXTRA_FLAGS)
# Linker
LD=g++
LIBS=-lpthread -lhts -lzstd
LDFLAGS=-O3 $(EXTRA_FLAGS) -L $(HTSLIB_PATH) -L $(ZSTD_PATH)

# Project specific :
TARGET := xsqueezeit
SOURCE := xsqueezeit.cpp
OBJ := $(SOURCE:.cpp=.o)
CPP_SOURCES := $(wildcard *.cpp)
CPP_OBJS := $(CPP_SOURCES:.cpp=.o)
CPP_OBJS := $(CPP_OBJS:.c=.o)
OBJS := xcf.o bcf_traversal.o accessor.o $(OBJ)
DEPENDENCIES := $(CPP_SOURCES:.cpp=.d)
DEPENDENCIES := $(DEPENDENCIES:.c=.d)

# Rules
all : $(TARGET) $(DEPENDENCIES)

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