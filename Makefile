BENCHMARK_PATH := ../benchmark/benchmark
HTSLIB_PATH := ../htslib/htslib

# C++ Compiler
CXX=g++
INCLUDE_DIRS=-I . -I $(HTSLIB_PATH)/htslib -I $(BENCHMARK_PATH)/include -I $(BENCHMARK_PATH)/googletest/googletest/include
#INCLUDE_DIRS=-I . -I benchmark/include -I benchmark/googletest/googletest/include
#CXXFLAGS=-O3 -Wall -std=c++17 $(INCLUDE_DIRS) $(CXXEXTRAFLAGS)
CXXFLAGS=-O3 -g -Wall -std=c++17 $(INCLUDE_DIRS) $(CXXEXTRAFLAGS)
# Linker
LD=g++
LIBS=-lpthread -lhts
#LDFLAGS=-L benchmark/build/src/
LDFLAGS=-O3 -L $(HTSLIB_PATH)
# Debugger
XDB=lldb

# Project specific :
TARGET := console_app
SOURCES := console_app.cpp #io/genotype_reader.cpp #$(wildcard *.cpp) $(wildcard objects/*.cpp) $(wildcard containers/*.cpp)
OBJS := $(SOURCES:.cpp=.o)
CPP_SOURCES := $(wildcard *.cpp)
DEPENDENCIES := $($CPP_SOURCES:.cpp=.d)

BENCHMARK_EXECUTABLE := pbwt_bench
BENCHMARK_SOURCES := main_bench.cpp
BENCHMARK_OBJS := $(BENCHMARK_SOURCES:.cpp=.o)

TEST_EXECUTABLE := pbwt_test
TEST_SOURCES := main_test.cpp $(BENCHMARK_PATH)/googletest/googletest/src/gtest_main.cc
TEST_OBJS := $(TEST_SOURCES:.cpp=.o)
TEST_OBJS := $(TEST_OBJS:.cc=.o)

MAIN_XCF := main_xcf
MAIN_XCF_SOURCES := main_xcf.cpp
MAIN_XCF_OBJS := $(MAIN_XCF_SOURCES:.cpp=.o)

# Rules
all : $(TARGET)

# Make and run
run : $(TARGET)
	./$(TARGET)

# Make and run with time
runtime : $(TARGET)
	time ./$(TARGET)

# Make and open in debugger
# TODO : Handle the debug flags !
debug : $(TARGET)
	$(XDB) ./$(TARGET)

benchmark : $(BENCHMARK_EXECUTABLE)
	./$(BENCHMARK_EXECUTABLE)

$(BENCHMARK_EXECUTABLE) : $(BENCHMARK_OBJS)
	$(LD) $(LDFLAGS) -L $(BENCHMARK_PATH)/build/src/ $^ $(LIBS) -lbenchmark -o $@

test : $(TEST_EXECUTABLE)
	./$(TEST_EXECUTABLE)

test_debug : $(TEST_EXECUTABLE)
	$(XDB) ./$(TEST_EXECUTABLE)

$(TEST_EXECUTABLE) : $(TEST_OBJS)
	$(LD) $(LDFLAGS) -L $(BENCHMARK_PATH)/build/lib/ $^ $(LIBS) -lgtest -o $@

main_xcf.cpp : xcf.hpp pbwt_exp.hpp pbwt_big.hpp
	touch $@

$(MAIN_XCF) : $(MAIN_XCF_OBJS)
	$(LD) $(LDFLAGS) $^ $(LIBS) -o $@

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
	rm -f $(OBJS) $(TARGET) $(DEPENDENCIES) $(TEST_OBJS) $(TEST_EXECUTABLE) $(BENCHMARK_EXECUTABLE) $(BENCHMARK_OBJS)

# Rules that don't generate artifacts
.PHONY :
	all clean run runtime debug test benchmark