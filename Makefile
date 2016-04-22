# *************************************************************************** #
#                      MAKEFILE FOR BASIC C++ PROGRAMS                        #
# *************************************************************************** #

# ---------------- Compiler Options;
CXX = g++
RM  = rm -f
LDFLAGS =

DEBUG = 0

executable = SE276C_HW3
includes = -I /usr/include/eigen3/
cxxStd = -std=c++11
optLevel = -O3
object_dir = ./obj
executableFull = $(executable)

# ----------------
ifeq ($(DEBUG), 1)
	optLevel = -O0
	CXXFLAGS += -g -DDEBUG
	LDFLAGS += -g -DDEBUG
endif

CXXFLAGS += $(includes) $(optLevel) $(cxxStd)

# object_dir = ./obj
SRCS = $(wildcard *.cpp)
OBJS = $(patsubst %.cpp,$(object_dir)/%.o,$(SRCS))
DEPS = $(patsubst %.cpp,$(object_dir)/%.d,$(SRCS))

all: directories $(executableFull)

.PHONY: directories
directories: $(object_dir)

$(object_dir):
	mkdir -p $(object_dir)


# ---------------- Main Target;
$(executableFull) : $(OBJS)
	$(CXX) $(LDFLAGS) -o $@ $(OBJS)

# ---------------- Automatic dependencies;

df = $(object_dir)/$(*F)

-include $(DEPS)

$(object_dir)/%.o: %.cpp
	$(CXX) $(CXXFLAGS) -MM -MP -MT $(df).o -MT $(df).d $< > $(df).d
	$(CXX) -c $< $(CXXFLAGS) -o $@

# ---------------- Auxiliary tools;
clean:
	$(RM) $(OBJS)

dist-clean: clean
	$(RM) $(DEPS)

