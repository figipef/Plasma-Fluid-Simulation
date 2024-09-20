# Compiler
CXX := g++
CXXFLAGS := -std=c++14 -Wall -Wextra -Wpedantic -Wconversion -Wunused-parameter -Wunused-but-set-parameter -O2

# Directories
ROOTDIR := ..
BUILDDIR := $(ROOTDIR)/build
SRCDIR := $(ROOTDIR)/src
BINDIR := $(ROOTDIR)/build
INCDIR := $(ROOTDIR)/inc
OBJDIR := $(BUILDDIR)/obj

# Files
SOURCES := $(wildcard $(SRCDIR)/*.cpp) #$(wildcard $(INCDIR)\chi2/*.cpp)
OBJECTS := $(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o,$(SOURCES))

$(info $(SRCDIR))

# Program name
PROGRAM_NAME := test
ifdef NAME
	PROGRAM_NAME := $(NAME)
endif

# Target binary
TARGET := $(BINDIR)/$(PROGRAM_NAME)

# Default target
all: $(TARGET)

# Linking
$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $^ -o $@

# Compilation
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp | $(OBJDIR)
	$(CXX) $(CXXFLAGS) -I$(INCDIR) -c $< -o $@

$(OBJDIR)/%.o: $(INCDIR)/chi2/%.cpp | $(OBJDIR)
	$(CXX) $(CXXFLAGS) -I$(INCDIR) -c $< -o $@

# Create build directory
$(OBJDIR):
	mkdir obj

# Clean
clean:
	rm -rf $(OBJDIR) $(TARGET)

# Phony targets
.PHONY: all clean
