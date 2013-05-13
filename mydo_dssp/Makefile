# Generated automatically from Makefile.in by configure.
#
# This is a Gromacs 3.0 template makefile for your own utility programs.
#
# Copy this file to whatever directory you are using for your own
# software and add more targets like the template one below.
#
# If you are using gmake it is relatively straightforward to add
# an include based on environment variables (like previous Gromacs versions)
# to select compiler flags and stuff automatically, but below it is static:
#

TARGET = mydo_dssp
GROMACS_DIR = ${HOME}/exec/gromacs-4.0.5/exec
BIN = ../bin

# Variables set by the configuration script:
LIBS      = -lmd -lgmx -lnsl -lfftw3f   
LDFLAGS   = -L${HOME}/exec/gromacs-4.0.5/exec/lib -L${HOME}/exec/fftw-3.2.2/exec/lib  
CFLAGS	  = -O3 -tpp7 -axW -ip -w -funroll-all-loops -std=gnu99  -I${HOME}/exec/gromacs-4.0.5/exec/include -I${HOME}/exec/gromacs-4.0.5/exec/include/gromacs	

# LDFLAGS      = -L/home/zhuyxue12/exec/gromacs-4.0.5/exec/lib -L/software/libraries/FFTW-3.3/mvapich2-intel/lib  
# CFLAGS	     = -O3 -fomit-frame-pointer -finline-functions -Wall -Wno-unused -funroll-all-loops -std=gnu99  -I/home/zhuyxue12/exec/gromacs-4.0.5/exec/include -I/home/zhuyxue12/exec/gromacs-4.0.5/exec/include/gromacs	
CC           = icc 
LD           = $(CC)

# The real make targets - note that most make programs support
# the shortcut $^ instead of listing all object files a second
# time, but we cannot count on it...

$(TARGET):	$(TARGET).o
		$(LD) $(LDFLAGS) $(TARGET).o $(LIBS) -o $(BIN)/$@

.PHONY: clean install
install:
	cp $(BIN}/$(TARGET) $(GROMACS_DIR)/bin
clean:
	rm -fv $(TARGET)