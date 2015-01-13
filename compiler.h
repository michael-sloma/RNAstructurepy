###############################################################################
## Many different C/C++ compilers can be used, defined as the following:
## Option 1: g++ (GNU C++ compiler)
## Option 2: CC (SGI and Sun compiler)
## Option 3: icc (Intel C compiler)
## Option 4: xlC (xlc compiler, e.g. on AIX)
###############################################################################

CXX = g++
#CXX = CC
#CXX = icc
#CXX = xlC

###############################################################################
## The Java compiler, defined below, must be active if any RNAstructure executable with a Java component is being built.
###############################################################################

JXX = javac 
JXXFLAGS = -source 1.6

###############################################################################
## For compilation of the C++ component of RNAstructure, the exact flags depend on the operating system.
## There are flags for compiling into object files and for linking into a shared library, defined below for the three major systems.
###############################################################################

# Linux
CXXFLAGS = -O3 -fsched-spec-load -fPIC -D NDEBUG
MULTIFINDFLAG = -DMULTIFIND
LIBFLAGS = -shared
OPSYSTEM = Linux

# Mac
#CXXFLAGS = -O3 -fPIC -arch x86_64
#LIBFLAGS = -dynamiclib -arch x86_64 -framework JavaVM
#OPSYSTEM = Mac

# Cygwin (running on Windows)
#CXXFLAGS = -O3 -fsched-spec-load -mno-cygwin
#LIBFLAGS = -Wl,--add-stdcall-alias -mno-cygwin
#OPSYSTEM = Windows

###############################################################################
## For compilation of SMP programs using OpenMP, the OpenMP flags need
##	to be set. ###############################################################################

#g++ compiler:
CXXOPENMPFLAGS = -fopenmp -D SMP ${CXXFLAGS} 

#intel compiler:
#CXXOPENMPFLAGS = -openmp -D SMP ${CXXFLAGS} 

###############################################################################
## The INCLUDEPATH flags are required to find proper linking headers for RNAstructure executables with Java components.
## On all systems, these paths are the required native headers needed to interface Java with C++ code.
## These flags must be changed manually to coincide with the current machine before the Java components of RNAstructure can run successfully.
##
##
## INCLUDEPATH1 must hold the directory which contains "jni.h".
## INCLUDEPATH2 must hold the directory which contains "jni_md.h".
## Absolute paths must be used for both.
##
## For Fedora Linux, these two paths are normally similar to:
##     INCLUDEPATH1 = /usr/java/jdk1.x.x_x/include
##     INCLUDEPATH2 = /usr/java/jdk1.x.x_x/include/linux
## For Ubuntu Linux, these two paths are normally similar to:
##     INCLUDEPATH1 = /usr/lib/jvm/java-6-sun-1.x.x.x/include
##     INCLUDEPATH2 = /usr/lib/jvm/java-6-sun-1.x.x.x/include/linux
## For Mac, these two paths are usually the same, and similar to:
##     INCLUDEPATH1 = /System/Library/Frameworks/JavaVM.framework/Headers
##     For Mac, INCLUDEPATH2 should be defined as "", since it's not used.
## For Windows (Cygwin build), these two paths are normally similar to:
##     INCLUDEPATH1 = "C:\Program Files\Java\jdk1.x.x_x\include"
##     INCLUDEPATH2 = "C:\Program Files\Java\jdk1.x.x_x\include\win32"
##     For Windows, quotes may be needed if paths contain spaces.
##
## Note that the example paths given above are only meant to be a general guide to the paths' forms on different systems.
## The actual paths on any given system may differ from these.
###############################################################################

#INCLUDEPATH1 = "/usr/java/jdk1.6.0_13/include"
#INCLUDEPATH2 = "/usr/java/jdk1.6.0_13/include/linux"

##################################################################
##If making Multifind, INCLUDESVM must indicate the location of svm.o,
##  from the libsvm package. 
####################

  INCLUDESVM = /usr/local/libsvm-3.17
  
###############################################################################
## The rest of this file is taken up by common compilation and linking commands.
## Nothing below this point needs to be, or should be, changed.
## It may, however, be useful to look at the commands below for reference.
##
## Note that for some of these commands to work, the ROOTPATH variable, which
## should hold the path to the RNAstructure root directory, must be defined.
##
## As a general rule of thumb, ROOTPATH should always be defined in a Makefile
## before compiler.h is included. ROOTPATH may be defined as an absolute path
## or relative to the current Makefile. All the Makefiles contained within
## RNAstructure follow this rule, and define ROOTPATH with a relative path.
###############################################################################

# Define operating system specific compiling and linking rules for C++ files that proxy with Java.
# Note that the linking rules don't pull in the target name as the executable name; that must be specified immediately after the linking rule, without spaces.
ifeq (${OPSYSTEM},Linux)
	COMPILE_PROXY_FOR_JAVA = ${CXX} ${CXXFLAGS} -I${INCLUDEPATH1} -I${INCLUDEPATH2} -O -D_JAVA_GUI -c -o $@ $<
	LINK_PROXY_FOR_JAVA = ${CXX} ${LIBFLAGS} -o ${ROOTPATH}/exe/
else ifeq (${OPSYSTEM},Mac)
	COMPILE_PROXY_FOR_JAVA = ${CXX} ${CXXFLAGS} -I${INCLUDEPATH1} -O -D_JAVA_GUI -c -o $@ $<
	LINK_PROXY_FOR_JAVA = ${CXX} ${LIBFLAGS} -o ${ROOTPATH}/exe/
else
	COMPILE_PROXY_FOR_JAVA = ${CXX} ${CXXFLAGS} -I${INCLUDEPATH1} -I${INCLUDEPATH2} -O -D_JAVA_GUI -c -o $@ $<
	LINK_PROXY_FOR_JAVA = ${CXX} ${LIBFLAGS} -shared -I${INCLUDEPATH1} -I${INCLUDEPATH2} -o ${ROOTPATH}/exe/
endif

# Define suffixes for compilation.
.SUFFIXES:      .cpp .cxx
.cpp.o:
	${CXX} ${CXXFLAGS} -c -o $@ $<
.cxx.o:
	${COMPILE_PROXY_FOR_JAVA}

# Define a Dynalign II compiling rule.
COMPILE_DYNALIGN_II = ${CXX} ${CXXFLAGS} -DDYNALIGN_II -c -o $@

# Define a Dynalign II SMP compiling rule.
COMPILE_DYNALIGN_II_SMP = ${CXX} ${CXXFLAGS} -DDYNALIGN_II -DCOMPILE_SMP -c -o $@

# Define an MULTIFIND compiling rule.
COMPILE_MULTIFIND = ${CXX} ${CXXFLAGS} -DMULTIFIND -c -o $@

# Define an SMP compiling rule.
COMPILE_SMP = ${CXX} ${CXXFLAGS} -DCOMPILE_SMP -c -o $@

# Define an SVM compiling rule.
COMPILE_SVM = ${CXX} ${CXXFLAGS} -I${INCLUDESVM} -c -o $@

# Define an SVM/SMP compiling rule.
COMPILE_SVM_SMP = ${CXX} ${CXXFLAGS} -I${INCLUDESVM} -DCOMPILE_SMP -c -o $@

# Define an INSTRUMENTED compiling rule.
COMPILE_INSTRUMENTED = ${CXX} ${CXXFLAGS} -DINSTRUMENTED -c -o $@

# Define a general linking rule.
LINK = ${CXX} ${ARCHITECTURE} ${CXXFLAGS} -o $@

# Define the SMP linking rule, for use with OpenMP programs.
LINKSMP = ${CXX} ${ARCHITECTURE} ${CXXOPENMPFLAGS} -o $@

# Define a secondary linking rule that creates an executable in the exe directory
# by appending the name of the target on the exe directory. For this rule to work,
# the name of the target must not include directory paths.
LINK_IN_EXE_DIR = ${CXX} ${ARCHITECTURE} ${CXXFLAGS} -o ${ROOTPATH}/exe/$@

