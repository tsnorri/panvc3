GCC_ROOT = CONDA_PREFIX
CLANG_ROOT = CONDA_PREFIX

BOOST_INCLUDE	=

GCC_CPPFLAGS	=	--sysroot CONDA_PREFIX/x86_64-conda-linux-gnu/sysroot \
					-isystem =/../../include \
					-DLIBBIO_NO_DISPATCH -DLIBBIO_NO_SAM_READER

LLVM_CPPFLAGS	=	--sysroot CONDA_PREFIX/x86_64-conda-linux-gnu/sysroot \
					-isystem CONDA_PREFIX/include \
					-DLIBBIO_NO_SAM_READER

SYSTEM_CXXFLAGS	=
LLVM_CFLAGS		= -fblocks
LLVM_CXXFLAGS	= -fblocks

GCC_LDFLAGS		= -L CONDA_PREFIX/lib -pthread -lz -lbz2 -lm -lc -ldl
LLVM_LDFLAGS	= -L CONDA_PREFIX/lib -pthread -lz -lbz2 -lm -lc -ldl -fblocks

LIBDISPATCH_CFLAGS      = 
LIBDISPATCH_CXXFLAGS    = 
LIBDISPATCH_LDFLAGS     = 

BOOST_LIBS = -L CONDA_PREFIX/lib -lboost_iostreams
