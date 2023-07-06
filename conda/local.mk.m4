GCC_ROOT = CONDA_PREFIX
CLANG_ROOT = CONDA_PREFIX

BOOST_INCLUDE	=

GCC_CPPFLAGS	=	-nostdinc \
					-isystem CONDA_PREFIX/lib/gcc/x86_64-conda-linux-gnu/12.3.0/include \
					-isystem CONDA_PREFIX/lib/gcc/x86_64-conda-linux-gnu/12.3.0/include-fixed \
					-isystem CONDA_PREFIX/x86_64-conda-linux-gnu/include/c++/12.3.0 \
					-isystem CONDA_PREFIX/x86_64-conda-linux-gnu/include/c++/12.3.0/x86_64-conda-linux-gnu \
					-isystem CONDA_PREFIX/x86_64-conda-linux-gnu/sysroot/usr/include \
					-isystem CONDA_PREFIX/include \
					-DBOOST_STACKTRACE_USE_NOOP -DLIBBIO_NO_DISPATCH
LLVM_CPPFLAGS	=	-nostdinc \
					-isystem CONDA_PREFIX/lib/clang/16/include \
					-isystem CONDA_PREFIX/x86_64-conda-linux-gnu/include/c++/12.3.0 \
					-isystem CONDA_PREFIX/x86_64-conda-linux-gnu/include/c++/12.3.0/x86_64-conda-linux-gnu \
					-isystem CONDA_PREFIX/x86_64-conda-linux-gnu/sysroot/usr/include \
					-isystem CONDA_PREFIX/include

SYSTEM_CXXFLAGS	=
LLVM_CFLAGS		= -fblocks
LLVM_CXXFLAGS	= -fblocks -nostdinc++

GCC_LDFLAGS		= -pthread -lz -lbz2 -lm -lc -ldl
LLVM_LDFLAGS	= -fblocks -pthread -lz -lbz2 -lm -lc -ldl

LIBDISPATCH_CFLAGS      = 
LIBDISPATCH_CXXFLAGS    = 
LIBDISPATCH_LDFLAGS     = 

BOOST_LIBS = -lboost_iostreams
