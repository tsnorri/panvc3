GCC_ROOT = CONDA_PREFIX
CLANG_ROOT = CONDA_PREFIX

BOOST_INCLUDE	=

GCC_CPPFLAGS	=	-nostdinc \
					-isystem CONDA_PREFIX/lib/gcc/x86_64-conda-linux-gnu/12.3.0/include \
					-isystem CONDA_PREFIX/lib/gcc/x86_64-conda-linux-gnu/12.3.0/include-fixed \
					-isystem CONDA_PREFIX/x86_64-conda-linux-gnu/include/c++/12.3.0 \
					-isystem CONDA_PREFIX/x86_64-conda-linux-gnu/include/c++/12.3.0/x86_64-conda-linux-gnu \
					-isystem CONDA_PREFIX/x86_64-conda-linux-gnu/sysroot/usr/include \
					-isystem CONDA_PREFIX/x86_64-conda-linux-gnu/sysroot/usr/include/linux \
					-isystem CONDA_PREFIX/include \
					-DBOOST_STACKTRACE_USE_NOOP -DLIBBIO_NO_DISPATCH -DLIBBIO_NO_SAM_READER
LLVM_CPPFLAGS	=	-nostdinc \
					-isystem CONDA_PREFIX/lib/clang/17/include \
					-isystem CONDA_PREFIX/x86_64-conda-linux-gnu/include/c++/12.3.0 \
					-isystem CONDA_PREFIX/x86_64-conda-linux-gnu/include/c++/12.3.0/x86_64-conda-linux-gnu \
					-isystem CONDA_PREFIX/x86_64-conda-linux-gnu/sysroot/usr/include \
					-isystem CONDA_PREFIX/x86_64-conda-linux-gnu/sysroot/usr/include/linux \
					-isystem CONDA_PREFIX/include \
					-DLIBBIO_NO_SAM_READER

SYSTEM_CXXFLAGS	=
LLVM_CFLAGS		= -fblocks
LLVM_CXXFLAGS	= -fblocks -nostdinc++

GCC_LDFLAGS		= -L CONDA_PREFIX/lib -pthread -lz -lbz2 -lm -lc -ldl
LLVM_LDFLAGS	= -L CONDA_PREFIX/lib -pthread -lz -lbz2 -lm -lc -ldl -fblocks

LIBDISPATCH_CFLAGS      = 
LIBDISPATCH_CXXFLAGS    = 
LIBDISPATCH_LDFLAGS     = 

BOOST_LIBS = -L CONDA_PREFIX/lib -lboost_iostreams
