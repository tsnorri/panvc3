GCC_ROOT		= CONDA_PREFIX

BOOST_INCLUDE	=

CPPFLAGS		=	--sysroot CONDA_PREFIX/x86_64-conda-linux-gnu/sysroot \
					-isystem =/../../include \
					-DLIBBIO_NO_DISPATCH

LDFLAGS			=	-L CONDA_PREFIX/lib -pthread -lz -lbz2 -lm -lc -ldl

BOOST_LIBS		=	-L CONDA_PREFIX/lib -lboost_iostreams
