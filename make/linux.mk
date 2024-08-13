OPT_FLAGS			?= -O2 -ggdb

CPPFLAGS			+= -D_LIBCPP_DISABLE_AVAILABILITY -DBOOST_STACKTRACE_USE_BACKTRACE -DBOOST_STACKTRACE_BACKTRACE_INCLUDE_FILE="</usr/lib/gcc/x86_64-linux-gnu/12/include/backtrace.h>" -DLIBBIO_NO_DISPATCH -I/usr/include/kqueue
#CFLAGS				+= -fsanitize=thread -fsanitize=undefined
#CXXFLAGS			+= -fsanitize=thread -fsanitize=undefined
#LDFLAGS			+= -fsanitize=thread -fsanitize=undefined
WARNING_CXXFLAGS_	= -Wno-interference-size

#BOOST_LIBS			?= $(BOOST_LIBDIR)/libboost_iostreams.a
BOOST_LIBS			?= -lboost_iostreams -lboost_stacktrace_backtrace
