OPT_FLAGS			?= -O2 -ggdb

CPPFLAGS			+= -D_LIBCPP_DISABLE_AVAILABILITY -DBOOST_STACKTRACE_USE_NOOP -DLIBBIO_NO_DISPATCH -I/usr/include/kqueue
#CFLAGS				+= -fsanitize=thread -fsanitize=undefined
#CXXFLAGS			+= -fsanitize=thread -fsanitize=undefined
#LDFLAGS			+= -fsanitize=thread -fsanitize=undefined
WARNING_CXXFLAGS_	= -Wno-interference-size

#BOOST_LIBS			?= $(BOOST_LIBDIR)/libboost_iostreams.a
BOOST_LIBS			?= -lboost_iostreams
