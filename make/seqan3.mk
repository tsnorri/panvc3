CPPFLAGS	+= -I../lib/seqan3/include -DSEQAN3_HAS_ZLIB=1 -DSEQAN3_HAS_BZIP2=1
LDFLAGS		+= -lz -lbz2
