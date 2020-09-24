ACLOCAL_AMFLAGS= -I m4
# Build the library in the hand subdirectory first.
LIB_TLDEVELDIR = @LIB_TLDEVELDIR@

@VALGRIND_CHECK_RULES@

CLEANFILES =

AM_CFLAGS = @TLDEVEL_CPPFLAGS@

AM_LDFLAGS = @DEAD_STRIP@
AM_LDFLAGS += -Wno-undef
AM_LDFLAGS += -static

LIBS = @TLDEVEL_LIB@ -lm

bin_PROGRAMS = cmsa mumsa

mumsa_SOURCES = \
global.h \
run_mumsa.c \
rwalign.c \
alphabet.h \
alphabet.c \
mumsa_data.h \
mumsa_data.c \
msa.h \
msa_ops.h \
msa_ops.c \
msa_info.h \
msa_info.c \
process_msa.h \
process_msa.c \
reporting.h \
reporting.c

cmsa_SOURCES = \
cmsa.h \
adler.h \
adler.c \
boot.h \
boot.c \
main.c \
misc.c \
print_reliable_alignment.h \
print_reliable_alignment_new.c \
score_sets.c \
upgma.c \
best.c \
input.h \
input.c \
make_sets.c \
output.h \
output.c \
r_output.c \
stat.c \
bit_array.h \
bit_array.c \
interface.c \
mem.c \
overlap.c \
rand.c \
string_matching.h \
string_matching.c

check_PROGRAMS =
TESTS =
