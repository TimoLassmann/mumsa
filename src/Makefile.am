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
run_mumsa.c \
rwalign.c

cmsa_SOURCES = \
adler.c \
boot.c \
main.c \
misc.c \
print_reliable_alignment_new.c \
score_sets.c \
upgma.c \
best.c \
input.c \
make_sets.c \
output.c \
r_output.c \
stat.c \
bit_array.c \
interface.c \
mem.c \
overlap.c \
rand.c \
string_matching.c

check_PROGRAMS =
TESTS =
