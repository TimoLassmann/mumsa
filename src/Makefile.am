# Build the library in the hand subdirectory first.

ACLOCAL_AMFLAGS = -I m4

@VALGRIND_CHECK_RULES@

AM_CFLAGS =  @TLDEVEL_CFLAGS@

CLEANFILES =

# top_builddir = @top_builddir@
AM_LDFLAGS = @DEAD_STRIP@
AM_LDFLAGS += -Wno-undef

bin_PROGRAMS = cmsa

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
