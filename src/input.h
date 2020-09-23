#include <unistd.h>
#include <string.h>
#include <ctype.h>

#define SEEK_START 0
#define SEEK_END 2

static struct alignment* read_alignment(struct alignment* aln,char* string);
static struct alignment* read_alignment_from_swissprot(struct alignment* aln,char* string);
static struct alignment* read_alignment_uniprot_xml(struct alignment* aln,char* string);
static struct alignment* read_alignment_macsim_xml(struct alignment* aln,char* string);
static struct feature* read_ft(struct feature* ft,char* p);
static struct alignment* read_alignment_clustal(struct alignment* aln,char* string);
static struct alignment* read_alignment_stockholm(struct alignment* aln,char* string);

static char* get_input_into_string(char* string,char* infile);



static int count_sequences_macsim(char* string);
static int count_sequences_swissprot(char* string);
static int count_sequences_uniprot(char* string);
static int count_sequences_stockholm(char* string);
static int count_sequences_clustalw(char* string);
static int count_sequences_fasta(char* string);
