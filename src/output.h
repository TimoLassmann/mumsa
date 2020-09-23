#include <ctype.h>

FILE* fout = 0;
static void msf_output(struct alignment* aln,char* outfile);
static void fasta_output(struct alignment* aln,char* outfile);
static void clustal_output(struct alignment* aln,char* outfile);
static void macsim_output(struct alignment* aln,char* outfile,char* infile);

