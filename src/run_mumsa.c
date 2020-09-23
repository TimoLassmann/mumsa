#include "tldevel.h"

#include <getopt.h>
#include <string.h>
#include "msa.h"

int print_kalign_help(char * argv[])
{
        const char usage[] = " -i <seq file> -o <out aln> ";
        fprintf(stdout,"\nUsage: %s %s\n\n",basename(argv[0]) ,usage);
        fprintf(stdout,"Options:\n\n");

        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--format","Output format." ,"[Fasta]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--reformat","Reformat existing alignment." ,"[NA]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--gpo","Gap open penalty." ,"[5.5]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--gpe","Gap extension penalty." ,"[2.0]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--tgpe","Terminal gap extension penalty." ,"[1.0]"  );

        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--version (-V/-v)","Prints version." ,"[NA]"  );

        fprintf(stdout,"\nExamples:\n\n");

        fprintf(stdout,"Passing sequences via stdin:\n\n   cat input.fa | kalign -f fasta > out.afa\n\n");
        fprintf(stdout,"Combining multiple input files:\n\n   kalign seqsA.fa seqsB.fa seqsC.fa -f fasta > combined.afa\n\n");


        return OK;
}


int main(int argc, char *argv[])
{

        return EXIT_SUCCESS;
ERROR:
        return EXIT_FAILURE;
}



int check_for_sequences(struct msa* msa)

{
        if(!msa){
                ERROR_MSG("No sequences were found in the input files or standard input.");
        }
        if(msa->numseq < 2){
                if(msa->numseq == 0){
                        ERROR_MSG("No sequences were found in the input files or standard input.");
                }else if (msa->numseq == 1){
                        ERROR_MSG("Only 1 sequence was found in the input files or standard input");
                }
        }
        return OK;
ERROR:
        return FAIL;
}
