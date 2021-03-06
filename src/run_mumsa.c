/*
  run_mumsa.c

  Released under GPL - see the 'COPYING' file

  Copyright (C) 2020 Timo Lassmann <timolassmann@gmail.com>

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

  Please send bug reports, comments etc. to:
  timolassmann@gmail.com
*/

#include "tldevel.h"
#include "tlmisc.h"

#include <getopt.h>
#include <string.h>
#include "global.h"
#include "mumsa_data.h"
#include "msa.h"
#include "msa_ops.h"
#include "msa_info.h"

#include "process_msa.h"

#include "reporting.h"

static int run_mumsa(struct parameters* param);
static int check_for_sequences(struct msa* msa);


static int print_mumsa_header(void);
static int print_mumsa_warranty(void);
static int print_mumsa_help(char * argv[]);
static void free_parameters(struct parameters* param);

#define OPT_MODE 1

int main(int argc, char *argv[])
{
        struct parameters* param = NULL;
        char* mode = NULL;
        int version = 0;
        int showw = 0;

        int c;
        MMALLOC(param, sizeof(struct parameters));
        param->help_flag = 0;
        param->num_infiles = 0;
        param->score_mode = MUMSA_SCORE_REF;
        param->infile = NULL;

        while (1){
                static struct option long_options[] ={
                        {"mode",  required_argument, 0, OPT_MODE},
                        {"help",   no_argument,0,'h'},
                        {"version",   no_argument,0,'v'},
                        {0, 0, 0, 0}
                };

                int option_index = 0;

                c = getopt_long_only (argc, argv,"hvV",long_options, &option_index);

                /* Detect the end of the options. */
                if (c == -1){
                        break;
                }
                switch(c) {
                case OPT_MODE:
                        mode = optarg;
                        break;
                case 'h':
                        param->help_flag = 1;
                        break;
                case 'v':
                case 'V':
                        version = 1;
                        break;
                case '?':
                        free_parameters(param);
                        return EXIT_FAILURE;
                        break;
                default:
                        abort ();
                }
        }

        if(version){
                fprintf(stdout,"%s %s\n",PACKAGE_NAME, PACKAGE_VERSION);
                free_parameters(param);
                return EXIT_SUCCESS;
        }
        print_mumsa_header();
        if(showw){
                print_mumsa_warranty();
                free_parameters(param);
                return EXIT_SUCCESS;
        }
        if(param->help_flag){
                RUN(print_mumsa_help(argv));
                free_parameters(param);
                return EXIT_SUCCESS;
        }
        if(mode != NULL){
                if(!strncmp(mode, "ref",16)){
                        param->score_mode = MUMSA_SCORE_REF;
                }else if(!strncmp(mode, "test",16)){
                        param->score_mode = MUMSA_SCORE_TEST;
                }else if(!strncmp(mode, "both",16)){
                        param->score_mode = MUMSA_SCORE_REF_TEST;
                }else{
                        ERROR_MSG("Mode: \n\n\"%s\"\n\n is not a valid option. Choose between ref,test or both.",mode);
                }

        }
        /* for good measure */
        param->num_infiles = 0;
        if (optind < argc){
                param->num_infiles += argc-optind;
                MMALLOC(param->infile, sizeof(char*) * param->num_infiles);
                c = 0;
                while (optind < argc){
                        param->infile[c] =  argv[optind++];
                        c++;
                }
        }
        if(param->num_infiles ==0){
                RUN(print_mumsa_help(argv));
                LOG_MSG("No input files");
                free_parameters(param);
                return EXIT_SUCCESS;
        }else{/* check if files exist  */
                for(c = 0; c < param->num_infiles;c++){
                        if(!my_file_exists(param->infile[c])){
                                ERROR_MSG("INPUT file: %s does not exist!", param->infile[c]);
                        }
                }
        }

        RUN(run_mumsa(param));

        free_parameters(param);
        return EXIT_SUCCESS;
ERROR:
        free_parameters(param);
        return EXIT_FAILURE;
}

int run_mumsa(struct parameters* param)
{

        struct mumsa_data* m_data = NULL;
        struct msa** msa = NULL;
        struct msa_info** msai = NULL;
        int i;
        int j;
        int c;

        RUN(alloc_mumsa_data(&m_data, param->num_infiles));

        /* Step one read in all alignments  */
        for(i = 0; i < param->num_infiles;i++){
                RUN(read_input(param->infile[i],&m_data->msa[i]));

                LOG_MSG("Aligned : %d", m_data->msa[i]->aligned);
        }
        for(i = 0; i < param->num_infiles;i++){
                if(m_data->msa[i]->aligned == ALN_STATUS_UNALIGNED){
                        ERROR_MSG("Sequences in file %s are not aligned.", param->infile[i]);
                }

        }
        /* do sanity checks AND SORT !!! */
        for(i = 0; i < param->num_infiles;i++){
                RUN(check_for_sequences(m_data->msa[i]));
                LOG_MSG("Detected: %d sequences.", m_data->msa[i]->numseq);
                RUN(sort_msa_by_seqname(m_data->msa[i]));
        }
        /* fill data structures needed for mumsa */
        for(i = 0; i < param->num_infiles;i++){
                /* LOG_MSG("reading: %s", param->infile[i]); */
                RUN(read_msa_into_msai(m_data->msa[i],&m_data->msai[i]));

        }
        RUN(alistat(m_data));

        for(i = 0; i < param->num_infiles;i++){
                RUN(msai_char_to_pos(m_data->msai[i]));

        }
        RUN(sanity_check_input(m_data));


        /* Ok I am ready to go */
        m_data->num_seq = m_data->msa[0]->numseq;
        RUN(galloc(&m_data->sim, m_data->num_seq,m_data->num_seq));

        RUN(calc_sim_pairs(m_data));
        RUN(calc_overlap(m_data));


        print_scores(m_data, param);
        free_mumsa_data(m_data);
        return OK;

ERROR:
        free_mumsa_data(m_data);
        return FAIL;
}



int print_mumsa_help(char * argv[])
{
        const char usage[] = " <refaln.msf> <testaln1.msf> <testaln2.msf>";
        fprintf(stdout,"\nUsage: %s %s\n\n",basename(argv[0]) ,usage);
        fprintf(stdout,"Options:\n\n");

        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--mode","specifies how to calculate scores." ,"[ref]"  );
        fprintf(stdout,"%*s%-*s: %s %s\n",3,"",MESSAGE_MARGIN-3,"--version (-V/-v)","Prints version." ,"[NA]"  );

        fprintf(stdout,"\nExamples:\n\n");

        fprintf(stdout,"mumsa aln.msf aln2.msf etc.... \n\n");

        fprintf(stdout,"Mumsa calculates the number of aligned residues shared between the reference and test alignment divided by:\n \"ref\" - the total number of aligned residues in the reference,\n\"test\" - the total number of aligned residues in the test alignment or:\n\"both\" - average number of aligned residues in test and ref alignment.\n\n");
        fprintf(stdout,"NOTE: mumsa always treats the first alignment as the reference.\n\n");



        return OK;
}

int print_mumsa_warranty(void)
{
        fprintf(stdout,"Here is the Disclaimer of Warranty section of the GNU General Public License (GPL):\n");
        fprintf(stdout,"\n");
        fprintf(stdout,"15. Disclaimer of Warranty.\n");
        fprintf(stdout,"THERE IS NO WARRANTY FOR THE PROGRAM, TO THE EXTENT PERMITTED BY\n");
        fprintf(stdout,"APPLICABLE LAW.  EXCEPT WHEN OTHERWISE STATED IN WRITING THE COPYRIGHT\n");
        fprintf(stdout,"HOLDERS AND/OR OTHER PARTIES PROVIDE THE PROGRAM \"AS IS\" WITHOUT WARRANTY\n");
        fprintf(stdout,"OF ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO,\n");
        fprintf(stdout,"THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR\n");
        fprintf(stdout,"PURPOSE.  THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE PROGRAM\n");
        fprintf(stdout,"IS WITH YOU.  SHOULD THE PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF\n");
        fprintf(stdout,"ALL NECESSARY SERVICING, REPAIR OR CORRECTION.\n");
        fprintf(stdout,"\n");
        fprintf(stdout,"A complete copy of the GPL can be found in the \"COPYING\" file.\n");
        return OK;
}

int print_mumsa_header(void)
{
        fprintf(stdout,"\n");
        fprintf(stdout,"Mumsa (%s)\n", PACKAGE_VERSION);
        fprintf(stdout,"\n");
        fprintf(stdout,"Copyright (C) 2020 Timo Lassmann\n");
        fprintf(stdout,"\n");
        fprintf(stdout,"This program comes with ABSOLUTELY NO WARRANTY; for details type:\n");
        fprintf(stdout,"`mumsa -showw'.\n");
        fprintf(stdout,"This is free software, and you are welcome to redistribute it\n");
        fprintf(stdout,"under certain conditions; consult the COPYING file for details.\n");
        fprintf(stdout,"\n");
        fprintf(stdout,"Please cite:\n");

        /*        fprintf(stdout,"  Mumsa 3: multiple sequence alignment of large data sets
                  Timo Lassmann
                  Bioinformatics, btz795, https://doi.org/10.1093/bioinformatics/btz795
        */
        fprintf(stdout,"  Lassmann T, Sonnhammer EL.\n");
        fprintf(stdout,"  \"Automatic extraction of reliable regions from multiple sequence alignments.\"\n");
        fprintf(stdout,"  BMC bioinformatics. 2007 May 1;8(S5):S9.\n");
        fprintf(stdout,"  https://doi.org/10.1186/1471-2105-8-S5-S9\n");
        fprintf(stdout,"\n");

        return OK;
}




void free_parameters(struct parameters* param)
{
        if(param){
                if(param->infile){
                        MFREE(param->infile);
                }
                MFREE(param);
        }
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
