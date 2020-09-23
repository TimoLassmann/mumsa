/*
	interface.c
	
	Released under GPL - see the 'COPYING' file   
	
	Copyright (C) 2006 Timo Lassmann <timolassmann@gmail.com>
	
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

#include "cmsa.h"

struct parameters* interface(struct parameters* param,int argc,char **argv)
{
	int i,c;
	static char  usage[] ="CMSA version 1.0, Copyright (C) 2006 Timo Lassmann <timolassmann@gmail.com>\n\
http://msa.cgb.ki.se/\n\n\
Usage: overlap [options] [test_alignment.fasta] [test_alignment 1.fasta] [test alignment_2.fasta]...\n\n\
Option:\n\
	-a	prints 'best' alignment to stdout\n\
	-g	consider residues aligned to gaps in addition to pairs of aligned residues\n\
	-r	calculate overlap scores in respect to first 'reference' alignment\n\
	-s	server\n\
	";
	
	param = tmalloc(sizeof(struct parameters));
	param->print_alignment_flag = 0;
	param->server_flag = 0;
	param->ref_flag = 0;
	param->gap_flag = 0;
	param->quiet_flag = 0;
	param->bootstrap = 0;
	param->bootcutoff = 0.5f;
	param->best_flag = 0.0f;
	param->print_column_cutoff = 0.0f;
	param->idal = 0.0f;
	param->format = 0;
	param->outfile = 0;
	param->infile = 0;
	param->print_spacer = 0;
	param->x1 = 0.0;
	param->x2 = 0.0;
	param->x3 = 0.0;
	param->set_diff = 0;
	
	param->sort = 0;
	param->order = 0;
	
	while (1){	
		static struct option long_options[] ={
			{"order", required_argument, 0,0},
			{"sort", required_argument, 0,0},
			{"idal", required_argument, 0,0},
			{"x1", required_argument, 0,0},
			{"x2", required_argument, 0,0},
			{"x3", required_argument, 0,0},
			{"boot",required_argument, 0,0},
			{"boot_cutoff",required_argument, 0,0},
			{"column",required_argument, 0,0},
			{"print_column",required_argument, 0,0},
			{"column_cutoff",required_argument, 0,0},
			{"spacer",0, 0,0},
			{"set_diff",required_argument, 0,0},
		
			{"quiet",  0, 0, 'q'},
			{"gap",  0, 0, 'g'},
			{"reference",  required_argument, 0, 'r'},
			{"format",  required_argument, 0, 'f'},
			{"print_best_alignment",  0, 0, 'a'},
			{"alignment",  0, 0, 'a'},
			{"output",  required_argument, 0, 'o'},
			{"outfile",  required_argument, 0, 'o'},
			{"out",  required_argument, 0, 'o'},
			{"server",  0, 0, 's'},
			{0, 0, 0, 0}
		};
		
		int option_index = 0;
                c = getopt_long_only (argc, argv,"sar:gqb:f:",long_options, &option_index);
		if (c == -1){
                        break;
                }
		//fprintf(stderr,"GAGA	%d\n",c);
		switch(c) {
			case 0:
				if (long_options[option_index].flag != 0){
					break;
				}
				switch (option_index){
					case 0:
					case 1:
						param->sort = optarg;
						break;
					case 2:
						param->idal = atof(optarg);
						break;
					case 3:
						param->x1 = atof(optarg);
						break;
					case 4:
						param->x2 = atof(optarg);
						break;
					case 5:
						param->x3 = atof(optarg);
						break;
					case 6:
						param->bootstrap = atoi(optarg);
						break;
					case 7:
						param->bootcutoff = atof(optarg);
						break;
						
					case 8:
					case 9:
					case 10:	
						param->print_column_cutoff = atof(optarg);
						break;
					case 11:
						param->print_spacer = 1;
						break;
					case 12:
						param->set_diff  = atoi(optarg);
						break;
					default:
						break;
				}
				break;
			case 'o':
				param->outfile = optarg;
				break;
			case 'f':
				param->format = optarg;
				break;
			case 'a':
				param->print_alignment_flag = 1;
				break;
			case 's':
				param->server_flag = 1;
				break;
			case 'r':
				param->ref_flag = optarg;
				break;
			case 'g':
				param->gap_flag = 1;
				break;
			case 'q':
				param->quiet_flag = 1;
				break;
			case 'b':
				param->best_flag = atof(optarg);
				break;
			case '?':
				(void)fprintf(stderr, "unknown argument: %c\n",optopt);
				(void)exit(EXIT_FAILURE);
				break;
			default:
				fprintf(stderr,"default\n\n\n\n");
                                abort ();
		}
	}
	
	/*if(param->idal  > 1.0 || param->idal < 0.0){
		fprintf(stderr,"Idal has to be: 0 < idal < 1\n");
		free_param(param);
		(void)exit(EXIT_FAILURE);
	}*/
	
	/*float test = 0.0;
	test = param->x1 + param->x2 + param->x3;
	if(test == 0.0){
		fprintf(stderr,"Something wrong with parameters x1,x2,x3\n");
		free_param(param);
		(void)exit(EXIT_FAILURE);
	}
	param->x1 /= test;
	param->x2 /= test;
	param->x3 /= test;*/
	
	if(param->quiet_flag){
		fclose(stderr);
	}
	if(param->server_flag){
		fclose(stderr);
	}
	if (!param->server_flag)fprintf(stderr,"%s\n",usage);

	if (optind < argc){
		num_alignments = argc - optind;
		param->infile = tmalloc(sizeof(char**)*num_alignments);
		for (i = 0 ; i < num_alignments;i++){
			param->infile[i] =  0;
		}
		i = 0;
		while (optind < argc){
			param->infile[i] =  argv[optind++];
			i++;
		}
        }else{
		(void)fprintf(stderr, "No Alignments...\n");
		(void)exit(EXIT_FAILURE);
		free_param(param);
	}

	return param;
}





