/*
 	r_output.c
	
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


void r_statistics_output(char ** name,float* scores,float** overlap,int num)
{
	FILE *fp;
	int i,j,c;
	int pid;
	char filename[] = "/tmp/exampleXXXXXX";
	int fd;
	char **p = 0;
	p = tmalloc(sizeof(char*)*num);
	
	if((fd = mkstemp(filename)) < 0){
		perror("mkstemp failed");
		exit(1);
	}
	if ((fp = fdopen(fd, "w")) == NULL){
		fprintf(stderr,"can't open output	%s\n",filename);
	}else{
		if(num == num_alignments){
			fprintf(fp,"postscript(file=\"quality.eps\",width=6.0,height=12.0,horizontal=F,paper=\"special\",family=\"Helvetica\",onefile = FALSE,pointsize=18,bg =\"white\")\n");
		}else{
			fprintf(fp,"postscript(file=\"sequences.eps\",width=6.0,height=12.0,horizontal=F,paper=\"special\",family=\"Helvetica\",onefile = FALSE,pointsize=4,bg =\"white\")\n");
		}
		fprintf(fp,"x <- matrix(ncol = %d,nrow = %d)\n",num,num);
		for (i = 0;i < num;i++){
			p[i] = name[i];
			j =0;
			c = 0;
			while(name[i][j]){
				if(name[i][j] == 47){//47 = '/'
					c = j;
				}
				j++;
			}
			if(c){
			p[i] += c+1;
			}
		}
		fprintf(fp,"colnames(x) = c(");
		for (i = 0;i < num-1;i++){
			fprintf(fp,"\"%s:	%0.2f\",",p[i],scores[i]);
		}
		fprintf(fp,"\"%s:	%0.2f\")\n",p[num-1],scores[num-1]);
		
		fprintf(fp,"rownames(x) = c(");
		for (i = 0;i < num-1;i++){
			fprintf(fp,"\"%s:	%0.2f\",",p[i],scores[i]);
		}
		fprintf(fp,"\"%s:	%0.2f\")\n",p[num-1],scores[num-1]);
		
		for (i = 0;i < num;i++){
			for ( j = 0; j < num;j++){
				fprintf(fp,"x[%d,%d] = %f\n",i+1,j+1,overlap[i][j]);
			}
		}
		fprintf(fp,"dst.mat<-as.dist(1-x)\n");
		fprintf(fp,"hcl.average<-hclust(dst.mat,method=\"average\")\n");
		if(num == num_alignments){
			fprintf(fp,"plclust(hcl.average,frame.plot = TRUE,ann = TRUE,hang = -1,axes=TRUE,unit=FALSE,xlab = \"Alignment Programs\",sub = \"\",ylab = \"Average Distance\")\n");
		}else{
			fprintf(fp,"plclust(hcl.average,frame.plot = TRUE,ann = TRUE,hang = -1,axes=TRUE,unit=FALSE,xlab = \"\",sub = \"\",ylab = \"\")\n");
		}
		fprintf(fp,"dev.off()\n");
		fclose(fp);
	}
	
	
	free(p);
	
	//printf("Forking process	%d\n",getpid());
	pid = fork();
	//printf("The process id is %d and return value is %d\n",	getpid(), return_value);
	if ( pid == 0 ){ /* Child process */ 
		execl("/usr/bin/R","R","CMD","BATCH",filename,"logfile.out",NULL);
	}else{ /* Parent process pid is child's pid */
		wait(0);
		if ((unlink(filename)) == -1){
			perror("Removing temporary file failed");
			exit(1);
		}
	}
}

