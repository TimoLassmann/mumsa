#include "cmsa.h"
#include "adler.h"

struct alignment* calculate_adler(struct alignment* alignment)
{
        char tmp_name[6];
        int i;
        int j;
        int len = 0;
        for (i = 0; i < alignment->numseq;i++){
                len = (5 < alignment->lsn[i]) ?  5 : alignment->lsn[i];
                for (j = 0;j < len;j++){
                        tmp_name[j] = alignment->sn[i][j];
                }
                tmp_name[len] = 0;

                alignment->adler[i] = adler(alignment->seq[i],alignment->len)+ adler(tmp_name,len);//+i if sequence occurs more than once...
                //fprintf(stderr,"%d->%u\n",i,alignment->adler[i]);
        }

        //fprintf(stderr,"\n\n\n");


        return alignment;
}

/*

  Source code for the adler-32 checksum algorithm was adopted from the code shown at:

  http://en.wikipedia.org/wiki/Adler-32


*/

unsigned int adler(char* seq,int len)
{
        unsigned int a = 1;
        unsigned int b = 0;
        while (len) {
                unsigned tlen = len > 5550 ? 5550 : len;
                len -= tlen;
                do {
                        if(isalpha ((int)*seq)){
                                a += (int)toupper(*seq);
                                b += a;
                        }
                        seq++;
                } while (--tlen);
                a = (a & 0xffff) + (a >> 16) * (65536-MOD_ADLER);
                b = (b & 0xffff) + (b >> 16) * (65536-MOD_ADLER);
        }
        if (a >= MOD_ADLER){
                a -= MOD_ADLER;
        }
        b = (b & 0xffff) + (b >> 16) * (65536-MOD_ADLER);
        if (b >= MOD_ADLER){
                b -= MOD_ADLER;
        }
        return (b << 16) | a;
}
