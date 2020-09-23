

struct aln_space* make_boot_sets(struct aln_space* aln_space, struct alignment** alignments, int gap_flag);
struct node** feed_boot_hash(struct node** hash, struct alignment* ap, unsigned int a ,unsigned int b, unsigned int aln, int gap_flag);
