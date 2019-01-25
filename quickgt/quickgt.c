#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <assert.h>
#include <ctype.h>
#include <zlib.h>
#include <bam.h>
#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <htslib/kfunc.h>
#include <htslib/khash.h>
#include <htslib/kstring.h>

#define MAX_DEP 0xFFFFF
#define basename(str) (strrchr(str, '/') ? strrchr(str, '/') + 1 : str)

static const unsigned BAM_FILTER = (BAM_FQCFAIL | BAM_FUNMAP | BAM_FSECONDARY | BAM_FDUP | BAM_FSUPPLEMENTARY);
static const char header[] ="CHR\tPOS\tREF\tDEP\tA\tC\tG\tT\tN\tINS\tDEL";

typedef struct
{
	bamFile fp;
	bam_hdr_t *hdr;
	bam_iter_t iter;
} aux_t;

int read_bam(void *data, bam1_t *b)
{
	aux_t *aux = (aux_t*)data;
	int ret;
	while (1)
	{
		ret = bam_iter_read(aux->fp, aux->iter, b);
		if (ret < 0) break;
		if (b->core.flag & BAM_FILTER) continue;
		break;
	}
	return ret;
}

void pileup(const char *bam, const char *ref, const char *reg)
{
	int tid = 0, pos = 0, beg = 0, end = 0, len = 0, n_plp = 0;
	aux_t *data = malloc(sizeof(aux_t)); assert(data);
	memset(data, 0, sizeof(aux_t));
	data->fp = bam_open(bam, "r"); assert(data->fp);
	data->hdr = bam_hdr_read(data->fp); assert(data->hdr);
	bam_hdr_t *h = data->hdr;
	hts_idx_t *idx = bam_index_load(bam); assert(idx);
	if (!hts_parse_reg(reg, &beg, &end))
	{
		fprintf(stderr, "Failed parsing region: %s\n", reg);
		exit(EXIT_FAILURE);
	}
	data->iter = bam_itr_querys(idx, data->hdr, reg); assert(data->iter);
	hts_idx_destroy(idx);
	beg = data->iter->beg;
	end = data->iter->end;
	const bam_pileup1_t *plp = calloc(1, sizeof(bam_pileup1_t)); assert(plp);
	bam_mplp_t mplp = bam_mplp_init(1, read_bam, (void **)&data);
	bam_mplp_set_maxcnt(mplp, MAX_DEP);
	faidx_t *fai = fai_load(ref); assert(fai);
	while (bam_mplp_auto(mplp, &tid, &pos, &n_plp, &plp) > 0) // pos: 0-based
	{
		if (pos < beg || pos >= end) continue;
		if (tid >= h->n_targets) continue;
		// reference sequence
		unsigned nt[9] = {0}, idl[7] = {0}, gap = 0;
		char *rseq = faidx_fetch_seq(fai, h->target_name[tid], pos, pos, &len);

		for (int j = 0; j < n_plp; ++j) // walk down through the pile
		{
			const bam_pileup1_t *p = plp + j;
			if (p->indel)
				++idl[(p->indel > 0 ? 0 : 1) + ((p->b->core.flag & BAM_FREVERSE) ? 5 : 0)];
			if (p->is_del || p->is_refskip)
			{
				++gap;
				continue;
			}
			char alt = bam_nt16_rev_table[bam1_seqi((char *)bam1_seq(p->b), p->qpos)];
			++nt[seq_nt16_int[seq_nt16_table[(int)alt]] + (((p->b->core.flag & BAM_FREVERSE) && (toupper(alt) != 'N')) ? 5 : 0)];
		}
		// output info
		printf("%s\t", h->target_name[tid]);
		printf("%d\t", pos + 1);
		printf("%s\t", rseq);
		printf("%d\t", n_plp - gap);
		printf("%d;%d\t", nt[0], nt[5]);
		printf("%d;%d\t", nt[1], nt[6]);
		printf("%d;%d\t", nt[2], nt[7]);
		printf("%d;%d\t", nt[3], nt[8]);
		printf("%d\t", nt[4]);
		printf("%d;%d\t", idl[0], idl[5]);
		printf("%d;%d\n", idl[1], idl[6]);
		free(rseq);
	}
	bam_mplp_destroy(mplp);
	bam_hdr_destroy(data->hdr);
	bam_itr_destroy(data->iter);
	fai_destroy(fai);
	bam_close(data->fp);
}

int main(int argc, char **argv)
{
	if (argc != 4)
	{
		fputc('\n', stderr);
		fputs("Quickly get genotype info from bam for given region\n", stderr);
		fputc('\n', stderr);
		fputs("\033[1mUsage\033[0m:\n", stderr);
		fprintf(stderr, "  \e[1;31m%s\e[0;0m <in.bam> <ref.fa> <chr:beg-end>\n", basename(argv[0]));
		fputc('\n', stderr);
		fputs("\033[1mNotes\033[0m:\n", stderr);
		fputs("  Input BAM index required\n", stderr);
		fputc('\n', stderr);
		fputs("  Output format (goes to stdout):\n", stderr);
		fputs("      CHR   contig name\n", stderr);
		fputs("      POS   1-based genomic coordinate\n", stderr);
		fputs("      REF   base on reference\n", stderr);
		fputs("      DEP   depth of coverage\n", stderr);
		fputs("      A     base A (on read mapped fwd;rev)\n", stderr);
		fputs("      C     base C (on read mapped fwd;rev)\n", stderr);
		fputs("      G     base G (on read mapped fwd;rev)\n", stderr);
		fputs("      T     base T (on read mapped fwd;rev)\n", stderr);
		fputs("      N     base N (on read mapped)\n", stderr);
		fputs("      INS   insertions (length insensitive, fwd;rev)\n", stderr);
		fputs("      DEL   deletions (length insensitive, fwd;rev)\n", stderr);
		fputc('\n', stderr);
		fputs("  Maximum supported depth: 65536\n", stderr);
		fputc('\n', stderr);
		fputs("  BAM filter criteria:\n", stderr);
		fputs("      BAM_FQCFAIL\n", stderr);
		fputs("      BAM_FUNMAP\n", stderr);
		fputs("      BAM_FSECONDARY\n", stderr);
		fputs("      BAM_FDUP\n", stderr);
		fputs("      BAM_FSUPPLEMENTARY\n", stderr);
		fputc('\n', stderr);
		exit(EXIT_FAILURE);
	}
	if (!strchr(argv[3], ':') || !strchr(argv[3], '-'))
	{
		fprintf(stderr, "[Error] Invalid region: %s\n", argv[3]);
		fputs("Valid region format: chr:beg-end\n", stderr);
		fputs("Both bed and end coordinates are 1-based\n", stderr);
		exit(EXIT_FAILURE);
	}
	fputs(header, stdout);
	fputc('\n', stdout);
	pileup(argv[1], argv[2], argv[3]);
}
