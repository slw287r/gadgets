# quickly get genotypes for one given region from bam

## How to install
* This tools requires <a href="https://github.com/samtools/samtools.git" target="_blank">samtools</a> and <a href="https://zlib.net" target="_blank">zlib</a> to compile. Please follow the link and install samtools first.

* Successful installation of samtools results in the following files (presumably in `/usr/local`):
    ```
    /usr/local/include/zlib.h
    /usr/local/include/bam.h
    /usr/local/include/htslib/*.h
    /usr/local/lib/libz.a
    /usr/local/lib/libbam.a
    /usr/local/lib/libhts.a
    ```

* Compile quickgt
    ```
    gcc -o quickgt quickgt.c -lhts -lbam -lz
    ```
    or
    ```
    gcc -I/usr/local/include -L/usr/local/lib -o quickgt quickgt.c -lhts -lbam -lz
    ```
    incase `gcc` can't find the headers or libs

## How to use it

* Usage:
    ```
    quickgt <in.bam> <ref.fa> <chr:beg-end>
    ```
* Notes:

    * Input BAM index required
    * Maximum supported depth: 65536

    * BAM filter criteria:

        * BAM_FQCFAIL
        * BAM_FUNMAP
        * BAM_FSECONDARY
        * BAM_FDUP
        * BAM_FSUPPLEMENTARY

* Output format (goes to stdout):

Column | Description
----|--------------------------
CHR | contig name
POS | 1-based genomic coordinate
REF | base on reference
DEP | depth of coverage
A   | base A (on reads mapped fwd;rev)
C   | base C (on reads mapped fwd;rev)
G   | base G (on reads mapped fwd;rev)
T   | base T (on reads mapped fwd;rev)
N   | base N (on reads mapped)
INS | insertions (length insensitive, fwd;rev)
DEL | deletions (length insensitive, fwd;rev)


