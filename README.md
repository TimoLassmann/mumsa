# Mumsa

Mumsa is a program to compare multiple sequence alignments. The current version only supports comparing alignments with each other. Given a reference alignment (first argument) and a number of other 'test' alignments, mumsa will:
1) compute the set of aligned pairs of residues in the reference alignment
2) determine the fraction of correctly aligned pairs of residues in test alignments.


# Installation

## Release Tarball

Download tarball from [releases](https://github.com/TimoLassmann/mumsa/releases). Then:

``` bash
tar -zxvf mumsa-<version>.tar.gz
cd mumsa-<version>
./autogen.sh
./configure
make
make check
make install
```

## Developer version

Mumsa relies on the library tldevel. Install system-wide first or extract a tldevel release in the mumsa directory, then run autogen.

``` bash
git clone https://github.com/TimoLassmann/mumsa.git
cd mumsa-<version>
./autogen.sh
./configure
make
make check
make install
```

# Usage


``` bash
Usage: mumsa  <msa1> <msa2> ...

Options:

   --version          : Print version and exit
```

Mumsa can read alignments in msf, clu and aligned fasta formats.

# Example

```
mumsa aln.msf aln2.msf etc....

```
NOTE: mumsa always treats the first alignment as the reference.

# Please cite:
1. Lassmann T, Sonnhammer EL.
  "Automatic extraction of reliable regions from multiple sequence alignments."
  BMC bioinformatics. 2007 May 1;8(S5):S9.
  https://doi.org/10.1186/1471-2105-8-S5-S9
