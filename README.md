# ibd-ends

The **ibd-ends** program estimates a probability distribution for each endpoint
of an identity-by-descent (IBD) segment. The input IBD segments for
the **ibd-ends** program are output IBD segments from the
[**hap-ibd**](https://github.com/browning-lab/hap-ibd) program.

If you use **ibd-ends** in a published analysis, please report the program
version printed in the first line of the output log file, and
please cite the article that describes the **ibd-ends** method:

> [Sharon R Browning](https://sites.uw.edu/sguy/),
[Brian L Browning](http://faculty.washington.edu/browning/). (2020)
Probabilistic estimation of identity by descent segment endpoints and
detection of recent selection. American Journal of Human Genetics 107(5):895-910.
[https://doi.org/10.1016/j.ajhg.2020.09.010](https://doi.org/10.1016/j.ajhg.2020.09.010)

Last updated May 20, 2022

## Contents

* [Installation](#installation)
* [Running ibd-ends](#running-ibd-ends)
  * [Required parameters](#required-parameters)
  * [Optional parameters](#optional-parameters)
* [Output files](#output-files)
* [License](#license)

## Installation

You can download the most recent version of the
[ibd-ends](https://faculty.washington.edu/browning/ibd-ends.jar) program
with the command:

    wget https://faculty.washington.edu/browning/ibd-ends.jar

or you can download the source files and create the executable file
with the commands:

    git clone https://github.com/browning-lab/ibd-ends.git
    javac -cp ibd-ends/src/ ibd-ends/src/ibdends/IbdEndsMain.java
    jar cfe ibd-ends.jar ibdends/IbdEndsMain -C ibd-ends/src/ ./
    jar -i ibd-ends.jar

[Back to Contents](#contents)

## Running ibd-ends

The **ibd-ends** program requires Java version 1.8 (or a later version). Use of an
earlier Java version will produce an "Unsupported Class Version" error.

The command:

    java -jar ibd-ends.jar

prints a summary of the command line arguments.

To run **ibd-ends**, enter the following command:

    java -Xmx[GB]g -jar ibd-ends.jar [arguments]

where [GB] is the maximum number of gigabytes of memory to use, and
[arguments] is a space-separated list of parameter values, each expressed as
**parameter=value**.

The shell script
[run.ibd-ends.test](https://raw.githubusercontent.com/browning-lab/ibd-ends/master/test/run.ibd-ends.test)
will run a test **ibd-ends** analysis.

[Back to Contents](#contents)

### Required Parameters

The **ibd-ends** program has four required parameters.  Any input file
with name ending in ".gz" is assumed to be gzip-compressed.

* **gt=[file]** where **[file]** is a
[Variant Call Format](https://faculty.washington.edu/browning/intro-to-vcf.html)
(VCF) file containing a GT FORMAT subfield.  All genotypes must be phased, have
the phased allele separator ('|'), and have no missing alleles. If your data
is unphased, you can phase your data using the
[Beagle](https://faculty.washington.edu/browning/beagle/beagle.html) program.
A VCF record may have multiple ALT alleles. Only one chromosome is processed
per analysis.  If the VCF file contains more than one chromosome,
the **chrom** parameter determines the chromosome that is analyzed. The
analysis relies on an accurate genetic map.  Any input VCF records that are
outside the genetic map are excluded from the analysis.

* **ibd=[file]** where **[file]** is the IBD output file from a
[hap-ibd](https://github.com/browning-lab/hap-ibd) analysis of the VCF file
specified with the **gt** parameter. An IBD segment will be ignored
if it involves a sample that is not present in the input VCF file or if the
segment lies completely outside the boundaries of the genetic map. If an
input IBD segment endpoint extends beyond a terminal genetic map position,
the input segment endpoint will set equal to the terminal map position.

* **map=[file]** where **[file]** is a
[PLINK format genetic map](http://zzz.bwh.harvard.edu/plink/data.shtml#map)
with cM units for the chromosome being analyzed (see the **gt** parameter).
The **ibd-ends** program uses linear interpolation to estimate the genetic
position of base coordinates between map positions.
The chromosome identifier in the input VCF file and map file must be identical.
HapMap genetic maps in cM units are available for
[GRCh36, GRCh37, and GRCh38](https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/).

* **out=[string]** where **[string]** is the output filename prefix.

[Back to Contents](#contents)

### Optional Parameters

Optional parameters have sensible default values.  Some analyses
will require specifying one or more of the the first four optional parameters:
**chrom**, **quantiles**, **samples**, and **nthreads**.

* **chrom=[chrom]:[start]‑[end]** specifies the chromosome or chromosome interval
to be analyzed: **[chrom]** is the CHROM identifier in the
input VCF and IBD files, **[start]** is the first base pair position, and
**[end]** is the last base pair position.
An entire chromosome, the beginning of a chromosome, or the end of a
chromosome may be specified with "**chrom=[chrom]**", "**chrom=[chrom]:‑[end]**",
and "**chrom=[chrom]:[start]‑**" respectively. If a **chrom** parameter is not
specified, the **chrom** parameter will be set to the CHROM field in the
first VCF record.

* **quantiles=[comma-separated numbers: 0.0 < number < 1.0]**
specifies the quantiles of the probability distribution of the IBD segment
endpoint that will be printed to the [IBD output file](#output-files)
(**default: quantiles=0.5**). For example, **quantiles=0.05,0.5,0.95** will
print the 0.05, 0.50, and 0.95 quantiles of each IBD segment endpoint.

* **nsamples=[integer ≥ 0]** specifies the number of random samples taken from
the probability distribution for each IBD segment endpoint
(**default: nsamples=0**). If **nsamples>0**, the sampled start and end
positions are printed to the [IBD output file](#output-files).

* **nthreads=[integer ≥ 1]** specifies the number of computational threads to
use. The default **nthreads** parameter is the number of CPU cores.

* **err=[0.0 < number < 1.0]** specifies the probability that corresponding
alleles on two IBD haplotypes are discordant (**default: err=0.0005**).

* **estimate-err=[true/false]**. If true, the output **log**
file will include an estimate of the **err** parameter
(**default: estimate-err=true**). If the estimated error differs
from the default **err** parameter by more than a factor of 3, we recommend
re-running the analysis using the estimated error parameter.

* **gc-err=[0.0 < number < 1.0]** specifies the probability that corresponding
alleles on two IBD haplotypes will be discordant within a gene conversion tract
(**default: gc-err=0.1**).

* **gc-bp=[integer ≥ 0]** specifies the maximal length of a gene conversion
tract (**default: gc-bp=1000**).

* **min-maf=[0.0 < number ≤ 0.5]** specifies the minimum minor allele frequency.
Markers with minor allele frequency less than the specified minimum
frequency will be excluded from the analysis (**default: min-maf=0.001**). For
multi-allelic markers, the minor allele frequency is the second-largest
allele frequency.

* **seed=[integer]** specifies a seed for generating random numbers
when sampling IBD endpoints (**default: seed=-99999**)

* **excludesamples=[file]** where **[file]** is a text file containing samples
(one sample per line) to be excluded from the analysis. Any input IBD segments
or VCF data involving an excluded sample will be ignored.

[Back to Contents](#contents)

## Output files
The **ibd-ends** program produces two output files: a **log** file, and
an **ibd** file.  The **log** file (.log) contains a summary of the analysis.
The **ibd** file (.ibd.gz) is tab delimited and gzip compressed.  The
**ibd** file contains a header line and a line for each input IBD segment that
is analyzed.

The first 8 columns report the input IBD segments.  The 9th column is the focus
position, which is the position from which the distances to the segment start
and end positions are estimated.  The header line labels and descriptions
for the first nine columns are:

Column | Label | Description
:---:  | :---  | :---
1 | ID1   | First sample identifier
2 | HAP1  | First sample haplotype (1 or 2)
3 | ID2   | Second sample identifier
4 | HAP2  | Second sample haplotype (1 or 2)
5 | CHROM | Chromosome identifier
6 | STA   | Segment start position
7 | END   | Segment end position
8 | CM    | cM length
9 | FOCUS | Focus position

After the first nine columns, the next set of columns report the segment start
position, end position, and cM length for each endpoint quantile in the order
specified by the **quantiles** parameter.  If **nsamples>0**, there is a final
set of columns that report the segment start position, end position, and cM
length for each pair of sampled endpoints.

Header line labels for columns reporting the start position,
end position, and cM length begin with "STA", "END", and CM" respectively.
For example, columns reporting the 0.05-th quantile of the endpoint
distribution are labeled "STA.05", "END.05", and "CM.05". Similarly,
columns reporting the 3rd set of sampled endpoints are labeled "STA-3",
"END-3", and "CM-3".  Note that "." is used in column labels for quantiles,
and "-" is used in column labels for sampled endpoints.


[Back to Contents](#contents)

## License
The **ibd-ends** program is licensed under the Apache License, Version 2.0 (the License).
You may obtain a copy of the License from http://www.apache.org/licenses/LICENSE-2.0

[Back to Contents](#contents)
