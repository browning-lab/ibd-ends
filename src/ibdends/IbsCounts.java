/*
 * Copyright 2020 Brian L. Browning
 *
 * This file is part of the ibd-ends program.
 *
 * Licensed under the Apache License, Version 2.0 (the License);
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an AS IS BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package ibdends;

import blbutil.Utilities;
import ints.CharArray;
import ints.IntArray;
import ints.IntList;
import ints.UnsignedByteArray;
import ints.WrappedIntArray;
import java.util.Arrays;
import java.util.Random;
import java.util.stream.IntStream;
import vcf.GT;
import vcf.RefGT;

/**
 * <p>Class {@code IbsCounts} counts the number of haplotype pairs
 * that are identical by state (IBS) on intervals between a start marker
 * and an end marker.</p>
 *
 * <p>Instances of class {@code IbsCounts} are immutable.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class IbsCounts {

    private final int nHaps;
    private final int[][] counts;

    /**
     * Constructs a new {@code IbsCounts} for the specified data.
     * @param par the analysis parameters
     * @param gt the phased, non-missing genotype data
     * @throws IllegalArgumentException if {@code gt.nHaps() < 2}
     * @throws IllegalArgumentException if
     * {@code ((long) n*(n-1)) >= Integer.MAX_VALUE} where
     * {@code n = Math.min(par.local_haps(), gt.nHaps()}
     * @throws NullPointerException if any argument is
     * {@code par == null || gt == null}
     */
    public IbsCounts(IbdEndsPar par, RefGT gt) {
        if (gt.nHaps()<2) {
            throw new IllegalArgumentException(String.valueOf(gt.nHaps()));
        }
        int[] hapList = hapList(par, gt);
        long n = hapList.length;
        if (n*(n-1)>=Integer.MAX_VALUE) {
            throw new IllegalArgumentException(String.valueOf(n));
        }
        IntArray[] alleles = alleles(gt, hapList);
        Boolean[] isMonomorphic = isMonomorphic(alleles);
        int maxAlleleCnt = maxAlleles(gt);

        double maxLocalCDF = par.max_local_cdf();
        long minIbsPairs = (long) Math.ceil(((1.0-maxLocalCDF)*n)*(n-1));
        this.nHaps = hapList.length;
        this.counts = IntStream.range(0, gt.nMarkers())
                .parallel()
                .mapToObj(m -> counts(gt, alleles, maxAlleleCnt, isMonomorphic,
                        m, minIbsPairs))
                .toArray(int[][]::new);
    }

    /* Constructor that reverses markers in specified {@code IbsCounts) instance */
    private IbsCounts(IbsCounts fwdCnts) {
        int nMarkers = fwdCnts.nMarkers();
        this.nHaps = fwdCnts.nHaps;
        this.counts = IntStream.range(0, nMarkers)
                .parallel()
                .mapToObj(m -> revCounts(fwdCnts, m))
                .toArray(int[][]::new);
    }

    private static int[] revCounts(IbsCounts fwdCnts, int revStart) {
        IntList revCnts = new IntList(1<<8);
        int inclEnd = fwdCnts.nMarkers() - 1 - revStart;
        for (int start=inclEnd; start>=0 && inclEnd<fwdCnts.end(start); --start) {
            revCnts.add(fwdCnts.counts(start, inclEnd));
        }
        return revCnts.toArray();
    }

    private static int[] hapList(IbdEndsPar par, GT gt) {
        int maxLocalHaps = par.local_haps();
        if (gt.nHaps()<=maxLocalHaps) {
            return IntStream.range(0, gt.nHaps()).toArray();
        }
        else {
            Random rand = new Random(par.seed());
            int[] allHaps = IntStream.range(0, gt.nHaps()).parallel().toArray();
            Utilities.shuffle(allHaps, maxLocalHaps, rand);
            Arrays.sort(allHaps, 0, maxLocalHaps);
            return Arrays.copyOf(allHaps, maxLocalHaps);
        }
    }


    private static IntArray[] alleles(GT gt, int[] hapList) {
        return IntStream.range(0, gt.nMarkers())
                .parallel()
                .mapToObj(m -> alleles(gt, m, hapList))
                .toArray(IntArray[]::new);
    }

    private static IntArray alleles(GT gt, int marker, int[] hapList) {
        int[] alleles = new int[hapList.length];
        for (int j=0; j<alleles.length; ++j) {
            alleles[j] = gt.allele(marker, hapList[j]);
        }
        int nAlleles = gt.marker(marker).nAlleles();
        if (nAlleles<=256) {
            return new UnsignedByteArray(alleles);
        }
        else if (nAlleles<=65536) {
            return new CharArray(alleles);
        }
        else {
            return new WrappedIntArray(alleles);
        }
    }

    private static Boolean[] isMonomorphic(IntArray[] alleles) {
        return Arrays.stream(alleles)
                .parallel()
                .map(ia -> isMonomorphic(ia))
                .toArray(Boolean[]::new);
    }

    private static Boolean isMonomorphic(IntArray alleles) {
        for (int j=1, n=alleles.size(); j<n; ++j) {
            if (alleles.get(j)!=alleles.get(j-1)) {
                return false;
            }
        }
        return true;
    }

    private static int maxAlleles(GT gt) {
        return IntStream.range(0, gt.nMarkers())
                .parallel()
                .map(m -> gt.marker(m).nAlleles())
                .max()
                .orElse(0);
    }

    private static int[] counts(GT gt, IntArray[] alleles, int maxAlleleCnt,
            Boolean[] isMonomorphic, int start, long minIbsPairs) {
        int nMarkers = alleles.length;
        int nHaps = alleles[start].size();
        IntList cnts = new IntList(1<<8);
        int[] hap2Seq = new int[nHaps];
        int[] seq2Cnt = new int[nHaps];
        int[] seqAlMap = new int[maxAlleleCnt*nHaps];
        seq2Cnt[0] = nHaps;
        int nSeq = 1;
        int ibsPairs = sumIbsPairs(seq2Cnt, nSeq);
        for (int m=start; m<nMarkers && ibsPairs>=minIbsPairs; ++m) {
            if (isMonomorphic[m]) {
                cnts.add(ibsPairs);
            }
            else {
                int nAlleles = gt.marker(m).nAlleles();
                Arrays.fill(seqAlMap, 0, nAlleles*nSeq, -1);
                Arrays.fill(seq2Cnt, 0, nSeq, 0);
                nSeq = 0;
                for (int j=0; j<nHaps; ++j) {
                    int seqAlIndex = hap2Seq[j]*nAlleles + alleles[m].get(j);
                    int newSeq = seqAlMap[seqAlIndex];
                    if (newSeq<0) {
                        newSeq = nSeq++;
                        seqAlMap[seqAlIndex] = newSeq;
                    }
                    hap2Seq[j] = newSeq;
                    ++seq2Cnt[newSeq];
                }
                ibsPairs = sumIbsPairs(seq2Cnt, nSeq);
                if (ibsPairs>=minIbsPairs) {
                    cnts.add(ibsPairs);
                }
            }
        }
        return cnts.toArray();
    }

    private static int sumIbsPairs(int[] seqCnts, int nSeq) {
        int sum = 0;
        for (int j=0; j<nSeq; ++j) {
            sum += seqCnts[j]*(seqCnts[j]-1);
        }
        return sum;
    }

    /**
     * Returns the number of markers.
     * @return the number of markers
     */
    public int nMarkers() {
        return counts.length;
    }

    /**
     * Returns the number of random haplotypes used to generate IBS segment
     * counts.
     * @return the number of random haplotypes used to generate IBS segment
     * counts
     */
    public int nHaps() {
        return nHaps;
    }

    /**
     * Returns the number of haplotype pairs from a random set of
     * {@code this.nHaps()} haplotypes that are IBS on the interval
     * from the specified start marker index (inclusive) to the specified
     * end marker index (inclusive).
     * @param start the start marker (inclusive)
     * @param inclEnd the end marker (inclusive)
     * @return the number of haplotype pairs from a random set of
     * {@code this.nHaps()} haplotypes that are IBS on the interval
     * from the specified start marker index (inclusive) to the specified
     * end marker index (inclusive)
     * @throws IndexOutOfBoundsException if
     * {@code start < 0 || start >= this.nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code inclEnd < start || inclEnd >= this.end(start)}
     */
    public int counts(int start, int inclEnd) {
        return counts[start][inclEnd - start];
    }

    /**
     * Returns the exclusive end of the set of end markers for which the
     * interval from the specified start marker to end marker (inclusive)
     * have stored IBS segments counts.
     * @param start a marker index
     * @return the exclusive end of the set of end markers for which the
     * interval from the specified start marker to end marker (inclusive)
     * have stored IBS segments counts
     * @throws IndexOutOfBoundsException if
     * {@code start < 0 || start >= this.nMarkers()}
     */
    public int end(int start) {
        return start + counts[start].length;
    }

    /**
     * Returns an {@code IbsCounts} instance obtained by reversing the marker
     * order of {@code this}.
     * @return an {@code IbsCounts} instance obtained by reversing the marker
     * order of {@code this}
     */
    public IbsCounts reverse() {
        return new IbsCounts(this);
    }
}
