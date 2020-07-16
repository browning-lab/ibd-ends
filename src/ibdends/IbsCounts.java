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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Objects;
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

        double maxLocalCDF = par.max_local_cdf();
        long minIbsPairs = (long) Math.ceil(((1.0-maxLocalCDF)*n)*(n-1));
        this.nHaps = hapList.length;
        this.counts = IntStream.range(0, gt.nMarkers())
                .parallel()
                .mapToObj(m -> counts(gt, alleles, isMonomorphic, m, minIbsPairs))
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

    private static int[] counts(GT gt, IntArray[] alleles, Boolean[] isMonomorphic,
            int start, long minIbsPairs) {
        int nMarkers = alleles.length;
        int n = alleles[start].size();
        IntList cnts = new IntList(1<<8);
        List<int[]> hapLists = new ArrayList<>(1);
        hapLists.add(IntStream.range(0, n).toArray());
        int lastIbsPairs = sumIbsPairs(hapLists);
        for (int m=start; m<nMarkers && lastIbsPairs>=minIbsPairs; ++m) {
            if (Objects.equals(isMonomorphic[m], Boolean.FALSE)) {
                hapLists = partitionLists(gt, alleles[m], m, hapLists);
                int ibsPairs = sumIbsPairs(hapLists);
                if (ibsPairs>=minIbsPairs) {
                    cnts.add(ibsPairs);
                }
                lastIbsPairs = ibsPairs;
            }
            else {
                cnts.add(lastIbsPairs);
            }
        }
        return cnts.toArray();
    }

    /**
     * Returns the probability that two distinct haplotypes in the same list
     * have different alleles at the specified marker.  The contract for this
     * method is undefined if any haplotype is present more than once in the
     * specified lists of haplotypes.
     * @param gt the genotype data
     * @param alleles the haplotype alleles
     * @param m the marker
     * @param hapLists lists of haplotypes
     * @return the probability that two distinct haplotypes in the same list
     * have different alleles at the specified marker
     */
    private static List<int[]> partitionLists(GT gt, IntArray alleles, int m,
            List<int[]> hapLists) {
        List<int[]> nextHapLists = new ArrayList<>(hapLists.size());
        IntList[] scratchLists = IntStream.range(0, gt.marker(m).nAlleles())
                .mapToObj(j -> new IntList())
                .toArray(IntList[]::new);
        for (int i=0, n=hapLists.size(); i<n; ++i) {
            int[] hapList = hapLists.get(i);
            for (int h : hapList) {
                scratchLists[alleles.get(h)].add(h);
            }
            for (IntList scratchList : scratchLists) {
                if (scratchList.size()>1) {
                    nextHapLists.add(scratchList.toArray());
                }
                scratchList.clear();
            }
        }
        hapLists.clear();
        return nextHapLists;
    }

    private static int sumIbsPairs(List<int[]> ibsLists) {
        int sum = 0;
        for (int j=0, n=ibsLists.size(); j<n; ++j) {
            long len = ibsLists.get(j).length;
            sum += len*(len-1);
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
