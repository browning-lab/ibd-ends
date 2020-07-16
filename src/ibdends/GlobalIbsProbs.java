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

import blbutil.DoubleArray;
import java.util.Arrays;
import java.util.Random;
import java.util.stream.IntStream;
import vcf.GT;
import vcf.RefGT;

/**
 * <p>Class {@code GlobalIbsProbs} estimates the one-sided global IBS
 * length distribution.  A one-sided IBS length is the distance in Morgans
 * in a specified direction from a focus position to the first position
 * for which two haplotypes have discordant alleles, or to the last marker
 * on the chromosome if there are no discordant alleles in that direction.</p>
 *
 * <p>Instances of class {@code GlobalIbsProbs} are immutable</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class GlobalIbsProbs {

    private final double[] lengths;
    private final double reciprocalSize;

    /**
     * Constructs a new {@code GlobalIbsProbs} for the specified data.  The
     * input phased genotype data is required to contain at least two haplotypes
     * because pairs of distinct haplotypes are sampled to estimate the global
     * one-directional IBS segment length distribution.
     * @param gt the phased, non-missing genotype data
     * @param morganPos an array containing the position of each marker in Morgans
     * @param par the analysis parameters
     * @throws IllegalArgumentException if {@code gt.isPhased() == false}
     * @throws IllegalArgumentException if {@code gt.nHaps() < 2}
     * @throws IllegalArgumentException if
     * {@code gt.nMarkers() != morganPos().size()}
     * @throws NullPointerException if
     * {@code gt == null || morganPos == null || par == null}
     */
    public GlobalIbsProbs(RefGT gt, DoubleArray morganPos, IbdEndsPar par) {
        checkArgs(gt, morganPos);
        double[][] lengths0 = IntStream.range(0, par.global_pos())
                .parallel()
                .mapToObj(i -> sampleIbsLengths(gt, morganPos,
                        par.global_segments(), (par.seed() + i)))
                .toArray(double[][]::new);

        lengths0 = filter(lengths0, par.global_quantile(), par.global_factor());

        this.lengths = Arrays.stream(lengths0)
                .parallel()
                .flatMapToDouble(da -> Arrays.stream(da))
                .sorted()
                .toArray();
        this.reciprocalSize = 1.0/this.lengths.length;
    }

    private static void checkArgs(RefGT gt, DoubleArray morganPos) {
        if (gt.isPhased()==false) {
            throw new IllegalArgumentException("gt.isPhased()=false");
        }
        if (gt.nMarkers()!=morganPos.size()) {
            throw new IllegalArgumentException("inconsistent data");
        }
        if (gt.nHaps()<2) {
            throw new IllegalArgumentException(String.valueOf(gt.nHaps()));
        }
    }

    private static double[] sampleIbsLengths(GT gt, DoubleArray genPos, int nSegs,
            long seed) {
        Random rand = new Random(seed);
        double pos = randomGenPos(rand, genPos);
        double midPos = 0.5*(genPos.get(0) + genPos.get(genPos.size()-1));
        return IntStream.range(0, nSegs)
                .mapToDouble(i -> sampleIbsLength(gt, genPos, midPos, pos, rand))
                .sorted()
                .toArray();
    }

    private static double randomGenPos(Random rand, DoubleArray genPos) {
        double startMorgans = genPos.get(0);
        double endMorgans = genPos.get(genPos.size()-1);
        double pos = startMorgans + rand.nextDouble()*(endMorgans - startMorgans);
        if (pos>=endMorgans) {
            pos = Math.nextDown(pos);
        }
        return pos;
    }

    private static double sampleIbsLength(GT gt, DoubleArray genPos,
            double midPos, double pos, Random rand) {
        int nHaps = gt.nHaps();
        assert nHaps>1;
        int h1 = rand.nextInt(nHaps);
        int h2 = rand.nextInt(nHaps);
        while (h1==h2) {
            h2 = rand.nextInt(nHaps);
        }
        if (pos < midPos) {
            return fwdLength(gt, genPos, pos, h1, h2);
        }
        else {
            return bwdLength(gt, genPos, pos, h1, h2);
        }
    }

    private static double fwdLength(GT gt, DoubleArray genPos, double pos,
            int h1, int h2) {
        int m = genPos.binarySearch(pos);
        if (m<0) {
            m = -m - 1;
        }
        while (m<genPos.size() && gt.allele(m, h1)==gt.allele(m, h2)) {
            ++m;
        }
        if (m==genPos.size()) {
            --m;
        }
        return (genPos.get(m) - pos);
    }

    private static double bwdLength(GT gt, DoubleArray genPos, double pos,
            int h1, int h2) {
        int m = genPos.binarySearch(pos);
        if (m<0) {
            m = -m - 2;
        }
        assert genPos.get(m) <= pos;
        while (m>=0 && gt.allele(m, h1)==gt.allele(m,h2)) {
            --m;
        }
        if (m<0) {
            ++m;
        }
        return (pos - genPos.get(m));
    }

    private static double[][] filter(double[][] lengths, float quantile,
            double factor) {
        int index = (int) Math.floor(quantile*lengths[0].length);
        double[] sortedVals = Arrays.stream(lengths)
                .parallel()
                .mapToDouble(da -> da[index])
                .sorted()
                .toArray();

        int n = lengths.length;
        double median = 0.5*(sortedVals[(n-1)>>1] + sortedVals[n>>1]);
        double maxValue = factor*median;
        return Arrays.stream(lengths)
                .parallel()
                .filter(da -> da[index]<=maxValue)
                .toArray(double[][]::new);
    }

    /**
     * Returns the number of filtered, sampled segments lengths.
     * @return the number of filtered, sampled segments lengths
     */
    public int nLengths() {
        return lengths.length;
    }

    /**
     * Returns the proportion of sampled, filtered one-sided discord distances
     * that are less than or equal to the specified length. Returns
     * {@code 1.0/this.nLengths()} if the the value is less than all filtered,
     * sampled discord distances and returns
     * {@code (this.nLengths() - 1.0)/this.nLengths()} if the the value is
     * greater than or equal to all filtered, sampled discord distances.
     * @param morgans a length in Morgans
     * @return the proportion of filtered, sampled one-sided discord distances
     * that are less than or equal to the specified length
     * @throws IllegalArgumentException if {@code Double.isNaN(morgans) == true}
     */
    public double cdf(double morgans) {
        if (Double.isNaN(morgans)) {
            throw new IllegalArgumentException(String.valueOf(morgans));
        }
        int index = Arrays.binarySearch(lengths, morgans);
        if (index>=0) {
            while ((index+1)<lengths.length && lengths[index]==lengths[index+1]) {
                ++index;
            }
            ++index;
        }
        else {
            index = -index - 1;
        }
        if (index==0) {
            ++index;
        }
        if (index==lengths.length) {
            --index;
        }
        return index*reciprocalSize;
    }
}
