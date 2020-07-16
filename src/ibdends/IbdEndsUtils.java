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
import ints.WrappedIntArray;
import vcf.GT;

/**
 * <p>Class {@code IbdEndsUtils} contains static utility methods
 * used when estimating IBD segment end points.</p>
 *
 * <p>The static methods in class {@code IbdEndsUtils} are thread-safe.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class IbdEndsUtils {

    private static final int BASE_POS_BACKOFF = 5_000_000;

    private IbdEndsUtils() {
        // private constructor to prevent instantiation
    }

    /**
     * Returns the estimated Morgan position of the specified base pair
     * position.  The Morgan position is estimated by linear interpolation.
     * @param basePos a map from marker index to base pair position
     * @param morganPos a map from marker index to Morgan position
     * @param inputBasePos a base pair position
     * @return the estimated Morgan position of the specified base pair
     * position.
     * @throws IllegalArgumentException if
     * {@code basePos.size() != morganPos.size()}
     * @throws NullPointerException if
     * {@code basePos == null || morganPos == null}
     */
    public static double morganPos(WrappedIntArray basePos,
            DoubleArray morganPos, int inputBasePos) {
        if (basePos.size()!=morganPos.size()) {
            throw new IllegalArgumentException("inconsistent data");
        }
        int index = basePos.binarySearch(inputBasePos);
        if (index>=0) {
            return morganPos.get(index);
        }
        else {
            int mapSizeM1 = basePos.size() - 1;
            int insPt = -index - 1;
            int aIndex = insPt - 1;
            int bIndex = insPt;
            if (aIndex==mapSizeM1) {
                insPt = basePos.binarySearch(basePos.get(mapSizeM1) - BASE_POS_BACKOFF);
                if (insPt<0) {
                    insPt = -insPt - 2;
                }
                assert insPt<mapSizeM1;
                aIndex = Math.max(insPt, 0);
                bIndex = mapSizeM1;
            }
            else if (bIndex==0) {
                insPt = basePos.binarySearch(basePos.get(0) + BASE_POS_BACKOFF);
                if (insPt<0) {
                    insPt = -insPt - 1;
                }
                assert insPt>0;
                aIndex = 0;
                bIndex = Math.min(insPt, mapSizeM1);
            }
            int x = inputBasePos;
            int a = basePos.get(aIndex);
            int b = basePos.get(bIndex);
            double fa = morganPos.get(aIndex);
            double fb = morganPos.get(bIndex);
            return fa + ( ( (double) (x-a)/(b-a)) * (fb-fa) );
        }
    }

    /**
     * Returns the first marker index greater than or equal to the specified
     * marker index such that the two specified haplotypes have discordant
     * alleles. Returns {@code this.gt.nMarkers()} if no such index exists.
     * @param gt the genotype data
     * @param hap1 a haplotype index
     * @param hap2 a haplotype index
     * @param marker a marker index
     * @return the first marker index greater than or equal to the specified
     * marker index such that the two specified haplotypes have discordant
     * alleles
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= gt.nMarkers()}
     * @throws IndexOutOfBoundsException of {@code hap1 < 0 || hap1 >= gt.nHaps()}
     * @throws IndexOutOfBoundsException of {@code hap2 < 0 || hap2 >= gt.nHaps()}
     * @throws NullPointerException if {@code gt == null}
     */
    public static int fwdDiscord(GT gt, int hap1, int hap2, int marker) {
        int end = gt.nMarkers();
        int m = marker;
        while (m<end && gt.allele(m, hap1)==gt.allele(m, hap2)) {
            ++m;
        }
        return m;
    }

    /**
     * Returns the probability that an IBD segment containing a focal point
     * has its right end point less than {@code y} Morgans from the focal point.
     * @param y a distance in Morgans
     * @param ne the constant effective population size
     * @return the probability that an IBD segment containing a focal point
     * has its right end point less than {@code y} Morgans from the focal point
     * @throws IllegalArgumentException if {@code y <= 0 || Double.isNaN(y)}
     * @throws IllegalArgumentException if
     * {@code ne <= 0.0 || Double.isFinite(ne) == false}
     */
    public static double F(double y, double ne) {
        if (y<=0 || Double.isNaN(y)) {
            throw new IllegalArgumentException(String.valueOf(y));
        }
        if (ne<=0.0 || Double.isFinite(ne)==false) {
            throw new IllegalArgumentException(String.valueOf(ne));
        }
        assert ne>0 && Double.isFinite(ne);
        double den = 2*ne*Math.expm1(2*y) + 1d;
        return 1d - 1d/den;
    }

    /**
     * Returns a value {@code y} such that {@code IbdEndsUtils.F(y, ne)} is
     * approximately equal to {@code p}.
     * @param p a probability satisfying {@code 0 < p && p < 1}
     * @param ne the (constant) effective population size
     * @return a value {@code y} such that {@code F(y, ne)} is approximately
     * equal to {@code p}
     * @throws IllegalArgumentException if
     * {@code p <= 0.0 || p >= 1.0 || Double.isNaN(p)}
     * @throws IllegalArgumentException if
     * {@code ne <= 0.0 || Double.isFinite(ne) == false}
     */
    public static double invF(double p, double ne) {
        if (p<=0d || p>=1d || Double.isNaN(p)) {
            throw new IllegalArgumentException(String.valueOf(p));
        }
        if (ne<=0.0 || Double.isFinite(ne)==false) {
            throw new IllegalArgumentException(String.valueOf(ne));
        }
        double d = 2*ne*(1 - p);
        return 0.5*Math.log((p + d)/d);
    }
}
