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
import java.util.Arrays;
import vcf.GT;

/**
 * <p>Class {@code QuantileEstimator} estimates the quantiles of an IBD segment
 * end point distribution.</p>
 *
 * <p>Instances of class {@code QuantileEstimator} are not thread-safe.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class QuantileEstimator {

    private static final double MIN_RATIO = 0.001;

    private final IbdEndsData data;
    private final Data fwdData;
    private final Data revData;
    private final double ne;
    private final double err;
    private final double gc_err;
    private final int gc_bp;

    // scratch space for storing CDF
    private final double[] cdf;
    private int cdfStart;
    private int cdfEnd;

    /**
     * Constructs a new {@code QuantileEstimator} instance for the specified data.
     * @param data the analysis input data
     * @throws NullPointerException if {@code data == null}
     */
    public QuantileEstimator(IbdEndsData data) {
        IbdEndsPar par = data.par();
        this.data = data;
        this.fwdData = new Data(data, true);
        this.revData = new Data(data, false);
        this.ne = par.ne();
        this.err = par.err();
        this.gc_bp = par.gc_bp();
        this.gc_err = par.gc_err();
        this.cdf = new double[data.fwdGT().nMarkers()];
    }

    /**
     * Returns the input data.
     * @return the input data
     */
    public IbdEndsData data() {
        return data;
    }

    /**
     * Estimates the specified quantiles of the forward endpoint of
     * the specified IBD segment.
     * @param hap1 the first IBD haplotype
     * @param hap2 the second IBD haplotype
     * @param ibdStartMorgans an estimated genetic position in Morgans of the
     * start of the IBD segment
     * @param focusPos the base pair coordinate of the focus position in
     * the IBD segment
     * @param probs a list of probabilities
     * @param quantiles an array which will contain the end-point quantiles
     * corresponding to the specified probabilities
     * @throws IndexOutOfBoundsException if
     * {@code hap1 < 0 || hap1 >= this.data().fwdGT().nHaps()}
     * @throws IndexOutOfBoundsException if
     * {@code hap2 < 0 || hap2 >= this.data().fwdGT().nHaps()}
     * @throws IllegalArgumentException if
     * {@code Double.isNaN(ibdStartMorgans) == true}
     * @throws IllegalArgumentException if the focus position is less than or
     * equal to the estimated base pair position of the specified start of the
     * IBD segment
     * @throws IllegalArgumentException if any element {@code p} of the
     * specified {@code probs} array satisfies
     * {@code (p <= 0.0 || p >= 1.0 || Double.isNaN(p))}
     * @throws IndexOutOfBoundsException if {@code probs.length < quantiles.length}
     * @throws NullPointerException if {@code probs == null || quantiles == null}
     */
    public void fwdQuantiles(int hap1, int hap2, double ibdStartMorgans,
            int focusPos, double[] probs, int[] quantiles) {

        double focusMorgans = fwdData.posToMorgans(focusPos);
        setCDF(fwdData, hap1, hap2, ibdStartMorgans, focusPos, focusMorgans);
        for (int j=0; j<quantiles.length; ++j) {
            quantiles[j] = quantile(fwdData, ibdStartMorgans, focusPos,
                    focusMorgans, probs[j]);
        }
    }

    /**
     * Estimates the specified quantiles of the backward endpoint of the
     * specified IBD segment.
     * @param hap1 the first IBD haplotype
     * @param hap2 the second IBD haplotype
     * @param ibdEndMorgans an estimated genetic position in Morgans of the
     * end of the IBD segment
     * @param focusPos the base pair coordinate of the focus position in
     * the IBD segment
     * @param probs a list of probabilities
     * @param quantiles an array which will contain the end-point quantiles
     * corresponding to the specified probabilities
     * @throws IllegalArgumentException if
     * {@code Double.isNaN(ibdEndMorgans) == true}
     * @throws IllegalArgumentException if the focus position is greater than or
     * equal to the estimated base pair position of the specified end of the IBD
     * segment
     * @throws IllegalArgumentException if any element {@code p} of the
     * specified {@code probs} array satisfies
     * {@code (p <= 0.0 || p >= 1.0 || Double.isNaN(p))}
     * @throws IndexOutOfBoundsException if
     * {@code hap1 < 0 || hap1 >= this.data().fwdGT().nHaps()}
     * @throws IndexOutOfBoundsException if
     * {@code hap2 < 0 || hap2 >= this.data().fwdGT().nHaps()}
     * @throws IndexOutOfBoundsException if {@code probs.length > quantiles.length}
     * @throws NullPointerException if {@code probs == null || quantiles == null}
     */
    public void bwdQuantiles(int hap1, int hap2, int focusPos,
            double ibdEndMorgans, double[] probs, int[] quantiles) {

        focusPos = -focusPos;
        ibdEndMorgans = -ibdEndMorgans;

        double focusMorgans = revData.posToMorgans(focusPos);
        setCDF(revData, hap1, hap2, ibdEndMorgans, focusPos, focusMorgans);
        for (int j=0; j<quantiles.length; ++j) {
            quantiles[j] = -quantile(revData, ibdEndMorgans, focusPos,
                    focusMorgans, probs[j]);
        }
    }

    /*
     * Stores the probability distribution of the position of the end
     * of an IBD segment. For markers {@code m} in the interval bounded by
     * {@code this.cdfStart+1} (inclusive) and {@code this.cdfEnd} (exclusive),
     * the probability that the IBD segment end is between markers
     * {@code (m-1)} and {@code m} is stored in the {@code m}-th element
     * of the specified {@code cdf} array. The probability that an
     * IBD segment end is between the focus position and {@code this.cdfStart}
     * is stored in the {@code this.cdf[cdfStart]}.
     */
    private void setCDF(Data data, int h1, int h2, double ibdStartMorgans,
            int focusPos, double focusMorgans) {
        this.cdfStart = data.nextMarker(focusPos);
        cdf[cdfStart-1] = 0.0;
        double constant = 1.0;
        double F1 = IbdEndsUtils.F((focusMorgans - ibdStartMorgans), ne);
        int start = cdfStart;
        int nextDiscord = IbdEndsUtils.fwdDiscord(data.gt, h1, h2, start);
        int minNextDiscordPos = data.pos(nextDiscord) + gc_bp;
        while (true) {
            cdfEnd = Math.min(nextDiscord+1, data.gt.nMarkers());
            for (int m=start; m<cdfEnd; ++m) {
                double F2 = IbdEndsUtils.F((data.morgans.get(m) - ibdStartMorgans), ne);
                cdf[m] = cdf[m-1] + (F2-F1)*data.ibsProbs.fwdProb(m, nextDiscord)*constant;
                F1 = F2;
            }
            if (finished(start)) {
                scale(cdf, cdfStart, cdfEnd, (1.0/cdf[cdfEnd-1]));
                return;
            }
            if (cdf[cdfEnd-1]>1e50) {
                double factor = (1.0/cdf[cdfEnd-1]);
                scale(cdf, cdfStart, cdfEnd, factor);
                constant *= factor;
            }
            start = cdfEnd;
            nextDiscord = IbdEndsUtils.fwdDiscord(data.gt, h1, h2, start);
            int discordPos = data.pos(nextDiscord);
            double num = gc_err;
            if (discordPos >= minNextDiscordPos) {
                num = err;
                minNextDiscordPos = discordPos + gc_bp;
            }
            constant *= (num/data.ibsProbs.fwdProb(start, nextDiscord));
        }
    }

    private boolean finished(int lastEnd) {
        if (cdfEnd==cdf.length) {
            return true;
        }
        return (cdf[cdfEnd-1] - cdf[lastEnd-1])<(MIN_RATIO*cdf[cdfEnd-1]);
    }

    private static void scale(double[] da, int start, int end, double factor) {
        for (int j=start; j<end; ++j) {
            da[j] *= factor;
        }
    }

    private int quantile(Data data, double ibdStartMorgans, int focusPos,
            double focusMorgans, double p) {
        if (p<=0d || p>=1d || Double.isNaN(p)==true) {
            throw new IllegalArgumentException(String.valueOf(p));
        }
        int index = Arrays.binarySearch(cdf, cdfStart, cdfEnd, p);
        if (index<0) {
            index = -index - 1;
        }
        double p1 = cdf[index-1];   // assert cdf[cdfStart-1]==0.0;
        double p2 = cdf[index];
        assert p1<=p && p<=p2;

        double x1 = (index==cdfStart) ? focusMorgans : data.morgans.get(index-1);
        double x2 = data.morgans.get(index);

        double F1 = IbdEndsUtils.F((x1 - ibdStartMorgans), ne);
        double F2 = IbdEndsUtils.F((x2 - ibdStartMorgans), ne);
        double pp = F1 + ((p-p1)/(p2-p1))*(F2-F1);
        double x = ibdStartMorgans + IbdEndsUtils.invF(pp, ne);
        assert x1<=x && x<=x2;
        double delta = (x-x1)/(x2-x1);

        // minimum quantile needs to be focusPos+1 to avoid division by 0
        int y1 = index==cdfStart ? focusPos+1 : data.pos.get(index-1);
        int y2 = data.pos.get(index);
        int y = (int) Math.rint(y1 + delta*(y2-y1));
        assert y1<=y && y<=y2;
        return y;
    }

    private static class Data {

        private final GT gt;
        private final WrappedIntArray pos;
        private final DoubleArray morgans;
        private final IbsLengthProbs ibsProbs;

        private Data(IbdEndsData data, boolean fwd) {
            if (fwd) {
                this.gt = data.fwdGT();
                this.pos = data.fwdPos();
                this.morgans = data.fwdMorgans();
                this.ibsProbs = data.fwdIbsProbs();
            }
            else {
                this.gt = data.revGT();
                this.pos = data.revPos();
                this.morgans = data.revMorgans();
                this.ibsProbs = data.revIbsProbs();
            }
        }

        private int pos(int marker) {
            return marker==gt.nMarkers() ? Integer.MAX_VALUE : pos.get(marker);
        }

        private double posToMorgans(int position) {
            return IbdEndsUtils.morganPos(pos, morgans, position);
        }

        private int nextMarker(int position) {
            int insPt = pos.binarySearch(position);
            return insPt<0 ? -insPt-1 : insPt+1;
        }
    }
}
