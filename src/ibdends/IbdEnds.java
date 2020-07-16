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

/**
 * <p>Class {@code IbdEnds} estimates IBD end-points.</p>
 *
 * <p>Instances of class {@code IbdEnds} are not thread-safe.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class IbdEnds {

    private final IbdEndsData data;
    private final QuantileEstimator quantEstimator;
    private final int maxIts;
    private final int maxItsM2;
    private final double maxRelDiff;
    private final boolean fixFocus;

    private int h1;
    private int h2;
    private int startPos;
    private int endPos;
    private int focusPos;
    private double startMorgans;
    private double endMorgans;
    private double focusMorgans;

    /**
     * Constructs an {@code IbdEnds} instance for the specified data.
     * @param data the analysis input data
     * @throws NullPointerException if {@code data == null}
     */
    public IbdEnds(IbdEndsData data) {
        IbdEndsPar par = data.par();
        this.data = data;
        this.quantEstimator = new QuantileEstimator(data);
        this.maxIts = par.max_its()<<1; // doubled since there are two ends
        this.maxItsM2 = this.maxIts - 2;
        this.maxRelDiff = par.max_diff();
        this.fixFocus = par.fix_focus();
    }

    /**
     * Returns the input data.
     * @return the input data
     */
    public IbdEndsData data() {
        return data;
    }

    /**
     * Calculates and stores the specified quantiles for the forward and
     * backward IBD segment end point distributions. The method returns the
     * focus base position from which the distances to the end points
     * are estimated. The focus base position is obtained from the specified
     * IBD segment and from the estimated first quantiles of the
     * {@code fwdProbs} and {@code bwdProbs} arrays.
     * @param seg an IBD segment
     * @param fwdProbs an array of probability values specifying quantiles
     * of the forward IBD segment end point distribution
     * @param fwdQuantiles an array in which quantiles of the forward IBD
     * segment end point distribution will be stored
     * @param bwdProbs an array of probability values specifying quantiles
     * of the backward IBD segment end point distribution
     * @param bwdQuantiles an array in which quantiles of the backward IBD
     * segment end point distribution will be stored
     * @return the focal base position from which distances to the segment
     * end points are estimated
     *
     * @throws IndexOutOfBoundsException if
     * {@code seg.hap1() < 0 || seg.hap1() >= this.data().fwdGT().nHaps()}
     * @throws IndexOutOfBoundsException if
     * {@code seg.hap2() < 0 || seg.hap2() >= this.data().fwdGT().nHaps()}
     * @throws IllegalArgumentException if any element {@code p} of the
     * specified {@code fwdProbs} array satisfies
     * {@code (p <= 0.0 || p >= 1.0 || Double.isNaN(p))}
     * @throws IllegalArgumentException if any element {@code p} of the
     * specified {@code bwdProbs} array satisfies
     * {@code (p <= 0.0 || p >= 1.0 || Double.isNaN(p))}
     * @throws IndexOutOfBoundsException if {@code fwdProbs.length < 1}
     * @throws IndexOutOfBoundsException if {@code bwdProbs.length < 1}
     * @throws IndexOutOfBoundsException if
     * {@code fwdProbs.length > fwdQuantiles.length}
     * @throws IndexOutOfBoundsException if
     * {@code bwdProbs.length > bwdQuantiles.length}
     * @throws NullPointerException if any parameter is {@code null}
     */
    public int getEnds(SharedSegment seg, double[] fwdProbs, int[] fwdQuantiles,
            double[] bwdProbs, int[] bwdQuantiles) {
        checkLengths(fwdProbs, bwdProbs);
        setFields(seg);
        int noEndUpdateCnt = 0; // 2 if fwd and bwd quantiles computed from same focus
        boolean updated;
        for (int it=0; noEndUpdateCnt<2; ++it) {
            if ((it & 1)==0) {
                quantEstimator.fwdQuantiles(h1, h2, startMorgans, focusPos,
                        fwdProbs, fwdQuantiles);
                updated = updateEndPos(seg, it, fwdQuantiles[0]);
            }
            else {
                quantEstimator.bwdQuantiles(h1, h2, focusPos, endMorgans,
                        bwdProbs, bwdQuantiles);
                updated = updateStartPos(seg, it, bwdQuantiles[0]);
            }
            noEndUpdateCnt = updated ? 0 : noEndUpdateCnt + 1;
        }
        return focusPos;
    }

    private void checkLengths(double[] fwdProbs, double[] bwdProbs) {
        if (fwdProbs.length<1) {
            throw new IllegalArgumentException(String.valueOf(fwdProbs.length));
        }
        if (bwdProbs.length<1) {
            throw new IllegalArgumentException(String.valueOf(bwdProbs.length));
        }
    }

    private void setFields(SharedSegment ss) {
        h1 = ss.hap1();
        h2 = ss.hap2();
        startPos = ss.start();
        endPos = ss.inclEnd();
        focusPos = (ss.start() + ss.inclEnd()) >>> 1;
        startMorgans = data.baseToMorgans(ss.start());
        endMorgans = data.baseToMorgans(ss.inclEnd());
        focusMorgans = data.baseToMorgans(focusPos);
    }

    private boolean updateEndPos(SharedSegment seg, int it, int newEndPos) {
        if (it>=maxItsM2) {
            return false;
        }
        if (newEndPos > seg.inclEnd()) {
            newEndPos = seg.inclEnd();
        }
        double newEndMorgans = data.baseToMorgans(newEndPos);
        if (noEndPointChange(focusMorgans, endMorgans, newEndMorgans)) {
            return false;
        }
        else {
            endPos = newEndPos;
            endMorgans = newEndMorgans;
            if (fixFocus==false) {
                focusPos = (startPos + endPos) >>> 1;
                focusMorgans = data.baseToMorgans(focusPos);
            }
            return true;
        }
    }

    private boolean updateStartPos(SharedSegment seg, int it, int newStartPos) {
        if (it>=maxItsM2) {
            return false;
        }
        if (newStartPos < seg.start()) {
            newStartPos = seg.start();
        }
        double newStartMorgans = data.baseToMorgans(newStartPos);
        if (noEndPointChange(focusMorgans, startMorgans, newStartMorgans)) {
            return false;
        }
        else {
            startPos = newStartPos;
            startMorgans = newStartMorgans;
            if (fixFocus==false) {
                focusPos = (startPos + endPos) >>> 1;
                focusMorgans = data.baseToMorgans(focusPos);
            }
            return true;
        }
    }

    private boolean noEndPointChange(double focusMorgans, double oldMorgans,
            double newMorgans) {
        double beforeLength = oldMorgans - focusMorgans;
        double afterLength = newMorgans - focusMorgans;
        if (beforeLength==0.0) {
            return false;
        }
        else {
            double absRelDiff = Math.abs((afterLength - beforeLength)/beforeLength);
            return absRelDiff < maxRelDiff;
        }
    }
}
