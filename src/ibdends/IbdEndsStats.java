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

import java.util.concurrent.atomic.LongAdder;

/**
 * <p>Class {@code IbdEndsStats} stores statistics from an IBD segment
 * endpoint analysis.</p>
 *
 * <p>Instances of class {@code IbdEndsStats} are thread-safe.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class IbdEndsStats {

    private final LongAdder nMarkers;
    private final LongAdder nSamples;
    private final LongAdder ibdSegmentCnt;
    private final LongAdder discordCnt;
    private final LongAdder totalCnt;

    /**
     * Constructs a new {@code IbdEndStats} instance.
     */
    public IbdEndsStats(){
        this.nMarkers = new LongAdder();
        this.nSamples = new LongAdder();
        this.ibdSegmentCnt = new LongAdder();
        this.discordCnt = new LongAdder();
        this.totalCnt = new LongAdder();
    }

   /**
     * Adds the specified number to the marker count. Invocation in the absence
     * of concurrent updates returns an accurate result, but concurrent updates
     * that occur while the count is being calculated might not be incorporated.
     * @param cnt the number to be added
     */
    public void addMarkers(int cnt) {
        nMarkers.add(cnt);
    }

    /**
     * Returns the marker count.
     * @return the marker count
     */
    public long nMarkers() {
        return nMarkers.sum();
    }

   /**
     * Adds the specified number to the sample count.
     * @param cnt the number to be added
     */
    public void addSamples(int cnt) {
        nSamples.add(cnt);
    }

    /**
     * Returns the sample count. Invocation in the absence of concurrent
     * updates returns an accurate result, but concurrent updates that occur
     * while the count is being calculated might not be incorporated.
     * @return the sample count
     */
    public long nSamples() {
        return nSamples.sum();
    }

    /**
     * Adds the specified number to the IBD segment count.
     * @param cnt the number to be added
     */
    public void addIbdSegmentCnt(long cnt) {
        ibdSegmentCnt.add(cnt);
    }

    /**
     * Increment the IBD segment count.
     */
    public void incrementIbdSegmentCnt() {
        ibdSegmentCnt.increment();
    }

    /**
     * Returns the IBD segment count. Invocation in the absence of concurrent
     * updates returns an accurate result, but concurrent updates that occur
     * while the count is being calculated might not be incorporated.
     * @return the IBD segment count
     */
    public long ibdSegmentCnt() {
        return ibdSegmentCnt.sum();
    }

    /**
     * Update the IBD segment allele discordant rate.
     * @param discordant the number of markers with discordant alleles
     * @param total the total number of markers examined
     * @throws IllegalArgumentException if {@code discordant > total}
     */
    public void updateDiscordRate(int discordant, int total) {
        if (discordant > total) {
            throw new IllegalArgumentException(discordant + ">" + total);
        }
        discordCnt.add(discordant);
        totalCnt.add(total);
    }

    /**
     * Returns the IBD segment allele discord rate. Invocation in the
     * absence of concurrent updates returns an accurate result, but
     * concurrent updates that occur while the count is being calculated might
     * not be incorporated.
     * @return the IBD sgment allele discord rate
     */
    public float discordRate() {
        long num = discordCnt.sum();
        long den = totalCnt.sum();
        return (float) ((double) num / den);
    }
}
