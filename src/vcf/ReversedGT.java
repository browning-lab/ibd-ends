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
package vcf;

/**
 * <p>Class {@code ReversedGT} is a wrapper for a {@code GT}
 * instance that reversed the marker order.
 * </p>
 * <p>Instances of class {@code ReversedGT} are immutable.
 * </p>
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class ReversedGT implements GT {

    private final GT gt;
    private final int lastMarker;

    /**
     * Constructs a new {@code ReversedGT} instance from the specified
     * data
     * @param gt the genotypes to be wrapped
     * @throws IllegalArgumentException if {@code gt.isReversed() == true}
     * @throws NullPointerException if {@code gt == null}
     */
    public ReversedGT(GT gt) {
        if (gt.isReversed()) {
            throw new IllegalArgumentException("gt.isReversed()");
        }
        this.lastMarker = gt.nMarkers() - 1;
        this.gt = gt;
    }

    @Override
    public boolean isReversed() {
        return true;
    }

    @Override
    public int nMarkers() {
        return gt.nMarkers();
    }

    @Override
    public Marker marker(int marker) {
        return gt.marker(lastMarker - marker);
    }

    @Override
    public Markers markers() {
        return gt.markers();
    }

    @Override
    public int nHaps() {
        return gt.nHaps();
    }

    @Override
    public int nSamples() {
        return gt.nSamples();
    }

    @Override
    public Samples samples() {
        return gt.samples();
    }

    @Override
    public boolean isPhased() {
        return gt.isPhased();
    }

    @Override
    public int allele1(int marker, int sample) {
        return gt.allele1(lastMarker - marker, sample);
    }

    @Override
    public int allele2(int marker, int sample) {
        return gt.allele2(lastMarker - marker, sample);
    }


    @Override
    public int allele(int marker, int hap) {
        return gt.allele(lastMarker - marker, hap);
    }

    @Override
    public GT restrict(Markers markers, int[] indices) {
        throw new UnsupportedOperationException();
    }

    @Override
    public String toString() {
        return ReversedGT.class.toString() + " : " + gt.toString();
    }
}
