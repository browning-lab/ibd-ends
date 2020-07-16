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

import blbutil.Const;

/**
 * <p>Class {@code SharedSegment} represents a shared chromosome segment for
 * a pair of haplotypes.</p>
 *
 * <p>Instances of class {@code SharedSegment} are immutable.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class SharedSegment implements Comparable<SharedSegment> {

    private final int hap1;
    private final int hap2;
    private final int start;
    private final int inclEnd;

    /**
     * Constructs a new {@code SharedSegment} instance for the specified data
     * @param hap1 the first haplotype index in the list of all sample haplotypes
     * @param hap2 the second haplotype index in the list of all sample haplotypes
     * @param start the starting coordinate for the shared segment (inclusive)
     * @param inclEnd the ending marker coordinate for the shared segment (inclusive)
     *
     * @throws IllegalArgumentException if
     * {@code  hap1 < 0 || hap2 < 0 || start > inclEnd}
     */
    public SharedSegment(int hap1, int hap2, int start, int inclEnd) {
        checkArguments(hap1, hap2, start, inclEnd);
        this.hap1 = hap1;
        this.hap2 = hap2;
        this.start = start;
        this.inclEnd = inclEnd;
    }

    private void checkArguments(int hap1, int hap2, int start, int end) {
        if (hap1<0) {
            throw new IllegalArgumentException(String.valueOf(hap1));
        }
        if (hap2<0) {
            throw new IllegalArgumentException(String.valueOf(hap2));
        }
        else if (end<start){
            String s = "start > end: start=" + start + " end=" + end;
            throw new IllegalArgumentException(s);
        }
    }

    /**
     * Compares the specified object with this {@code SharedSegment} for
     * equality.  Returns {@code  true} if the specified object is a
     * {@code SharedSegment} instance whose {@code hap1()}, {@code hap2()},
     * {@code start()}, and {@code end()} methods return the same values
     * as the corresponding methods for {@code  this}, and returns {@code  false}
     * otherwise.
     * @param obj the object to be compared for equality with this
     * {@code SharedSegment}.
     *
     * @return {@code  true} if the specified object is a {@code SharedSegment}
     * instance whose {@code hap1()}, {@code hap2()}, {@code start()}, and
     * {@code inclEnd()} methods return the same values as the corresponding
     * methods for {@code  this}
     */
    @Override
    public boolean equals(Object obj) {
        if (this==obj) {
            return true;
        }
        if ((obj instanceof SharedSegment)==false) {
            return false;
        }
        SharedSegment other = (SharedSegment) obj;
        if (this.hap1!=other.hap1) {
            return false;
        }
        if (this.hap2!=other.hap2) {
            return false;
        }
        if (this.start!=other.start) {
            return false;
        }
        return this.inclEnd == other.inclEnd;
    }

     /**
     * Returns a hash code value for this object.
     * @return a hash code value for this object
     */
    @Override
    public int hashCode() {
        int hash = 3;
        hash = 43 * hash + this.hap1;
        hash = 43 * hash + this.hap2;
        hash = 43 * hash + this.start;
        hash = 43 * hash + this.inclEnd;
        return hash;
    }

    /**
     * Returns -1, 0, or 1 depending on whether this {@code SharedSegment} is
     * less than, equal to, or greater than the specified {@code SharedSegment}.
     * Two shared segments are ordered by their first haplotype index,
     * then by their second haplotype index, then by their starting
     * marker index, and finally by their ending marker index. This ordering
     * is consistent with equals.
     *
     * @param o the IBD segment to be compared.
     * @return -1, 0, or 1 depending on whether this {@code SharedSegment} is
     * less than, equal to, or greater than the specified {@code SharedSegment}.
     */
    @Override
    public int compareTo(SharedSegment o) {
        if (this.hap1 != o.hap1) {
            return (this.hap1 < o.hap1) ? -1 : 1;
        }
        if (this.hap2 != o.hap2) {
            return (this.hap2 < o.hap2) ? -1 : 1;
        }
        if (this.start != o.start) {
            return (this.start < o.start) ? -1 : 1;
        }
        if (this.inclEnd != o.inclEnd) {
            return (this.inclEnd < o.inclEnd) ? -1 : 1;
        }
        return 0;
    }

    /**
     * Returns the first haplotype index in the list of all sample haplotypes
     * @return the first haplotype index in the list of all sample haplotypes
     */
    public int hap1() {
        return hap1;
    }

    /**
     * Returns the second haplotype index in the list of all sample haplotypes
     * @return the second haplotype index in the list of all sample haplotypes
     */
    public int hap2() {
        return hap2;
    }

    /**
     * Returns the starting coordinate (inclusive)
     * @return the starting coordinate (inclusive)
     */
    public int start() {
        return start;
    }

    /**
     * Returns the ending coordinate (inclusive)
     * @return the ending coordinate (inclusive)
     */
    public int inclEnd() {
        return inclEnd;
    }

    /**
     * Returns a string representation of {@code this}.  The exact
     * details of the representation are unspecified and subject to change.
     * @return a string representation of {@code this}
     */
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append(hap1);
        sb.append(Const.tab);
        sb.append(hap2);
        sb.append(Const.tab);
        sb.append(start);
        sb.append(Const.tab);
        sb.append(inclEnd);
        return sb.toString();
    }
}
