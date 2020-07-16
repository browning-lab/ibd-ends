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

import blbutil.StringUtil;
import java.util.HashMap;
import java.util.Map;
import vcf.GT;
import vcf.Samples;

/**
 * <p>Class {@code HapIbdParser} parses string representation of IBD segments
 * output by the hap-ibd program.</p>
 *
 * <p>Instances of class {@code HapIbdParser} are immutable.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class HapIbdParser {

    /**
     * The segment returned by the {@code this.parse()} method if the
     * chromosome or sample identifier in a parsed IBD segment
     * is not present in {@code this.gt()}.
     */
    public static final SharedSegment NIL = new SharedSegment(0, 0, 0, 0);

    private final GT gt;
    private final String chrom;
    private final int minStart;
    private final int maxInclEnd;
    private final Map<String, Integer> sampleMap;

    /**
     * Constructs a new {@code HapIbdParser} from the specified data.
     * @param gt the phase, non-missing genotype data
     * @throws NullPointerException if {@code gt==null}
     */
    public HapIbdParser(GT gt) {
        this.gt = gt;
        this.chrom = gt.markers().marker(0).chrom();
        this.minStart = gt.marker(0).pos();
        this.maxInclEnd = gt.marker(gt.nMarkers()-1).pos();
        this.sampleMap = sampleMap(gt.samples());
    }

    private Map<String, Integer> sampleMap(Samples samples) {
        String[] ids = samples.ids();
        Map<String, Integer> map = new HashMap<>();
        for (int j=0, n=samples.nSamples(); j<n; ++j) {
            map.put(ids[j], j);
        }
        return map;
    }

    /**
     * Returns a {@code SharedSegment} object corresponding to the specified
     * IBD segment string, or {@code HapIbdParser.NIL} if a sample or chromosome
     * identifier in the IBD segment string is not present in {@code this.gt()}.
     * The specified IBD segment string must contain at least
     * 7 white-space delimited fields.  The first 7 fields must be:<pre>
     * 1) the first sample identifier
     * 2) the haplotype index (1 or 2) for the first sample
     * 3) the second sample identifier
     * 4) the haplotype index (1 or 2) for the second sample
     * 5) the chromosome
     * 6) the starting marker position (inclusive)
     * 7) the ending marker position (inclusive)</pre>
     * If the starting marker position is less than the first marker position
     * in {@code this.gt()}, the starting marker position is set to the
     * first marker position in {@code this.gt()}.  Similarly, if the
     * ending marker position is greater than the last marker position
     * in {@code this.gt()}, the ending marker position is set to the
     * last marker position in {@code this.gt()}.
     *
     * @param line the string describing an IBD segment that will be parsed
     *
     * @return a {@code SharedSegment} object corresponding to the specified
     * IBD segment string, or {@code HapIbdParser.NIL} if a sample or chromosome
     * identifier in the IBD segment string is not present in {@code this.gt()}
     *
     * @throws NullPointerException if {@code line==null}
     * @throws IllegalArgumentException if the specified line does
     * not contain at least 7 white-space delimited fields
     * @throws IllegalArgumentException if a haplotype index is not "1" or "2"
     * @throws IllegalArgumentException if the starting marker position
     * is greater than the ending marker position
     * @throws IllegalArgumentException if the segment starting position
     * is greater than or equal to the last position in {@code this.gt()}
     * @throws IllegalArgumentException if the segment ending position
     * is less than or equal to the first position in {@code this.gt()}
     * @throws NumberFormatException if the segment starting position is not
     * a parsable integer
     * @throws NumberFormatException if the segment ending position is not
     * a parsable integer
     */
    public SharedSegment parse(String line) {
        String[] fields = StringUtil.getFields(line, 8);
        if (fields.length < 7) {
            String s = "IBD segment does not have at least 7"
                    + " white-space delimited fields [" + line + "]";
            throw new IllegalArgumentException(s);
        }
        int s1 = sampleIndex(fields[0]);
        int s2 = sampleIndex(fields[2]);
        if (s1 == -1 || s2 == -1 || chrom.equals(fields[4])==false) {
            return NIL;
        }
        else {
            int hap1 = (s1<<1) + parseHap(line, fields[1]) - 1;
            int hap2 = (s2<<1) + parseHap(line, fields[3]) - 1;
            int start = Integer.parseInt(fields[5]);
            int inclEnd = Integer.parseInt(fields[6]);
            if (start>inclEnd) {
                String s = "start > end. start=" + start + " end=" + inclEnd
                    + " [" + line + "]";
                throw new IllegalArgumentException(s);
            }
            if (inclEnd<=minStart || start>=maxInclEnd) {
                return NIL;
            }
            if (start<minStart) {
                start = minStart;
            }
            if (inclEnd>maxInclEnd) {
                inclEnd = maxInclEnd;
            }
            return new SharedSegment(hap1, hap2, start, inclEnd);
        }
    }

    private int sampleIndex(String sample) {
        Integer sIndex = sampleMap.get(sample);
        return sIndex==null ? -1 : sIndex;
    }

    private int parseHap(String line, String hapString) {
        int offset = hapString.charAt(0) - '0';
        if (hapString.length()!=1 || (offset!=1 && offset!=2)) {
            throw new IllegalArgumentException("haplotype index (" + hapString
                    + ") is not 1 or 2 [" + line + "]");
        }
        return offset;
    }

    /**
     * Returns the genotype data.
     * @return the genotype data
     */
    public GT gt() {
        return gt;
    }

    /**
     * Returns a string representation of {@code this}.  The details of
     * this representation are unspecified and subject to change.
     * @return a string representation of {@code this}
     */
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append("[");
        sb.append(HapIbdParser.class.toString());
        sb.append(": ");
        sb.append(sampleMap.size());
        sb.append(" samples]");
        return sb.toString();
    }
}
