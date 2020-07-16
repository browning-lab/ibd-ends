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

import beagleutil.ChromInterval;
import blbutil.Const;
import blbutil.DoubleArray;
import blbutil.FileIt;
import blbutil.Filter;
import blbutil.InputIt;
import blbutil.SampleFileIt;
import blbutil.Utilities;
import bref.Bref3It;
import ints.WrappedIntArray;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;
import vcf.FilterUtil;
import vcf.GT;
import vcf.GTRec;
import vcf.GeneticMap;
import vcf.IntervalVcfIt;
import vcf.Marker;
import vcf.PlinkGenMap;
import vcf.RefGT;
import vcf.RefGTRec;
import vcf.RefIt;
import vcf.ReversedGT;

/**
 * <p>Class {@code IbdEndsData} represents the immutable data for an ibd-ends
 * analysis.</p>
 *
 * <p>Instances of class {@code IbdEndsData} are immutable.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class IbdEndsData {

    private final double MIN_CM_DIST = 1e-6;

    private final IbdEndsPar par;
    private final PlinkGenMap genMap;
    private final RefGT fwdGT;
    private final GT revGT;
    private final WrappedIntArray fwdBasePos;
    private final WrappedIntArray revBasePos;
    private final DoubleArray fwdMorganPos;
    private final DoubleArray revMorganPos;
    private final HapIbdParser hapIbdParser;

    private final IbsLengthProbs fwdIbsProbs;
    private final IbsLengthProbs revIbsProbs;

    /**
     * Constructs a new {@code IbdEndData} instance from the specified data.
     * The constructor will terminate the Java Virtual Machine with
     * an error message if an I/O error or file format error is detected.
     *
     * @param par the analysis parameters
     *
     * @throws IllegalArgumentException if a format error is detected in a
     * line of an input data file
     * @throws NumberFormatException if an input integer field is not
     * a parsable integer
     * @throws NullPointerException if {@code par == null}
     */
    public IbdEndsData(IbdEndsPar par) {
        this.par = par;
        this.genMap = PlinkGenMap.fromPlinkMapFile(par.map());
        this.fwdGT = readChrom(par, genMap);
        this.revGT = new ReversedGT(fwdGT);

        this.fwdBasePos = IbdEndsData.fwdBasePos(fwdGT);
        this.revBasePos = revBasePos(fwdGT);
        double[] cmPos = GeneticMap.genPos(genMap, MIN_CM_DIST, fwdGT.markers());
        this.fwdMorganPos = fwdMorganPos(cmPos);
        this.revMorganPos = revMorganPos(cmPos);
        this.hapIbdParser = new HapIbdParser(fwdGT);

        GlobalIbsProbs gip = new GlobalIbsProbs(fwdGT, fwdMorganPos, par);
        IbsCounts fwdIbsCnts = new IbsCounts(par, fwdGT);
        IbsCounts revIbsCnts = fwdIbsCnts.reverse();
        this.fwdIbsProbs = new IbsLengthProbs(par, fwdMorganPos, gip, fwdIbsCnts);
        this.revIbsProbs = new IbsLengthProbs(par, revMorganPos, gip, revIbsCnts);
    }

    private static RefGT readChrom(IbdEndsPar par, PlinkGenMap map) {
        List<RefGTRec> list = new ArrayList<>(8192);
        try (SampleFileIt<RefGTRec> it = refIt(par)) {
            int nHaps = it.samples().nSamples()<<1;
            int minMac = Math.max(1, (int) Math.ceil(nHaps*par.min_maf()));
            RefGTRec next = advanceToFirstMapRec(it, map);
            int chromIndex = next.marker().chromIndex();
            int nMapPos = map.nMapPositions(chromIndex);
            int lastMapPos = map.index2BasePos(chromIndex, nMapPos-1);
            while (next!=null && next.marker().chromIndex()==chromIndex
                    && next.marker().pos()<=lastMapPos) {
                list.add(next);
                next = it.hasNext() ? it.next() : null;
            }
            applyMacFilter(list, minMac);
        } catch (Throwable t) {
            Utilities.exit(t);
        }
        return new RefGT(list.toArray(new RefGTRec[0]));
    }

    private static RefGTRec advanceToFirstMapRec(SampleFileIt<RefGTRec> it,
            PlinkGenMap map) {
        if (it.hasNext()==false) {
            Utilities.exit(Const.nl + "No VCF records found in " + it.file());
        }
        RefGTRec rec = it.next();
        int chromIndex = rec.marker().chromIndex();
        int firstMapPos = map.index2BasePos(chromIndex, 0);
        while (rec.marker().pos()<firstMapPos) {
            if (it.hasNext()==false) {
                Utilities.exit(Const.nl
                        + "No VCF records within chromosome interval and genetic map. Exiting program.");
            }
            rec = it.next();
        }
        return rec;
    }

    private static SampleFileIt<RefGTRec> refIt(IbdEndsPar par) {
        SampleFileIt<RefGTRec> refIt;
        if (par.gt().toString().endsWith(".bref3")) {
            if (par.excludesamples()!=null) {
                Utilities.exit("ERROR: The \"excludesamples\" parameter cannot "
                        + "be used if the \"gt\" parameter is in bref3 format");
            }
            refIt = new Bref3It(par.gt());
        }
        else {
            Filter<String> sFilter = FilterUtil.sampleFilter(par.excludesamples());
            Filter<Marker> mFilter = Filter.acceptAllFilter();
            FileIt<String> it0 = InputIt.fromGzipFile(par.gt());
            refIt = RefIt.create(it0, sFilter, mFilter);
        }
        ChromInterval chromInt = par.chromInt();
        if (chromInt != null) {
            refIt = new IntervalVcfIt<>(refIt, chromInt);
        }
        return refIt;
    }

    private static void applyMacFilter(List<RefGTRec> list, int minMac) {
        if (minMac>0) {
            RefGTRec[] filteredArray = list.stream()
                    .parallel()
                    .filter(rec -> mac(rec)>=minMac)
                    .toArray(RefGTRec[]::new);
            if (filteredArray.length<list.size()) {
                list.clear();
                list.addAll(Arrays.asList(filteredArray));
            }
        }
    }

    private static int mac(RefGTRec rec) {
        int[] alCnts = GTRec.alleleCounts(rec);
        Arrays.sort(alCnts);
        return alCnts.length<=1 ? 0 : alCnts[alCnts.length-2];
    }

    private static WrappedIntArray fwdBasePos(GT gt) {
        return new WrappedIntArray(
                IntStream.range(0, gt.nMarkers())
                .parallel()
                .map(j -> gt.marker(j).pos())
        );
    }

    private static WrappedIntArray revBasePos(GT gt) {
        int nMarkersM1 = gt.nMarkers() - 1;
        return new WrappedIntArray(
                IntStream.range(0, gt.nMarkers())
                .parallel()
                .map(j -> -gt.marker(nMarkersM1 - j).pos())
        );

    }

    private static DoubleArray fwdMorganPos(double[] cmPos) {
        DoubleStream ds = Arrays.stream(cmPos).parallel().map(d -> 0.01*d);
        return new DoubleArray(ds);
    }

    private static DoubleArray revMorganPos(double[] cmPos) {
        int sizeM1 = cmPos.length - 1;
        DoubleStream ds = IntStream.range(0, cmPos.length)
                .parallel()
                .mapToDouble(j -> -0.01*cmPos[sizeM1 - j]);
        return new DoubleArray(ds);

    }

    /**
     * Returns the genetic Morgan position of the specified base position.
     * The Morgan and base positions are in chromosome order.
     * @param pos a base position
     * @return the genetic Morgan position of the specified base position
     */
    public double baseToMorgans(int pos) {
        return IbdEndsUtils.morganPos(fwdBasePos, fwdMorganPos, pos);
    }

    /**
     * Returns the analysis parameters.
     * @return the analysis parameters
     */
    public IbdEndsPar par() {
        return par;
    }

    /**
     * Returns the genetic map.
     * @return the genMap
     */
    public GeneticMap genMap() {
        return genMap;
    }

    /**
     * Returns the input genotype data.
     * @return the input genotype data
     */
    public RefGT fwdGT() {
        return fwdGT;
    }

    /**
     * Returns the input genotype data.
     * @return the input genotype data
     */
    public GT revGT() {
        return revGT;
    }

    /**
     * Returns the base position of each marker in chromosome order.
     * @return the base position of each marker in chromosome order
     */
    public WrappedIntArray fwdPos() {
        return fwdBasePos;
    }

    /**
     * Returns the negated base position of each marker in reverse
     * chromosome order.
     *
     * @return the negated base position of each marker in reverse
     * chromosome order
     */
    public WrappedIntArray revPos() {
        return revBasePos;
    }

    /**
     * Returns the Morgan position of each marker in chromosome order.
     * @return the Morgan position of each marker in chromosome order
     */
    public DoubleArray fwdMorgans() {
        return fwdMorganPos;
    }

    /**
     * Returns the negated Morgan position of each marker in reverse
     * chromosome order.
     * @return the negated Morgan position of each marker in reverse
     * chromosome order
     */
    public DoubleArray revMorgans() {
        return revMorganPos;
    }

    /**
     * Returns a parser for shared haplotype segments reported by
     * the hap-ibd program.
     * @return a parser for shared haplotype segments reported by
     * the hap-ibd program
     */
    public HapIbdParser hapIbdParser() {
        return hapIbdParser;
    }

    /**
     * Returns the one-sided forward IBS length probabilities.
     * @return the one-sided forward IBS length probabilities
     */
    public IbsLengthProbs fwdIbsProbs() {
        return fwdIbsProbs;
    }

    /**
     * Returns the one-sided reverse IBS length probabilities.
     * @return the one-sided reverse IBS length probabilities
     */
    public IbsLengthProbs revIbsProbs() {
        return revIbsProbs;
    }
}
