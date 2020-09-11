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
import blbutil.FloatArray;
import blbutil.StringUtil;
import blbutil.Utilities;
import blbutil.Validate;
import java.io.File;
import java.util.Map;

/**
 * <p>Class {@code IbdEndsPar} represents the analysis parameters
 * for an ibd-ends analysis.</p>
 *
 * <p>Class {@code IbdEndsPar} is immutable.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class IbdEndsPar {

    private final String[] args;

    // data input/output parameters
    private final File gt;
    private final ChromInterval chromInt;
    private final File ibd;
    private final File map;
    private final String out;
    private final File excludesamples;

    // algorithm parameters
    private final FloatArray quantiles;
    private final int nsamples;
    private final int nthreads;
    private final float err;
    private final boolean estimate_err;
    private final float gc_err;
    private final int gc_bp;
    private final float min_maf;
    private final long seed;

    private static final String DEF_QUANTILES = "0.5";
    private static final int DEF_NSAMPLES = 0;
    private static final int DEF_NTHREADS = Runtime.getRuntime().availableProcessors();
    private static final float DEF_ERR = 0.0005f;
    private static final boolean DEF_ESTIMATE_ERR = true;
    private static final float DEF_GC_ERR = 0.1f;
    private static final int DEF_GC_BP = 1000;
    private static final float DEF_MIN_MAF = 0.001f;
    private static final int DEF_SEED = -99999;

    // undocumented parameters
    private final float ne;              // effective population size
    private final int local_haps;        // number haplotypes used for local IBS destribution
    private final int global_pos;        // number of sampled positions for global IBS length distribution
    private final int global_segments;   // number of segments sampled per position for global IBS distribution
    private final float global_quantile; // quantile of global IBS destribution used for filtering positions
    private final float global_factor;   // factor to multiply by median global-quantile value
    private final float max_local_cdf;   // max CDF probability for using local IBS length distribution
    private final int max_its;           // max number of iterative updates to end points
    private final float end_morgans;     // morgans from first/last marker to hypothetical next discordance

    private static final float DEF_NE = 10000f;
    private static final int DEF_LOCAL_HAPS = 10000;
    private static final int DEF_GLOBAL_POS = 1000;
    private static final int DEF_GLOBAL_SEGMENTS = 2000;
    private static final float DEF_GLOBAL_QUANTILE = 0.9f;
    private static final float DEF_GLOBAL_FACTOR = 3.0f;
    private static final float DEF_MAX_LOCAL_CDF = 0.999f;
    private static final int DEF_MAX_ITS = 10;
    private static final float DEF_END_MORGANS = 1.0f;

    private final boolean fix_focus;     // if true, do not iteratively update focus
    private final float length_quantile; // quantile used to estimate pre-focus IBD length
    private final float max_diff;        // max permitted relative difference between IBD length before focus
                                         //   and distance from length-quantile of start distribution to focus

    private static final boolean DEF_FIX_FOCUS = false;
    private static final float DEF_LENGTH_QUANTILE = 0.05f;
    private static final float DEF_MAX_DIFF = 0.1f;

    private final boolean debug;         // if true, print debug output

    private static final boolean DEF_DEBUG = false;

    /**
     * Constructs an {@code IbdEndPar} object that represents the
     * analysis parameters for an ibd-ends analysis.  See the
     * {@code usage()} method for a description of the command line
     * parameters.
     *
     * @param args the command line arguments
     * @throws IllegalArgumentException if the command line arguments
     * are incorrectly specified
    */
    public IbdEndsPar(String[] args) {
        float MIN_PROP = Float.MIN_VALUE;
        float MAX_PROP = Math.nextDown(1.0f);
        int IMAX = Integer.MAX_VALUE;
        int MAX_LOCAL_HAPS = 40000;     // prevents 4-byte signed int overflow
        float FMIN = Float.MIN_VALUE;
        float FMAX = Float.MAX_VALUE;
        long LMIN = Long.MIN_VALUE;
        long LMAX = Long.MAX_VALUE;

        Map<String, String> argsMap = Validate.argsToMap(args, '=');
        this.args = args.clone();

        // data input/output parameters
        gt = Validate.getFile(Validate.stringArg("gt", argsMap, true, null,
                null));
        chromInt = parseChromInt(Validate.stringArg("chrom", argsMap, false, null, null));
        map = Validate.getFile(Validate.stringArg("map", argsMap, true, null,
                null));
        ibd = Validate.getFile(Validate.stringArg("ibd", argsMap, true, null,
                null));
        out = Validate.stringArg("out", argsMap, true, null, null);
        excludesamples = Validate.getFile(Validate.stringArg("excludesamples",
                argsMap, false, null, null));

        // algorithm parameters
        quantiles = setQuantiles(Validate.stringArg("quantiles", argsMap, false,
                DEF_QUANTILES, null));
        nsamples = Validate.intArg("nsamples", argsMap, false, DEF_NSAMPLES, 0, IMAX);
        nthreads = Validate.intArg("nthreads", argsMap, false, DEF_NTHREADS, 1, IMAX);
        err = Validate.floatArg("err", argsMap, false, DEF_ERR, MIN_PROP, MAX_PROP);
        estimate_err = Validate.booleanArg("estimate-err", argsMap, false, DEF_ESTIMATE_ERR);
        gc_err = Validate.floatArg("gc-err", argsMap, false, DEF_GC_ERR, MIN_PROP, MAX_PROP);
        gc_bp = Validate.intArg("gc-bp", argsMap, false, DEF_GC_BP, 0, IMAX);
        min_maf = Validate.floatArg("min-maf", argsMap, false, DEF_MIN_MAF, FMIN, 0.5f);
        seed = Validate.longArg("seed", argsMap, false, DEF_SEED, LMIN, LMAX);

        // undocumented parameters
        ne = Validate.floatArg("ne", argsMap, false, DEF_NE, 1, LMAX);
        local_haps = Validate.intArg("local-haps", argsMap, false, DEF_LOCAL_HAPS, 1, MAX_LOCAL_HAPS);
        global_pos = Validate.intArg("global-pos", argsMap, false, DEF_GLOBAL_POS, 1, IMAX);
        global_segments = Validate.intArg("global-segments", argsMap, false, DEF_GLOBAL_SEGMENTS, 1, IMAX);
        global_quantile = Validate.floatArg("global-quantile", argsMap, false,
                DEF_GLOBAL_QUANTILE, MIN_PROP, MAX_PROP);
        global_factor = Validate.floatArg("global-factor", argsMap, false,
                DEF_GLOBAL_FACTOR, FMIN, FMAX);
        max_local_cdf = Validate.floatArg("max-local-cdf", argsMap, false,
                DEF_MAX_LOCAL_CDF, MIN_PROP, MAX_PROP);
        max_its = Validate.intArg("max-its", argsMap, false, DEF_MAX_ITS, 1, IMAX);
        end_morgans = Validate.floatArg("end-morgans", argsMap, false,
                DEF_END_MORGANS, 0.0f, Float.MAX_VALUE);

        fix_focus = Validate.booleanArg("fix-focus", argsMap, false, DEF_FIX_FOCUS);
        length_quantile = Validate.floatArg("length-quantile", argsMap, false,
                DEF_LENGTH_QUANTILE, FMIN, FMAX);
        max_diff = Validate.floatArg("max-diff", argsMap, false, DEF_MAX_DIFF,
                MIN_PROP, MAX_PROP);
        debug = Validate.booleanArg("debug", argsMap, false, DEF_DEBUG);

        Validate.confirmEmptyMap(argsMap);
    }

    private static ChromInterval parseChromInt(String str) {
        ChromInterval chromInt = ChromInterval.parse(str);
        if (str!=null && str.length()>0 && chromInt==null) {
            throw new IllegalArgumentException("Invalid chrom parameter: " + str);
        }
        return chromInt;
    }

    private FloatArray setQuantiles(String stringQuantiles) {
        String[] sa = StringUtil.getFields(stringQuantiles, Const.comma);
        float[] fa = new float[sa.length];
        for (int j=0; j<sa.length; ++j) {
            try {
                fa[j] = Float.parseFloat(sa[j]);
                if (Float.isFinite(fa[j])==false || fa[j]<=0f || fa[j]>=1f) {
                    quantileError(sa[j]);
                }
            }
            catch(NumberFormatException e) {
                quantileError(sa[j]);
            }
        }
        return new FloatArray(fa);
    }

    private void quantileError(String s) {
        String msg = "ERROR: invalid quantile: " + s;
        throwErr(msg);
    }

    private void throwErr(String msg) {
        System.out.println(usage());
        Utilities.exit(new IllegalArgumentException(msg));
    }

    /**
     * Returns the command line arguments.
     * @return the command line arguments
     */
    public String[] args() {
        return args.clone();
    }

    /**
     * Returns a string describing the command line arguments.
     * The format of the returned string is unspecified and subject to change.
     *
     * @return a string describing the command line arguments.
     */
    public static String usage() {
        String nl = Const.nl;
        return "Usage: " + IbdEndsMain.COMMAND + " [arguments]" + nl
                + nl
                + "Data Parameters: " + nl
                + "  gt=<VCF file with GT field>                           (required)" + nl
                + "  ibd=<hap-ibd output file from specified VCF file>     (required)" + nl
                + "  map=<PLINK map file with cM units>                    (required)" + nl
                + "  chrom=<chromosome to be analyzed>                     (default: 1st VCF chrom)" + nl
                + "  out=<output file prefix>                              (required)" + nl
                + "  excludesamples=<excluded samples file>                (optional)" + nl + nl

                + "Algorithm Parameters: " + nl
                + "  quantiles=<comma-separated list of quantiles>         (default: " + DEF_QUANTILES + ")" + nl
                + "  nsamples=<number of sampled endpoints>                (default: " + DEF_NSAMPLES + ")" + nl
                + "  nthreads=<number of computational threads>            (default: all CPU cores)" + nl
                + "  err=<IBD allele mismatch probability>                 (default: " + DEF_ERR + ")" + nl
                + "  estimate-err=<true/false>                             (default: " + DEF_ESTIMATE_ERR + ")" + nl
                + "  gc-err=<gene conversion allele mismatch probability>  (default: " + DEF_GC_ERR + ")" + nl
                + "  gc-bp=<max gene conversion base pair length>          (default: " + DEF_GC_BP + ")" + nl
                + "  min-maf=<min permitted minor allele frequency>        (default: " + DEF_MIN_MAF + ")" + nl
                + "  seed=<random seed>                                    (default=" + DEF_SEED + ")" + nl + nl;
    }

    // data input/output parameters

    /**
     * Returns the gt parameter.
     * @return the gt parameter
     */
    public File gt() {
        return gt;
    }

    /**
     * Returns the chromosome interval or {@code null} if no chrom
     * parameter was specified.
     *
     * @return the chromosome interval or {@code null} if no chrom
     * parameter was specified.
     */
    public ChromInterval chromInt() {
        return chromInt;
    }

    /**
     * Returns the map parameter.
     * @return the map parameter
     */
    public File map() {
        return map;
    }

    /**
     * Returns the ibd parameter.
     * @return the ibd parameter
     */
    public File ibd() {
        return ibd;
    }

    /**
     * Returns the out parameter.
     * @return the out parameter
     */
    public String out() {
        return out;
    }

    /**
     * Returns the excludesamples parameter or {@code null}
     * if no excludesamples parameter was specified.
     *
     * @return the excludesamples parameter or {@code null}
     * if no excludesamples parameter was specified
     */
    public File excludesamples() {
        return excludesamples;
    }

    // algorithm parameterst

    /**
     * Returns the quantiles parameter
     * @return the quantiles parameter
     */
    public FloatArray quantiles() {
        return quantiles;
    }

    /**
     * Returns the nsamples parameter.
     * @return the nsamples parameter
     */
    public int nsamples() {
        return nsamples;
    }

    /**
     * Returns the nthreads parameter.
     * @return the nthreads parameter
     */
    public int nthreads() {
        return nthreads;
    }

    /**
     * Returns the err parameter.
     * @return the err parameter
     */
    public float err() {
        return err;
    }

    /**
     * Returns the estimate-err parameter.
     * @return the estimate-err parameter
     */
    public boolean estimate_err() {
        return estimate_err;
    }

    /**
     * Returns the gc-err parameter.
     * @return the gc-err parameter
     */
    public float gc_err() {
        return gc_err;
    }

    /**
     * Returns the gc-bp parameter.
     * @return the gc-bp parameter
     */
    public int gc_bp() {
        return gc_bp;
    }

    /**
     * Returns the min-maf parameter.
     * @return the min-maf parameter
     */
    public float min_maf() {
        return min_maf;
    }

    /**
     * Returns the seed parameter.
     * @return the seed parameter
     */
    public long seed() {
        return seed;
    }

    // undocumented parameters


    /**
     * Returns the ne parameter.
     * @return the ne parameter
     */
    public float ne() {
        return ne;
    }

    /**
     * Returns the local-haps parameter.
     * @return the local-haps parameter
     */
    public int local_haps() {
        return local_haps;
    }

    /**
     * Returns the global-pos parameter.
     * @return the global-pos parameter
     */
    public int global_pos() {
        return global_pos;
    }

    /**
     * Returns the global-segments parameter.
     * @return the global-segments parameter
     */
    public int global_segments() {
        return global_segments;
    }

    /**
     * Returns the global-quantile parameter.
     * @return the global-quantile parameter
     */
    public float global_quantile() {
        return global_quantile;
    }

    /**
     * Returns the global-factor parameter.
     * @return the global-factor parameter
     */
    public float global_factor() {
        return global_factor;
    }

    /**
     * Returns the max-local-cdf parameter.
     * @return the max-local-cdf parameter
     */
    public float max_local_cdf() {
        return max_local_cdf;
    }

    /**
     * Returns the max-its parameter.
     * @return the max-its parameter
     */
    public int max_its() {
        return max_its;
    }

    /**
     * Returns the end-morgans parameter.
     * @return the end-morgans parameter
     */
    public float end_morgans() {
        return end_morgans;
    }

    /**
     * Returns the fix-focus parameter.
     * @return the fix-focus parameter
     */
    public boolean fix_focus() {
        return fix_focus;
    }

    /**
     * Returns the length-quantile parameter.
     * @return the length-quantile parameter
     */
    public float length_quantile() {
        return length_quantile;
    }

    /**
     * Returns the max_diff parameter.
     * @return the max_diff parameter
     */
    public float max_diff() {
        return max_diff;
    }

    /**
     * Returns the debug parameter.
     * @return the debug parameter
     */
    public boolean debug() {
        return debug;
    }
}