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

import blbutil.SynchFileOutputStream;
import blbutil.BGZIPOutputStream;
import blbutil.BlockLineReader;
import blbutil.Const;
import blbutil.FileIt;
import blbutil.FileUtil;
import blbutil.FloatArray;
import blbutil.InputIt;
import blbutil.MultiThreadUtils;
import blbutil.Utilities;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.math.BigDecimal;
import java.math.MathContext;
import java.text.DecimalFormat;
import java.util.Locale;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

/**
 * <p>Class {@code IbdEndsMain} contains the main() method for the ibd-ends
 * program.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class IbdEndsMain {

    private static final String EXECUTABLE = "ibd-ends.jar";
    private static final String PROGRAM = EXECUTABLE + "  [ version 1.0, __REV__ ]";
    private static final String COPYRIGHT = "Copyright (C) 2020 Brian L. Browning";

    /**
     * The java command to run this version of the program.
     */
    static final String COMMAND = "java -jar " + EXECUTABLE;

    private static final String HELP_MESSAGE = "Enter \"" + COMMAND
            + "\" to print usage instructions";

    private static final int BLOCK_SIZE = 10_000;
    private static final float MAX_ERR_RATIO = 3.0f;
    private static final DecimalFormat DF = new DecimalFormat("0.00E0");

    private IbdEndsMain() {
        // private constructor to prevent instantiation
    }

    /**
     * Entry point to the ibd-ends program.  See the {@code hapibd.IbdEndsPar}
     * class for details of the program arguments.
     *
     * @param args the command line arguments
     */
    public static void main(String[] args) {
	Locale.setDefault(Locale.US);
        if (args.length==0 || args[0].toLowerCase().startsWith("help")) {
            System.out.print("Program: ");
            System.out.println(PROGRAM);
            System.out.println(COPYRIGHT);
            System.out.println();
            System.out.println(IbdEndsPar.usage());
            System.exit(0);
        }
        IbdEndsPar par = new IbdEndsPar(args);
        System.setProperty("java.util.concurrent.ForkJoinPool.common.parallelism",
                String.valueOf(par.nthreads()));
        checkOutputPrefix(par);
        IbdEndsStats stats = new IbdEndsStats();

        try (PrintWriter log = FileUtil.printWriter(new File(par.out() + ".log"))) {
            long t0 = System.nanoTime();
            Utilities.duoPrintln(log, startInfo(par));
            runAnalysis(par, stats);
            Utilities.duoPrintln(log, statistics(par, stats));
            Utilities.duoPrintln(log, endInfo(t0));
        }
    }

    /**
     * Runs an IBD endpoint analysis.
     * @param par the command line parameters
     * @throws NullPointerException if {@code par == null}
     */
    private static void runAnalysis(IbdEndsPar par, IbdEndsStats stats) {
        File ibdOutFile = new File(par.out() + ".ibd.gz");
        int nBlocks = par.nthreads() << 1;
        try (SynchFileOutputStream os = new SynchFileOutputStream(ibdOutFile);
                FileIt<String> ibdIt = InputIt.fromGzipFile(par.ibd());
                BlockLineReader ibdReader = BlockLineReader.create(ibdIt, nBlocks, BLOCK_SIZE)) {
            IbdEndsData data = new IbdEndsData(par);
            writeHeaderLine(par, ibdOutFile, os);
            int nThreads = par.nthreads();
            ExecutorService es = Executors.newFixedThreadPool(nThreads);
            for (int j=0; j<nThreads; ++j) {
                IbdEndsRunnable runnable = new IbdEndsRunnable(data, ibdReader,
                        stats, os);
                es.submit(runnable);
            }
            MultiThreadUtils.shutdownExecService(es);
            stats.addMarkers(data.fwdGT().nMarkers());
            stats.addSamples(data.fwdGT().nSamples());
            os.writeEmptyBgzipBlock();
        } catch (IOException ex) {
            Utilities.exit(ex);
        }
    }

    private static void checkOutputPrefix(IbdEndsPar par) {
        File outPrefix = new File(par.out());
        if (outPrefix.isDirectory()) {
            String s = "ERROR: \"out\" parameter cannot be a directory: \""
                    + par.out() + "\"";
            Utilities.exit(IbdEndsPar.usage() + s);
        }
        checkOutputFilename(par, ".ibd.gz");
        checkOutputFilename(par, ".log");
    }

    private static void checkOutputFilename(IbdEndsPar par, String outSuffix) {
        File file = new File(par.out() + outSuffix);
        if (file.equals(par.gt())
                || file.equals(par.map())
                || file.equals(par.ibd())
                || file.equals(par.excludesamples())) {
            String s = "ERROR: Output file same as input file: " + file;
            Utilities.exit(IbdEndsPar.usage() + s);
        }
    }

    private static void writeHeaderLine(IbdEndsPar par, File ibdFile,
            SynchFileOutputStream os) {
        ByteArrayOutputStream baos = new ByteArrayOutputStream();
        try (PrintWriter out = new PrintWriter(new BGZIPOutputStream(baos, false))) {
            out.println(IbdEndsRunnable.ibdOutHeader(par));
        }
        try {
            os.write(baos.toByteArray());
        } catch (IOException e) {
            Utilities.exit("Error writing to file: " + ibdFile, e);
        }
    }

    private static String startInfo(IbdEndsPar par) {
        StringBuilder sb = new StringBuilder(300);
        sb.append(COPYRIGHT);
        sb.append(Const.nl);
        sb.append(HELP_MESSAGE);
        sb.append(Const.nl);
        sb.append(Const.nl);
        sb.append("Program            :  ");
        sb.append(PROGRAM);
        sb.append(Const.nl);
        sb.append("Start Time         :  ");
        sb.append(Utilities.timeStamp());
        sb.append(Const.nl);
        sb.append("Max Memory         :  ");
        long maxMemory = Runtime.getRuntime().maxMemory();
        if (maxMemory != Long.MAX_VALUE) {
            long maxMb = maxMemory / (1024*1024);
            sb.append(maxMb);
            sb.append(" MB");
        }
        else {
            sb.append("[no limit])");
        }
        sb.append(Const.nl);
        sb.append(Const.nl);
        sb.append(parameters(par));
        return sb.toString();
    }

    private static String parameters(IbdEndsPar par) {
        StringBuilder sb = new StringBuilder(150);
        sb.append("Parameters");
        sb.append(Const.nl);
        sb.append("  gt               :  ");
        sb.append(par.gt());
        sb.append(Const.nl);
        if (par.chromInt()!=null) {
            sb.append("  chrom            :  ");
            sb.append(par.chromInt().toString());
            sb.append(Const.nl);
        }
        sb.append("  ibd              :  ");
        sb.append(par.ibd());
        sb.append(Const.nl);
        sb.append("  map              :  ");
        sb.append(par.map());
        sb.append(Const.nl);
        sb.append("  out              :  ");
        sb.append(par.out());
        if (par.excludesamples()!=null) {
            sb.append(Const.nl);
            sb.append("  excludesamples   :  ");
            sb.append(par.excludesamples());
        }
        sb.append(Const.nl);
        sb.append("  quantiles        :  ");
        FloatArray q = par.quantiles();
        if (q!=null) {
            for (int j=0, n=q.size(); j<n; ++j) {
                if (j>0) {
                    sb.append(Const.comma);
                }
                sb.append(q.get(j));
            }
        }
        else {
            sb.append("none specified");
        }
        sb.append(Const.nl);
        sb.append("  nsamples         :  ");
        sb.append(par.nsamples());
        sb.append(Const.nl);
        sb.append("  err              :  ");
        sb.append(par.err());
        sb.append(Const.nl);
        sb.append("  estimate-err     :  ");
        sb.append(par.estimate_err());
        sb.append(Const.nl);
        sb.append("  gc-err           :  ");
        sb.append(par.gc_err());
        sb.append(Const.nl);
        sb.append("  gc-bp            :  ");
        sb.append(par.gc_bp());
        sb.append(Const.nl);
        sb.append("  min-maf          :  ");
        sb.append(par.min_maf());
        sb.append(Const.nl);
        sb.append("  seed             :  ");
        sb.append(par.seed());
        sb.append(Const.nl);
        sb.append("  nthreads         :  ");
        sb.append(par.nthreads());
        return sb.toString();
    }

    private static String statistics(IbdEndsPar par, IbdEndsStats stats) {
        StringBuilder sb = new StringBuilder(300);
        sb.append(Const.nl);
        sb.append("Analysis summary");
        sb.append(Const.nl);
        sb.append("  samples          :  ");
        sb.append(stats.nSamples());
        sb.append(Const.nl);
        sb.append("  markers          :  ");
        sb.append(stats.nMarkers());
        sb.append(Const.nl);
        sb.append("  segments         :  ");
        sb.append(stats.ibdSegmentCnt());
        sb.append(Const.nl);
        if (par.estimate_err()) {
            double estErr = stats.discordRate();
            sb.append("  estimated err    :  ");
            sb.append(DF.format(estErr));
            boolean warning = (estErr>=par.err() && estErr/par.err()>=MAX_ERR_RATIO)
                    || (par.err()>=estErr && estErr>0f && par.err()/estErr>=MAX_ERR_RATIO);
            if (warning) {
                sb.append("      Recommendation: reanalyze with err=");
                sb.append(DF.format(estErr));
            }
            sb.append(Const.nl);
         }
        return sb.toString();
    }

    private static String endInfo(long startNanoTime) {
        StringBuilder sb = new StringBuilder(300);
        long elapsedNanoTime = System.nanoTime() - startNanoTime;
        sb.append(Const.nl);
        sb.append("Wallclock Time:    :  ");
        sb.append(Utilities.elapsedNanos(elapsedNanoTime));
        sb.append(Const.nl);
        sb.append("End Time           :  ");
        sb.append(Utilities.timeStamp());
        return sb.toString();
    }

    private static double round(double d, int significantDigits) {
        assert significantDigits>0;
        BigDecimal bd = new BigDecimal(d);
        bd.round(new MathContext(significantDigits));
        return bd.doubleValue();
    }
}
