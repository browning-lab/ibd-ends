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
import blbutil.BlockLineReader;
import blbutil.BGZIPOutputStream;
import blbutil.Const;
import blbutil.DoubleArray;
import blbutil.FloatArray;
import blbutil.Utilities;
import ints.WrappedIntArray;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.Random;
import vcf.RefGT;

/**
 * <p>Class {@code IbdEndsRunnable} performs IBD segment end point estimation.</p>
 *
 * <p>Instances of class {@code IbdEndsRunnable} are not thread-safe.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class IbdEndsRunnable implements Runnable {

    private static final int BAOS_THRESHOLD = 1<<18;
    private static final float MIN_ERR_ESTIMATION_MORGANS = 0.02f;

    private final IbdEndsPar par;
    private final IbdEndsData data;
    private final IbdEndsStats stats;
    private final int nQuantilesP1;
    private final BlockLineReader reader;
    private final HapIbdParser parser;
    private final IbdEnds ibdEnds;
    private final boolean estErr;

    private final Random rand;
    private final double[] fwdQuantiles;
    private final double[] bwdQuantiles;
    private final int[] fwdEnds;
    private final int[] bwdEnds;

    private final ByteArrayOutputStream baos;
    private final SynchFileOutputStream outStream;

    /**
     * Constructs a new instance of {@code IbdEndsRunnable} from the specified
     * data.
     * @param data the analysis input data
     * @param reader an iterator that string arrays containing IBD segment
     * in hap-ibd format
     * @param stats an object which stores analysis statistics
     * @param outStream the output stream to which bgzip-compressed output
     * will be written
     * @throws NullPointerException if any parameter is {@code null}
     */
    public IbdEndsRunnable(IbdEndsData data, BlockLineReader reader,
            IbdEndsStats stats, SynchFileOutputStream outStream) {
        if (reader==null) {
            throw new NullPointerException(BlockLineReader.class.toString());
        }
        if (stats==null) {
            throw new NullPointerException(IbdEndsStats.class.toString());
        }
        if (outStream==null) {
            throw new NullPointerException(OutputStream.class.toString());
        }
        this.par = data.par();
        this.data = data;
        this.stats = stats;
        this.nQuantilesP1 = par.quantiles().size() + 1;
        this.reader = reader;
        this.parser = data.hapIbdParser();
        this.ibdEnds = new IbdEnds(data);
        this.estErr = par.estimate_err();

        this.rand = (par.nsamples()>0 ? new Random(par.seed()) : null);
        this.fwdQuantiles = quantilesArray(par);
        this.bwdQuantiles = fwdQuantiles.clone();
        this.fwdEnds = new int[fwdQuantiles.length];
        this.bwdEnds = new int[fwdQuantiles.length];
        this.baos = new ByteArrayOutputStream(3*BAOS_THRESHOLD/2 + 1);
        this.outStream = outStream;
    }

    /**
     * Estimates IBD segment endpoints and writes a bgzip-compressed string
     * representation of the IBD segments and estimated endpoints to an
     * output file.
     */
    @Override
    public void run() {
        try {
            try (PrintWriter out = new PrintWriter(new BGZIPOutputStream(baos, false))) {
                assert reader.hasNext();
                String[] lines = reader.next();
                while (lines!=BlockLineReader.SENTINAL) {
                    Arrays.stream(lines)
                            .forEach(line -> processSegment(line, out));
                    flushBuffer(BAOS_THRESHOLD, out);
                    assert reader.hasNext();
                    lines = reader.next();
                }
                assert lines==BlockLineReader.SENTINAL;
                flushBuffer(0, out);
            }
        }
        catch (Throwable t) {
            Utilities.exit(t);
        }
    }

    private void processSegment(String line, PrintWriter out) {
        SharedSegment ss = parser.parse(line);
        if (ss!=HapIbdParser.NIL) {
            stats.incrementIbdSegmentCnt();
            if (rand!=null) {
                rand.setSeed(par.seed() + ss.hashCode());
                for (int j=nQuantilesP1; j<fwdQuantiles.length; ++j) {
                    fwdQuantiles[j] = rand.nextDouble();
                    bwdQuantiles[j] = rand.nextDouble();
                }
            }
            int focusPos = ibdEnds.getEnds(ss, fwdQuantiles, fwdEnds,
                    bwdQuantiles, bwdEnds);
            if (estErr) {
                updateDiscordRate(ss.hap1(), ss.hap2(), bwdEnds[0], fwdEnds[0]);
            }
            out.print(line);
            out.print(Const.tab);
            out.print(focusPos);
            for (int j=1; j<bwdEnds.length; ++j) {
                // 0-th element is from "length-quantile" parameter
                double startMorgans = data.baseToMorgans(bwdEnds[j]);
                double endMorgans = data.baseToMorgans(fwdEnds[j]);
                double cm = 100*(endMorgans - startMorgans);
                out.print(Const.tab);
                out.print(bwdEnds[j]);
                out.print(Const.tab);
                out.print(fwdEnds[j]);
                out.print(Const.tab);
                print3(cm, out);
            }
            out.println();
        }
    }

    public void updateDiscordRate(int hap1, int hap2, int startPos, int endPos) {
        assert startPos <= endPos;
        WrappedIntArray basePos = data.fwdPos();
        DoubleArray morganPos = data.fwdMorgans();
        RefGT fwdGT = data.fwdGT();
        int startMarker = basePos.binarySearch(startPos);
        if (startMarker<0) {
            startMarker = -startMarker-1;
        }
        int endMarker = basePos.binarySearch(startMarker, basePos.size(), endPos);
        if (endMarker<0) {
            endMarker = -endMarker-2; // -2 so that marker is within interval
        }
        if (startMarker<=endMarker ) {
            double length = morganPos.get(endMarker) - morganPos.get(startMarker);
            if (length >= MIN_ERR_ESTIMATION_MORGANS) {
                int discordCnt = 0;
                for (int m=startMarker; m<=endMarker; ++m) {
                    if (fwdGT.allele(m, hap1)!=fwdGT.allele(m, hap2)) {
                        ++discordCnt;
                    }
                }
                int totalChecked = endMarker - startMarker + 1;
                stats.updateDiscordRate(discordCnt, totalChecked);
            }
        }
    }

    /* Prints the specified double with 3 decimal places to the specified
       PrintWriter
    */
    private static void print3(double d, PrintWriter out) {
        if (d<0) {
            out.print('-');
            d = -d;
        }
        d += 5e-4;
        long integerPart = (long) d;
        double fraction = Math.floor(1000*(d - integerPart));
        out.print(integerPart);
        out.print('.');
        if (fraction<100) {
            out.print('0');
            if (fraction<10) {
                out.print('0');
            }
        }
        out.print((int) fraction);
    }

    /**
     * Returns the header line for the ibd-ends analysis output.
     * @param par the analysis parameters
     * @return the header line for the ibd-ends analysis output
     * @throws NullPointerException if {@code par == null}
     */
    public static String ibdOutHeader(IbdEndsPar par) {
        // Need DecimalFormat.format() must remain within float's accuracy
        // so that values like 0.1 printed are not printed as 0.10000001
        DecimalFormat df = new DecimalFormat("0.00#####");
        StringBuilder sb = new StringBuilder(80);
        FloatArray quantiles = par.quantiles();
        int nSamples = par.nsamples();
        sb.append("ID1");
        sb.append(Const.tab);
        sb.append("HAP1");
        sb.append(Const.tab);
        sb.append("ID2");
        sb.append(Const.tab);
        sb.append("HAP2");
        sb.append(Const.tab);
        sb.append("CHROM");
        sb.append(Const.tab);
        sb.append("START");
        sb.append(Const.tab);
        sb.append("END");
        sb.append(Const.tab);
        sb.append("CM");
        sb.append(Const.tab);
        sb.append("FOCUS");
        for (int j=0, n=quantiles.size(); j<n; ++j) {
            String quantile = df.format(quantiles.get(j));
            int dotIndex = quantile.indexOf('.');
            assert dotIndex>=0;
            quantile = quantile.substring(dotIndex);
            sb.append(Const.tab);
            sb.append("STA");
            sb.append(quantile);
            sb.append(Const.tab);
            sb.append("END");
            sb.append(quantile);
            sb.append(Const.tab);
            sb.append("CM");
            sb.append(quantile);
        }
        for (int j=1; j<=nSamples; ++j) {
            String index = String.valueOf(j);
            sb.append(Const.tab);
            sb.append("STA-");
            sb.append(index);
            sb.append(Const.tab);
            sb.append("END-");
            sb.append(index);
            sb.append(Const.tab);
            sb.append("CM-");
            sb.append(index);
        }
        return sb.toString();
    }

    private static double[] quantilesArray(IbdEndsPar par) {
        FloatArray quants = par.quantiles();
        double[] extQuants = new double[quants.size() + par.nsamples() + 1];
        extQuants[0] = par.length_quantile();
        for (int j=0, n=quants.size(); j<n; ++j) {
            extQuants[j+1] = quants.get(j);
        }
        return extQuants;
    }

    private void flushBuffer(int byteThreshold, PrintWriter out) {
        if (baos.size() >= byteThreshold) {
            try {
                out.flush();
                outStream.write(baos.toByteArray());
                baos.reset();
            } catch (IOException ex) {
                Utilities.exit("ERROR: ", ex);
            }
        }
    }
}
