/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2011-2021 Broad Institute, Aiden Lab, Rice University, Baylor College of Medicine
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 *  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 *  FITNESS FOR A PARTICULAR PURPOSE AND NON-INFRINGEMENT. IN NO EVENT SHALL THE
 *  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 *  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 *  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 *  THE SOFTWARE.
 */


package javastraw.reader.expected;


import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.datastructures.ListOfDoubleArrays;
import javastraw.reader.datastructures.ListOfFloatArrays;
import javastraw.reader.iterators.IteratorContainer;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationType;

import java.util.*;

/**
 * Computes an "expected" density vector.  Essentially there are 3 steps to using this class
 * <p/>
 * (1) instantiate it with a collection of Chromosomes (representing a genome) and a grid size
 * (2) loop through the pair data,  calling addDistance for each pair, to accumulate all counts
 * (3) when data loop is complete, call computeDensity to do the calculation
 * <p/>
 * <p/>
 * Methods are provided to save the result of the calculation to a binary file, and restore it.  See the
 * DensityUtil class for example usage.
 */
public class ExpectedValueCalculation {

    // used for sparsity cutoff
    private static final int MIN_VALS_NEEDED = 400;
    private final int binSize;

    private final int numberOfBins;
    /**
     * chromosome index -> "normalization factor", essentially a fudge factor to make
     * the "expected total"  == observed total
     */
    private final double[] chrScaleFactors;
    private final NormalizationType type;
    /**
     * chromosome index -> total count for that chromosome
     */
    private final double[] chromosomeCounts;
    /**
     * Genome wide count of binned reads at a given distance
     */
    private final List<List<Double>> actualDistances;
    /**
     * Chromosome in this genome, needed for normalizations
     */
    private final Chromosome[] chromosomesMap;
    /**
     * Expected count at a given binned distance from diagonal
     */
    private ListOfDoubleArrays densityAvg;

    private int indexWhereSparsityStarts;

    /**
     * Instantiate a DensityCalculation.  This constructor is used to compute the "expected" density from pair data.
     *
     * @param chromosomeHandler Handler for list of chromosomesMap, mainly used for size
     * @param binSize           Grid size, used for binning appropriately
     * @param type              Identifies the observed matrix type,  either NONE (observed), VC, or KR.
     */
    public ExpectedValueCalculation(ChromosomeHandler chromosomeHandler, int binSize, NormalizationType type) {

        this.type = type;
        this.binSize = binSize;
        this.indexWhereSparsityStarts = 100000000 / binSize;
        long maxLen = 0;

        int maxNumChromosomes = chromosomeHandler.getMaxChromIndex() + 1;
        chromosomesMap = new Chromosome[maxNumChromosomes];
        Arrays.fill(chromosomesMap, null);
        chromosomeCounts = new double[maxNumChromosomes];
        chrScaleFactors = new double[maxNumChromosomes];

        for (Chromosome chromosome : chromosomeHandler.getChromosomeArrayWithoutAllByAll()) {
            if (chromosome != null) {
                chromosomesMap[chromosome.getIndex()] = chromosome;
                try {
                    maxLen = Math.max(maxLen, chromosome.getLength());
                } catch (NullPointerException error) {
                    System.err.println("Problem with creating fragment-delimited maps, NullPointerException.\n" +
                            "This could be due to a null fragment map or to a mismatch in the chromosome name in " +
                            "the fragment map vis-a-vis the input file or chrom.sizes file.\n" +
                            "Exiting.");
                    System.exit(63);
                } catch (ArrayIndexOutOfBoundsException error) {
                    System.err.println("Problem with creating fragment-delimited maps, ArrayIndexOutOfBoundsException.\n" +
                            "This could be due to a null fragment map or to a mismatch in the chromosome name in " +
                            "the fragment map vis-a-vis the input file or chrom.sizes file.\n" +
                            "Exiting.");
                    System.exit(22);
                }
            }
        }

        numberOfBins = (int) (maxLen / binSize) + 1;
        actualDistances = new ArrayList<List<Double>>(numberOfBins);
        for (int z = 0; z < numberOfBins; z++) {
            actualDistances.add(new ArrayList<>());
        }
    }

    public static boolean isValidNormValue(float v) {
        return !Float.isNaN(v) && v > 0;
    }

    public int getBinSize() {
        return binSize;
    }

    public void addDistance(int chrIdx, ContactRecord cr) {
        addDistance(chrIdx, cr.getBinX(), cr.getBinY(), cr.getCounts());
    }

    /**
     * Add an observed distance.  This is called for each pair in the data set
     *
     * @param chrIdx index of chromosome where observed, so can increment count
     * @param bin1   Position1 observed in units of "bins"
     * @param bin2   Position2 observed in units of "bins"
     */
    public synchronized void addDistance(int chrIdx, int bin1, int bin2, double weight) {

        // Ignore NaN values
        if (Double.isNaN(weight)) return;
        if (chromosomesMap[chrIdx] == null) return;

        chromosomeCounts[chrIdx] += weight;
        int dist = Math.abs(bin1 - bin2);
        actualDistances.get(dist).add(weight);
    }

    public void merge(ExpectedValueCalculation otherEVCalc) {
        for (int z = 0; z < chromosomeCounts.length; z++) {
            chromosomeCounts[z] += otherEVCalc.chromosomeCounts[z];
        }

        for (int i = 0; i < actualDistances.size(); i++) {
            if (otherEVCalc.actualDistances.get(i).size() > 0) {
                actualDistances.get(i).addAll(otherEVCalc.actualDistances.get(i));
            }
        }
    }

    public boolean hasData() {
        for (double val : chromosomeCounts) {
            if (val > 0) return true;
        }
        return false;
    }

    /**
     * Compute the "density" -- port of python function getDensityControls().
     * The density is a measure of the average distribution of counts genome-wide for a ligated molecule.
     * The density will decrease as distance from the center diagonal increases.
     * First compute "possible distances" for each bin.
     * "possible distances" provides a way to normalize the counts. Basically it's the number of
     * slots available in the diagonal.  The sum along the diagonal will then be the count at that distance,
     * an "expected" or average uniform density.
     */
    public synchronized void computeDensity() {

        long maxNumBins = 0;

        //System.err.println("# of bins=" + numberOfBins);
        /**
         * Genome wide binned possible distances
         */
        double[] possibleDistances = new double[numberOfBins];

        long[] nChrBins = new long[chromosomesMap.length];
        for (int z = 0; z < chromosomesMap.length; z++) {
            Chromosome chromosome = chromosomesMap[z];
            if (chromosome == null || chromosomeCounts[z] < 1) continue;
            long len = chromosome.getLength();
            nChrBins[z] = len / binSize;
        }


        for (int z = 0; z < chromosomesMap.length; z++) {
            Chromosome chromosome = chromosomesMap[z];
            if (chromosome == null || chromosomeCounts[z] < 1) continue;

            maxNumBins = Math.max(maxNumBins, nChrBins[z]);

            for (int i = 0; i < nChrBins[z]; i++) {
                possibleDistances[i] += (nChrBins[z] - i);
            }
        }
        //System.err.println("max # bins " + maxNumBins);
        densityAvg = new ListOfDoubleArrays(maxNumBins);

        // Smoothing.  Keep pointers to window size.  When read counts drops below 400 (= 5% shot noise), smooth

        ListOfDoubleArrays numEntriesForDist = new ListOfDoubleArrays(maxNumBins);
        for (int q = 0; q < maxNumBins; q++) {
            List<Double> values = new ArrayList<>(actualDistances.get(q));
            if (values.size() < MIN_VALS_NEEDED && q < indexWhereSparsityStarts) {
                indexWhereSparsityStarts = q;
            }

            int window = 0;
            while (values.size() < MIN_VALS_NEEDED) {
                window++;
                if (q - window > -1) {
                    values.addAll(actualDistances.get(q - window));
                }
                if (q + window < actualDistances.size()) {
                    values.addAll(actualDistances.get(q + window));
                }
            }

            numEntriesForDist.set(q, values.size());
            densityAvg.set(q, QuickMedian.fastMedian(values));
        }
        actualDistances.clear();
        // previously mean window 400
        // now median window, 400

        // rolling mean from 1st sparse index
        densityAvg.doRollingMeanFromIndex(indexWhereSparsityStarts, 100);

        // Compute fudge factors for each chromosome so the total "expected" count for that chromosome == the observed

        for (int z = 0; z < chromosomesMap.length; z++) {
            Chromosome chromosome = chromosomesMap[z];
            if (chromosome == null || chromosomeCounts[z] < 1) continue;

            double expectedCount = 0;
            for (long n = 0; n < nChrBins[z]; n++) {
                if (n < maxNumBins) {
                    final double v = densityAvg.get(n);
                    // this is the sum of the diagonal for this particular chromosome.
                    // the value in each bin is multiplied by the length of the diagonal to get expected count
                    // the total at the end should be the sum of the expected matrix for this chromosome
                    // i.e., for each chromosome, we calculate sum (genome-wide actual)/(genome-wide possible) == v
                    // then multiply it by the chromosome-wide possible == nChrBins - n.
                    expectedCount += (nChrBins[z] - n) * v;
                }
            }

            double observedCount = chromosomeCounts[chromosome.getIndex()];
            chrScaleFactors[z] = expectedCount / observedCount;
        }
    }

    /**
     * Accessor for the normalization factors
     *
     * @return The normalization factors
     */
    public Map<Integer, Double> getChrScaleFactors() {
        Map<Integer, Double> scaleFactorsMap = new HashMap<>();
        for (int z = 0; z < chrScaleFactors.length; z++) {
            if (chrScaleFactors[z] > 0) {
                scaleFactorsMap.put(z, chrScaleFactors[z]);
            }
        }
        return scaleFactorsMap;
    }

    /**
     * Accessor for the densities
     *
     * @return The densities
     */
    public ListOfDoubleArrays getDensityAvg() {
        return densityAvg;
    }

    /**
     * Accessor for the normalization type
     *
     * @return The normalization type
     */
    public NormalizationType getType() {
        return type;
    }

    public ExpectedValueFunctionImpl getExpectedValueFunction() {
        computeDensity();
        return new ExpectedValueFunctionImpl(type, HiCZoom.HiCUnit.BP,
                binSize, densityAvg, getChrScaleFactors());
    }

    // TODO: this is often inefficient, we have all of the contact records when we leave norm calculations, should do this there if possible
    public void addDistancesFromIterator(int chrIndx, IteratorContainer ic, ListOfFloatArrays vector) {
        Iterator<ContactRecord> iterator = ic.getNewContactRecordIterator();
        while (iterator.hasNext()) {
            ContactRecord cr = iterator.next();
            int x = cr.getBinX();
            int y = cr.getBinY();
            final float counts = cr.getCounts();
            float xVal = vector.get(x);
            float yVal = vector.get(y);
            if (isValidNormValue(xVal) & isValidNormValue(yVal)) {
                double value = counts / (xVal * yVal);
                addDistance(chrIndx, x, y, value);
            }
        }
    }
}
