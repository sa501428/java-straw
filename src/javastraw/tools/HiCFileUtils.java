/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2011-2020 Broad Institute, Aiden Lab, Rice University, Baylor College of Medicine
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
 *  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 *  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 *  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 *  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 *  THE SOFTWARE.
 */

package javastraw.tools;

import javastraw.reader.Dataset;
import javastraw.reader.DatasetReaderV2;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.expected.ExpectedValueFunction;
import javastraw.reader.expected.ExpectedValueFunctionImpl;
import javastraw.reader.norm.NormalizationVector;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationHandler;
import javastraw.reader.type.NormalizationType;

import java.io.IOException;
import java.util.Map;

/**
 * Utility functions to dump various bits of a hic file to stdout or file(s)
 */

class HiCFileUtils {

    private final Dataset dataset;

    private HiCFileUtils(String hicfile, boolean useCache, boolean useDynamicBlockIndex) throws IOException {
        DatasetReaderV2 reader = new DatasetReaderV2(hicfile, useCache, useDynamicBlockIndex);
        dataset = reader.read();
    }

    public static void main(String[] args) throws IOException {
        HiCFileUtils utils = new HiCFileUtils(args[0], false, true);
        utils.dumpNormalizationVectors(NormalizationHandler.KR, "1", HiCZoom.HiCUnit.BP, 250000);
        utils.dumpExpectedVectors(NormalizationHandler.KR, HiCZoom.HiCUnit.BP, 1000000);
    }

    private void dumpNormalizationVectors(NormalizationType normType, String chrName, HiCZoom.HiCUnit unit, int binSize) {
        Chromosome chromosome = findChromosome(chrName);
        HiCZoom zoom = new HiCZoom(unit, binSize);
        NormalizationVector nv = dataset.getNormalizationVector(chromosome.getIndex(), zoom, normType);
        String label = "Normalization vector: type = " + normType.getLabel() + " chr = " + chrName +
                " resolution = " + binSize + " " + unit;
        System.out.println(label);
        /*
        for(int i=0; i<data.length; i++) {
            System.out.println(data[i]);
        }
        */
        for (double[] array : nv.getData().getValues()) {
            for (double datum : array) {
                System.out.println(datum);
            }
        }
    }

    private void dumpExpectedVectors(NormalizationType normType, HiCZoom.HiCUnit unit, int binSize) {

        Map<String, ExpectedValueFunction> expValFunMap = dataset.getExpectedValueFunctionMap();
        for (Map.Entry<String, ExpectedValueFunction> entry : expValFunMap.entrySet()) {

            ExpectedValueFunctionImpl ev = (ExpectedValueFunctionImpl) entry.getValue();

            if (ev.getUnit().equals(unit) && ev.getBinSize() == binSize && ev.getNormalizationType().equals(normType)) {
                String label = ev.getNormalizationType() + "\t" + ev.getUnit() + "\t" + ev.getBinSize();

                System.out.println("Norm factors: " + label);
                for (Map.Entry<Integer, Double> nf : ev.getNormFactors().entrySet()) {
                    System.out.println(nf.getKey() + "\t" + nf.getValue());
                }
    
                System.out.println("Expected values: " + label);
                /*
                for (int i = 0; i < values.length; i++) {
                    System.out.println(values[i]);
                }
                */
                for (double[] values : ev.getExpectedValuesNoNormalization().getValues()) {
                    for (double datum : values) {
                        System.out.println(datum);
                    }
                }
    
                System.out.println("End expected values: " + label);
                System.out.println();
            }
        }
    }

    private Chromosome findChromosome(String name) {
        for (Chromosome chr : dataset.getChromosomeHandler().getChromosomeArray()) {
            if (chr.getName().equalsIgnoreCase(name)) return chr;
        }
        return null;
    }
}
