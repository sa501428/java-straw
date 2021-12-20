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

package javastraw.reader.pearsons;

import javastraw.matrices.BasicMatrix;
import javastraw.matrices.InMemoryMatrix;
import javastraw.tools.ParallelizedJuicerTools;

import java.util.BitSet;
import java.util.concurrent.atomic.AtomicInteger;

public class Pearsons {

    private static double getVectorMean(double[] vector) {
        double sum = 0;
        int count = 0;
        for (double aVector : vector) {
            if (!Double.isNaN(aVector)) {
                sum += aVector;
                count++;
            }
        }
        return count == 0 ? 0 : sum / count;
    }

    public static BasicMatrix computeParallelizedPearsons(double[][] matrix, int dim, BitSet bitSet) {
        // Subtract row means
        inPlaceParallelSubtractRowMeans(matrix, dim, bitSet);

        BasicMatrix pearsons = new InMemoryMatrix(dim);

        // Set diagonal to 1, set centromere to NaN
        AtomicInteger dCounter = new AtomicInteger(0);
        ParallelizedJuicerTools.launchParallelizedCode(() -> {
            int i = dCounter.getAndIncrement();
            while (i < dim) {
                if (bitSet.get(i)) {
                    for (int j = i + 1; j < dim; j++) {
                        float corr;
                        if (bitSet.get(j)) {
                            corr = (float) PearsonCorrelationMetric.corr(matrix[i], matrix[j]);
                        } else {
                            corr = Float.NaN;
                        }
                        synchronized (pearsons) {
                            pearsons.setEntry(i, j, corr);
                            pearsons.setEntry(j, i, corr);
                        }
                    }
                } else {
                    synchronized (pearsons) {
                        for (int j = i + 1; j < dim; j++) {
                            pearsons.setEntry(i, j, Float.NaN);
                            pearsons.setEntry(j, i, Float.NaN);
                        }
                    }
                }

                i = dCounter.getAndIncrement();
            }
        });

        for (int i = 0; i < dim; i++) {
            if (bitSet.get(i)) {
                pearsons.setEntry(i, i, 1.0f);
            } else {
                pearsons.setEntry(i, i, Float.NaN);
            }
        }
        return pearsons;
    }

    private static void inPlaceParallelSubtractRowMeans(double[][] matrix, int dim, BitSet bitSet) {
        double[] rowMeans = new double[dim];
        AtomicInteger index = new AtomicInteger(0);
        ParallelizedJuicerTools.launchParallelizedCode(() -> {
            int i = index.getAndIncrement();
            while (i < dim) {
                if (bitSet.get(i)) {
                    rowMeans[i] = getVectorMean(matrix[i]);
                }
                i = index.getAndIncrement();
            }
        });

        AtomicInteger index2 = new AtomicInteger(0);
        ParallelizedJuicerTools.launchParallelizedCode(() -> {
            int i = index2.getAndIncrement();
            while (i < dim) {
                if (bitSet.get(i)) {
                    for (int j = 0; j < dim; j++) {
                        matrix[i][j] -= rowMeans[j];
                    }
                }
                i = index2.getAndIncrement();
            }
        });
    }
}
