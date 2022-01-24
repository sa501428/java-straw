/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2011-2020 Rice University, Baylor College of Medicine, Aiden Lab
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

import javastraw.reader.block.Block;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.expected.ExpectedValueFunction;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.NormalizationType;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class ExtractingOEDataUtils {

    private static final double e = Math.exp(1);

    public static float[][] simpleLog(float[][] matrix, float pseudocount) {
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                float val = matrix[i][j];
                if (!Float.isNaN(val)) { // && val > 0
                    matrix[i][j] = (float) Math.log(val + pseudocount);
                }
            }
        }
        return matrix;
    }

    public static float[][] simpleLogWithCleanup(float[][] matrix, float pseudocount) {
        matrix = simpleLog(matrix, pseudocount);
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[0].length; j++) {
                if (Float.isInfinite(matrix[i][j]) || Math.abs(matrix[i][j]) < 1E-10) {
                    matrix[i][j] = Float.NaN;
                }
            }
        }
        return matrix;
    }

    public static float[][] logOEP1(float[][] matrix, double averageCount) {
        double denom = Math.log(averageCount + e);
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                matrix[i][j] = (float) (Math.log(matrix[i][j] + e) / denom);
            }
        }
        return matrix;
    }

    public static float[][] extractObsOverExpBoundedRegion(MatrixZoomData zd, int binXStart, int binXEnd,
                                                           int binYStart, int binYEnd, int numRows, int numCols,
                                                           NormalizationType normalizationType,
                                                           ExpectedValueFunction df, int chrIndex, double threshold,
                                                           boolean isIntra, boolean fillUnderDiagonal, ThresholdType thresholdType,
                                                           float pseudocount, float invalidReplacement) throws IOException {
        if (isIntra && df == null) {
            System.err.println("DF is null");
            return null;
        }

        List<ContactRecord> records = getOEBlocks(zd, binXStart, binXEnd, binYStart, binYEnd,
                normalizationType, isIntra, fillUnderDiagonal, df, chrIndex, pseudocount, thresholdType,
                threshold, invalidReplacement);

        // numRows/numCols is just to ensure a set size in case bounds are approximate
        // left upper corner is reference for 0,0
        float[][] data = new float[numRows][numCols];
        for (ContactRecord rec : records) {
            HiCFileTools.fillInMatrixWithRecords(binXStart, binYStart, numRows, numCols,
                    fillUnderDiagonal, data, rec);
        }
        // force cleanup
        records = null;
        return data;
    }

    public static List<ContactRecord> getOEBlocks(MatrixZoomData zd, int binXStart, int binXEnd,
                                                  int binYStart, int binYEnd, NormalizationType normalizationType,
                                                  boolean isIntra, boolean fillUnderDiagonal, ExpectedValueFunction df,
                                                  int chrIndex, float pseudocount, ThresholdType thresholdType,
                                                  double threshold, float invalidReplacement) throws IOException {
        List<Block> blocks = HiCFileTools.getAllRegionBlocks(zd, binXStart, binXEnd, binYStart, binYEnd,
                normalizationType, fillUnderDiagonal);
        List<ContactRecord> records = new ArrayList<>();

        for (Block b : blocks) {
            if (b != null) {
                for (ContactRecord rec : b.getContactRecords()) {
                    double expected = getExpected(rec, df, chrIndex, isIntra, zd.getAverageCount());
                    double val = rec.getCounts();

                    double observed = val + pseudocount;
                    expected = expected + pseudocount;
                    double answer = Double.NaN;

                    if (thresholdType.equals(ThresholdType.LOG_BASE_EXP_OF_OBS)) {

                        answer = (Math.log(observed) / Math.log(expected));

                    } else if (thresholdType.equals(ThresholdType.TRUE_OE_LOG)
                            || thresholdType.equals(ThresholdType.LOG_OE_BOUNDED)
                            || thresholdType.equals(ThresholdType.LOG_OE_BOUNDED_SCALED_BTWN_ZERO_ONE)) {

                        answer = Math.log(observed / expected);
                        if (thresholdType.equals(ThresholdType.LOG_OE_BOUNDED)) {
                            answer = Math.min(Math.max(-threshold, answer), threshold);
                        } else if (thresholdType.equals(ThresholdType.LOG_OE_BOUNDED_SCALED_BTWN_ZERO_ONE)) {
                            answer = Math.min(Math.max(-threshold, answer), threshold);
                            answer = (answer + threshold) / (2 * threshold);
                        }
                    } else if (thresholdType.equals(ThresholdType.TRUE_OE)
                            || thresholdType.equals(ThresholdType.TRUE_OE_LINEAR)) {
                        answer = observed / expected;
                        if (thresholdType.equals(ThresholdType.TRUE_OE_LINEAR)) {
                            if (answer < 1) {
                                answer = 1 - 1 / answer;
                            } else {
                                answer -= 1;
                            }
                        }
                    }

                    float floatAnswer = (float) answer;
                    if (Float.isNaN(floatAnswer) || Float.isInfinite(floatAnswer)) {
                        floatAnswer = invalidReplacement;
                    }

                    records.add(new ContactRecord(rec.getBinX(), rec.getBinY(), floatAnswer));

                }
            }
        }
        blocks = null;
        return records;
    }

    private static double getExpected(ContactRecord rec, ExpectedValueFunction df, int chrIndex,
                                      boolean isIntra, double averageCount) {
        int x = rec.getBinX();
        int y = rec.getBinY();
        double expected;
        if (isIntra) {
            int dist = Math.abs(x - y);
            expected = df.getExpectedValue(chrIndex, dist);
        } else {
            expected = (averageCount > 0 ? averageCount : 1);
        }
        return expected;
    }

    public enum ThresholdType {
        LOG_BASE_EXP_OF_OBS, // log(obs+e)/log(exp+e)
        LOG_OE_BOUNDED, // log((obs+1)/(exp+1)) with thres
        TRUE_OE, // (obs+1)/(exp+1)
        TRUE_OE_LOG, // log((obs+1)/(exp+1))
        TRUE_OE_LINEAR,
        LOG_OE_BOUNDED_SCALED_BTWN_ZERO_ONE
    }
}
