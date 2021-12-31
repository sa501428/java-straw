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

package javastraw.reader.expected;

import javastraw.reader.datastructures.ListOfDoubleArrays;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationType;

public interface ExpectedValueFunction {

    int DEFAULT_CORRECTION_DISTANCE = 5000000;

    double getExpectedValue(int chrIdx, long distance);

    long getLength();

    NormalizationType getNormalizationType();

    HiCZoom.HiCUnit getUnit();

    int getBinSize();

    ListOfDoubleArrays getExpectedValuesNoNormalization();

    ListOfDoubleArrays getExpectedValuesWithNormalization(int chrIdx);

    static String getKey(HiCZoom zoom, NormalizationType normType, boolean isCorrected, int window) {
        if (isCorrected) {
            return zoom.getKey() + "_" + normType + "_" + window;
        } else {
            return zoom.getKey() + "_" + normType;
        }
    }

    static String getKey(HiCZoom.HiCUnit unit, int binSize, NormalizationType normType) {
        return getKey(new HiCZoom(unit, binSize), normType, false, 0);
    }

    ExpectedValueFunction getCorrectedVersion(int window);
}
