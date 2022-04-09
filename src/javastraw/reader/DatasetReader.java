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

package javastraw.reader;

import javastraw.reader.block.Block;
import javastraw.reader.datastructures.ListOfDoubleArrays;
import javastraw.reader.norm.NormalizationVector;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationType;

import java.io.IOException;
import java.util.List;

public interface DatasetReader {

    boolean isActive();

    void setActive(boolean status);

    int getVersion();

    Dataset read() throws IOException;

    Matrix readMatrix(String key, int resolution) throws IOException;

    Block readNormalizedBlock(int blockNumber, String zdKey, NormalizationType no,
                              int chr1Index, int chr2Index, HiCZoom zoom) throws IOException;

    List<Integer> getBlockNumbers(String zdKey);

    NormalizationVector readNormalizationVector(NormalizationType type, int chrIdx, HiCZoom.HiCUnit unit, int binSize) throws IOException;

    NormalizationVector readNormalizationVectorPart(NormalizationType type, int chrIdx, HiCZoom.HiCUnit unit, int binSize, int bound1, int bound2) throws IOException;

    ListOfDoubleArrays readExpectedVectorPart(long position, long nVals) throws IOException;

    String getPath();

    String readStats() throws IOException;

    NormalizationVector getNormalizationVector(int chr1Idx, HiCZoom zoom, NormalizationType normalizationType);

    int getDepthBase();
}
