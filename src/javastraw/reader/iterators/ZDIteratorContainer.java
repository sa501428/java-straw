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

package javastraw.reader.iterators;

import javastraw.reader.DatasetReader;
import javastraw.reader.block.Block;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.datastructures.ListOfDoubleArrays;
import javastraw.reader.datastructures.ListOfFloatArrays;
import javastraw.reader.mzd.MatrixZoomData;
import org.broad.igv.util.collections.LRUCache;

import java.util.Iterator;

public class ZDIteratorContainer extends IteratorContainer {

    private final LRUCache<String, Block> blockCache;
    private final DatasetReader reader;
    private final MatrixZoomData zd;
    private final boolean useCache;

    public ZDIteratorContainer(DatasetReader reader, MatrixZoomData zd, LRUCache<String, Block> blockCache,
                               boolean useCache) {
        super(zd.getMatrixSize());
        this.reader = reader;
        this.zd = zd;
        this.blockCache = blockCache;
        this.useCache = useCache;
    }

    public static ListOfFloatArrays matrixVectorMultiplyOnIterator(Iterator<ContactRecord> iterator,
                                                                   ListOfFloatArrays vector, long vectorLength) {
        ListOfDoubleArrays sumVector = new ListOfDoubleArrays(vectorLength);
        while (iterator.hasNext()) {
            ContactRecord cr = iterator.next();
            ListIteratorContainer.matrixVectorMult(vector, sumVector, cr);
        }
        return sumVector.convertToFloats();
    }

    @Override
    public Iterator<ContactRecord> getNewContactRecordIterator() {
        return new ContactRecordIterator(reader, zd, blockCache, useCache);
    }

    @Override
    public ListOfFloatArrays sparseMultiply(ListOfFloatArrays vector, long vectorLength) {
        return matrixVectorMultiplyOnIterator(getNewContactRecordIterator(), vector, vectorLength);
    }

    @Override
    public void clear() {
        //blockCache.clear();
    }
}
