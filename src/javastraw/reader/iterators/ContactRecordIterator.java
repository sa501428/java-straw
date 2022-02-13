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
import javastraw.reader.mzd.BlockCache;
import javastraw.reader.mzd.BlockLoader;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationHandler;

import java.io.IOException;
import java.util.Iterator;
import java.util.List;

/**
 * Class for iterating over the contact records
 */
public class ContactRecordIterator implements Iterator<ContactRecord> {

    private final List<Integer> blockNumbers;
    private int blockIdx;
    private Iterator<ContactRecord> currentBlockIterator;
    private final DatasetReader reader;
    private final String zdKey;
    private final BlockCache blockCache;
    private final int chr1Idx, chr2Idx;
    private final HiCZoom zoom;

    /**
     * Initializes the iterator
     */
    public ContactRecordIterator(DatasetReader reader, String zdKey, BlockCache blockCache,
                                 int chr1Idx, int chr2Idx, HiCZoom zoom) {
        this.reader = reader;
        this.zdKey = zdKey;
        this.chr1Idx = chr1Idx;
        this.chr2Idx = chr2Idx;
        this.zoom = zoom;
        this.blockCache = blockCache;
        this.blockIdx = -1;
        this.blockNumbers = reader.getBlockNumbers(zdKey);
    }

    /**
     * Indicates whether or not there is another block waiting; checks current block
     * iterator and creates a new one if need be
     *
     * @return true if there is another block to be read
     */
    @Override
    public boolean hasNext() {

        if (currentBlockIterator != null && currentBlockIterator.hasNext()) {
            return true;
        } else {
            blockIdx++;
            if (blockIdx < blockNumbers.size()) {
                try {
                    int blockNumber = blockNumbers.get(blockIdx);

                    // Optionally check the cache
                    String key = BlockLoader.getBlockKey(zdKey, blockNumber, NormalizationHandler.NONE);
                    Block nextBlock;
                    if (blockCache.containsKey(key)) {
                        nextBlock = blockCache.get(key);
                    } else {
                        nextBlock = reader.readNormalizedBlock(blockNumber, zdKey, NormalizationHandler.NONE,
                                chr1Idx, chr2Idx, zoom);
                    }
                    currentBlockIterator = nextBlock.getContactRecords().iterator();
                    return true;
                } catch (IOException e) {
                    System.err.println("Error fetching block " + e.getMessage());
                    return false;
                }
            }
        }

        return false;
    }

    /**
     * Returns the next contact record
     *
     * @return The next contact record
     */
    @Override
    public ContactRecord next() {
        return currentBlockIterator == null ? null : currentBlockIterator.next();
    }

    /**
     * Not supported
     */
    @Override
    public void remove() {
        //Not supported
        throw new RuntimeException("remove() is not supported");
    }
}
