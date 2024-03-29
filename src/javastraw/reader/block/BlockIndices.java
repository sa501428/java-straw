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

package javastraw.reader.block;

import htsjdk.tribble.util.LittleEndianInputStream;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class BlockIndices {
    protected final Map<Integer, IndexEntry> blockIndices;
    protected final int numBlocks;

    public BlockIndices(int nBlocks) {
        numBlocks = nBlocks;
        blockIndices = new HashMap<>(nBlocks);
    }

    public void populateBlocks(LittleEndianInputStream dis) throws IOException {
        for (int b = 0; b < numBlocks; b++) {
            int blockNumber = dis.readInt();
            long filePosition = dis.readLong();
            int blockSizeInBytes = dis.readInt();
            blockIndices.put(blockNumber, new IndexEntry(filePosition, blockSizeInBytes));
        }
    }

    public List<Integer> getBlockNumbers() {
        return new ArrayList<>(blockIndices.keySet());
    }

    public Integer getBlockSize(int num) {
        if (blockIndices.containsKey(num)) {
            return blockIndices.get(num).size;
        } else {
            return null;
        }
    }

    public IndexEntry getBlock(int blockNumber) {
        return blockIndices.get(blockNumber);
    }

    public void clearCache() {
        blockIndices.clear();
    }
}
