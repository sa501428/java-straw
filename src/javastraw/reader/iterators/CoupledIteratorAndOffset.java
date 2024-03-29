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

import javastraw.reader.block.ContactRecord;

import java.util.Iterator;

public class CoupledIteratorAndOffset implements Iterator<ContactRecord> {

    private final Iterator<ContactRecord> internalIterator;
    private final int xOffset, yOffset;
    private final boolean isIntra;

    public CoupledIteratorAndOffset(Iterator<ContactRecord> iterator, int xOffset, int yOffset,
                                    boolean isIntra) {
        internalIterator = iterator;
        this.xOffset = xOffset;
        this.yOffset = yOffset;
        this.isIntra = isIntra;
    }

    public boolean getIsIntra() {
        return isIntra;
    }

    @Override
    public boolean hasNext() {
        return internalIterator.hasNext();
    }

    @Override
    public ContactRecord next() {
        ContactRecord cr = internalIterator.next();
        int binX = cr.getBinX() + xOffset;
        int binY = cr.getBinY() + yOffset;
        return new ContactRecord(binX, binY, cr.getCounts());
    }
}
