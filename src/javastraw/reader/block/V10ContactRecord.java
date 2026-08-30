/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2011-2024 Broad Institute, Aiden Lab, Rice University, Baylor College of Medicine
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

/**
 * A contact record that also carries its exact integer count.
 * <p>
 * V10 stores raw contacts as unsigned integers and forbids converting them to
 * floating point merely because they are large, so a count above 2^24 can no
 * longer be represented exactly by {@link ContactRecord}'s float. Records with
 * such counts are emitted as this subclass, which keeps the float value for
 * every existing caller while exposing the exact count through
 * {@link #getExactCount()}. Counts that a float represents exactly still use
 * the plain {@link ContactRecord}, so nothing changes for ordinary data.
 */
public class V10ContactRecord extends ContactRecord {

    private final long exactCount;

    private V10ContactRecord(int binX, int binY, long exactCount) {
        super(binX, binY, (float) exactCount);
        this.exactCount = exactCount;
    }

    /**
     * Returns a plain {@link ContactRecord} when the count is exactly
     * representable as a float, and an exact-count record otherwise.
     */
    public static ContactRecord create(int binX, int binY, long count) {
        if ((long) (float) count == count) {
            return new ContactRecord(binX, binY, (float) count);
        }
        return new V10ContactRecord(binX, binY, count);
    }

    @Override
    public long getExactCount() {
        return exactCount;
    }
}
