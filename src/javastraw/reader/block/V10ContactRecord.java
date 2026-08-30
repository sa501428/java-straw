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

import java.math.BigInteger;

/**
 * A contact record that also carries its exact integer count.
 * <p>
 * V10 stores raw contacts as unsigned integers and forbids converting them to
 * floating point merely because they are large, so a count above 2^24 can no
 * longer be represented exactly by {@link ContactRecord}'s float. Records with
 * such counts are emitted as this subclass, which keeps the float value for
 * every existing caller while exposing the exact count through
 * {@link #getExactCount()}. Every V10 integer record uses this representation
 * so the exact value is preserved consistently through all code paths.
 */
public class V10ContactRecord extends ContactRecord {

    private final BigInteger exactCount;

    private V10ContactRecord(int binX, int binY, BigInteger exactCount) {
        super(binX, binY, exactCount.floatValue());
        this.exactCount = exactCount;
    }

    /** Create a record from the raw two's-complement bits of a V10 uint64. */
    public static ContactRecord create(int binX, int binY, long unsignedBits) {
        return create(binX, binY, javastraw.reader.v10.V10.unsignedLong(unsignedBits));
    }

    public static ContactRecord create(int binX, int binY, BigInteger count) {
        if (count.signum() <= 0 || count.bitLength() > 64) {
            throw new IllegalArgumentException("V10 count is outside uint64");
        }
        return new V10ContactRecord(binX, binY, count);
    }

    @Override
    public BigInteger getExactCount() {
        return exactCount;
    }

    @Override
    public ContactRecord transpose() {
        return new V10ContactRecord(getBinY(), getBinX(), exactCount);
    }
}
