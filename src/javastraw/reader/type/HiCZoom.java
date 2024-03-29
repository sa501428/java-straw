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


package javastraw.reader.type;

import java.util.Objects;

public class HiCZoom implements Comparable<HiCZoom> {

    private final HiCUnit unit;
    private final Integer binSize;

    public HiCZoom(HiCUnit unit, int binSize) {
        this.unit = unit;
        this.binSize = binSize;
    }

    // can assume this default given FRAG resolution is no longer used
    public HiCZoom(int binSize) {
        this.unit = HiCUnit.BP;
        this.binSize = binSize;
    }

    @SuppressWarnings("MethodDoesntCallSuperMethod")
    public HiCZoom clone() {
        return new HiCZoom(unit, binSize);
    }

    public static HiCUnit valueOfUnit(String unit) {
        if (unit.equalsIgnoreCase(HiCUnit.BP.toString())) {
            return HiCUnit.BP;
        } else if (unit.equalsIgnoreCase(HiCUnit.FRAG.toString())) {
            return HiCUnit.FRAG;
        }
        return null;
    }

    public int getBinSize() {
        return binSize;
    }

    public String getKey() {
        return unit.toString() + "_" + binSize;
    }

    public String toString() {
        return getKey();
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o instanceof HiCZoom) {
            HiCZoom hiCZoom = (HiCZoom) o;
            return (binSize.equals(hiCZoom.binSize)) && (unit == hiCZoom.unit);
        }
        return false;
    }

    @Override
    public int hashCode() {
        return Objects.hash(unit.hashCode(), binSize);
    }

    @Override
    public int compareTo(HiCZoom o) {
        return binSize.compareTo(o.binSize);
    }

    public HiCUnit getUnit() {
        return unit;
    }

    public enum HiCUnit {BP, FRAG}
}
