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

import java.util.ArrayList;
import java.util.List;
import java.util.Objects;

public class Block implements Comparable<Block> {

    private final int number;
    private final String uniqueRegionID;
    private final List<ContactRecord> records;

    public Block(int number, String regionID) {
        this.number = number;
        records = new ArrayList<>();
        uniqueRegionID = regionID + "_" + number;
    }

    public Block(int number, List<ContactRecord> records, String regionID) {
        this.number = number;
        this.records = records;
        this.uniqueRegionID = regionID + "_" + number;
    }

    public int getNumber() {
        return number;
    }

    public String getUniqueRegionID() {
        return uniqueRegionID;
    }

    public List<ContactRecord> getContactRecords() {
        return records;
    }

    @Override
    public int compareTo(Block o) {
        if (this == o) return 0;
        int[] comparisons = new int[]{Integer.compare(number, o.number),
                uniqueRegionID.compareTo(o.uniqueRegionID),
                Integer.compare(records.size(), o.records.size())
        };
        for (int val : comparisons) {
            if (val != 0) {
                return val;
            }
        }
        return 0;
    }

    @Override
    public boolean equals(Object obj) {
        return obj instanceof Block
                && ((Block) obj).number == number
                && ((Block) obj).uniqueRegionID.equals(uniqueRegionID)
                && ((Block) obj).records.size() == records.size();
    }

    @Override
    public int hashCode() {
        return Objects.hash(number, uniqueRegionID, records.size());
    }
}
