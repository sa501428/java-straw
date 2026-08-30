package javastraw.reader.v10;

import java.util.ArrayList;
import java.util.List;

import static javastraw.reader.v10.V10.require;

/**
 * One 76-byte matrix resolution descriptor (Section E.2), plus the parser for
 * the enclosing matrix metadata record (Section E.1).
 */
public class V10Zoom {

    public int unit;
    public int storageMode;
    public int aggregation;
    public int valueType;
    public int resolutionIndex;
    public int binSize;
    public long sourceResolutionIndex;
    public int gridType;
    /**
     * Raw 8 bytes: an exact uint64 count sum, or the bits of an f64 score sum.
     */
    public long sumBits;
    public long occupiedCellCount;
    public float stdDev;
    public float percent95;
    /**
     * Logical block scale in bins.
     */
    public int blockBinCount;
    public int blockColumnCount;
    public V10Locator pageIndex;
    public int pageCount;
    public int logicalBlockCount;

    public boolean isDerived() {
        return storageMode == V10.DERIVED;
    }

    public boolean isScore() {
        return valueType == V10.SCORE_FLOAT32;
    }

    public boolean isRotatedCis() {
        return gridType == V10.GRID_ROTATED_CIS;
    }

    /**
     * Exact contact-count sum for COUNT_UINT matrices, or the stored f64 score
     * sum for SCORE_FLOAT32 matrices.
     */
    public double getSum() {
        return isScore() ? Double.longBitsToDouble(sumBits) : (double) sumBits;
    }

    public long getCountSum() {
        require(!isScore(), "score matrices have no integer count sum");
        return sumBits;
    }

    /**
     * Parses a complete matrix metadata record for the given chromosome pair.
     */
    public static List<V10Zoom> parseMatrixRecord(byte[] bytes, int chr1, int chr2, V10Header header,
                                                  V10Source source) {
        V10Cursor c = new V10Cursor(bytes);
        c.magic("H10M");
        require(c.word() == 1, "unknown matrix version");
        require(c.word() == chr1 && c.word() == chr2, "matrix key mismatch");
        long n = c.word();
        c.zero(4);
        int bpCount = header.resolutions[V10.UNIT_BP].size();
        require(n == header.totalResolutionCount() && n * 76 == c.left(),
                "matrix descriptor count/length mismatch");

        List<V10Zoom> result = new ArrayList<>((int) n);
        for (int i = 0; i < (int) n; i++) {
            V10Zoom z = new V10Zoom();
            z.unit = c.byteValue();
            z.storageMode = c.byteValue();
            z.aggregation = c.byteValue();
            z.valueType = c.byteValue();
            z.resolutionIndex = c.wordAsInt("resolution index");
            z.binSize = c.wordAsInt("bin size");
            z.sourceResolutionIndex = c.word();
            z.gridType = c.byteValue();
            c.zero(3);
            z.sumBits = c.integer(8);
            z.occupiedCellCount = c.wide();
            z.stdDev = V10.bitsToFloat(c.word());
            z.percent95 = V10.bitsToFloat(c.word());
            z.blockBinCount = c.wordAsInt("block bin count");
            z.blockColumnCount = c.wordAsInt("block column count");
            z.pageIndex = V10Locator.read(c);
            z.pageCount = c.wordAsInt("page count");
            z.logicalBlockCount = c.wordAsInt("logical block count");

            int unit = i < bpCount ? V10.UNIT_BP : V10.UNIT_FRAG;
            int ri = i - (unit == V10.UNIT_FRAG ? bpCount : 0);
            V10Resolution r = header.resolutions[unit].get(ri);
            require(z.unit == unit && z.resolutionIndex == ri && z.binSize == r.binSize
                            && z.storageMode == r.storageMode && z.aggregation == r.aggregation
                            && z.sourceResolutionIndex == r.sourceResolutionIndex,
                    "resolution descriptor mismatch");
            require(z.valueType <= 1 && z.gridType == (chr1 == chr2 ? V10.GRID_ROTATED_CIS : V10.GRID_RECTANGULAR)
                            && z.blockBinCount > 0 && z.blockColumnCount > 0,
                    "invalid matrix geometry/type");
            long expectedColumns = (header.bins(chr1, unit, ri) + z.blockBinCount - 1) / z.blockBinCount;
            require(z.blockColumnCount == expectedColumns, "invalid block column count");
            if (z.isDerived()) {
                require(z.pageIndex.position == 0 && z.pageCount == 0 && z.logicalBlockCount == 0,
                        "derived resolution has storage");
            } else if (z.occupiedCellCount > 0) {
                require(z.pageIndex.position != 0 && z.pageCount > 0 && z.logicalBlockCount > 0,
                        "missing matrix pages");
            } else {
                require(z.pageCount == 0 && z.logicalBlockCount == 0, "empty matrix has pages");
            }
            if (z.pageIndex.isPresent()) {
                source.interval(z.pageIndex.position, z.pageIndex.length);
            }
            result.add(z);
        }
        c.done();

        for (V10Zoom z : result) {
            if (z.isDerived()) {
                int base = z.unit == V10.UNIT_FRAG ? bpCount : 0;
                require(z.valueType == result.get(base + (int) z.sourceResolutionIndex).valueType,
                        "derived value type mismatch");
            }
        }
        return result;
    }
}
