package javastraw.reader.v10;

import static javastraw.reader.v10.V10.require;

/**
 * Decodes one V10 logical block (Section I).
 * <p>
 * All three representations (sparse-delta, bitmap, dense) and all three value
 * modes (all-default, default-with-exceptions, direct) are supported. Every
 * reconstructed cell is validated against the chromosome bin counts, the
 * canonical cis triangle, and the block number under the declared grid.
 */
public final class V10BlockDecoder {

    private V10BlockDecoder() {
    }

    /**
     * Reads and validates only the fixed 40-byte block header, leaving the
     * cursor positioned at the start of the position stream.
     */
    public static V10BlockHeader readHeader(V10Cursor c, V10Zoom zoom, long nBins1, long nBins2) {
        require(c.byteValue() == 1, "unknown block version");
        int representation = c.byteValue();
        int valueMode = c.byteValue();
        int valueType = c.byteValue();
        int flags = c.byteValue();
        c.zero(3);
        long binColumnOffset = c.word();
        long binRowOffset = c.word();
        long width = c.word();
        long height = c.word();
        long occupied = c.wide();
        long positionStreamBytes = c.word();
        long valueStreamBytes = c.word();

        require(representation <= 2 && valueMode <= 2 && valueType == zoom.valueType && flags <= 1
                        && width > 0 && height > 0 && occupied > 0,
                "invalid block header");
        long cells = V10.multiply(width, height);
        long slots = representation == V10.DENSE ? cells : occupied;
        require(occupied <= cells && slots <= V10.ALLOCATION_LIMIT / 8
                        && V10.add(positionStreamBytes, valueStreamBytes) == c.left(),
                "invalid block stream lengths/size");
        require(binColumnOffset < nBins1 && binRowOffset < nBins2, "block offsets exceed chromosome");
        return new V10BlockHeader(representation, valueMode, valueType, flags, binColumnOffset,
                binRowOffset, width, height, occupied, positionStreamBytes, valueStreamBytes);
    }

    /**
     * @param c           cursor positioned at the start of the block, sized to exactly this block
     * @param blockNumber logical block number this payload was indexed under
     * @param zoom        the resolution descriptor the block belongs to
     * @param nBins1      bin count of chromosome 1 (columns) at this resolution
     * @param nBins2      bin count of chromosome 2 (rows) at this resolution
     * @param isIntra     whether the matrix is cis, in which case binRow >= binColumn
     */
    public static void decode(V10Cursor c, int blockNumber, V10Zoom zoom, long nBins1, long nBins2,
                              boolean isIntra, V10RecordHandler handler) {
        V10BlockHeader h = readHeader(c, zoom, nBins1, nBins2);
        int representation = h.representation;
        int valueMode = h.valueMode;
        int valueType = h.valueType;
        int flags = h.flags;
        long binColumnOffset = h.binColumnOffset;
        long binRowOffset = h.binRowOffset;
        long width = h.blockWidth;
        long occupied = h.occupiedCellCount;
        long positionStreamBytes = h.positionStreamBytes;
        long valueStreamBytes = h.valueStreamBytes;
        long cells = V10.multiply(width, h.blockHeight);
        long slots = representation == V10.DENSE ? cells : occupied;

        V10Cursor positions = c.take(positionStreamBytes);
        V10Cursor values = c.take(valueStreamBytes);

        long[] occupiedPositions = null;
        if (representation == V10.SPARSE_DELTA) {
            require(flags == 0 && occupied <= positionStreamBytes, "invalid sparse position stream");
            occupiedPositions = new long[(int) occupied];
            long previous = 0;
            for (int i = 0; i < (int) occupied; i++) {
                long delta = positions.varint();
                require(i == 0 || delta > 0, "duplicate sparse cell");
                long p = i == 0 ? delta : V10.add(previous, delta);
                require(p < cells, "sparse position out of bounds");
                occupiedPositions[i] = p;
                previous = p;
            }
        } else if (representation == V10.BITMAP || valueType == V10.SCORE_FLOAT32) {
            // Bitmap blocks, and dense score blocks, both carry an explicit
            // presence bitmap so absence stays distinguishable from a zero or NaN.
            require(flags == 1 && positionStreamBytes == (cells + 7) / 8, "invalid presence bitmap");
            if (cells % 8 != 0) {
                require((positions.peekAt(positionStreamBytes - 1) >>> (cells % 8)) == 0,
                        "nonzero bitmap padding");
            }
            occupiedPositions = new long[(int) occupied];
            int found = 0;
            for (long i = 0; i < cells; i++) {
                if ((positions.peekAt(i / 8) & (1 << (int) (i % 8))) != 0) {
                    require(found < occupied, "bitmap population mismatch");
                    occupiedPositions[found++] = i;
                }
            }
            require(found == occupied, "bitmap population mismatch");
            positions.skipToEnd();
        } else {
            require(flags == 0 && positionStreamBytes == 0, "dense counts have no presence stream");
        }
        positions.done();
        require(representation != V10.DENSE || valueMode == V10.DIRECT, "dense values must be direct");

        long[] decoded = decodeValues(values, valueMode, valueType, slots);
        values.done();

        long emitted = 0;
        int nextOccupied = 0;
        for (long i = 0; i < slots; i++) {
            boolean present = true;
            long position = representation == V10.DENSE ? i : occupiedPositions[(int) i];
            long value = decoded[(int) i];
            if (representation == V10.DENSE) {
                if (valueType == V10.COUNT_UINT) {
                    present = value != 0;
                } else {
                    present = nextOccupied < occupiedPositions.length && occupiedPositions[nextOccupied] == i;
                    if (present) {
                        nextOccupied++;
                    } else {
                        require(value == 0, "absent dense score must be positive zero");
                    }
                }
            } else if (valueType == V10.COUNT_UINT) {
                require(value > 0, "sparse/bitmap count must be positive");
            }
            if (!present) continue;

            long binColumn = binColumnOffset + position % width;
            long binRow = binRowOffset + position / width;
            require(binColumn < nBins1 && binRow < nBins2, "cell violates block geometry");
            require(!isIntra || binRow >= binColumn, "cell violates the canonical cis triangle");
            int x = V10.toInt(binColumn, "bin column");
            int y = V10.toInt(binRow, "bin row");
            require(V10Grid.blockNumber(x, y, zoom) == blockNumber, "cell violates block geometry");

            if (valueType == V10.SCORE_FLOAT32) {
                handler.record(x, y, 0, V10.bitsToFloat(value), true);
            } else {
                handler.record(x, y, value, 0f, false);
            }
            emitted++;
        }
        require(emitted == occupied, "occupied cell count mismatch");
    }

    private static long[] decodeValues(V10Cursor values, int valueMode, int valueType, long slots) {
        long[] decoded = new long[(int) slots];
        if (valueMode == V10.ALL_DEFAULT) {
            require(slots > 0, "all-default mode requires at least one value slot");
            long defaultValue = scalar(values, valueType);
            for (int i = 0; i < decoded.length; i++) {
                decoded[i] = defaultValue;
            }
        } else if (valueMode == V10.DEFAULT_EXCEPTIONS) {
            long defaultValue = scalar(values, valueType);
            long exceptionCount = values.varint();
            require(exceptionCount > 0 && exceptionCount < slots && exceptionCount <= values.left(),
                    "invalid exception count");
            long[] ordinals = new long[(int) exceptionCount];
            long previous = 0;
            for (int i = 0; i < ordinals.length; i++) {
                long delta = values.varint();
                require(i == 0 || delta > 0, "duplicate exception ordinal");
                long ordinal = i == 0 ? delta : V10.add(previous, delta);
                require(ordinal < slots, "exception out of range");
                ordinals[i] = ordinal;
                previous = ordinal;
            }
            for (int i = 0; i < decoded.length; i++) {
                decoded[i] = defaultValue;
            }
            for (long ordinal : ordinals) {
                long exception = scalar(values, valueType);
                require(exception != defaultValue, "exception equals default");
                decoded[(int) ordinal] = exception;
            }
        } else {
            require(slots <= values.left() / (valueType == V10.SCORE_FLOAT32 ? 4 : 1), "truncated values");
            for (int i = 0; i < decoded.length; i++) {
                decoded[i] = scalar(values, valueType);
            }
        }
        return decoded;
    }

    /**
     * A count scalar is canonical ULEB128; a score scalar is the four raw
     * little-endian bytes of its f32 bit pattern (Section I.5).
     */
    private static long scalar(V10Cursor values, int valueType) {
        return valueType == V10.SCORE_FLOAT32 ? values.word() : values.varint();
    }
}
