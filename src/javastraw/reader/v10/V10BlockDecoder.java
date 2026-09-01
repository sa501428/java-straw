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

        require(representation != V10.DENSE || valueMode == V10.DIRECT, "dense values must be direct");
        ValueReader valueReader = new ValueReader(values, valueMode, valueType, slots);

        if (representation == V10.SPARSE_DELTA) {
            require(flags == 0 && occupied <= positionStreamBytes, "invalid sparse position stream");
            long previous = 0;
            for (int i = 0; i < (int) occupied; i++) {
                long delta = positions.varint();
                require(i == 0 || delta > 0, "duplicate sparse cell");
                long p = i == 0 ? delta : V10.add(previous, delta);
                require(p < cells, "sparse position out of bounds");
                long value = valueReader.value(i);
                require(valueType != V10.COUNT_UINT || value != 0,
                        "sparse/bitmap count must be positive");
                emit(p, value, h, blockNumber, zoom, nBins1, nBins2, isIntra, handler);
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
            long found = 0;
            for (long i = 0; i < cells; i++) {
                boolean present = (positions.peekAt(i / 8) & (1 << (int) (i % 8))) != 0;
                if (representation == V10.DENSE) {
                    long value = valueReader.value(i);
                    if (present) {
                        require(found < occupied, "bitmap population mismatch");
                        found++;
                        emit(i, value, h, blockNumber, zoom, nBins1, nBins2, isIntra, handler);
                    } else {
                        require(value == 0, "absent dense score must be positive zero");
                    }
                } else if (present) {
                    require(found < occupied, "bitmap population mismatch");
                    long value = valueReader.value(found++);
                    require(valueType != V10.COUNT_UINT || value != 0,
                            "sparse/bitmap count must be positive");
                    emit(i, value, h, blockNumber, zoom, nBins1, nBins2, isIntra, handler);
                }
            }
            require(found == occupied, "bitmap population mismatch");
            positions.skipToEnd();
        } else {
            require(flags == 0 && positionStreamBytes == 0, "dense counts have no presence stream");
            long emitted = 0;
            for (long i = 0; i < slots; i++) {
                long value = valueReader.value(i);
                if (value != 0) {
                    emit(i, value, h, blockNumber, zoom, nBins1, nBins2, isIntra, handler);
                    emitted++;
                }
            }
            require(emitted == occupied, "occupied cell count mismatch");
        }
        positions.done();
        valueReader.done();
    }

    private static void emit(long position, long value, V10BlockHeader h, int blockNumber,
                             V10Zoom zoom, long nBins1, long nBins2, boolean isIntra,
                             V10RecordHandler handler) {
        long binColumn = h.binColumnOffset + position % h.blockWidth;
        long binRow = h.binRowOffset + position / h.blockWidth;
        require(binColumn < nBins1 && binRow < nBins2, "cell violates block geometry");
        require(!isIntra || binRow >= binColumn, "cell violates the canonical cis triangle");
        int x = V10.toInt(binColumn, "bin column");
        int y = V10.toInt(binRow, "bin row");
        require(V10Grid.blockNumber(x, y, zoom) == blockNumber, "cell violates block geometry");
        if (h.valueType == V10.SCORE_FLOAT32) {
            handler.record(x, y, 0, V10.bitsToFloat(value), true);
        } else {
            handler.record(x, y, value, 0f, false);
        }
    }

    /** Streaming value decoder; only exception ordinals require temporary storage. */
    private static final class ValueReader {
        private final V10Cursor values;
        private final int mode;
        private final int type;
        private final long defaultValue;
        private final long[] exceptionOrdinals;
        private int nextException;

        ValueReader(V10Cursor values, int mode, int type, long slots) {
            this.values = values;
            this.mode = mode;
            this.type = type;
            if (mode == V10.ALL_DEFAULT) {
                require(slots > 0, "all-default mode requires at least one value slot");
                defaultValue = scalar(values, type);
                exceptionOrdinals = null;
            } else if (mode == V10.DEFAULT_EXCEPTIONS) {
                defaultValue = scalar(values, type);
                long count = values.varint();
                require(count > 0 && count < slots && count <= values.left(), "invalid exception count");
                exceptionOrdinals = new long[(int) count];
                long previous = 0;
                for (int i = 0; i < exceptionOrdinals.length; i++) {
                    long delta = values.varint();
                    require(i == 0 || delta > 0, "duplicate exception ordinal");
                    long ordinal = i == 0 ? delta : V10.add(previous, delta);
                    require(ordinal < slots, "exception out of range");
                    exceptionOrdinals[i] = ordinal;
                    previous = ordinal;
                }
            } else {
                defaultValue = 0;
                exceptionOrdinals = null;
                require(slots <= values.left() / (type == V10.SCORE_FLOAT32 ? 4 : 1),
                        "truncated values");
            }
        }

        long value(long ordinal) {
            if (mode == V10.DIRECT) return scalar(values, type);
            if (mode == V10.DEFAULT_EXCEPTIONS && nextException < exceptionOrdinals.length
                    && exceptionOrdinals[nextException] == ordinal) {
                nextException++;
                long exception = scalar(values, type);
                require(exception != defaultValue, "exception equals default");
                return exception;
            }
            return defaultValue;
        }

        void done() {
            require(exceptionOrdinals == null || nextException == exceptionOrdinals.length,
                    "unconsumed exception ordinal");
            values.done();
        }
    }

    /**
     * A count scalar is canonical ULEB128; a score scalar is the four raw
     * little-endian bytes of its f32 bit pattern (Section I.5).
     */
    private static long scalar(V10Cursor values, int valueType) {
        return valueType == V10.SCORE_FLOAT32 ? values.word() : values.unsignedVarintBits();
    }
}
