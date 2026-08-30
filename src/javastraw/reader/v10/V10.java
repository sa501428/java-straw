package javastraw.reader.v10;

/**
 * Enumerations, limits and checked arithmetic for the V10 wire format.
 * See hic-format/HiCFormatV10.md Section A.
 */
public final class V10 {

    public static final int VERSION = 10;

    // Unit
    public static final int UNIT_BP = 0;
    public static final int UNIT_FRAG = 1;

    // StorageMode
    public static final int MATERIALIZED = 0;
    public static final int DERIVED = 1;

    // Aggregation
    public static final int AGGREGATION_NONE = 0;
    public static final int AGGREGATION_SUM = 1;

    // ValueType
    public static final int COUNT_UINT = 0;
    public static final int SCORE_FLOAT32 = 1;

    // GridType
    public static final int GRID_RECTANGULAR = 0;
    public static final int GRID_ROTATED_CIS = 1;

    // BlockRepresentation
    public static final int SPARSE_DELTA = 0;
    public static final int BITMAP = 1;
    public static final int DENSE = 2;

    // ValueMode
    public static final int ALL_DEFAULT = 0;
    public static final int DEFAULT_EXCEPTIONS = 1;
    public static final int DIRECT = 2;

    // VectorTransform
    public static final int TRANSFORM_RAW = 0;
    public static final int TRANSFORM_BYTE_SHUFFLE = 1;
    public static final int TRANSFORM_XOR32 = 2;

    // PageCodec
    public static final int CODEC_ZSTD = 1;

    public static final int NO_SOURCE_RESOLUTION = 0xFFFFFFFF;

    /**
     * Per-record allocation ceiling; independent of the size of the .hic file.
     */
    public static final long ALLOCATION_LIMIT = 512L * 1024 * 1024;

    private V10() {
    }

    public static void require(boolean ok, String message) {
        if (!ok) {
            throw new V10FormatException(message);
        }
    }

    /**
     * Checked addition of two non-negative values.
     */
    public static long add(long a, long b) {
        require(a >= 0 && b >= 0 && b <= Long.MAX_VALUE - a, "integer overflow");
        return a + b;
    }

    /**
     * Checked multiplication of two non-negative values.
     */
    public static long multiply(long a, long b) {
        require(a >= 0 && b >= 0, "integer overflow");
        require(a == 0 || b <= Long.MAX_VALUE / a, "integer overflow");
        return a * b;
    }

    public static long unsigned32(long value) {
        require(value >= 0 && value <= 0xFFFFFFFFL, "value exceeds uint32");
        return value;
    }

    /**
     * Java exposes bin coordinates and block numbers as signed ints. The format
     * permits the full uint32 range, so files beyond that are rejected cleanly
     * rather than silently wrapping.
     */
    public static int toInt(long value, String what) {
        require(value >= 0 && value <= Integer.MAX_VALUE,
                what + " (" + value + ") exceeds the signed 32-bit range supported by this API");
        return (int) value;
    }

    public static float bitsToFloat(long bits) {
        require(bits >= 0 && bits <= 0xFFFFFFFFL, "value exceeds uint32");
        return Float.intBitsToFloat((int) bits);
    }

    public static String unitName(int unit) {
        return unit == UNIT_FRAG ? "FRAG" : "BP";
    }

    public static int unitId(String unit) {
        if ("BP".equalsIgnoreCase(unit)) return UNIT_BP;
        if ("FRAG".equalsIgnoreCase(unit)) return UNIT_FRAG;
        throw new V10FormatException("unit must be BP or FRAG");
    }
}
