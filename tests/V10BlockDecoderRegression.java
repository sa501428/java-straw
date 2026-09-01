import javastraw.reader.v10.V10;
import javastraw.reader.v10.V10BlockDecoder;
import javastraw.reader.v10.V10Cursor;
import javastraw.reader.v10.V10Grid;
import javastraw.reader.v10.V10RecordHandler;
import javastraw.reader.v10.V10Zoom;

import java.io.ByteArrayOutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/** Allocation-free streaming decoder coverage for every representation/value mode. */
public final class V10BlockDecoderRegression {
    private static final int[] POSITIONS = {0, 3, 6};

    public static void main(String[] args) {
        checkExactFarCisRange();
        for (int mode = V10.ALL_DEFAULT; mode <= V10.DIRECT; mode++) {
            check(V10.SPARSE_DELTA, mode, V10.COUNT_UINT);
            check(V10.BITMAP, mode, V10.COUNT_UINT);
        }
        check(V10.DENSE, V10.DIRECT, V10.COUNT_UINT);
        check(V10.DENSE, V10.DIRECT, V10.SCORE_FLOAT32);
        System.out.println("V10 streaming block decoder: all representations and value modes passed");
    }

    private static void checkExactFarCisRange() {
        V10Zoom zoom = new V10Zoom();
        zoom.gridType = V10.GRID_ROTATED_CIS;
        zoom.blockBinCount = 4;
        zoom.blockColumnCount = 20;
        List<V10Grid.BlockRange> ranges = V10Grid.blockRangesForRegion(2, 20, 2, 20,
                zoom, 80, 80);
        if (ranges.size() != 1 || ranges.get(0).first != 42 || ranges.get(0).last != 42) {
            throw new AssertionError("far-cis query selected non-intersecting distance bands");
        }
    }

    private static void check(int representation, int mode, int type) {
        final List<String> actual = new ArrayList<>();
        V10Zoom zoom = new V10Zoom();
        zoom.valueType = type;
        zoom.gridType = V10.GRID_RECTANGULAR;
        zoom.blockBinCount = 4;
        zoom.blockColumnCount = 2;
        V10BlockDecoder.decode(new V10Cursor(block(representation, mode, type)), 0, zoom,
                8, 8, false, new V10RecordHandler() {
                    @Override
                    public void record(int x, int y, long count, float score, boolean isScore) {
                        actual.add(x + "," + y + "=" + (isScore ? Float.toString(score) : Long.toString(count)));
                    }
                });

        String middle = mode == V10.ALL_DEFAULT ? "1" : "5";
        String last = mode == V10.DIRECT ? "7" : "1";
        if (type == V10.SCORE_FLOAT32) {
            middle = "5.0";
            last = "7.0";
        }
        String first = type == V10.SCORE_FLOAT32 ? "1.0" : "1";
        List<String> expected = Arrays.asList("0,0=" + first, "3,0=" + middle, "2,1=" + last);
        if (!expected.equals(actual)) {
            throw new AssertionError("decoder mismatch for representation=" + representation
                    + " mode=" + mode + " type=" + type + ": " + actual);
        }
    }

    private static byte[] block(int representation, int mode, int type) {
        Buffer positions = new Buffer();
        if (representation == V10.SPARSE_DELTA) {
            positions.varint(0);
            positions.varint(3);
            positions.varint(3);
        } else if (representation == V10.BITMAP || type == V10.SCORE_FLOAT32) {
            positions.u8(0x49);
        }

        Buffer values = new Buffer();
        if (representation == V10.DENSE) {
            long[] dense = {1, 0, 0, 5, 0, 0, 7, 0};
            for (long value : dense) scalar(values, type, value);
        } else if (mode == V10.ALL_DEFAULT) {
            scalar(values, type, 1);
        } else if (mode == V10.DEFAULT_EXCEPTIONS) {
            scalar(values, type, 1);
            values.varint(1);
            values.varint(1);
            scalar(values, type, 5);
        } else {
            scalar(values, type, 1);
            scalar(values, type, 5);
            scalar(values, type, 7);
        }

        byte[] positionBytes = positions.toArray();
        byte[] valueBytes = values.toArray();
        Buffer out = new Buffer();
        out.u8(1);
        out.u8(representation);
        out.u8(mode);
        out.u8(type);
        out.u8(representation == V10.BITMAP || type == V10.SCORE_FLOAT32 ? 1 : 0);
        out.raw(new byte[3]);
        out.u32(0);
        out.u32(0);
        out.u32(4);
        out.u32(2);
        out.u64(POSITIONS.length);
        out.u32(positionBytes.length);
        out.u32(valueBytes.length);
        out.raw(positionBytes);
        out.raw(valueBytes);
        return out.toArray();
    }

    private static void scalar(Buffer out, int type, long value) {
        if (type == V10.SCORE_FLOAT32) {
            out.u32(Float.floatToRawIntBits((float) value) & 0xFFFFFFFFL);
        } else {
            out.varint(value);
        }
    }

    private static final class Buffer {
        private final ByteArrayOutputStream out = new ByteArrayOutputStream();

        void u8(int value) { out.write(value); }
        void u32(long value) { little(value, 4); }
        void u64(long value) { little(value, 8); }
        void raw(byte[] bytes) { out.write(bytes, 0, bytes.length); }

        void varint(long value) {
            do {
                int b = (int) (value & 0x7F);
                value >>>= 7;
                out.write(value == 0 ? b : b | 0x80);
            } while (value != 0);
        }

        private void little(long value, int bytes) {
            for (int i = 0; i < bytes; i++) out.write((int) (value >>> (8 * i)) & 0xFF);
        }

        byte[] toArray() { return out.toByteArray(); }
    }
}
