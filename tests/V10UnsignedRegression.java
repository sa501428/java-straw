import javastraw.reader.block.ContactRecord;
import javastraw.reader.block.V10ContactRecord;
import javastraw.reader.v10.V10;
import javastraw.reader.v10.V10BlockDecoder;
import javastraw.reader.v10.V10Cursor;
import javastraw.reader.v10.V10Zoom;

import java.math.BigInteger;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;

public final class V10UnsignedRegression {
    public static void main(String[] args) {
        byte[] maximumUleb128 = {
            (byte) 0xff, (byte) 0xff, (byte) 0xff, (byte) 0xff, (byte) 0xff,
            (byte) 0xff, (byte) 0xff, (byte) 0xff, (byte) 0xff, 0x01
        };
        long bits = new V10Cursor(maximumUleb128).unsignedVarintBits();
        assertEquals(V10.UINT64_MAX, V10.unsignedLong(bits));

        ContactRecord record = V10ContactRecord.create(3, 7, bits);
        assertEquals(V10.UINT64_MAX, record.getExactCount());
        assertEquals(V10.UINT64_MAX, record.transpose().getExactCount());

        ByteBuffer block = ByteBuffer.allocate(51).order(ByteOrder.LITTLE_ENDIAN);
        block.put((byte) 1);                 // block version
        block.put((byte) V10.SPARSE_DELTA);
        block.put((byte) V10.DIRECT);
        block.put((byte) V10.COUNT_UINT);
        block.put((byte) 0).put(new byte[3]);
        block.putInt(0).putInt(0);           // bin offsets
        block.putInt(1).putInt(1);           // width, height
        block.putLong(1);                    // one occupied cell
        block.putInt(1).putInt(10);          // stream lengths
        block.put((byte) 0);                 // sparse position zero
        block.put(maximumUleb128);
        V10Zoom zoom = new V10Zoom();
        zoom.valueType = V10.COUNT_UINT;
        zoom.gridType = V10.GRID_RECTANGULAR;
        zoom.blockBinCount = 1;
        zoom.blockColumnCount = 1;
        final long[] decodedBits = new long[1];
        V10BlockDecoder.decode(new V10Cursor(block.array()), 0, zoom, 1, 1, false,
            (x, y, count, score, isScore) -> decodedBits[0] = count);
        assertEquals(V10.UINT64_MAX, V10.unsignedLong(decodedBits[0]));

        BigInteger highBit = BigInteger.ONE.shiftLeft(63);
        assertEquals(highBit, V10.addUnsigned(BigInteger.ZERO, Long.MIN_VALUE));
        assertEquals(V10.UINT64_MAX,
                     V10.addUnsigned(highBit.subtract(BigInteger.ONE), Long.MIN_VALUE));

        boolean overflowRejected = false;
        try {
            V10.addUnsigned(BigInteger.ONE, -1L);
        } catch (RuntimeException expected) {
            overflowRejected = true;
        }
        if (!overflowRejected) throw new AssertionError("uint64 aggregation overflow was accepted");
    }

    private static void assertEquals(Object expected, Object actual) {
        if (!expected.equals(actual)) {
            throw new AssertionError("expected " + expected + ", got " + actual);
        }
    }
}
