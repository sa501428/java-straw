import javastraw.reader.v10.V10;
import javastraw.reader.v10.V10Header;
import javastraw.reader.v10.V10Resolution;

import java.io.ByteArrayOutputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.charset.StandardCharsets;

/** Regression coverage for user-selected V10 materialization layouts. */
public final class V10ResolutionPolicyRegression {

    private static final int[] BINS = {50, 100, 150, 200, 300, 500, 500000};

    public static void main(String[] args) {
        V10Header header = V10Header.parse(customHeader());
        assertEquals("resolution count", BINS.length, header.resolutions[V10.UNIT_BP].size());
        for (int i = 0; i < BINS.length; i++) {
            V10Resolution resolution = header.resolutions[V10.UNIT_BP].get(i);
            assertEquals("bin size " + i, BINS[i], resolution.binSize);
            assertEquals("storage mode " + i, i == 0 ? V10.MATERIALIZED : V10.DERIVED,
                    resolution.storageMode);
            assertEquals("source " + i, i == 0 ? 0xFFFFFFFFL : 0L,
                    resolution.sourceResolutionIndex);
        }
        System.out.println("V10 resolution policy: custom materialized base and direct derivations passed");
    }

    private static byte[] customHeader() {
        Buffer out = new Buffer();
        out.raw(new byte[88]);
        out.string("testGenome");
        out.u32(0);                         // attributes
        out.u32(1);                         // chromosomes
        out.string("chr1");
        out.u64(1_000_000);
        out.u32(BINS.length);               // BP resolutions
        for (int i = 0; i < BINS.length; i++) {
            out.u32(BINS[i]);
            out.u8(i == 0 ? V10.MATERIALIZED : V10.DERIVED);
            out.u8(V10.AGGREGATION_SUM);
            out.u16(0);
            out.u32(i == 0 ? 0xFFFFFFFFL : 0); // every derived level reads 50 bp directly
        }
        out.u32(0);                         // FRAG resolutions
        out.u32(0);                         // normalizations

        byte[] bytes = out.toArray();
        ByteBuffer prefix = ByteBuffer.wrap(bytes).order(ByteOrder.LITTLE_ENDIAN);
        prefix.put("HIC\0".getBytes(StandardCharsets.UTF_8));
        prefix.putInt(V10.VERSION);
        prefix.putLong(bytes.length);
        prefix.putLong(bytes.length);       // footer begins immediately after this header
        prefix.putLong(24);                 // non-empty footer locator
        return bytes;
    }

    private static void assertEquals(String what, long expected, long actual) {
        if (expected != actual) {
            throw new AssertionError(what + ": expected " + expected + ", got " + actual);
        }
    }

    private static final class Buffer {
        private final ByteArrayOutputStream out = new ByteArrayOutputStream();

        void raw(byte[] value) {
            out.write(value, 0, value.length);
        }

        void u8(int value) {
            out.write(value & 0xFF);
        }

        void u16(int value) {
            u8(value);
            u8(value >>> 8);
        }

        void u32(long value) {
            for (int i = 0; i < 4; i++) u8((int) (value >>> (8 * i)));
        }

        void u64(long value) {
            for (int i = 0; i < 8; i++) u8((int) (value >>> (8 * i)));
        }

        void string(String value) {
            raw(value.getBytes(StandardCharsets.UTF_8));
            u8(0);
        }

        byte[] toArray() {
            return out.toByteArray();
        }
    }
}
