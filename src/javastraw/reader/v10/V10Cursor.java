package javastraw.reader.v10;

import java.io.UnsupportedEncodingException;

import static javastraw.reader.v10.V10.require;

/**
 * Bounds-checked little-endian reader over an in-memory V10 record.
 * Every accessor validates that the requested bytes remain inside the record,
 * so a truncated or hostile file cannot drive an out-of-range read.
 * See hic-format/HiCFormatV10.md Sections A.1-A.3 and N.
 */
public final class V10Cursor {

    private static final int MAX_STRING_LENGTH = 1024 * 1024;

    private final byte[] data;
    private final int start;
    private final int end;
    private int at;

    public V10Cursor(byte[] data) {
        this(data, 0, data.length);
    }

    public V10Cursor(byte[] data, int start, int length) {
        require(start >= 0 && length >= 0 && start + length <= data.length && start + length >= start,
                "invalid cursor window");
        this.data = data;
        this.start = start;
        this.end = start + length;
        this.at = start;
    }

    public int size() {
        return end - start;
    }

    public int left() {
        return end - at;
    }

    /**
     * Number of bytes consumed since the start of this window.
     */
    public int consumed() {
        return at - start;
    }

    /**
     * Absolute index of the cursor inside the backing array.
     */
    public int rawPosition() {
        return at;
    }

    public byte[] rawArray() {
        return data;
    }

    public void need(long n) {
        require(n >= 0 && n <= left(), "truncated record");
    }

    public long integer(int nBytes) {
        need(nBytes);
        long value = 0;
        for (int i = 0; i < nBytes; i++) {
            value |= (long) (data[at++] & 0xFF) << (8 * i);
        }
        return value;
    }

    public int byteValue() {
        return (int) integer(1);
    }

    /**
     * u32 widened into a long; always non-negative.
     */
    public long word() {
        return integer(4);
    }

    public int wordAsInt(String what) {
        return V10.toInt(word(), what);
    }

    /**
     * u64. Values above 2^63-1 are rejected because Java has no unsigned long
     * and the format already recommends that signed-language readers do so.
     */
    public long wide() {
        long value = integer(8);
        require(value >= 0, "value exceeds the signed 64-bit range supported by this API");
        return value;
    }

    /**
     * Canonical unsigned LEB128 (Section A.2).
     */
    public long varint() {
        long value = varintBits();
        require(value >= 0, "ULEB128 value exceeds the signed 64-bit range supported by this API");
        return value;
    }

    /** Read a canonical uint64 ULEB128 and return its raw long bits. */
    public long unsignedVarintBits() {
        return varintBits();
    }

    private long varintBits() {
        long value = 0;
        for (int i = 0; i < 10; i++) {
            int b = byteValue();
            require(i < 9 || b <= 1, "overflowing ULEB128");
            value |= (long) (b & 0x7F) << (i * 7);
            if ((b & 0x80) == 0) {
                require(i == 0 || b != 0, "non-canonical ULEB128");
                return value;
            }
        }
        throw new V10FormatException("unterminated ULEB128");
    }

    public void zero(int n) {
        for (int i = 0; i < n; i++) {
            require(byteValue() == 0, "nonzero reserved field");
        }
    }

    public void magic(String expected) {
        need(expected.length());
        for (int i = 0; i < expected.length(); i++) {
            require((data[at + i] & 0xFF) == (expected.charAt(i) & 0xFF), "bad record magic");
        }
        at += expected.length();
    }

    /**
     * Zero-terminated UTF-8 string; overlong forms, surrogates and values above
     * U+10FFFF are rejected.
     */
    public String str() {
        int begin = at;
        while (at < end && data[at] != 0) {
            at++;
        }
        require(at < end && at - begin <= MAX_STRING_LENGTH, "invalid string length or terminator");
        int length = at - begin;
        at++;
        validateUtf8(data, begin, length);
        try {
            return new String(data, begin, length, "UTF-8");
        } catch (UnsupportedEncodingException e) {
            throw new V10FormatException("UTF-8 is unavailable", e);
        }
    }

    private static void validateUtf8(byte[] bytes, int offset, int length) {
        int i = offset;
        int limit = offset + length;
        while (i < limit) {
            int c = bytes[i++] & 0xFF;
            if (c < 0x80) continue;
            int n;
            if (c >= 0xC2 && c <= 0xDF) n = 1;
            else if (c >= 0xE0 && c <= 0xEF) n = 2;
            else if (c >= 0xF0 && c <= 0xF4) n = 3;
            else n = 0;
            require(n > 0 && i + n <= limit, "invalid UTF-8");
            int v = c & ((1 << (6 - n)) - 1);
            for (int j = 0; j < n; j++) {
                int d = bytes[i++] & 0xFF;
                require((d & 0xC0) == 0x80, "invalid UTF-8");
                v = (v << 6) | (d & 0x3F);
            }
            int minimum = n == 1 ? 0x80 : n == 2 ? 0x800 : 0x10000;
            require(v >= minimum && v <= 0x10FFFF && !(v >= 0xD800 && v <= 0xDFFF), "invalid UTF-8");
        }
    }

    public V10Cursor take(long n) {
        need(n);
        V10Cursor sub = new V10Cursor(data, at, (int) n);
        at += (int) n;
        return sub;
    }

    /**
     * Reads a raw byte without advancing, relative to the start of this window.
     */
    public int peekAt(long index) {
        require(index >= 0 && index < size(), "bitmap index out of bounds");
        return data[start + (int) index] & 0xFF;
    }

    public void skipToEnd() {
        at = end;
    }

    public void done() {
        require(at == end, "trailing bytes in record");
    }
}
