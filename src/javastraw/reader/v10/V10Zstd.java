package javastraw.reader.v10;

import io.airlift.compress.zstd.ZstdDecompressor;

import static javastraw.reader.v10.V10.require;

/**
 * Zstandard support for V10 stored blocks and vector chunks (Sections H and J.2).
 * <p>
 * The stored payload must be exactly one ordinary Zstandard data frame as
 * defined by RFC 8878: not a skippable frame, not a concatenation of frames,
 * and never using a preset dictionary. Those constraints are checked here by
 * walking the frame header and block headers before any decompression happens,
 * which also yields the exact compressed frame length.
 */
public final class V10Zstd {

    private static final int MAGIC = 0xFD2FB528;
    private static final ThreadLocal<ZstdDecompressor> DECOMPRESSOR =
            new ThreadLocal<ZstdDecompressor>() {
                @Override
                protected ZstdDecompressor initialValue() {
                    return new ZstdDecompressor();
                }
            };

    private V10Zstd() {
    }

    /**
     * Decompresses the remainder of {@code cursor}, which must be exactly one
     * Zstandard frame decompressing to {@code expectedBytes} bytes.
     */
    public static byte[] decompress(V10Cursor cursor, long expectedBytes) {
        require(expectedBytes > 0 && expectedBytes <= V10.ALLOCATION_LIMIT,
                "decompressed record exceeds allocation limit");
        byte[] input = cursor.rawArray();
        int offset = cursor.rawPosition();
        int length = cursor.left();
        long frameLength = frameCompressedSize(input, offset, length);
        require(frameLength == length, "invalid or concatenated Zstandard frame");

        byte[] output = new byte[(int) expectedBytes];
        int produced;
        try {
            produced = DECOMPRESSOR.get().decompress(input, offset, length, output, 0, output.length);
        } catch (RuntimeException e) {
            throw new V10FormatException("Zstandard decompression failure: " + e.getMessage(), e);
        }
        require(produced == expectedBytes, "Zstandard decompressed length mismatch");
        cursor.skipToEnd();
        return output;
    }

    /**
     * Walks the frame header and block headers of a single Zstandard frame and
     * returns its exact stored length. Rejects dictionaries, skippable frames
     * and reserved encodings.
     */
    static long frameCompressedSize(byte[] data, int offset, int length) {
        require(length >= 6, "not a Zstandard data frame");
        long magic = readLE(data, offset, 4);
        require(magic == (MAGIC & 0xFFFFFFFFL), "not a Zstandard data frame");

        int at = offset + 4;
        int limit = offset + length;
        require(at < limit, "truncated Zstandard frame header");
        int descriptor = data[at++] & 0xFF;
        int contentSizeFlag = (descriptor >>> 6) & 3;
        boolean singleSegment = ((descriptor >>> 5) & 1) != 0;
        boolean hasChecksum = ((descriptor >>> 2) & 1) != 0;
        int dictionaryIdFlag = descriptor & 3;
        require(((descriptor >>> 3) & 3) == 0, "reserved Zstandard frame header bits are set");

        if (!singleSegment) {
            require(at < limit, "truncated Zstandard frame header");
            at++; // window descriptor
        }
        int dictionaryIdBytes = dictionaryIdFlag == 0 ? 0 : dictionaryIdFlag == 3 ? 4 : dictionaryIdFlag;
        require(at + dictionaryIdBytes <= limit, "truncated Zstandard frame header");
        if (dictionaryIdBytes > 0) {
            require(readLE(data, at, dictionaryIdBytes) == 0, "Zstandard dictionaries are forbidden");
            at += dictionaryIdBytes;
        }
        int contentSizeBytes = contentSizeFlag == 0 ? (singleSegment ? 1 : 0)
                : contentSizeFlag == 1 ? 2 : contentSizeFlag == 2 ? 4 : 8;
        require(at + contentSizeBytes <= limit, "truncated Zstandard frame header");
        at += contentSizeBytes;

        while (true) {
            require(at + 3 <= limit, "truncated Zstandard block header");
            int header = (int) readLE(data, at, 3);
            at += 3;
            boolean last = (header & 1) != 0;
            int blockType = (header >>> 1) & 3;
            int blockSize = header >>> 3;
            require(blockType != 3, "reserved Zstandard block type");
            int stored = blockType == 1 ? 1 : blockSize; // RLE blocks store a single byte
            require(at + stored <= limit && stored >= 0, "truncated Zstandard block");
            at += stored;
            if (last) break;
        }
        if (hasChecksum) {
            require(at + 4 <= limit, "truncated Zstandard frame checksum");
            at += 4;
        }
        return at - offset;
    }

    private static long readLE(byte[] data, int offset, int nBytes) {
        long value = 0;
        for (int i = 0; i < nBytes; i++) {
            value |= (long) (data[offset + i] & 0xFF) << (8 * i);
        }
        return value;
    }
}
