package javastraw.reader.v10;

import static javastraw.reader.v10.V10.require;

/**
 * One decompressed logical block: the payload of an independently compressed
 * H10B stored-block record (Section H).
 * <p>
 * The Zstandard frame of a stored record holds exactly one logical block, so
 * decoding a block never requires reading or decompressing any other block.
 */
public class V10StoredBlock {

    private final int blockNumber;
    private final byte[] payload;

    private V10StoredBlock(int blockNumber, byte[] payload) {
        this.blockNumber = blockNumber;
        this.payload = payload;
    }

    public int getBlockNumber() {
        return blockNumber;
    }

    /**
     * A cursor over the whole logical block, positioned at its 40-byte header.
     */
    public V10Cursor cursor() {
        return new V10Cursor(payload);
    }

    /**
     * Approximate retained size, used to bound the decompressed-block cache.
     */
    public int byteSize() {
        return payload.length;
    }

    public static V10StoredBlock parse(byte[] stored, V10BlockIndexEntry entry) {
        return parse(stored, 0, stored.length, entry);
    }

    /**
     * Validates the stored-record header of the block at
     * {@code [offset, offset + length)} and decompresses its single frame.
     * The window may be a slice of one coalesced range read.
     */
    public static V10StoredBlock parse(byte[] stored, int offset, int length, V10BlockIndexEntry entry) {
        V10Cursor c = new V10Cursor(stored, offset, length);
        c.magic("H10B");
        require(c.byteValue() == V10.CODEC_ZSTD, "unknown stored block codec");
        require(c.byteValue() == 1, "unknown stored block version");
        c.zero(2);
        long uncompressedBytes = c.word();
        require(c.word() == entry.blockNumber, "block record/index number mismatch");
        require(uncompressedBytes >= 40, "invalid uncompressed block length");
        byte[] payload = V10Zstd.decompress(c, uncompressedBytes);
        c.done();
        return new V10StoredBlock(entry.blockNumber, payload);
    }
}
