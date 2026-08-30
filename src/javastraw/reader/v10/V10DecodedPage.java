package javastraw.reader.v10;

import static javastraw.reader.v10.V10.require;

/**
 * A decompressed matrix page: its block directory plus the payload bytes
 * (Section H). Once a page has been decompressed, nearby block queries need no
 * further I/O or decompression, so instances are worth caching.
 */
public class V10DecodedPage {

    private final byte[] payload;
    private final int[] blockNumbers;
    private final int[] offsets;
    private final int[] lengths;

    private V10DecodedPage(byte[] payload, int[] blockNumbers, int[] offsets, int[] lengths) {
        this.payload = payload;
        this.blockNumbers = blockNumbers;
        this.offsets = offsets;
        this.lengths = lengths;
    }

    public int blockCount() {
        return blockNumbers.length;
    }

    public int blockNumber(int i) {
        return blockNumbers[i];
    }

    /**
     * Index of {@code blockNumber} within this page, or -1 if the page does not
     * store it. Block numbers inside a page range may legitimately be absent.
     */
    public int indexOf(int blockNumber) {
        int lo = 0;
        int hi = blockNumbers.length - 1;
        while (lo <= hi) {
            int mid = (lo + hi) >>> 1;
            if (blockNumber < blockNumbers[mid]) {
                hi = mid - 1;
            } else if (blockNumber > blockNumbers[mid]) {
                lo = mid + 1;
            } else {
                return mid;
            }
        }
        return -1;
    }

    public V10Cursor cursorFor(int index) {
        return new V10Cursor(payload, offsets[index], lengths[index]);
    }

    /**
     * Approximate retained size, used to bound the page cache.
     */
    public int byteSize() {
        return payload.length + 12 * blockNumbers.length;
    }

    /**
     * Validates the stored page header, decompresses the Zstandard frame and
     * reconstructs the block directory.
     */
    public static V10DecodedPage parse(byte[] stored, V10Page page) {
        V10Cursor c = new V10Cursor(stored);
        c.magic("H10P");
        require(c.byteValue() == V10.CODEC_ZSTD, "unknown page codec");
        require(c.byteValue() == 1, "unknown page version");
        c.zero(2);
        require(c.word() == page.uncompressedBytes, "page size mismatch");
        long blockCount = c.word();
        require(blockCount > 0 && blockCount <= page.uncompressedBytes / 42, "invalid page block count");

        byte[] payload = V10Zstd.decompress(c, page.uncompressedBytes);
        c.done();

        V10Cursor body = new V10Cursor(payload);
        long directoryLength = body.word();
        V10Cursor directory = body.take(directoryLength);

        int n = (int) blockCount;
        int[] numbers = new int[n];
        int[] lengths = new int[n];
        int[] offsets = new int[n];
        long previous = 0;
        long totalLength = 0;
        for (int i = 0; i < n; i++) {
            long delta = directory.varint();
            require(i == 0 || delta > 0, "duplicate block number");
            long number = i == 0 ? delta : V10.add(previous, delta);
            V10.unsigned32(number);
            long length = directory.varint();
            require(length >= 40, "invalid block length");
            numbers[i] = V10.toInt(number, "block number");
            lengths[i] = V10.toInt(length, "block length");
            totalLength = V10.add(totalLength, length);
            previous = number;
        }
        directory.done();
        require(numbers[0] == page.firstBlock && numbers[n - 1] == page.lastBlock, "page range mismatch");
        require(totalLength == body.left(), "block payload length mismatch");

        int start = body.rawPosition();
        for (int i = 0; i < n; i++) {
            offsets[i] = start;
            start += lengths[i];
        }
        return new V10DecodedPage(payload, numbers, offsets, lengths);
    }
}
