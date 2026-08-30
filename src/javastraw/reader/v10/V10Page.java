package javastraw.reader.v10;

import java.util.ArrayList;
import java.util.List;

import static javastraw.reader.v10.V10.require;

/**
 * One stored matrix page and the parser for the checkpointed resolution page
 * index (Sections G and H).
 */
public class V10Page {

    /**
     * Inclusive first logical block number in this page.
     */
    public final int firstBlock;
    /**
     * Inclusive last logical block number in this page.
     */
    public final int lastBlock;
    /**
     * Decompressed payload length.
     */
    public final int uncompressedBytes;
    public final long position;
    /**
     * Stored length, including the 16-byte page header.
     */
    public final long storedByteLength;

    V10Page(int firstBlock, int lastBlock, int uncompressedBytes, long position, long storedByteLength) {
        this.firstBlock = firstBlock;
        this.lastBlock = lastBlock;
        this.uncompressedBytes = uncompressedBytes;
        this.position = position;
        this.storedByteLength = storedByteLength;
    }

    public boolean contains(int blockNumber) {
        return blockNumber >= firstBlock && blockNumber <= lastBlock;
    }

    /**
     * Parses the uncompressed page index of a materialized resolution. Pages are
     * returned in ordinal order with non-overlapping, increasing block ranges.
     */
    public static List<V10Page> parseIndex(byte[] bytes, V10Zoom zoom, V10Source source) {
        V10Cursor c = new V10Cursor(bytes);
        c.magic("H10I");
        require(c.word() == 1, "page index version mismatch");
        require(c.word() == zoom.pageCount, "page index count mismatch");
        long interval = c.word();
        long groups = c.word();
        c.zero(4);
        long blobLength = c.wide();
        require(interval > 0 && groups == ((long) zoom.pageCount + interval - 1) / interval
                        && groups * 32 <= c.left(),
                "invalid checkpoints");

        long[] firstPageOrdinal = new long[(int) groups];
        long[] pagesInGroup = new long[(int) groups];
        long[] firstBlockNumber = new long[(int) groups];
        long[] firstPagePosition = new long[(int) groups];
        long[] descriptorOffset = new long[(int) groups];
        for (int i = 0; i < (int) groups; i++) {
            firstPageOrdinal[i] = c.word();
            pagesInGroup[i] = c.word();
            firstBlockNumber[i] = c.word();
            c.zero(4);
            firstPagePosition[i] = c.wide();
            descriptorOffset[i] = c.wide();
        }
        require(blobLength == c.left(), "page descriptor length mismatch");
        V10Cursor blob = c.take(blobLength);

        List<V10Page> out = new ArrayList<>(zoom.pageCount);
        for (int i = 0; i < (int) groups; i++) {
            // Descriptors are consumed strictly sequentially, so each group must
            // begin exactly where the previous group ended.
            require(firstPageOrdinal[i] == out.size() && pagesInGroup[i] > 0 && pagesInGroup[i] <= interval
                            && descriptorOffset[i] == blob.consumed(),
                    "invalid checkpoint coverage");
            require(i != 0 || descriptorOffset[i] == 0, "first descriptor offset must be zero");
            long position = firstPagePosition[i];
            long first = firstBlockNumber[i];
            if (!out.isEmpty()) {
                V10Page previous = out.get(out.size() - 1);
                require(position == V10.add(previous.position, previous.storedByteLength)
                                && first > previous.lastBlock,
                        "pages are not contiguous/ordered");
            }
            for (long j = 0; j < pagesInGroup[i]; j++) {
                if (j > 0) {
                    V10Page previous = out.get(out.size() - 1);
                    first = V10.add(V10.add(previous.lastBlock, 1), blob.varint());
                }
                long last = V10.add(first, blob.varint());
                long storedLength = blob.varint();
                long raw = blob.varint();
                require(storedLength > 16 && raw >= 4 && raw <= V10.ALLOCATION_LIMIT, "invalid page length");
                V10.unsigned32(last);
                source.interval(position, storedLength);
                out.add(new V10Page(V10.toInt(first, "block number"), V10.toInt(last, "block number"),
                        (int) raw, position, storedLength));
                position = V10.add(position, storedLength);
            }
        }
        blob.done();
        require(out.size() == zoom.pageCount, "page count mismatch");
        return out;
    }

    /**
     * Index of the page whose block range contains {@code blockNumber}, or -1.
     * Page ranges are ordered and non-overlapping, so a binary search suffices.
     */
    public static int findPage(List<V10Page> pages, int blockNumber) {
        int lo = 0;
        int hi = pages.size() - 1;
        while (lo <= hi) {
            int mid = (lo + hi) >>> 1;
            V10Page page = pages.get(mid);
            if (blockNumber < page.firstBlock) {
                hi = mid - 1;
            } else if (blockNumber > page.lastBlock) {
                lo = mid + 1;
            } else {
                return mid;
            }
        }
        return -1;
    }
}
