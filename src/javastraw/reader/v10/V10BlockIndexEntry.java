package javastraw.reader.v10;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static javastraw.reader.v10.V10.require;

/**
 * One locator of the exact resolution block index, plus the parser for the index
 * itself (Section G).
 * <p>
 * Every non-empty logical block is stored independently and has exactly one
 * entry here, so a block is addressed by an absolute file position rather than
 * reconstructed from cumulative sizes. Entries are sorted by strictly increasing
 * block number, which makes lookup a binary search and leaves absent numbers
 * meaning "empty logical block".
 */
public class V10BlockIndexEntry {

    public final int blockNumber;
    /**
     * Stored-block header plus one Zstandard frame; always greater than 16.
     */
    public final long storedByteLength;
    /**
     * Absolute file position of the H10B stored-block record.
     */
    public final long blockPosition;

    V10BlockIndexEntry(int blockNumber, long storedByteLength, long blockPosition) {
        this.blockNumber = blockNumber;
        this.storedByteLength = storedByteLength;
        this.blockPosition = blockPosition;
    }

    /**
     * Parses and validates the block index of a materialized resolution. The
     * index is exactly {@code 24 + 16 * blockCount} bytes.
     */
    public static List<V10BlockIndexEntry> parseIndex(byte[] bytes, V10Zoom zoom, V10Source source) {
        V10Cursor c = new V10Cursor(bytes);
        c.magic("H10I");
        long version = c.word();
        // Pre-final V10 drafts grouped blocks into compression pages behind a
        // version-1 index. That layout is not guessed at; it is rejected.
        require(version != 1, "this file uses the withdrawn page-based V10 draft "
                + "(block index version 1) and must be rewritten");
        require(version == 2, "unknown block index version " + version);
        require(c.wide() == zoom.blockIndex.length, "block index length mismatch");
        long count = c.word();
        c.zero(4);
        require(count == zoom.logicalBlockCount, "block index count mismatch");
        require(zoom.blockIndex.length == V10.add(24, V10.multiply(count, 16)),
                "invalid block index length");

        List<V10BlockIndexEntry> out = new ArrayList<>((int) count);
        // Writers should lay blocks out contiguously in block-number order. When
        // they do, non-overlap follows from a running comparison and the sorted
        // sweep below is not needed.
        boolean ascending = true;
        for (int i = 0; i < (int) count; i++) {
            long number = c.word();
            long storedLength = c.word();
            long position = c.wide();
            require(storedLength > 16, "invalid stored block length");
            source.interval(position, storedLength);
            if (i > 0) {
                V10BlockIndexEntry previous = out.get(i - 1);
                require(number > previous.blockNumber, "unordered block index");
                ascending &= position >= V10.add(previous.blockPosition, previous.storedByteLength);
            }
            out.add(new V10BlockIndexEntry(V10.toInt(number, "block number"), storedLength, position));
        }
        c.done();
        if (!ascending) {
            requireNonOverlapping(out);
        }
        return out;
    }

    /**
     * Indexed intervals must not overlap. Sorting the interval starts and ends
     * separately detects any overlap: a second record begins before an earlier
     * one ends exactly when {@code ends[i] > starts[i + 1]} somewhere.
     */
    private static void requireNonOverlapping(List<V10BlockIndexEntry> entries) {
        int n = entries.size();
        long[] starts = new long[n];
        long[] ends = new long[n];
        for (int i = 0; i < n; i++) {
            V10BlockIndexEntry entry = entries.get(i);
            starts[i] = entry.blockPosition;
            ends[i] = V10.add(entry.blockPosition, entry.storedByteLength);
        }
        Arrays.sort(starts);
        Arrays.sort(ends);
        for (int i = 0; i + 1 < n; i++) {
            require(ends[i] <= starts[i + 1], "overlapping stored block intervals");
        }
    }

    /**
     * Index of {@code blockNumber} within the block index, or -1 when the file
     * stores no such block. Entries are ordered by block number, so this is a
     * binary search.
     */
    public static int find(List<V10BlockIndexEntry> index, int blockNumber) {
        int lo = lowerBound(index, blockNumber);
        return lo < index.size() && index.get(lo).blockNumber == blockNumber ? lo : -1;
    }

    /** First entry whose block number is at least {@code blockNumber}. */
    public static int lowerBound(List<V10BlockIndexEntry> index, int blockNumber) {
        int lo = 0;
        int hi = index.size();
        while (lo < hi) {
            int mid = (lo + hi) >>> 1;
            int number = index.get(mid).blockNumber;
            if (number < blockNumber) {
                lo = mid + 1;
            } else {
                hi = mid;
            }
        }
        return lo;
    }
}
