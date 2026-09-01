package javastraw.reader.v10;

import java.util.ArrayList;
import java.util.List;

/**
 * Logical block geometry (Section F): the rectangular trans grid and the
 * rotated cis grid, plus enumeration of the blocks a query rectangle can touch.
 */
public final class V10Grid {

    /** Inclusive, ordered logical block-number interval. */
    public static final class BlockRange {
        public final int first;
        public final int last;

        BlockRange(int first, int last) {
            this.first = first;
            this.last = last;
        }
    }

    /**
     * Largest long whose square still fits a signed 64-bit integer.
     */
    private static final long SQRT_LONG_MAX = 3037000499L;

    private V10Grid() {
    }

    /**
     * Anti-diagonal band index: the largest integer a >= 0 such that
     * d^2 >= 2 * B^2 * (2^a - 1)^2. This is the normative integer form of
     * floor(log2(1 + d / (sqrt(2) * B))) and avoids platform-dependent
     * floating-point behavior near a band boundary (Section F.3).
     */
    public static int alongAntiDiagonal(long binColumn, long binRow, long blockBinCount) {
        return depthForDistance(Math.abs(binRow - binColumn), blockBinCount);
    }

    private static int depthForDistance(long d, long blockBinCount) {
        long dd = V10.multiply(d, d);
        long b = V10.multiply(2, V10.multiply(blockBinCount, blockBinCount));
        long quotient = dd / b;
        int depth = 0;
        while (depth < 32) {
            long t = (1L << (depth + 1)) - 1;
            if (t > SQRT_LONG_MAX || t * t > quotient) break;
            depth++;
        }
        return depth;
    }

    /**
     * Exact candidate block-number ranges for an inclusive bin rectangle.
     * One range is returned per trans row or cis distance band. This is the
     * representation consumed by the block index so sparse queries never walk
     * the empty logical block numbers between stored entries.
     */
    public static List<BlockRange> blockRangesForRegion(long binX1, long binY1,
                                                         long binX2, long binY2,
                                                         V10Zoom zoom, long nBins1, long nBins2) {
        List<BlockRange> ranges = new ArrayList<>();
        long b = zoom.blockBinCount;
        long x1 = Math.max(0, Math.min(binX1, binX2));
        long x2 = Math.min(nBins1 - 1, Math.max(binX1, binX2));
        long y1 = Math.max(0, Math.min(binY1, binY2));
        long y2 = Math.min(nBins2 - 1, Math.max(binY1, binY2));
        if (x1 > x2 || y1 > y2) return ranges;

        if (zoom.isRotatedCis()) {
            long firstPad = V10.add(x1, y1) / (2 * b);
            long lastPad = V10.add(x2, y2) / (2 * b);
            long nearest = x2 < y1 ? y1 - x2 : y2 < x1 ? x1 - y2 : 0;
            long farthest = Math.max(Math.abs(y2 - x1), Math.abs(x2 - y1));
            int firstDepth = depthForDistance(nearest, b);
            int lastDepth = depthForDistance(farthest, b);
            for (int depth = firstDepth; depth <= lastDepth; depth++) {
                addRange(ranges, depth * (long) zoom.blockColumnCount + firstPad,
                        depth * (long) zoom.blockColumnCount + lastPad);
            }
        } else {
            long firstColumn = x1 / b;
            long lastColumn = x2 / b;
            for (long row = y1 / b; row <= y2 / b; row++) {
                addRange(ranges, row * zoom.blockColumnCount + firstColumn,
                        row * zoom.blockColumnCount + lastColumn);
            }
        }
        return ranges;
    }

    /**
     * Logical block number of a canonical cell under the descriptor's grid.
     */
    public static int blockNumber(long binColumn, long binRow, V10Zoom zoom) {
        long b = zoom.blockBinCount;
        if (!zoom.isRotatedCis()) {
            long number = V10.add(V10.multiply(binRow / b, zoom.blockColumnCount), binColumn / b);
            return V10.toInt(V10.unsigned32(number), "block number");
        }
        long depth = alongAntiDiagonal(binColumn, binRow, b);
        long alongDiagonal = (binColumn + binRow) / (2 * b);
        long number = V10.add(V10.multiply(depth, zoom.blockColumnCount), alongDiagonal);
        return V10.toInt(V10.unsigned32(number), "block number");
    }

    /**
     * Every logical block number a query rectangle could intersect. Bin bounds
     * are inclusive, matching the rest of java-straw. The result is a permissive
     * superset: numbers with no stored block are simply absent from the file.
     */
    public static List<Integer> blockNumbersForRegion(long binX1, long binY1, long binX2, long binY2,
                                                      V10Zoom zoom, long nBins1, long nBins2) {
        List<Integer> blocks = new ArrayList<>();
        for (BlockRange range : blockRangesForRegion(binX1, binY1, binX2, binY2,
                zoom, nBins1, nBins2)) {
            for (int number = range.first; ; number++) {
                blocks.add(number);
                if (number == range.last) break;
            }
        }
        return blocks;
    }

    private static void addRange(List<BlockRange> ranges, long first, long last) {
        ranges.add(new BlockRange(V10.toInt(V10.unsigned32(first), "block number"),
                V10.toInt(V10.unsigned32(last), "block number")));
    }
}
