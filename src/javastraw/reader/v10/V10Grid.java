package javastraw.reader.v10;

import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

/**
 * Logical block geometry (Section F): the rectangular trans grid and the
 * rotated cis grid, plus enumeration of the blocks a query rectangle can touch.
 */
public final class V10Grid {

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
        long d = Math.abs(binRow - binColumn);
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
        long b = zoom.blockBinCount;
        Set<Integer> blocks = new LinkedHashSet<>();

        long x1 = Math.max(0, Math.min(binX1, binX2));
        long x2 = Math.min(nBins1 - 1, Math.max(binX1, binX2));
        long y1 = Math.max(0, Math.min(binY1, binY2));
        long y2 = Math.min(nBins2 - 1, Math.max(binY1, binY2));
        if (x1 > x2 || y1 > y2) return new ArrayList<>(blocks);

        if (zoom.isRotatedCis()) {
            long lowerPAD = (x1 + y1) / 2 / b;
            long higherPAD = (x2 + y2) / 2 / b + 1;
            int nearerDepth = alongAntiDiagonal(x1, y2, b);
            int furtherDepth = alongAntiDiagonal(x2, y1, b);
            int nearest = Math.min(nearerDepth, furtherDepth);
            // A rectangle straddling the diagonal reaches distance zero, whatever
            // its corner distances are.
            if (x1 <= y2 && y1 <= x2) {
                nearest = 0;
            }
            int furthest = Math.max(nearerDepth, furtherDepth) + 1;
            for (long depth = nearest; depth <= furthest; depth++) {
                for (long pad = lowerPAD; pad <= higherPAD; pad++) {
                    addIfValid(blocks, depth * zoom.blockColumnCount + pad);
                }
            }
        } else {
            long col1 = x1 / b;
            long col2 = x2 / b;
            long row1 = y1 / b;
            long row2 = y2 / b;
            for (long r = row1; r <= row2; r++) {
                for (long c = col1; c <= col2; c++) {
                    addIfValid(blocks, r * zoom.blockColumnCount + c);
                }
            }
        }
        return new ArrayList<>(blocks);
    }

    private static void addIfValid(Set<Integer> blocks, long number) {
        if (number >= 0 && number <= Integer.MAX_VALUE) {
            blocks.add((int) number);
        }
    }
}
