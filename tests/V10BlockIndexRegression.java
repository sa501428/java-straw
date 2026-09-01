import io.airlift.compress.zstd.ZstdCompressor;
import javastraw.reader.Dataset;
import javastraw.reader.DatasetReaderV10;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.block.Block;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.mzd.Matrix;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.mzd.V10MatrixZoomData;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationHandler;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

/**
 * End-to-end regression for the final V10 layout, in which every logical block
 * is independently compressed and located through an exact H10I version-2 block
 * index. The withdrawn page-based draft is not readable and must be rejected.
 * <p>
 * The fixture is written here byte by byte from the normative tables in
 * hic-format/HiCFormatV10.md, so it exercises the reader rather than agreeing
 * with it.
 */
public final class V10BlockIndexRegression {

    private static final Charset UTF8 = Charset.forName("UTF-8");

    /**
     * Cells of the materialized 10 bp chrB x chrB matrix: {binColumn, binRow, count}.
     * With blockBinCount 4 and 20 block columns these fall into seven distinct
     * logical blocks, including three distance bands of the rotated cis grid.
     */
    private static final long[][] CELLS = {
            {0, 0, 1}, {1, 2, 1}, {3, 3, 5},        // block 0
            {5, 5, 2}, {6, 7, 3},                   // block 1
            {10, 11, 4},                            // block 2
            {20, 21, 6},                            // block 5
            {0, 10, 7},                             // block 21
            {2, 20, 8},                             // block 42
            {30, 79, 9}                             // block 73
    };

    private static final int BLOCK_BIN_COUNT = 4;
    private static final int BLOCK_COLUMN_COUNT = 20;   // ceil(80 bins / 4)
    private static final int CHR_B = 1;

    /**
     * @param args optional directory to write fixtures into; the system
     *             temporary directory is used when absent
     */
    public static void main(String[] args) throws Exception {
        Path dir = args.length > 0
                ? Files.createTempDirectory(new File(args[0]).toPath(), "java-straw-v10")
                : Files.createTempDirectory("java-straw-v10");
        try {
            readsIndependentBlocks(dir, false);
            // Physical contiguity is an optimization only; padded and reversed
            // layouts must decode identically through the exact locators.
            readsIndependentBlocks(dir, true);
            readsBlocksStoredOutOfOrder(dir);
            rejectsCorruptFiles(dir);
            System.out.println("V10 block index: exact locators, derived resolution, "
                    + "coalesced and split reads, draft rejection passed");
        } finally {
            deleteRecursively(dir.toFile());
        }
    }

    // ------------------------------------------------------------- assertions

    private static void readsIndependentBlocks(Path dir, boolean padBetweenBlocks) throws IOException {
        Fixture fixture = new Fixture(padBetweenBlocks, false);
        Path path = dir.resolve(padBetweenBlocks ? "padded.hic" : "packed.hic");
        Files.write(path, fixture.bytes);

        DatasetReaderV10 reader = new DatasetReaderV10(path.toString(), true);
        Dataset dataset = reader.read();
        assertEquals("chromosome count", 3, dataset.getChromosomeHandler().size());

        Chromosome chrB = dataset.getChromosomeHandler().getChromosomeFromIndex(CHR_B);
        Matrix matrix = dataset.getMatrix(chrB, chrB);
        assertTrue("matrix present", matrix != null);

        V10MatrixZoomData raw = (V10MatrixZoomData) matrix.getZoomData(new HiCZoom(HiCZoom.HiCUnit.BP, 10));
        assertTrue("10 bp is materialized", !raw.isDerivedResolution());
        assertEquals("occupied cells", (long) CELLS.length, raw.getOccupiedCellCount());

        // Every stored block, walked in index order.
        assertEquals("full iteration", expectedRaw(), collect(raw.getDirectIterator()));

        // One block only: the region touches block 73 alone.
        assertEquals("single-block region", cellsOf(new long[][]{{30, 79, 9}}),
                collect(raw, 28, 76, 31, 79));

        // A near-diagonal window spanning several adjacent blocks. Block reads
        // are permissive supersets, so the result is filtered to the region.
        Map<String, Long> window = collect(raw, 0, 0, 11, 11);
        for (long[] cell : CELLS) {
            boolean inside = cell[0] <= 11 && cell[1] <= 11;
            assertEquals("window cell " + cell[0] + "," + cell[1],
                    inside ? Long.valueOf(cell[2]) : null, window.get(key(cell[0], cell[1])));
        }

        // Off-diagonal bands must be reachable without walking every band.
        assertEquals("far band", cellsOf(new long[][]{{2, 20, 8}}), collect(raw, 2, 20, 2, 20));

        // The block index carries the stored size, so no read is needed for it.
        assertEquals("stored size of block 73", Integer.valueOf(fixture.storedLength(73)),
                raw.getBlockSize(73));
        assertEquals("absent block has no size", null, raw.getBlockSize(4));

        // The derived resolution stores nothing and is summed from the source.
        V10MatrixZoomData derived = (V10MatrixZoomData) matrix.getZoomData(new HiCZoom(HiCZoom.HiCUnit.BP, 20));
        assertTrue("20 bp is derived", derived.isDerivedResolution());
        assertEquals("derived iteration", expectedDerived(), collect(derived.getDirectIterator()));
        assertEquals("derived region", cellsOf(new long[][]{{15, 39, 9}}), collect(derived, 15, 39, 15, 39));
    }

    /**
     * Writers should lay blocks out in increasing block-number order, but the
     * index locators are authoritative. A file whose blocks are stored in
     * reverse physical order is still valid and must read identically.
     */
    private static void readsBlocksStoredOutOfOrder(Path dir) throws IOException {
        Path path = dir.resolve("reversed.hic");
        Files.write(path, new Fixture(false, true).bytes);
        DatasetReaderV10 reader = new DatasetReaderV10(path.toString(), true);
        Dataset dataset = reader.read();
        Chromosome chrB = dataset.getChromosomeHandler().getChromosomeFromIndex(CHR_B);
        MatrixZoomData raw = dataset.getMatrix(chrB, chrB).getZoomData(new HiCZoom(HiCZoom.HiCUnit.BP, 10));
        assertEquals("reversed physical order", expectedRaw(), collect(raw.getDirectIterator()));
    }

    private static void rejectsCorruptFiles(Path dir) throws IOException {
        // The withdrawn draft grouped blocks into pages behind a version-1 index.
        rejects(dir, "draft-index", "withdrawn page-based V10 draft", new Patch() {
            public void apply(Fixture f) {
                f.patchWord(f.indexPosition + 4, 1);
            }
        });
        rejects(dir, "reserved1", "nonzero reserved field", new Patch() {
            public void apply(Fixture f) {
                f.patchWord(f.descriptorPosition + 72, 1);
            }
        });
        rejects(dir, "index-count", "block index count mismatch", new Patch() {
            public void apply(Fixture f) {
                f.patchWord(f.indexPosition + 16, CELLS.length);
            }
        });
        rejects(dir, "unordered-index", "unordered block index", new Patch() {
            public void apply(Fixture f) {
                // Second entry's block number drops below the first.
                f.patchWord(f.indexPosition + 24 + 16, 0);
            }
        });
        rejects(dir, "overlapping-blocks", "overlapping stored block intervals", new Patch() {
            public void apply(Fixture f) {
                // Stretch the first record over the second without moving it.
                f.patchWord(f.indexPosition + 24 + 4, f.storedLength(0) + 8);
            }
        });
        rejects(dir, "block-number-mismatch", "block record/index number mismatch", new Patch() {
            public void apply(Fixture f) {
                f.patchWord(f.blockPosition(0) + 12, 99);
            }
        });
        rejects(dir, "block-codec", "unknown stored block codec", new Patch() {
            public void apply(Fixture f) {
                f.bytes[f.blockPosition(0) + 4] = 2;
            }
        });
    }

    private static void rejects(Path dir, String name, String reason, Patch patch) throws IOException {
        Fixture fixture = new Fixture(false, false);
        patch.apply(fixture);
        Path path = dir.resolve(name + ".hic");
        Files.write(path, fixture.bytes);
        try {
            // Read through the reader rather than Dataset, which logs and
            // swallows a malformed matrix record instead of propagating it.
            DatasetReaderV10 reader = new DatasetReaderV10(path.toString(), false);
            reader.read();
            Matrix matrix = reader.readMatrix(CHR_B + "_" + CHR_B, -1);
            if (matrix != null) {
                MatrixZoomData zd = matrix.getZoomData(new HiCZoom(HiCZoom.HiCUnit.BP, 10));
                collect(zd.getDirectIterator());
            }
        } catch (RuntimeException | IOException expected) {
            assertTrue(name + " was rejected for the wrong reason: " + expected.getMessage(),
                    expected.getMessage() != null && expected.getMessage().contains(reason));
            return;
        }
        throw new AssertionError("accepted a malformed V10 file: " + name);
    }

    private interface Patch {
        void apply(Fixture fixture);
    }

    // ------------------------------------------------------------- expected

    private static Map<String, Long> expectedRaw() {
        return cellsOf(CELLS);
    }

    private static Map<String, Long> expectedDerived() {
        Map<String, Long> sums = new TreeMap<>();
        for (long[] cell : CELLS) {
            String k = key(cell[0] / 2, cell[1] / 2);
            Long previous = sums.get(k);
            sums.put(k, (previous == null ? 0L : previous) + cell[2]);
        }
        return sums;
    }

    private static Map<String, Long> cellsOf(long[][] cells) {
        Map<String, Long> map = new TreeMap<>();
        for (long[] cell : cells) {
            map.put(key(cell[0], cell[1]), cell[2]);
        }
        return map;
    }

    private static String key(long binColumn, long binRow) {
        return binColumn + "," + binRow;
    }

    private static Map<String, Long> collect(Iterator<ContactRecord> it) {
        Map<String, Long> map = new TreeMap<>();
        while (it.hasNext()) {
            ContactRecord record = it.next();
            map.put(key(record.getBinX(), record.getBinY()), (long) record.getCounts());
        }
        return map;
    }

    /**
     * Records of the blocks overlapping a region, filtered to it exactly the way
     * a caller of the V6-V9 block API would.
     */
    private static Map<String, Long> collect(MatrixZoomData zd, long x1, long y1, long x2, long y2) {
        Map<String, Long> map = new TreeMap<>();
        List<Block> blocks = zd.getNormalizedBlocksOverlapping(x1, y1, x2, y2,
                NormalizationHandler.NONE, false);
        for (Block block : blocks) {
            for (ContactRecord record : block.getContactRecords()) {
                if (record.getBinX() < x1 || record.getBinX() > x2
                        || record.getBinY() < y1 || record.getBinY() > y2) {
                    continue;
                }
                map.put(key(record.getBinX(), record.getBinY()), (long) record.getCounts());
            }
        }
        return map;
    }

    // -------------------------------------------------------------- fixture

    /**
     * A complete V10 file: header, one cis matrix record with a materialized and
     * a derived resolution, independently compressed blocks, and their exact
     * block index.
     */
    private static final class Fixture {

        byte[] bytes;
        int descriptorPosition;
        int indexPosition;
        private final List<int[]> storedBlocks = new ArrayList<>();   // {blockNumber, position, length}

        Fixture(boolean padBetweenBlocks, boolean reversePhysicalOrder) {
            Buffer out = new Buffer();
            out.raw(new byte[88]);                       // backpatched prefix

            out.string("testGenome");
            out.u32(1);
            out.string("software");
            out.string("java-straw-test");

            out.u32(3);
            out.string("chrA");
            out.u64(500);
            out.string("chrB");
            out.u64(800);
            out.string("chrC");
            out.u64(700);

            out.u32(2);                                   // BP resolutions
            out.u32(10);
            out.u8(0);                                    // MATERIALIZED
            out.u8(1);                                    // SUM
            out.u16(0);
            out.u32(0xFFFFFFFFL);                         // no source
            out.u32(20);
            out.u8(1);                                    // DERIVED
            out.u8(1);
            out.u16(0);
            out.u32(0);                                   // source is index 0
            out.u32(0);                                   // no FRAG resolutions
            out.u32(0);                                   // no normalizations
            int headerLength = out.size();

            int matrixPosition = out.size();
            out.raw("H10M".getBytes(UTF8));
            out.u32(1);
            out.u32(CHR_B);
            out.u32(CHR_B);
            out.u32(2);
            out.u32(0);
            descriptorPosition = out.size();
            out.raw(new byte[76 * 2]);                    // backpatched descriptors

            List<byte[]> records = new ArrayList<>();
            List<Integer> numbers = new ArrayList<>();
            for (Map.Entry<Integer, List<long[]>> entry : groupByBlock().entrySet()) {
                numbers.add(entry.getKey());
                records.add(storedRecord(entry.getKey(), entry.getValue()));
            }
            int[] order = new int[records.size()];
            for (int i = 0; i < order.length; i++) {
                order[i] = reversePhysicalOrder ? order.length - 1 - i : i;
            }
            int[][] placed = new int[records.size()][];
            for (int i : order) {
                if (padBetweenBlocks && !storedBlocks.isEmpty()) {
                    out.raw(new byte[7]);                 // unindexed filler
                }
                placed[i] = new int[]{numbers.get(i), out.size(), records.get(i).length};
                out.raw(records.get(i));
            }
            for (int[] entry : placed) {
                storedBlocks.add(entry);                  // index order, by block number
            }

            indexPosition = out.size();
            long indexLength = 24L + 16L * storedBlocks.size();
            out.raw("H10I".getBytes(UTF8));
            out.u32(2);
            out.u64(indexLength);
            out.u32(storedBlocks.size());
            out.u32(0);
            for (int[] entry : storedBlocks) {
                out.u32(entry[0]);
                out.u32(entry[2]);
                out.u64(entry[1]);
            }

            int footerPosition = out.size();
            out.raw("H10F".getBytes(UTF8));
            out.u32(1);
            out.u64(48);
            out.u32(1);
            out.u32(0);
            out.u32(CHR_B);
            out.u32(CHR_B);
            out.u64(matrixPosition);
            out.u64(24 + 76 * 2);

            bytes = out.toArray();

            long sum = 0;
            for (long[] cell : CELLS) sum += cell[2];
            writeDescriptor(descriptorPosition, 0, 10, false, sum, CELLS.length,
                    BLOCK_COLUMN_COUNT, indexPosition, indexLength, storedBlocks.size());
            writeDescriptor(descriptorPosition + 76, 1, 20, true, sum, expectedDerived().size(),
                    10, 0, 0, 0);

            ByteBuffer prefix = ByteBuffer.wrap(bytes, 0, 88).order(ByteOrder.LITTLE_ENDIAN);
            prefix.put("HIC\0".getBytes(UTF8));
            prefix.putInt(10);
            prefix.putLong(headerLength);
            prefix.putLong(footerPosition);
            prefix.putLong(bytes.length - footerPosition);
            // No normalization, expected or normalized-expected vectors.
        }

        private void writeDescriptor(int at, int resolutionIndex, int binSize, boolean derived,
                                     long sum, long occupied, int blockColumnCount,
                                     long indexPos, long indexLen, int blockCount) {
            ByteBuffer b = ByteBuffer.wrap(bytes, at, 76).order(ByteOrder.LITTLE_ENDIAN);
            b.put((byte) 0);                              // BP
            b.put((byte) (derived ? 1 : 0));
            b.put((byte) 1);                              // SUM
            b.put((byte) 0);                              // COUNT_UINT
            b.putInt(resolutionIndex);
            b.putInt(binSize);
            b.putInt(derived ? 0 : 0xFFFFFFFF);
            b.put((byte) 1);                              // ROTATED_CIS
            b.put(new byte[3]);
            b.putLong(sum);
            b.putLong(occupied);
            b.putInt(0x7fc00000);                         // stdDev NaN
            b.putInt(0x7fc00000);                         // percent95 NaN
            b.putInt(BLOCK_BIN_COUNT);
            b.putInt(blockColumnCount);
            b.putLong(indexPos);
            b.putLong(indexLen);
            b.putInt(blockCount);
            b.putInt(0);                                  // reserved1
        }

        int storedLength(int blockNumber) {
            return entry(blockNumber)[2];
        }

        int blockPosition(int blockNumber) {
            return entry(blockNumber)[1];
        }

        private int[] entry(int blockNumber) {
            for (int[] entry : storedBlocks) {
                if (entry[0] == blockNumber) return entry;
            }
            throw new AssertionError("no such block " + blockNumber);
        }

        void patchWord(int at, long value) {
            ByteBuffer.wrap(bytes, at, 4).order(ByteOrder.LITTLE_ENDIAN).putInt((int) value);
        }

        /**
         * One H10B record wrapping a single independently compressed block.
         */
        private static byte[] storedRecord(int blockNumber, List<long[]> cells) {
            byte[] logical = logicalBlock(cells);
            byte[] frame = compress(logical);
            Buffer out = new Buffer();
            out.raw("H10B".getBytes(UTF8));
            out.u8(1);                                    // ZSTD
            out.u8(1);                                    // record version
            out.u16(0);
            out.u32(logical.length);
            out.u32(blockNumber);
            out.raw(frame);
            return out.toArray();
        }

        /**
         * A sparse-delta block with direct ULEB128 counts (Sections I.1-I.5).
         */
        private static byte[] logicalBlock(List<long[]> cells) {
            long minColumn = Long.MAX_VALUE;
            long minRow = Long.MAX_VALUE;
            long maxColumn = Long.MIN_VALUE;
            long maxRow = Long.MIN_VALUE;
            for (long[] cell : cells) {
                minColumn = Math.min(minColumn, cell[0]);
                maxColumn = Math.max(maxColumn, cell[0]);
                minRow = Math.min(minRow, cell[1]);
                maxRow = Math.max(maxRow, cell[1]);
            }
            long width = maxColumn - minColumn + 1;
            long height = maxRow - minRow + 1;

            Map<Long, Long> byPosition = new TreeMap<>();
            for (long[] cell : cells) {
                byPosition.put((cell[1] - minRow) * width + (cell[0] - minColumn), cell[2]);
            }

            Buffer positions = new Buffer();
            Buffer values = new Buffer();
            long previous = 0;
            for (Map.Entry<Long, Long> cell : byPosition.entrySet()) {
                positions.varint(cell.getKey() - previous);
                previous = cell.getKey();
                values.varint(cell.getValue());
            }
            byte[] positionStream = positions.toArray();
            byte[] valueStream = values.toArray();

            Buffer out = new Buffer();
            out.u8(1);                                    // block version
            out.u8(0);                                    // SPARSE_DELTA
            out.u8(2);                                    // DIRECT
            out.u8(0);                                    // COUNT_UINT
            out.u8(0);                                    // no EXPLICIT_PRESENCE
            out.raw(new byte[3]);
            out.u32(minColumn);
            out.u32(minRow);
            out.u32(width);
            out.u32(height);
            out.u64(byPosition.size());
            out.u32(positionStream.length);
            out.u32(valueStream.length);
            out.raw(positionStream);
            out.raw(valueStream);
            return out.toArray();
        }

        /**
         * Cells grouped by rotated-cis block number, computed here directly from
         * Section F.3 rather than through the reader's own grid code.
         */
        private static Map<Integer, List<long[]>> groupByBlock() {
            Map<Integer, List<long[]>> byBlock = new TreeMap<>();
            for (long[] cell : CELLS) {
                int number = rotatedCisBlockNumber(cell[0], cell[1]);
                List<long[]> list = byBlock.get(number);
                if (list == null) {
                    list = new ArrayList<>();
                    byBlock.put(number, list);
                }
                list.add(cell);
            }
            return byBlock;
        }

        private static int rotatedCisBlockNumber(long binColumn, long binRow) {
            long d = binRow - binColumn;
            int depth = 0;
            while (d * d >= 2L * BLOCK_BIN_COUNT * BLOCK_BIN_COUNT
                    * ((1L << (depth + 1)) - 1) * ((1L << (depth + 1)) - 1)) {
                depth++;
            }
            long alongDiagonal = (binColumn + binRow) / (2L * BLOCK_BIN_COUNT);
            return (int) (depth * BLOCK_COLUMN_COUNT + alongDiagonal);
        }

        private static byte[] compress(byte[] input) {
            ZstdCompressor compressor = new ZstdCompressor();
            byte[] output = new byte[compressor.maxCompressedLength(input.length)];
            int n = compressor.compress(input, 0, input.length, output, 0, output.length);
            return Arrays.copyOf(output, n);
        }
    }

    /**
     * Little-endian byte assembler.
     */
    private static final class Buffer {
        private final ByteArrayOutputStream out = new ByteArrayOutputStream();

        int size() {
            return out.size();
        }

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

        void varint(long value) {
            while (value >= 128) {
                u8((int) (value & 127) | 128);
                value >>>= 7;
            }
            u8((int) value);
        }

        void string(String value) {
            raw(value.getBytes(UTF8));
            u8(0);
        }

        byte[] toArray() {
            return out.toByteArray();
        }
    }

    // ------------------------------------------------------------- plumbing

    private static void assertEquals(String what, Object expected, Object actual) {
        if (expected == null ? actual != null : !expected.equals(actual)) {
            throw new AssertionError(what + ": expected " + expected + ", got " + actual);
        }
    }

    private static void assertTrue(String what, boolean ok) {
        if (!ok) throw new AssertionError(what);
    }

    private static void deleteRecursively(File file) {
        File[] children = file.listFiles();
        if (children != null) {
            for (File child : children) deleteRecursively(child);
        }
        // Best effort: the fixtures live under a temporary directory.
        file.delete();
    }

    private V10BlockIndexRegression() {
    }
}
