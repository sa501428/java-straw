/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2011-2024 Broad Institute, Aiden Lab, Rice University, Baylor College of Medicine
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 *  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 *  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 *  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 *  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 *  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 *  THE SOFTWARE.
 */

package javastraw.reader.mzd;

import javastraw.reader.DatasetReaderV10;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.block.Block;
import javastraw.reader.block.BlockModifier;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.block.V10ContactRecord;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationHandler;
import javastraw.reader.type.NormalizationType;
import javastraw.reader.v10.V10;
import javastraw.reader.v10.V10BlockDecoder;
import javastraw.reader.v10.V10BlockHeader;
import javastraw.reader.v10.V10BlockIndexEntry;
import javastraw.reader.v10.V10Cursor;
import javastraw.reader.v10.V10FormatException;
import javastraw.reader.v10.V10Grid;
import javastraw.reader.v10.V10RecordHandler;
import javastraw.reader.v10.V10StoredBlock;
import javastraw.reader.v10.V10Zoom;

import java.io.IOException;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.TreeMap;
import java.util.function.LongToDoubleFunction;

/**
 * A V10 matrix at one logical resolution.
 * <p>
 * Two storage modes are handled transparently. A <b>materialized</b> resolution
 * locates each logical block through its exact block index and decodes that
 * block's own Zstandard record.
 * A <b>derived</b> resolution stores no blocks: its raw values are reconstructed
 * by exact summation from the declared finer materialized source, always before
 * normalization or expected-value division (Section K). Either way, callers see
 * the same {@link Block} and {@link ContactRecord} objects as for V6-V9 files.
 */
public class V10MatrixZoomData extends MatrixZoomData {

    /**
     * How many stored blocks a bulk scan fetches at once. Adjacent records are
     * coalesced into one range read, so a batch bounds both the number of reads
     * and how many decompressed blocks are held at a time.
     */
    private static final int STREAM_BATCH_SIZE = 16;

    private final DatasetReaderV10 v10Reader;
    private final V10Zoom v10Zoom;
    /**
     * Non-null exactly when this resolution is derived.
     */
    private final V10MatrixZoomData sourceZD;
    private final long nBins1;
    private final long nBins2;
    private final int derivationFactor;

    private List<V10BlockIndexEntry> blockIndex;
    // View tiles have stable bin rectangles. Reuse their exact raw aggregates
    // across redraws and normalization changes without retaining source records.
    private final Map<String, List<ContactRecord>> derivedRawCache = new LinkedHashMap<String, List<ContactRecord>>(16, .75f, true);
    private long derivedRawBytes;
    private static final long DERIVED_CACHE_BYTES = 64L * 1024 * 1024;
    private final Map<Integer, RecordBlock> rawBlockCache = new LinkedHashMap<Integer, RecordBlock>(16, .75f, true);
    private long rawBlockBytes;
    private static final long RAW_BLOCK_CACHE_BYTES = 128L * 1024 * 1024;

    public V10MatrixZoomData(Chromosome chr1, Chromosome chr2, V10Zoom zoom, DatasetReaderV10 reader,
                             V10MatrixZoomData sourceZD, boolean useCache, long nBins1, long nBins2) {
        super(chr1, chr2, new HiCZoom(zoom.unit == V10.UNIT_FRAG ? HiCZoom.HiCUnit.FRAG : HiCZoom.HiCUnit.BP,
                        zoom.binSize), zoom.blockBinCount, zoom.blockColumnCount, null, null, reader, null,
                useCache, zoom.isScore() ? Double.NaN : zoom.getCountSum().doubleValue());
        this.v10Reader = reader;
        this.v10Zoom = zoom;
        this.sourceZD = sourceZD;
        this.nBins1 = nBins1;
        this.nBins2 = nBins2;
        this.derivationFactor = sourceZD == null ? 1 : zoom.binSize / sourceZD.v10Zoom.binSize;
    }

    // ------------------------------------------------------------- metadata

    public V10Zoom getZoomDescriptor() {
        return v10Zoom;
    }

    /**
     * True when this resolution stores no blocks and is reconstructed by exact
     * summation from a finer materialized resolution.
     */
    public boolean isDerivedResolution() {
        return v10Zoom.isDerived();
    }

    /**
     * The materialized resolution this one is derived from, or null.
     */
    public V10MatrixZoomData getSourceZoomData() {
        return sourceZD;
    }

    /**
     * Exact number of non-absent cells in the full logical matrix.
     */
    public long getOccupiedCellCount() {
        return v10Zoom.occupiedCellCount;
    }

    /**
     * Population standard deviation among occupied values, or NaN if the writer
     * did not compute it.
     */
    public float getStdDev() {
        return v10Zoom.stdDev;
    }

    /**
     * Estimated 95th percentile among occupied values, or NaN if not computed.
     */
    public float getPercent95() {
        return v10Zoom.percent95;
    }

    /**
     * True when values are genuine f32 scores rather than integer contact counts.
     */
    public boolean isScoreMatrix() {
        return v10Zoom.isScore();
    }

    public long getNumberOfBinsX() {
        return nBins1;
    }

    public long getNumberOfBinsY() {
        return nBins2;
    }

    @Override
    public double getAverageCount() {
        // V10 geometry includes partial BP bins and uses fragment counts for FRAG.
        // The raw whole-matrix mean is the trans expected denominator, including
        // for callers that normalize the observed numerator (legacy V9 behavior).
        if (nBins1 == 0 || nBins2 == 0) return Double.NaN;
        return (v10Zoom.getSum() / nBins1) / nBins2;
    }

    @Override
    public double getExpectedValueForTrans() {
        double average = getAverageCount();
        return average != 0 && Double.isFinite(average) ? average : Double.NaN;
    }

    /**
     * Stored size of one block, matching the V6-V9 contract: the bytes the block
     * occupies on disk, which the exact block index now gives without any read.
     */
    @Override
    public Integer getBlockSize(int blockNum) {
        try {
            List<V10BlockIndexEntry> index = blockIndex();
            int i = V10BlockIndexEntry.find(index, blockNum);
            if (i < 0) return null;
            return V10.toInt(index.get(i).storedByteLength, "stored block length");
        } catch (IOException e) {
            return null;
        }
    }

    @Override
    public void printFullDescription() {
        super.printFullDescription();
        System.out.println("storage: " + (isDerivedResolution()
                ? "derived from " + sourceZD.getBinSize() + " bp" : "materialized"));
        System.out.println("grid: " + (v10Zoom.isRotatedCis() ? "rotated cis" : "rectangular"));
        System.out.println("occupied cells: " + v10Zoom.occupiedCellCount);
        System.out.println();
    }

    // ---------------------------------------------------------- block index

    private synchronized List<V10BlockIndexEntry> blockIndex() throws IOException {
        if (blockIndex == null) {
            blockIndex = v10Reader.readBlockIndex(v10Zoom);
        }
        return blockIndex;
    }

    /**
     * Index entries for the block numbers a query touches, in index order.
     * Numbers with no entry denote empty logical blocks and are dropped.
     */
    private List<V10BlockIndexEntry> entriesFor(List<V10Grid.BlockRange> ranges) throws IOException {
        List<V10BlockIndexEntry> index = blockIndex();
        List<V10BlockIndexEntry> selected = new ArrayList<>();
        for (V10Grid.BlockRange range : ranges) {
            int i = V10BlockIndexEntry.lowerBound(index, range.first);
            while (i < index.size() && index.get(i).blockNumber <= range.last) {
                selected.add(index.get(i++));
            }
        }
        return selected;
    }

    @Override
    void clearCache() {
        super.clearCache();
        blockIndex = null;
        synchronized (derivedRawCache) {
            derivedRawCache.clear();
            derivedRawBytes = 0;
        }
        synchronized (rawBlockCache) {
            rawBlockCache.clear();
            rawBlockBytes = 0;
        }
    }

    // --------------------------------------------------------------- blocks

    @Override
    public List<Block> getNormalizedBlocksOverlapping(long binX1, long binY1, long binX2, long binY2,
                                                      final NormalizationType no, boolean fillUnderDiagonal,
                                                      BlockModifier modifier) {
        try {
            return isDerivedResolution()
                    ? derivedBlocks(binX1, binY1, binX2, binY2, no, modifier)
                    : materializedBlocks(binX1, binY1, binX2, binY2, no, modifier);
        } catch (IOException e) {
            System.err.println("Error reading V10 blocks for " + getKey() + ": " + e.getMessage());
            return new ArrayList<>();
        }
    }

    /**
     * Whole logical blocks overlapping the region, matching the V6-V9 contract:
     * records outside the requested rectangle may be included and are filtered
     * by the caller.
     */
    private List<Block> materializedBlocks(long binX1, long binY1, long binX2, long binY2,
                                           NormalizationType no, BlockModifier modifier) throws IOException {
        List<V10Grid.BlockRange> ranges = V10Grid.blockRangesForRegion(binX1, binY1, binX2, binY2,
                v10Zoom, nBins1, nBins2);
        if (ranges.isEmpty()) return new ArrayList<>();

        List<V10BlockIndexEntry> entries = entriesFor(ranges);
        List<Block> blockList = new ArrayList<>(entries.size());
        List<V10BlockIndexEntry> wanted = new ArrayList<>(entries.size());
        for (V10BlockIndexEntry entry : entries) {
            int blockNumber = entry.blockNumber;
            String key = BlockLoader.getBlockKey(getKey(), blockNumber, no);
            if (blockCache.containsKey(key)) {
                blockList.add(blockCache.get(key));
            } else {
                wanted.add(entry);
            }
        }
        if (wanted.isEmpty()) return blockList;

        LongToDoubleFunction[] norms = normalizationVectors(no);
        List<V10BlockIndexEntry> toRead = new ArrayList<>();
        for (V10BlockIndexEntry entry : wanted) {
            RecordBlock packed = cachedRecordBlock(entry.blockNumber);
            if (packed == null) { toRead.add(entry); continue; }
            String key = BlockLoader.getBlockKey(getKey(), entry.blockNumber, no);
            Block block = modifier.modify(buildPackedBlock(entry.blockNumber, packed, no, norms), key, getBinSize(), chr1, chr2);
            blockCache.put(key, block);
            blockList.add(block);
        }
        BlockScan scan = new BlockScan(v10Reader, toRead);
        V10StoredBlock stored;
        while ((stored = scan.next()) != null) {
            int blockNumber = stored.getBlockNumber();
            String key = BlockLoader.getBlockKey(getKey(), blockNumber, no);
            Block block = buildStoredBlock(stored, no, norms);
            block = modifier.modify(block, key, getBinSize(), chr1, chr2);
            blockCache.put(key, block);
            blockList.add(block);
        }
        return blockList;
    }

    /**
     * Raw source values summed into target bins, then grouped by the block
     * geometry of the target resolution. Unlike the materialized path this
     * returns only cells inside the requested rectangle, because a partial block
     * of a derived resolution would otherwise be indistinguishable from a
     * complete one; such blocks are therefore never cached.
     */
    private List<Block> derivedBlocks(long binX1, long binY1, long binX2, long binY2,
                                      NormalizationType no, BlockModifier modifier) throws IOException {
        String regionKey = binX1 + "_" + binY1 + "_" + binX2 + "_" + binY2;
        List<ContactRecord> raw;
        synchronized (derivedRawCache) { raw = derivedRawCache.get(regionKey); }
        if (raw == null) {
            Map<Long, Accumulator> aggregated = aggregateDerived(binX1, binY1, binX2, binY2, true);
            raw = new ArrayList<>(aggregated.size());
            for (Map.Entry<Long, Accumulator> entry : aggregated.entrySet()) {
                raw.add(entry.getValue().toRecord(columnOf(entry.getKey()), rowOf(entry.getKey()), v10Zoom.isScore()));
            }
            // Conservative allowance for each record, reference and map overhead.
            long bytes = 64L + raw.size() * 64L;
            if (blockCache.getUseCache() && bytes <= DERIVED_CACHE_BYTES) synchronized (derivedRawCache) {
                List<ContactRecord> old = derivedRawCache.remove(regionKey);
                if (old != null) derivedRawBytes -= 64L + old.size() * 64L;
                while (!derivedRawCache.isEmpty() && (derivedRawCache.size() >= 24 || derivedRawBytes + bytes > DERIVED_CACHE_BYTES)) {
                    String oldest = derivedRawCache.keySet().iterator().next();
                    derivedRawBytes -= 64L + derivedRawCache.remove(oldest).size() * 64L;
                }
                derivedRawCache.put(regionKey, raw);
                derivedRawBytes += bytes;
            }
        }
        Map<Integer, List<ContactRecord>> byBlock = new LinkedHashMap<>();
        for (ContactRecord record : raw) {
            int x = record.getBinX();
            int y = record.getBinY();
            int blockNumber = V10Grid.blockNumber(x, y, v10Zoom);
            List<ContactRecord> records = byBlock.get(blockNumber);
            if (records == null) {
                records = new ArrayList<>();
                byBlock.put(blockNumber, records);
            }
            records.add(record);
        }

        LongToDoubleFunction[] norms = normalizationVectors(no);
        List<Block> blockList = new ArrayList<>(byBlock.size());
        for (Map.Entry<Integer, List<ContactRecord>> entry : byBlock.entrySet()) {
            Block block = buildBlock(entry.getKey(), entry.getValue(), no, norms);
            blockList.add(modifier.modify(block, BlockLoader.getBlockKey(getKey(), entry.getKey(), no),
                    getBinSize(), chr1, chr2));
        }
        return blockList;
    }

    /** Exact physical-block access for adapters with a block-based API. */
    public List<Integer> getPhysicalBlockNumbers() throws IOException {
        List<Integer> numbers = new ArrayList<>();
        if (!isDerivedResolution()) for (V10BlockIndexEntry entry : blockIndex()) numbers.add(entry.blockNumber);
        return numbers;
    }

    public Block readNormalizedBlock(int blockNumber, NormalizationType normalization) throws IOException {
        if (isDerivedResolution()) throw new IOException("Derived v10 resolutions have no physical blocks; use regional queries or iterators");
        List<V10BlockIndexEntry> entries = blockIndex();
        int index = V10BlockIndexEntry.lowerBound(entries, blockNumber);
        if (index < entries.size() && entries.get(index).blockNumber == blockNumber) {
            V10BlockIndexEntry entry = entries.get(index);
            RecordBlock packed = cachedRecordBlock(blockNumber);
            LongToDoubleFunction[] norms = normalizationVectors(normalization);
            return packed == null ? buildStoredBlock(v10Reader.readBlock(entry), normalization, norms)
                    : buildPackedBlock(blockNumber, packed, normalization, norms);
        }
        return null;
    }

    private List<ContactRecord> decodeBlock(V10StoredBlock stored) {
        RecordBlock packed = recordBlock(stored);
        if (packed != null) return packed.records(v10Zoom.isScore());
        final List<ContactRecord> records = new ArrayList<>();
        V10BlockDecoder.decode(stored.cursor(), stored.getBlockNumber(), v10Zoom, nBins1, nBins2, isIntra,
                new V10RecordHandler() {
                    @Override
                    public void record(int binColumn, int binRow, long count, float score, boolean isScore) {
                        records.add(isScore ? new ContactRecord(binColumn, binRow, score)
                                : V10ContactRecord.create(binColumn, binRow, count));
                    }
                });
        return records;
    }

    /** Coordinates and exact value bits, decoded once and shared by every norm. */
    private RecordBlock cachedRecordBlock(int number) {
        if (!blockCache.getUseCache()) return null;
        synchronized (rawBlockCache) { return rawBlockCache.get(number); }
    }

    private RecordBlock recordBlock(V10StoredBlock stored) {
        if (!blockCache.getUseCache()) return null;
        RecordBlock cached = cachedRecordBlock(stored.getBlockNumber());
        if (cached != null) return cached;
        V10BlockHeader header = V10BlockDecoder.readHeader(stored.cursor(), v10Zoom, nBins1, nBins2);
        if (header.occupiedCellCount * 16 > RAW_BLOCK_CACHE_BYTES) return null;
        final RecordBlock packed = new RecordBlock((int) header.occupiedCellCount);
        V10BlockDecoder.decode(stored.cursor(), stored.getBlockNumber(), v10Zoom, nBins1, nBins2, isIntra,
                (x, y, count, score, scored) -> {
                    int i = packed.size++;
                    V10.require(i < packed.x.length, "occupied-cell count mismatch");
                    packed.x[i] = x;
                    packed.y[i] = y;
                    packed.values[i] = scored ? Float.floatToRawIntBits(score) : count;
                });
        synchronized (rawBlockCache) {
            RecordBlock previous = rawBlockCache.remove(stored.getBlockNumber());
            if (previous != null) rawBlockBytes -= previous.byteSize();
            while (!rawBlockCache.isEmpty() && (rawBlockCache.size() >= 16
                    || rawBlockBytes + packed.byteSize() > RAW_BLOCK_CACHE_BYTES)) {
                Integer oldest = rawBlockCache.keySet().iterator().next();
                rawBlockBytes -= rawBlockCache.remove(oldest).byteSize();
            }
            rawBlockCache.put(stored.getBlockNumber(), packed);
            rawBlockBytes += packed.byteSize();
        }
        return packed;
    }

    private Block buildStoredBlock(V10StoredBlock stored, NormalizationType no, LongToDoubleFunction[] norms) {
        RecordBlock packed = recordBlock(stored);
        if (packed == null) return buildBlock(stored.getBlockNumber(), decodeBlock(stored), no, norms);
        return buildPackedBlock(stored.getBlockNumber(), packed, no, norms);
    }

    private Block buildPackedBlock(int blockNumber, RecordBlock packed, NormalizationType no, LongToDoubleFunction[] norms) {
        if (no == null || no.equals(NormalizationHandler.NONE)) {
            return new Block(blockNumber, packed.records(v10Zoom.isScore()),
                    BlockLoader.getBlockKey(getKey(), blockNumber, no));
        }
        String key = BlockLoader.getBlockKey(getKey(), blockNumber, no);
        if (norms == null) return new Block(blockNumber, key);
        List<ContactRecord> records = new ArrayList<>(packed.size);
        for (int i = 0; i < packed.size; i++) {
            double denominator = norms[0].applyAsDouble(packed.x[i]) * norms[1].applyAsDouble(packed.y[i]);
            if (!Double.isFinite(denominator) || denominator == 0) continue;
            float value = (float) (packed.value(i, v10Zoom.isScore()) / denominator);
            if (!Float.isNaN(value)) records.add(new ContactRecord(packed.x[i], packed.y[i], value));
        }
        return new Block(blockNumber, records, key);
    }

    private static final class RecordBlock {
        final int[] x, y;
        final long[] values;
        int size;
        RecordBlock(int capacity) {
            x = new int[capacity]; y = new int[capacity]; values = new long[capacity];
        }
        long byteSize() { return 16L * x.length; }
        float value(int i, boolean scored) {
            return scored ? Float.intBitsToFloat((int) values[i]) : values[i] >= 0 ? (float) values[i]
                    : V10.unsignedLong(values[i]).floatValue();
        }
        List<ContactRecord> records(boolean scored) {
            List<ContactRecord> records = new ArrayList<>(size);
            for (int i = 0; i < size; i++) records.add(scored ? new ContactRecord(x[i], y[i], value(i, true))
                    : V10ContactRecord.create(x[i], y[i], values[i]));
            return records;
        }
        int rowStart(long row) {
            int lo = 0, hi = size;
            while (lo < hi) {
                int mid = (lo + hi) >>> 1;
                if (y[mid] < row) lo = mid + 1;
                else hi = mid;
            }
            return lo;
        }
        void visit(long x1, long y1, long x2, long y2, boolean filter, boolean transpose,
                   boolean scored, V10RecordHandler handler) {
            int first = filter ? rowStart(y1) : 0, last = filter ? rowStart(y2 + 1) : size;
            int secondFirst = filter && transpose ? rowStart(x1) : 0;
            int secondLast = filter && transpose ? rowStart(x2 + 1) : 0;
            if (secondFirst < first) {
                int start = first, end = last;
                first = secondFirst; last = secondLast; secondFirst = start; secondLast = end;
            }
            for (int pass = 0; pass < 2; pass++) {
                int from = pass == 0 ? first : Math.max(last, secondFirst);
                int to = pass == 0 ? last : secondLast;
                for (int i = from; i < to; i++) {
                    boolean inside = x[i] >= x1 && x[i] <= x2 && y[i] >= y1 && y[i] <= y2;
                    if (filter && !inside && !(transpose && y[i] >= x1 && y[i] <= x2 && x[i] >= y1 && x[i] <= y2)) continue;
                    handler.record(x[i], y[i], scored ? 0 : values[i], scored ? value(i, true) : 0f, scored);
                }
            }
        }
    }

    // -------------------------------------------------------- normalization

    private LongToDoubleFunction[] normalizationVectors(NormalizationType no) {
        if (no == null || no.equals(NormalizationHandler.NONE)) return null;
        LongToDoubleFunction nv1 = v10Reader.getNormalizationVectorForQuery(getChr1Idx(), getZoom(), no);
        LongToDoubleFunction nv2 = getChr1Idx() == getChr2Idx() ? nv1
                : v10Reader.getNormalizationVectorForQuery(getChr2Idx(), getZoom(), no);
        if (nv1 == null || nv2 == null) {
            System.err.println("Norm " + no + " missing for: " + getKey());
            return null;
        }
        return new LongToDoubleFunction[]{nv1, nv2};
    }

    private Block buildBlock(int blockNumber, List<ContactRecord> records, NormalizationType no,
                             LongToDoubleFunction[] norms) {
        String key = BlockLoader.getBlockKey(getKey(), blockNumber, no);
        if (no == null || no.equals(NormalizationHandler.NONE)) {
            return new Block(blockNumber, records, key);
        }
        if (norms == null) {
            return new Block(blockNumber, key);
        }
        LongToDoubleFunction nv1 = norms[0];
        LongToDoubleFunction nv2 = norms[1];
        List<ContactRecord> normalized = new ArrayList<>(records.size());
        for (ContactRecord record : records) {
            double denominator = nv1.applyAsDouble(record.getBinX()) * nv2.applyAsDouble(record.getBinY());
            if (!Double.isFinite(denominator) || denominator == 0) continue;
            float counts = (float) (record.getCounts() / denominator);
            if (!Float.isNaN(counts)) {
                normalized.add(new ContactRecord(record.getBinX(), record.getBinY(), counts));
            }
        }
        return new Block(blockNumber, normalized, key);
    }

    // ------------------------------------------------------ raw record access

    /**
     * Streams raw (unnormalized) cells of a materialized resolution.
     *
     * @param filterToRegion when true, only cells inside the inclusive bin
     *                       rectangle are emitted; otherwise whole blocks are
     */
    public void streamMaterializedRecords(long binX1, long binY1, long binX2, long binY2,
                                          final boolean filterToRegion, final V10RecordHandler handler)
            throws IOException {
        streamMaterializedRecords(binX1, binY1, binX2, binY2, filterToRegion, false, handler);
    }

    private void streamMaterializedRecords(long binX1, long binY1, long binX2, long binY2,
                                           boolean filterToRegion, boolean includeTransposed, V10RecordHandler handler)
            throws IOException {
        V10.require(!isDerivedResolution(), "resolution is derived; query its source instead");
        List<V10Grid.BlockRange> ranges = V10Grid.blockRangesForRegion(binX1, binY1, binX2, binY2,
                v10Zoom, nBins1, nBins2);
        final long x1 = binX1;
        final long x2 = binX2;
        final long y1 = binY1;
        final long y2 = binY2;

        List<V10BlockIndexEntry> selected = entriesFor(ranges);
        List<V10BlockIndexEntry> missing = new ArrayList<>();
        Map<Integer, RecordBlock> resident = new HashMap<>();
        for (V10BlockIndexEntry entry : selected) {
            RecordBlock packed = cachedRecordBlock(entry.blockNumber);
            if (packed == null) missing.add(entry);
            else resident.put(entry.blockNumber, packed);
        }
        BlockScan scan = new BlockScan(v10Reader, missing);
        for (V10BlockIndexEntry entry : selected) {
            RecordBlock packed = resident.get(entry.blockNumber);
            V10StoredBlock stored = null;
            if (packed == null) {
                stored = scan.next();
                V10.require(stored != null && stored.getBlockNumber() == entry.blockNumber, "missing source block");
                packed = recordBlock(stored);
            }
            if (packed != null) {
                packed.visit(x1, y1, x2, y2, filterToRegion, includeTransposed, v10Zoom.isScore(), handler);
                continue;
            }
            V10BlockDecoder.decode(stored.cursor(), stored.getBlockNumber(), v10Zoom, nBins1, nBins2,
                    isIntra, new V10RecordHandler() {
                        @Override
                        public void record(int binColumn, int binRow, long count, float score, boolean isScore) {
                            boolean inside = binColumn >= x1 && binColumn <= x2 && binRow >= y1 && binRow <= y2;
                            if (filterToRegion && !inside && !(includeTransposed && binRow >= x1 && binRow <= x2
                                    && binColumn >= y1 && binColumn <= y2)) {
                                return;
                            }
                            handler.record(binColumn, binRow, count, score, isScore);
                        }
                    });
        }
    }

    // ------------------------------------------------------------- derivation

    /**
     * Exact summation of raw source values into target bins (Section K).
     *
     * @param includeTransposed for cis matrices, also accept cells whose
     *                          transpose falls inside the requested rectangle
     */
    private Map<Long, Accumulator> aggregateDerived(final long binX1, final long binY1,
                                                    final long binX2, final long binY2,
                                                    final boolean includeTransposed) throws IOException {
        final int factor = derivationFactor;
        final long x1 = Math.max(0, binX1);
        final long x2 = Math.min(nBins1 - 1, binX2);
        final long y1 = Math.max(0, binY1);
        final long y2 = Math.min(nBins2 - 1, binY2);
        final Map<Long, Accumulator> sums = new TreeMap<>();
        if (x1 > x2 || y1 > y2) return sums;

        long sourceX1 = x1 * factor;
        long sourceX2 = Math.min((x2 + 1) * factor - 1, sourceZD.nBins1 - 1);
        long sourceY1 = y1 * factor;
        long sourceY2 = Math.min((y2 + 1) * factor - 1, sourceZD.nBins2 - 1);

        final boolean cis = isIntra;
        sourceZD.streamMaterializedRecords(sourceX1, sourceY1, sourceX2, sourceY2, true, includeTransposed && cis,
                new V10RecordHandler() {
                    @Override
                    public void record(int binColumn, int binRow, long count, float score, boolean isScore) {
                        int x = binColumn / factor;
                        int y = binRow / factor;
                        boolean inside = x >= x1 && x <= x2 && y >= y1 && y <= y2;
                        if (!inside && includeTransposed && cis) {
                            inside = y >= x1 && y <= x2 && x >= y1 && x <= y2;
                        }
                        if (!inside) return;
                        accumulate(sums, x, y, binColumn, binRow, count, score, isScore);
                    }
                });
        return sums;
    }

    private static void accumulate(Map<Long, Accumulator> sums, int x, int y, int sourceX, int sourceY, long count, float score,
                                   boolean isScore) {
        Long key = cellKey(x, y);
        Accumulator accumulator = sums.get(key);
        if (accumulator == null) {
            accumulator = new Accumulator();
            sums.put(key, accumulator);
        }
        if (isScore) {
            V10.require(!Float.isNaN(score) && !Float.isInfinite(score), "nonfinite derived source score");
            if (accumulator.scores == null) accumulator.scores = new TreeMap<>();
            V10.require(accumulator.scores.put(cellKey(sourceX, sourceY), score) == null,
                    "duplicate derived source cell");
        } else {
            accumulator.count = V10.addUnsigned(accumulator.count, count);
        }
    }

    private static Long cellKey(int binColumn, int binRow) {
        return ((long) binRow << 32) | (binColumn & 0xFFFFFFFFL);
    }

    private static int columnOf(long key) {
        return (int) (key & 0xFFFFFFFFL);
    }

    private static int rowOf(long key) {
        return (int) (key >>> 32);
    }

    /**
     * Counts accumulate exactly as unsigned 64-bit integers; scores accumulate
     * in binary64 and are rounded once to binary32, as Section K requires.
     */
    private static class Accumulator {
        BigInteger count = BigInteger.ZERO;
        TreeMap<Long, Float> scores;

        ContactRecord toRecord(int binColumn, int binRow, boolean isScore) {
            if (!isScore) return V10ContactRecord.create(binColumn, binRow, count);
            double score = 0;
            for (float value : scores.values()) {
                score += value;
                V10.require(Double.isFinite(score), "derived score overflow");
            }
            return new ContactRecord(binColumn, binRow, (float) score);
        }
    }

    // ------------------------------------------------------------ iterators

    @Override
    public Iterator<ContactRecord> getDirectIterator() {
        return getNormalizedIterator(NormalizationHandler.NONE);
    }

    @Override
    public Iterator<ContactRecord> getNormalizedIterator(NormalizationType normType) {
        try {
            Iterator<ContactRecord> raw = isDerivedResolution()
                    ? new DerivedIterator() : new MaterializedIterator();
            if (normType == null || normType.equals(NormalizationHandler.NONE)) return raw;
            LongToDoubleFunction[] norms = normalizationVectors(normType);
            if (norms == null) return Collections.<ContactRecord>emptyList().iterator();
            return new NormalizingIterator(raw, norms[0], norms[1]);
        } catch (IOException e) {
            throw new V10FormatException("could not iterate " + getKey(), e);
        }
    }

    /**
     * Walks every stored block of this resolution in block-number order,
     * fetching them in batches so adjacent records coalesce into one read.
     */
    private class MaterializedIterator implements Iterator<ContactRecord> {
        private final BlockScan scan;
        private Iterator<ContactRecord> current;

        MaterializedIterator() throws IOException {
            scan = new BlockScan(v10Reader, blockIndex());
        }

        @Override
        public boolean hasNext() {
            while (current == null || !current.hasNext()) {
                try {
                    V10StoredBlock stored = scan.next();
                    if (stored == null) return false;
                    current = decodeBlock(stored).iterator();
                } catch (IOException e) {
                    throw new V10FormatException("could not read block", e);
                }
            }
            return true;
        }

        @Override
        public ContactRecord next() {
            if (!hasNext()) throw new NoSuchElementException();
            return current.next();
        }

        @Override
        public void remove() {
            throw new UnsupportedOperationException("remove() is not supported");
        }
    }

    /**
     * Sequential cursor over a list of stored blocks. Blocks are fetched
     * {@link #STREAM_BATCH_SIZE} at a time so a broad query does not issue one
     * read per block, while only one batch is decompressed at any moment.
     */
    private static class BlockScan {
        private final DatasetReaderV10 reader;
        private final List<V10BlockIndexEntry> entries;
        private int at = 0;
        private Map<Integer, V10StoredBlock> batch;

        BlockScan(DatasetReaderV10 reader, List<V10BlockIndexEntry> entries) {
            this.reader = reader;
            this.entries = entries;
        }

        /**
         * The next stored block, or null once the entries are exhausted.
         */
        V10StoredBlock next() throws IOException {
            if (at >= entries.size()) return null;
            if (at % STREAM_BATCH_SIZE == 0) {
                batch = reader.readBlocks(
                        entries.subList(at, Math.min(at + STREAM_BATCH_SIZE, entries.size())));
            }
            V10StoredBlock stored = batch.get(entries.get(at++).blockNumber);
            V10.require(stored != null, "indexed block was not fetched");
            return stored;
        }
    }

    /**
     * Iterates a derived resolution without materializing the whole matrix.
     * <p>
     * Work is driven by the source blocks. For each source block the target
     * bins it touches are aggregated over a target-aligned neighbourhood, so
     * every contributing source cell is seen even when a target bin straddles a
     * source block boundary. A target bin is emitted by the single source block
     * that owns its smallest present source cell, which yields each derived cell
     * exactly once.
     */
    private class DerivedIterator implements Iterator<ContactRecord> {
        private final BlockScan scan;
        private Iterator<ContactRecord> current;

        DerivedIterator() throws IOException {
            scan = new BlockScan(sourceZD.v10Reader, sourceZD.blockIndex());
        }

        @Override
        public boolean hasNext() {
            while (current == null || !current.hasNext()) {
                try {
                    V10StoredBlock stored = scan.next();
                    if (stored == null) return false;
                    current = recordsOwnedBy(stored).iterator();
                } catch (IOException e) {
                    throw new V10FormatException("could not read block", e);
                }
            }
            return true;
        }

        private List<ContactRecord> recordsOwnedBy(V10StoredBlock sourceBlock) throws IOException {
            final int blockNumber = sourceBlock.getBlockNumber();
            V10Cursor cursor = sourceBlock.cursor();
            V10BlockHeader header = V10BlockDecoder.readHeader(cursor, sourceZD.v10Zoom,
                    sourceZD.nBins1, sourceZD.nBins2);

            final int factor = derivationFactor;
            long sourceX1 = (header.binColumnOffset / factor) * factor;
            long sourceX2 = Math.min((header.lastBinColumn() / factor + 1) * factor - 1, sourceZD.nBins1 - 1);
            long sourceY1 = (header.binRowOffset / factor) * factor;
            long sourceY2 = Math.min((header.lastBinRow() / factor + 1) * factor - 1, sourceZD.nBins2 - 1);

            final Map<Long, Accumulator> sums = new TreeMap<>();
            final Map<Long, Long> smallestSourceCell = new HashMap<>();
            sourceZD.streamMaterializedRecords(sourceX1, sourceY1, sourceX2, sourceY2, true,
                    new V10RecordHandler() {
                        @Override
                        public void record(int binColumn, int binRow, long count, float score, boolean isScore) {
                            int x = binColumn / factor;
                            int y = binRow / factor;
                            Long key = cellKey(x, y);
                            accumulate(sums, x, y, binColumn, binRow, count, score, isScore);
                            Long candidate = cellKey(binColumn, binRow);
                            Long existing = smallestSourceCell.get(key);
                            if (existing == null || candidate < existing) {
                                smallestSourceCell.put(key, candidate);
                            }
                        }
                    });

            List<ContactRecord> owned = new ArrayList<>();
            for (Map.Entry<Long, Accumulator> entry : sums.entrySet()) {
                Long smallest = smallestSourceCell.get(entry.getKey());
                if (V10Grid.blockNumber(columnOf(smallest), rowOf(smallest), sourceZD.v10Zoom) != blockNumber) {
                    continue;
                }
                owned.add(entry.getValue().toRecord(columnOf(entry.getKey()), rowOf(entry.getKey()),
                        v10Zoom.isScore()));
            }
            return owned;
        }

        @Override
        public ContactRecord next() {
            if (!hasNext()) throw new NoSuchElementException();
            return current.next();
        }

        @Override
        public void remove() {
            throw new UnsupportedOperationException("remove() is not supported");
        }
    }

    /**
     * Applies normalization vectors lazily, skipping cells whose normalized
     * value is not available (Section J.7).
     */
    private static class NormalizingIterator implements Iterator<ContactRecord> {
        private final Iterator<ContactRecord> source;
        private final LongToDoubleFunction nv1;
        private final LongToDoubleFunction nv2;
        private ContactRecord next;

        NormalizingIterator(Iterator<ContactRecord> source, LongToDoubleFunction nv1, LongToDoubleFunction nv2) {
            this.source = source;
            this.nv1 = nv1;
            this.nv2 = nv2;
        }

        @Override
        public boolean hasNext() {
            while (next == null && source.hasNext()) {
                ContactRecord record = source.next();
                double denominator = nv1.applyAsDouble(record.getBinX()) * nv2.applyAsDouble(record.getBinY());
                if (!Double.isFinite(denominator) || denominator == 0) continue;
                float counts = (float) (record.getCounts() / denominator);
                if (!Float.isNaN(counts)) {
                    next = new ContactRecord(record.getBinX(), record.getBinY(), counts);
                }
            }
            return next != null;
        }

        @Override
        public ContactRecord next() {
            if (!hasNext()) throw new NoSuchElementException();
            ContactRecord result = next;
            next = null;
            return result;
        }

        @Override
        public void remove() {
            throw new UnsupportedOperationException("remove() is not supported");
        }
    }
}
