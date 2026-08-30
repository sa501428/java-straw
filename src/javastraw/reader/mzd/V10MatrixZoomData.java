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
import javastraw.reader.datastructures.ListOfDoubleArrays;
import javastraw.reader.norm.NormalizationVector;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationHandler;
import javastraw.reader.type.NormalizationType;
import javastraw.reader.v10.V10;
import javastraw.reader.v10.V10BlockDecoder;
import javastraw.reader.v10.V10BlockHeader;
import javastraw.reader.v10.V10Cursor;
import javastraw.reader.v10.V10DecodedPage;
import javastraw.reader.v10.V10FormatException;
import javastraw.reader.v10.V10Grid;
import javastraw.reader.v10.V10Page;
import javastraw.reader.v10.V10RecordHandler;
import javastraw.reader.v10.V10Zoom;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.TreeMap;

/**
 * A V10 matrix at one logical resolution.
 * <p>
 * Two storage modes are handled transparently. A <b>materialized</b> resolution
 * reads Zstandard-compressed pages and decodes the logical blocks they contain.
 * A <b>derived</b> resolution stores no blocks: its raw values are reconstructed
 * by exact summation from the declared finer materialized source, always before
 * normalization or expected-value division (Section K). Either way, callers see
 * the same {@link Block} and {@link ContactRecord} objects as for V6-V9 files.
 */
public class V10MatrixZoomData extends MatrixZoomData {

    private final DatasetReaderV10 v10Reader;
    private final V10Zoom v10Zoom;
    /**
     * Non-null exactly when this resolution is derived.
     */
    private final V10MatrixZoomData sourceZD;
    private final long nBins1;
    private final long nBins2;
    private final int derivationFactor;

    private List<V10Page> pages;

    public V10MatrixZoomData(Chromosome chr1, Chromosome chr2, V10Zoom zoom, DatasetReaderV10 reader,
                             V10MatrixZoomData sourceZD, boolean useCache, long nBins1, long nBins2) {
        super(chr1, chr2, new HiCZoom(zoom.unit == V10.UNIT_FRAG ? HiCZoom.HiCUnit.FRAG : HiCZoom.HiCUnit.BP,
                        zoom.binSize), zoom.blockBinCount, zoom.blockColumnCount, null, null, reader, null,
                useCache, zoom.isScore() ? Double.NaN : (double) zoom.getCountSum());
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
        if (v10Zoom.isScore()) return Double.NaN;
        return super.getAverageCount();
    }

    @Override
    public Integer getBlockSize(int blockNum) {
        try {
            List<V10Page> pageList = pages();
            int pageIndex = V10Page.findPage(pageList, blockNum);
            if (pageIndex < 0) return null;
            V10DecodedPage page = v10Reader.readPage(pageList.get(pageIndex));
            int i = page.indexOf(blockNum);
            return i < 0 ? null : page.cursorFor(i).size();
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

    // ---------------------------------------------------------------- pages

    private synchronized List<V10Page> pages() throws IOException {
        if (pages == null) {
            pages = v10Reader.readPageIndex(v10Zoom);
        }
        return pages;
    }

    @Override
    void clearCache() {
        super.clearCache();
        pages = null;
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
        List<Integer> blockNumbers = V10Grid.blockNumbersForRegion(binX1, binY1, binX2, binY2,
                v10Zoom, nBins1, nBins2);
        if (blockNumbers.isEmpty()) return new ArrayList<>();

        List<Block> blockList = new ArrayList<>(blockNumbers.size());
        // Group the blocks still needed by the page that stores them so each page
        // is fetched and decompressed at most once.
        Map<Integer, List<Integer>> byPage = new LinkedHashMap<>();
        List<V10Page> pageList = pages();
        for (Integer blockNumber : blockNumbers) {
            String key = BlockLoader.getBlockKey(getKey(), blockNumber, no);
            if (blockCache.containsKey(key)) {
                blockList.add(blockCache.get(key));
                continue;
            }
            int pageIndex = V10Page.findPage(pageList, blockNumber);
            if (pageIndex < 0) continue;
            List<Integer> wanted = byPage.get(pageIndex);
            if (wanted == null) {
                wanted = new ArrayList<>();
                byPage.put(pageIndex, wanted);
            }
            wanted.add(blockNumber);
        }

        NormalizationVector[] norms = normalizationVectors(no);
        for (Map.Entry<Integer, List<Integer>> entry : byPage.entrySet()) {
            V10DecodedPage page = v10Reader.readPage(pageList.get(entry.getKey()));
            for (Integer blockNumber : entry.getValue()) {
                String key = BlockLoader.getBlockKey(getKey(), blockNumber, no);
                int indexInPage = page.indexOf(blockNumber);
                Block block;
                if (indexInPage < 0) {
                    // Numbers inside a page range may legitimately be absent.
                    block = new Block(blockNumber, key);
                } else {
                    block = buildBlock(blockNumber, decodeBlock(page, indexInPage, blockNumber), no, norms);
                }
                block = modifier.modify(block, key, getBinSize(), chr1, chr2);
                blockCache.put(key, block);
                blockList.add(block);
            }
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
        Map<Long, Accumulator> aggregated = aggregateDerived(binX1, binY1, binX2, binY2, true);
        Map<Integer, List<ContactRecord>> byBlock = new LinkedHashMap<>();
        for (Map.Entry<Long, Accumulator> entry : aggregated.entrySet()) {
            int x = columnOf(entry.getKey());
            int y = rowOf(entry.getKey());
            int blockNumber = V10Grid.blockNumber(x, y, v10Zoom);
            List<ContactRecord> records = byBlock.get(blockNumber);
            if (records == null) {
                records = new ArrayList<>();
                byBlock.put(blockNumber, records);
            }
            records.add(entry.getValue().toRecord(x, y, v10Zoom.isScore()));
        }

        NormalizationVector[] norms = normalizationVectors(no);
        List<Block> blockList = new ArrayList<>(byBlock.size());
        for (Map.Entry<Integer, List<ContactRecord>> entry : byBlock.entrySet()) {
            Block block = buildBlock(entry.getKey(), entry.getValue(), no, norms);
            blockList.add(modifier.modify(block, BlockLoader.getBlockKey(getKey(), entry.getKey(), no),
                    getBinSize(), chr1, chr2));
        }
        return blockList;
    }

    private List<ContactRecord> decodeBlock(V10DecodedPage page, int indexInPage, int blockNumber) {
        final List<ContactRecord> records = new ArrayList<>();
        V10BlockDecoder.decode(page.cursorFor(indexInPage), blockNumber, v10Zoom, nBins1, nBins2, isIntra,
                new V10RecordHandler() {
                    @Override
                    public void record(int binColumn, int binRow, long count, float score, boolean isScore) {
                        records.add(isScore ? new ContactRecord(binColumn, binRow, score)
                                : V10ContactRecord.create(binColumn, binRow, count));
                    }
                });
        return records;
    }

    // -------------------------------------------------------- normalization

    private NormalizationVector[] normalizationVectors(NormalizationType no) {
        if (no == null || no.equals(NormalizationHandler.NONE)) return null;
        NormalizationVector nv1 = v10Reader.getNormalizationVector(getChr1Idx(), getZoom(), no);
        NormalizationVector nv2 = v10Reader.getNormalizationVector(getChr2Idx(), getZoom(), no);
        if (nv1 == null || nv2 == null) {
            System.err.println("Norm " + no + " missing for: " + getKey());
            return null;
        }
        return new NormalizationVector[]{nv1, nv2};
    }

    private Block buildBlock(int blockNumber, List<ContactRecord> records, NormalizationType no,
                             NormalizationVector[] norms) {
        String key = BlockLoader.getBlockKey(getKey(), blockNumber, no);
        if (no == null || no.equals(NormalizationHandler.NONE)) {
            return new Block(blockNumber, records, key);
        }
        if (norms == null) {
            return new Block(blockNumber, key);
        }
        ListOfDoubleArrays nv1 = norms[0].getData();
        ListOfDoubleArrays nv2 = norms[1].getData();
        List<ContactRecord> normalized = new ArrayList<>(records.size());
        for (ContactRecord record : records) {
            double denominator = nv1.get(record.getBinX()) * nv2.get(record.getBinY());
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
        V10.require(!isDerivedResolution(), "resolution is derived; query its source instead");
        List<V10Page> pageList = pages();
        List<Integer> blockNumbers = V10Grid.blockNumbersForRegion(binX1, binY1, binX2, binY2,
                v10Zoom, nBins1, nBins2);
        final long x1 = binX1;
        final long x2 = binX2;
        final long y1 = binY1;
        final long y2 = binY2;

        Map<Integer, List<Integer>> byPage = new LinkedHashMap<>();
        for (Integer blockNumber : blockNumbers) {
            int pageIndex = V10Page.findPage(pageList, blockNumber);
            if (pageIndex < 0) continue;
            List<Integer> wanted = byPage.get(pageIndex);
            if (wanted == null) {
                wanted = new ArrayList<>();
                byPage.put(pageIndex, wanted);
            }
            wanted.add(blockNumber);
        }
        for (Map.Entry<Integer, List<Integer>> entry : byPage.entrySet()) {
            V10DecodedPage page = v10Reader.readPage(pageList.get(entry.getKey()));
            for (Integer blockNumber : entry.getValue()) {
                int indexInPage = page.indexOf(blockNumber);
                if (indexInPage < 0) continue;
                V10BlockDecoder.decode(page.cursorFor(indexInPage), blockNumber, v10Zoom, nBins1, nBins2, isIntra,
                        new V10RecordHandler() {
                            @Override
                            public void record(int binColumn, int binRow, long count, float score, boolean isScore) {
                                if (filterToRegion && (binColumn < x1 || binColumn > x2
                                        || binRow < y1 || binRow > y2)) {
                                    return;
                                }
                                handler.record(binColumn, binRow, count, score, isScore);
                            }
                        });
            }
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
        sourceZD.streamMaterializedRecords(sourceX1, sourceY1, sourceX2, sourceY2, false,
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
                        accumulate(sums, x, y, count, score, isScore);
                    }
                });
        return sums;
    }

    private static void accumulate(Map<Long, Accumulator> sums, int x, int y, long count, float score,
                                   boolean isScore) {
        Long key = cellKey(x, y);
        Accumulator accumulator = sums.get(key);
        if (accumulator == null) {
            accumulator = new Accumulator();
            sums.put(key, accumulator);
        }
        if (isScore) {
            V10.require(!Float.isNaN(score) && !Float.isInfinite(score), "nonfinite derived source score");
            accumulator.score += score;
            V10.require(!Double.isNaN(accumulator.score) && !Double.isInfinite(accumulator.score),
                    "derived score overflow");
        } else {
            accumulator.count = V10.add(accumulator.count, count);
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
     * Counts accumulate exactly in a long; scores accumulate in binary64 and are
     * rounded once to binary32, as Section K requires.
     */
    private static class Accumulator {
        long count;
        double score;

        ContactRecord toRecord(int binColumn, int binRow, boolean isScore) {
            return isScore ? new ContactRecord(binColumn, binRow, (float) score)
                    : V10ContactRecord.create(binColumn, binRow, count);
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
            NormalizationVector[] norms = normalizationVectors(normType);
            if (norms == null) return Collections.<ContactRecord>emptyList().iterator();
            return new NormalizingIterator(raw, norms[0].getData(), norms[1].getData());
        } catch (IOException e) {
            throw new V10FormatException("could not iterate " + getKey(), e);
        }
    }

    /**
     * Walks every stored page in order, decoding one logical block at a time.
     */
    private class MaterializedIterator implements Iterator<ContactRecord> {
        private final List<V10Page> pageList;
        private int pageIndex = 0;
        private int blockIndex = 0;
        private V10DecodedPage page;
        private Iterator<ContactRecord> current;

        MaterializedIterator() throws IOException {
            pageList = pages();
        }

        @Override
        public boolean hasNext() {
            while (current == null || !current.hasNext()) {
                try {
                    if (page == null || blockIndex >= page.blockCount()) {
                        if (pageIndex >= pageList.size()) return false;
                        page = v10Reader.readPage(pageList.get(pageIndex++));
                        blockIndex = 0;
                    }
                    int blockNumber = page.blockNumber(blockIndex);
                    current = decodeBlock(page, blockIndex, blockNumber).iterator();
                    blockIndex++;
                } catch (IOException e) {
                    throw new V10FormatException("could not read page", e);
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
        private final List<V10Page> pageList;
        private int pageIndex = 0;
        private int blockIndex = 0;
        private V10DecodedPage page;
        private Iterator<ContactRecord> current;

        DerivedIterator() throws IOException {
            pageList = sourceZD.pages();
        }

        @Override
        public boolean hasNext() {
            while (current == null || !current.hasNext()) {
                try {
                    if (page == null || blockIndex >= page.blockCount()) {
                        if (pageIndex >= pageList.size()) return false;
                        page = sourceZD.v10Reader.readPage(pageList.get(pageIndex++));
                        blockIndex = 0;
                    }
                    current = recordsOwnedBy(page, blockIndex).iterator();
                    blockIndex++;
                } catch (IOException e) {
                    throw new V10FormatException("could not read page", e);
                }
            }
            return true;
        }

        private List<ContactRecord> recordsOwnedBy(V10DecodedPage sourcePage, int indexInPage)
                throws IOException {
            final int blockNumber = sourcePage.blockNumber(indexInPage);
            V10Cursor cursor = sourcePage.cursorFor(indexInPage);
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
                            accumulate(sums, x, y, count, score, isScore);
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
        private final ListOfDoubleArrays nv1;
        private final ListOfDoubleArrays nv2;
        private ContactRecord next;

        NormalizingIterator(Iterator<ContactRecord> source, ListOfDoubleArrays nv1, ListOfDoubleArrays nv2) {
            this.source = source;
            this.nv1 = nv1;
            this.nv2 = nv2;
        }

        @Override
        public boolean hasNext() {
            while (next == null && source.hasNext()) {
                ContactRecord record = source.next();
                double denominator = nv1.get(record.getBinX()) * nv2.get(record.getBinY());
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
