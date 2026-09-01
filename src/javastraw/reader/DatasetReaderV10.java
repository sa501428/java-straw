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

package javastraw.reader;

import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.block.Block;
import javastraw.reader.block.IndexEntry;
import javastraw.reader.datastructures.ListOfDoubleArrays;
import javastraw.reader.expected.ExpectedValueFunction;
import javastraw.reader.mzd.Matrix;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.mzd.V10MatrixZoomData;
import javastraw.reader.norm.NormalizationVector;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationHandler;
import javastraw.reader.type.NormalizationType;
import javastraw.reader.v10.V10;
import javastraw.reader.v10.V10BlockIndexEntry;
import javastraw.reader.v10.V10Cursor;
import javastraw.reader.v10.V10ExpectedValueFunction;
import javastraw.reader.v10.V10FormatException;
import javastraw.reader.v10.V10Header;
import javastraw.reader.v10.V10Locator;
import javastraw.reader.v10.V10Source;
import javastraw.reader.v10.V10StoredBlock;
import javastraw.reader.v10.V10VectorEntry;
import javastraw.reader.v10.V10VectorIndex;
import javastraw.reader.v10.V10Zoom;
import org.broad.igv.util.collections.LRUCache;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import static javastraw.reader.v10.V10.require;

/**
 * Reader for the consolidated V10 .hic format (hic-format/HiCFormatV10.md).
 * <p>
 * V10 is a distinct wire format rather than an extension of V6-V9: numeric
 * chromosome-pair keys, independently Zstandard-compressed logical blocks
 * addressed by an exact binary block index, flat delta-coded cell positions with
 * separate value streams, explicit integer counts, chunked normalization and
 * expected-value vectors, and resolutions that may be derived by exact summation
 * from a finer materialized resolution.
 * <p>
 * The public java-straw API is unchanged: this reader produces the same
 * {@link Dataset}, {@link Matrix}, {@link MatrixZoomData} and
 * {@link javastraw.reader.block.ContactRecord} objects as
 * {@link DatasetReaderV2}, so existing tools work against V10 files without
 * modification.
 */
public class DatasetReaderV10 extends AbstractDatasetReader {

    /**
     * The logical block is the decompression and cache unit of the format
     * (Section H), so a small decompressed-block cache is kept even when
     * block caching is off. It is simply larger when caching is on.
     */
    private static final int BLOCK_CACHE_SIZE_WITH_CACHING = 256;
    private static final int BLOCK_CACHE_SIZE_MINIMAL = 16;

    /**
     * Upper bound on one coalesced range read. Only physically adjacent stored
     * records are ever merged, so no unrequested block is fetched; this cap just
     * keeps a broad query from turning into one enormous read.
     */
    private static final long MAX_COALESCED_READ = 8L * 1024 * 1024;

    private static final Comparator<V10BlockIndexEntry> BY_POSITION = new Comparator<V10BlockIndexEntry>() {
        @Override
        public int compare(V10BlockIndexEntry a, V10BlockIndexEntry b) {
            return Long.compare(a.blockPosition, b.blockPosition);
        }
    };

    private final Dataset dataset;
    private final boolean useCache;

    private V10Source source;
    private V10Header header;
    private final Map<String, V10Locator> matrixDirectory = new LinkedHashMap<>();
    private final Map<String, List<V10Zoom>> zoomCache = Collections.synchronizedMap(new HashMap<String, List<V10Zoom>>());
    private final LRUCache<Long, V10StoredBlock> decodedBlockCache;

    private V10VectorIndex normIndex;
    private V10VectorIndex expectedIndex;
    private V10VectorIndex normExpectedIndex;

    private boolean activeStatus = true;

    public DatasetReaderV10(String path, boolean useCache) {
        super(path);
        this.useCache = useCache;
        this.decodedBlockCache = new LRUCache<>(useCache ? BLOCK_CACHE_SIZE_WITH_CACHING : BLOCK_CACHE_SIZE_MINIMAL);
        this.dataset = new Dataset(this);
    }

    /**
     * True when the file at {@code path} declares version 10.
     */
    public static boolean isV10(String path) throws IOException {
        V10Source probe = new V10Source(path);
        V10Cursor c = new V10Cursor(probe.read(0, 8));
        c.magic("HIC\0");
        long version = c.word();
        require(version >= 6 && version <= V10.VERSION, "unsupported .hic version " + version);
        return version == V10.VERSION;
    }

    @Override
    public Dataset read() throws IOException {
        source = new V10Source(path);

        V10Cursor prefix = new V10Cursor(source.read(0, 88));
        prefix.magic("HIC\0");
        require(prefix.word() == V10.VERSION, "not a V10 file");
        long headerLength = prefix.wide();
        require(headerLength >= 88 && headerLength <= V10.ALLOCATION_LIMIT, "invalid header length");

        header = V10Header.parse(source.read(0, headerLength));
        for (V10Locator locator : new V10Locator[]{header.footer, header.normVectorIndex,
                header.expectedValueIndex, header.normExpectedValueIndex}) {
            if (locator.isPresent()) {
                source.interval(locator.position, locator.length);
            }
        }

        readFooter();

        dataset.setAttributes(new LinkedHashMap<>(header.attributes));
        // V10 stores exact chromosome names and lengths and the numeric index is
        // the matrix key, so entry zero is never rewritten into an all-by-all
        // pseudo-chromosome the way V6-V9 files assume.
        ChromosomeHandler chromosomeHandler = new ChromosomeHandler(
                new ArrayList<>(header.chromosomes), header.genomeId, false, false);
        dataset.setChromosomeHandler(chromosomeHandler);
        dataset.setGenomeId(chromosomeHandler.getGenomeID());
        dataset.setBpZooms(header.binSizes(V10.UNIT_BP));
        dataset.setFragZooms(header.binSizes(V10.UNIT_FRAG));
        if (!header.resolutions[V10.UNIT_FRAG].isEmpty()) {
            Map<String, Integer> fragmentCounts = new HashMap<>();
            for (int i = 0; i < header.chromosomes.size(); i++) {
                fragmentCounts.put(header.chromosomes.get(i).getName(), header.fragmentSiteCounts[i]);
            }
            dataset.setFragmentCounts(fragmentCounts);
        }

        readVectorIndexes();

        return dataset;
    }

    private void readFooter() throws IOException {
        V10Cursor c = new V10Cursor(source.read(header.footer.position, header.footer.length));
        c.magic("H10F");
        require(c.word() == 1, "unknown footer version");
        require(c.wide() == header.footer.length, "footer length mismatch");
        long entryCount = c.word();
        c.zero(4);
        require(entryCount * 24 == c.left(), "invalid footer count");

        long previous = -1;
        for (long i = 0; i < entryCount; i++) {
            long chr1 = c.word();
            long chr2 = c.word();
            V10Locator locator = V10Locator.read(c);
            require(chr1 <= chr2 && chr2 < header.chromosomes.size(), "invalid matrix directory key");
            long key = (chr1 << 32) | chr2;
            require(i == 0 || key > previous, "invalid matrix directory key");
            previous = key;
            source.interval(locator.position, locator.length);
            matrixDirectory.put(chr1 + "_" + chr2, locator);
        }
        c.done();
    }

    private void readVectorIndexes() throws IOException {
        normIndex = V10VectorIndex.parse(source, header.normVectorIndex, V10VectorIndex.KIND_NORM, header);
        expectedIndex = V10VectorIndex.parse(source, header.expectedValueIndex, V10VectorIndex.KIND_EXPECTED, header);
        normExpectedIndex = V10VectorIndex.parse(source, header.normExpectedValueIndex,
                V10VectorIndex.KIND_NORM_EXPECTED, header);

        // Index presence is the capability advertisement: only normalizations
        // that actually carry vectors are offered to callers.
        Set<Integer> available = new LinkedHashSet<>();
        for (V10VectorEntry entry : normIndex.getEntries()) {
            available.add(entry.normalizationTypeId);
        }
        for (Integer id : available) {
            dataset.addNormalizationType(
                    dataset.getNormalizationHandler().getNormTypeFromString(header.normalizations.get(id)));
        }

        Map<String, ExpectedValueFunction> expectedValues = new LinkedHashMap<>();
        for (V10VectorEntry entry : expectedIndex.getEntries()) {
            addExpectedValueFunction(expectedValues, entry, NormalizationHandler.NONE);
        }
        for (V10VectorEntry entry : normExpectedIndex.getEntries()) {
            NormalizationType norm = dataset.getNormalizationHandler()
                    .getNormTypeFromString(header.normalizations.get(entry.normalizationTypeId));
            addExpectedValueFunction(expectedValues, entry, norm);
        }
        dataset.setExpectedValueFunctionMap(expectedValues);
    }

    private void addExpectedValueFunction(Map<String, ExpectedValueFunction> map, V10VectorEntry entry,
                                          NormalizationType norm) {
        HiCZoom.HiCUnit unit = entry.unit == V10.UNIT_FRAG ? HiCZoom.HiCUnit.FRAG : HiCZoom.HiCUnit.BP;
        map.put(ExpectedValueFunction.getKey(unit, entry.binSize, norm),
                new V10ExpectedValueFunction(norm, unit, entry.binSize, entry, source));
    }

    @Override
    public Matrix readMatrix(String key, int specificResolution) throws IOException {
        V10Locator locator = matrixDirectory.get(key);
        if (locator == null) {
            // A chromosome pair with no directory entry is an all-zero matrix.
            return null;
        }
        String[] parts = key.split("_");
        int chr1Index = Integer.parseInt(parts[0]);
        int chr2Index = Integer.parseInt(parts[1]);
        Chromosome chr1 = dataset.getChromosomeHandler().getChromosomeFromIndex(chr1Index);
        Chromosome chr2 = dataset.getChromosomeHandler().getChromosomeFromIndex(chr2Index);

        List<V10Zoom> zooms = getZooms(key, chr1Index, chr2Index, locator);
        int bpCount = header.resolutions[V10.UNIT_BP].size();

        List<MatrixZoomData> zdList = new ArrayList<>(zooms.size());
        V10MatrixZoomData[] byDescriptor = new V10MatrixZoomData[zooms.size()];
        for (int i = 0; i < zooms.size(); i++) {
            V10Zoom zoom = zooms.get(i);
            V10MatrixZoomData sourceZD = null;
            if (zoom.isDerived()) {
                int base = zoom.unit == V10.UNIT_FRAG ? bpCount : 0;
                sourceZD = byDescriptor[base + (int) zoom.sourceResolutionIndex];
                require(sourceZD != null, "derived resolution precedes its source");
            }
            V10MatrixZoomData zd = new V10MatrixZoomData(chr1, chr2, zoom, this, sourceZD, useCache,
                    header.bins(chr1Index, zoom.unit, zoom.resolutionIndex),
                    header.bins(chr2Index, zoom.unit, zoom.resolutionIndex));
            byDescriptor[i] = zd;
            zdList.add(zd);
        }
        return new Matrix(chr1Index, chr2Index, zdList);
    }

    private List<V10Zoom> getZooms(String key, int chr1Index, int chr2Index, V10Locator locator) throws IOException {
        List<V10Zoom> cached = zoomCache.get(key);
        if (cached != null) return cached;
        List<V10Zoom> zooms = V10Zoom.parseMatrixRecord(source.read(locator.position, locator.length),
                chr1Index, chr2Index, header, source);
        zoomCache.put(key, zooms);
        return zooms;
    }

    /**
     * Resolution descriptors for a chromosome pair, or null when the pair has no
     * stored matrix. Exposes exact per-resolution statistics without decoding
     * any contact data.
     */
    public List<V10Zoom> getZoomDescriptors(int chr1Index, int chr2Index) throws IOException {
        int a = Math.min(chr1Index, chr2Index);
        int b = Math.max(chr1Index, chr2Index);
        String key = a + "_" + b;
        V10Locator locator = matrixDirectory.get(key);
        if (locator == null) return null;
        return getZooms(key, a, b, locator);
    }

    // --------------------------------------------------------------- blocks

    public V10Source getSource() {
        return source;
    }

    public V10Header getHeader() {
        return header;
    }

    /**
     * Reads and validates the exact block index of a materialized resolution.
     */
    public List<V10BlockIndexEntry> readBlockIndex(V10Zoom zoom) throws IOException {
        if (!zoom.blockIndex.isPresent()) return Collections.emptyList();
        return V10BlockIndexEntry.parseIndex(source.read(zoom.blockIndex.position, zoom.blockIndex.length),
                zoom, source);
    }

    /**
     * Fetches, validates and decompresses one stored block, reusing the cached
     * copy when caching is enabled.
     */
    public V10StoredBlock readBlock(V10BlockIndexEntry entry) throws IOException {
        V10StoredBlock cached = cachedBlock(entry);
        if (cached != null) return cached;
        V10StoredBlock decoded = V10StoredBlock.parse(
                source.read(entry.blockPosition, entry.storedByteLength), entry);
        cacheBlock(entry, decoded);
        return decoded;
    }

    /**
     * Fetches several stored blocks, keyed by block number.
     * <p>
     * Every block remains independently decompressed, but stored records that
     * are physically adjacent are fetched with one range read. That keeps a
     * broad query from issuing one disk seek or HTTP request per block without
     * ever reading bytes belonging to a block that was not requested.
     */
    public Map<Integer, V10StoredBlock> readBlocks(List<V10BlockIndexEntry> entries) throws IOException {
        Map<Integer, V10StoredBlock> result = new HashMap<>();
        List<V10BlockIndexEntry> missing = new ArrayList<>(entries.size());
        for (V10BlockIndexEntry entry : entries) {
            V10StoredBlock cached = cachedBlock(entry);
            if (cached != null) {
                result.put(entry.blockNumber, cached);
            } else {
                missing.add(entry);
            }
        }
        if (missing.isEmpty()) return result;

        Collections.sort(missing, BY_POSITION);
        int from = 0;
        while (from < missing.size()) {
            long start = missing.get(from).blockPosition;
            long end = V10.add(start, missing.get(from).storedByteLength);
            int to = from + 1;
            while (to < missing.size() && missing.get(to).blockPosition == end
                    && V10.add(end, missing.get(to).storedByteLength) - start <= MAX_COALESCED_READ) {
                end = V10.add(end, missing.get(to).storedByteLength);
                to++;
            }
            byte[] bytes = source.read(start, end - start);
            for (int i = from; i < to; i++) {
                V10BlockIndexEntry entry = missing.get(i);
                V10StoredBlock decoded = V10StoredBlock.parse(bytes,
                        (int) (entry.blockPosition - start), (int) entry.storedByteLength, entry);
                cacheBlock(entry, decoded);
                result.put(entry.blockNumber, decoded);
            }
            from = to;
        }
        return result;
    }

    private V10StoredBlock cachedBlock(V10BlockIndexEntry entry) {
        synchronized (decodedBlockCache) {
            return decodedBlockCache.get(entry.blockPosition);
        }
    }

    private void cacheBlock(V10BlockIndexEntry entry, V10StoredBlock block) {
        synchronized (decodedBlockCache) {
            decodedBlockCache.put(entry.blockPosition, block);
        }
    }

    public void clearBlockCache() {
        synchronized (decodedBlockCache) {
            decodedBlockCache.clear();
        }
    }

    // ------------------------------------------------------------- vectors

    private int normalizationId(NormalizationType type) {
        for (int i = 0; i < header.normalizations.size(); i++) {
            if (header.normalizations.get(i).equalsIgnoreCase(type.getLabel())) return i;
        }
        return -1;
    }

    @Override
    public NormalizationVector readNormalizationVector(NormalizationType type, int chrIdx,
                                                       HiCZoom.HiCUnit unit, int binSize) throws IOException {
        return readNormalizationVector(type, chrIdx, unit, binSize, 0, Long.MAX_VALUE);
    }

    @Override
    public NormalizationVector readNormalizationVectorPart(NormalizationType type, int chrIdx, HiCZoom.HiCUnit unit,
                                                           int binSize, int bound1, int bound2) throws IOException {
        return readNormalizationVector(type, chrIdx, unit, binSize, bound1, (long) bound2 + 1);
    }

    private NormalizationVector readNormalizationVector(NormalizationType type, int chrIdx, HiCZoom.HiCUnit unit,
                                                        int binSize, long begin, long end) throws IOException {
        if (normIndex == null) return null;
        int unitId = unit == HiCZoom.HiCUnit.FRAG ? V10.UNIT_FRAG : V10.UNIT_BP;
        int resolutionIndex = header.resolutionIndex(unitId, binSize);
        if (resolutionIndex < 0) return null;

        boolean useVCForVCSQRT = false;
        int normId = normalizationId(type);
        V10VectorEntry entry = normId < 0 ? null : normIndex.get(normId, chrIdx, unitId, resolutionIndex);
        if (entry == null && type.equals(NormalizationHandler.VC_SQRT)) {
            normId = normalizationId(NormalizationHandler.VC);
            entry = normId < 0 ? null : normIndex.get(normId, chrIdx, unitId, resolutionIndex);
            useVCForVCSQRT = true;
        }
        if (entry == null) return null;

        double[] raw = entry.read(source, begin, end);
        ListOfDoubleArrays values = new ListOfDoubleArrays(raw.length);
        boolean allNaN = true;
        for (int i = 0; i < raw.length; i++) {
            double value = useVCForVCSQRT ? Math.sqrt(raw[i]) : raw[i];
            values.set(i, value);
            if (!Double.isNaN(raw[i])) allNaN = false;
        }
        if (allNaN) return null;
        return new NormalizationVector(type, chrIdx, unit, binSize, values);
    }

    @Override
    public NormalizationVector getNormalizationVector(int chr1Idx, HiCZoom zoom, NormalizationType normalizationType) {
        return dataset.getNormalizationVector(chr1Idx, zoom, normalizationType);
    }

    @Override
    public ListOfDoubleArrays readExpectedVectorPart(long position, long nVals) throws IOException {
        // V10 expected vectors are chunked and addressed through EVI0/NEVI, not by
        // a raw file offset; V10ExpectedValueFunction serves them instead.
        throw new IOException("V10 expected vectors are not addressed by file position; "
                + "use Dataset.getExpectedValues(...)");
    }

    // -------------------------------------------------------------- blocks

    @Override
    public Block readNormalizedBlock(int blockNumber, String zdKey, NormalizationType no, int chr1Index,
                                     int chr2Index, HiCZoom zoom, IndexEntry idx) throws IOException {
        // V10 blocks are located through the resolution's own exact block index
        // and need its descriptor to be decoded, not a V6-V9 (position, size) entry.
        throw new V10FormatException("V10 blocks are addressed through the resolution block index; "
                + "use MatrixZoomData.getNormalizedBlocksOverlapping(...)");
    }

    // --------------------------------------------------------------- misc

    @Override
    public boolean isActive() {
        return activeStatus;
    }

    @Override
    public void setActive(boolean status) {
        activeStatus = status;
    }

    @Override
    public int getVersion() {
        return V10.VERSION;
    }

    @Override
    public int getDepthBase() {
        // The rotated cis grid is normatively log base 2 (Section F.3).
        return 2;
    }

    public Dataset getDataset() {
        return dataset;
    }
}
