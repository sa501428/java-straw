package javastraw.reader.v10;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import static javastraw.reader.v10.V10.require;

/**
 * A parsed normalization (NVI0), raw expected (EVI0) or normalized expected
 * (NEVI) index (Section J.3).
 * <p>
 * Index presence is the capability advertisement: a missing entry means the
 * file does not offer that normalization or expected vector at that resolution,
 * and a reader must never infer it from another resolution or normalization.
 */
public class V10VectorIndex {

    public static final int KIND_NORM = 0;
    public static final int KIND_EXPECTED = 1;
    public static final int KIND_NORM_EXPECTED = 2;

    private final int kind;
    private final Map<String, V10VectorEntry> entries = new LinkedHashMap<>();
    private final List<V10VectorEntry> ordered = new ArrayList<>();

    private V10VectorIndex(int kind) {
        this.kind = kind;
    }

    public int getKind() {
        return kind;
    }

    public List<V10VectorEntry> getEntries() {
        return ordered;
    }

    public V10VectorEntry get(int normalizationTypeId, int chrIndex, int unit, int resolutionIndex) {
        return entries.get(key(normalizationTypeId, chrIndex, unit, resolutionIndex));
    }

    private static String key(int normalizationTypeId, int chrIndex, int unit, int resolutionIndex) {
        return normalizationTypeId + "_" + chrIndex + "_" + unit + "_" + resolutionIndex;
    }

    /**
     * Reads and fully validates one vector index.
     */
    public static V10VectorIndex parse(V10Source source, V10Locator locator, int kind, V10Header header)
            throws IOException {
        V10VectorIndex index = new V10VectorIndex(kind);
        if (!locator.isPresent()) return index;

        byte[] bytes = source.read(locator.position, locator.length);
        V10Cursor c = new V10Cursor(bytes);
        c.magic(kind == KIND_NORM ? "NVI0" : kind == KIND_EXPECTED ? "EVI0" : "NEVI");
        require(c.word() == 1, "unknown vector index version");
        long entryCount = c.word();
        c.zero(4);
        require(entryCount <= c.left() / 40, "invalid vector entry count");

        long[] previousKey = null;
        for (long i = 0; i < entryCount; i++) {
            long entryLength = c.word();
            require(entryLength >= 40, "invalid vector entry size");
            V10Cursor e = c.take(entryLength - 4);

            int normId = kind == KIND_EXPECTED ? 0 : (int) e.word();
            int chrIndex = kind == KIND_NORM ? (int) e.word() : 0;
            int unit = e.byteValue();
            e.zero(3);
            int resolutionIndex = (int) e.word();
            long binSize = e.word();
            require(unit <= 1 && resolutionIndex >= 0 && resolutionIndex < header.resolutions[unit].size()
                            && binSize == header.resolutions[unit].get(resolutionIndex).binSize
                            && (kind == KIND_EXPECTED || (normId >= 0 && normId < header.normalizations.size()))
                            && (kind != KIND_NORM || (chrIndex >= 0 && chrIndex < header.chromosomes.size())),
                    "invalid vector key");

            long[] currentKey = {normId, chrIndex, unit, resolutionIndex};
            require(previousKey == null || compare(previousKey, currentKey) < 0,
                    "unsorted/duplicate vector key");
            previousKey = currentKey;

            long valueCount = e.wide();
            long nominal = e.word();
            long chunkCount = e.word();
            long required;
            if (kind == KIND_NORM) {
                required = header.bins(chrIndex, unit, resolutionIndex);
            } else {
                // Raw and normalized expected vectors cover every possible cis
                // distance, so the longest chromosome sets the length.
                required = 0;
                for (int chr = 0; chr < header.chromosomes.size(); chr++) {
                    required = Math.max(required, header.bins(chr, unit, resolutionIndex));
                }
            }
            require(valueCount == required && (valueCount == 0 || (nominal > 0 && chunkCount > 0))
                            && (valueCount > 0 || chunkCount == 0),
                    "invalid vector length/chunks");

            Map<Integer, Float> scaleFactors = new HashMap<>();
            if (kind != KIND_NORM) {
                long scaleFactorCount = e.word();
                e.zero(4);
                require(scaleFactorCount <= e.left() / 8, "invalid scale factor count");
                long previousChr = -1;
                for (long j = 0; j < scaleFactorCount; j++) {
                    long chr = e.word();
                    long bits = e.word();
                    require(chr < header.chromosomes.size() && (j == 0 || chr > previousChr),
                            "invalid scale factor key");
                    scaleFactors.put((int) chr, V10.bitsToFloat(bits));
                    previousChr = chr;
                }
            }
            require(chunkCount * 32 == e.left(), "vector descriptor length mismatch");

            List<V10Chunk> chunks = new ArrayList<>((int) chunkCount);
            long next = 0;
            for (long j = 0; j < chunkCount; j++) {
                long first = e.wide();
                long count = e.word();
                int transform = e.byteValue();
                int codec = e.byteValue();
                e.zero(2);
                long position = e.wide();
                long stored = e.word();
                long raw = e.word();
                require(first == next && count > 0 && count * 4 == raw && transform <= 2
                                && codec == V10.CODEC_ZSTD && stored > 16,
                        "invalid vector chunk descriptor");
                next = V10.add(next, count);
                require(next <= valueCount, "vector chunk exceeds length");
                source.interval(position, stored);
                chunks.add(new V10Chunk(first, (int) count, transform, codec, position, stored, (int) raw));
            }
            require(next == valueCount, "incomplete vector coverage");
            e.done();

            V10VectorEntry entry = new V10VectorEntry(normId, chrIndex, unit, resolutionIndex,
                    (int) binSize, valueCount, (int) nominal, chunks, scaleFactors);
            index.entries.put(key(normId, chrIndex, unit, resolutionIndex), entry);
            index.ordered.add(entry);
        }
        c.done();
        return index;
    }

    private static int compare(long[] a, long[] b) {
        for (int i = 0; i < a.length; i++) {
            if (a[i] != b[i]) return a[i] < b[i] ? -1 : 1;
        }
        return 0;
    }
}
