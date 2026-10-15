package javastraw.reader.v10;

import java.io.IOException;
import java.util.List;
import java.util.Map;

import static javastraw.reader.v10.V10.require;

/**
 * One entry of a normalization (NVI0), raw expected (EVI0) or normalized
 * expected (NEVI) index, together with the decoder for its chunked f32 payload
 * (Sections J.2 and J.4-J.6).
 */
public class V10VectorEntry {

    public final int normalizationTypeId;
    public final int chrIndex;
    public final int unit;
    public final int resolutionIndex;
    public final int binSize;
    public final long vectorValueCount;
    public final int nominalChunkValueCount;
    public final List<V10Chunk> chunks;
    /**
     * Chromosome scale factors; empty for normalization vectors. A missing
     * factor means 1.0.
     */
    public final Map<Integer, Float> scaleFactors;

    V10VectorEntry(int normalizationTypeId, int chrIndex, int unit, int resolutionIndex, int binSize,
                   long vectorValueCount, int nominalChunkValueCount, List<V10Chunk> chunks,
                   Map<Integer, Float> scaleFactors) {
        this.normalizationTypeId = normalizationTypeId;
        this.chrIndex = chrIndex;
        this.unit = unit;
        this.resolutionIndex = resolutionIndex;
        this.binSize = binSize;
        this.vectorValueCount = vectorValueCount;
        this.nominalChunkValueCount = nominalChunkValueCount;
        this.chunks = chunks;
        this.scaleFactors = scaleFactors;
    }

    /**
     * Chromosome scale factor for this expected vector; 1.0 when absent.
     */
    public double getScaleFactor(int chromosomeIndex) {
        Float factor = scaleFactors.get(chromosomeIndex);
        return factor == null ? 1.0 : factor;
    }

    /**
     * Decodes values with array indices in [begin, end), fetching only the
     * chunks that actually overlap the request.
     */
    public double[] read(V10Source source, long begin, long end) throws IOException {
        long last = Math.min(end, vectorValueCount);
        require(begin >= 0 && begin <= last && last - begin <= V10.ALLOCATION_LIMIT / 8,
                "vector range exceeds allocation limit");
        double[] result = new double[(int) (last - begin)];
        for (int i = chunkIndex(begin); i < chunks.size(); i++) {
            V10Chunk chunk = chunks.get(i);
            if (chunk.firstValueIndex >= last) break;
            float[] values = source.readVectorChunk(chunk);
            long from = Math.max(begin, chunk.firstValueIndex);
            long to = Math.min(last, chunk.firstValueIndex + chunk.valueCount);
            for (long index = from; index < to; index++) {
                result[(int) (index - begin)] = values[(int) (index - chunk.firstValueIndex)];
            }
        }
        return result;
    }
    /** First chunk whose half-open interval contains or follows the index. */
    public int chunkIndex(long index) {
        int lo = 0, hi = chunks.size();
        while (lo < hi) {
            int mid = (lo + hi) >>> 1;
            V10Chunk chunk = chunks.get(mid);
            if (chunk.firstValueIndex + chunk.valueCount <= index) lo = mid + 1;
            else hi = mid;
        }
        return lo;
    }

    static float[] decodeChunk(V10Source source, V10Chunk chunk) throws IOException {
        V10Cursor c = new V10Cursor(source.read(chunk.filePosition, chunk.storedByteLength));
        c.magic("H10V");
        require(c.byteValue() == chunk.codec && c.byteValue() == chunk.transform,
                "vector chunk codec/transform mismatch");
        c.zero(2);
        require(c.word() == chunk.uncompressedByteLength && c.word() == chunk.valueCount,
                "vector chunk size mismatch");
        byte[] data = V10Zstd.decompress(c, chunk.uncompressedByteLength);
        c.done();
        V10Cursor words = new V10Cursor(data);
        float[] values = new float[chunk.valueCount];
        long previous = 0;
        for (int i = 0; i < values.length; i++) {
            long bits;
            if (chunk.transform == V10.TRANSFORM_BYTE_SHUFFLE) {
                bits = (data[i] & 255L) | ((data[values.length + i] & 255L) << 8)
                        | ((data[2 * values.length + i] & 255L) << 16)
                        | ((data[3 * values.length + i] & 255L) << 24);
            } else {
                bits = words.word();
                if (chunk.transform == V10.TRANSFORM_XOR32) bits ^= previous;
            }
            previous = bits;
            values[i] = V10.bitsToFloat(bits);
        }
        return values;
    }

}
