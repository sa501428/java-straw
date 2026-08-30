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
        for (V10Chunk chunk : chunks) {
            long chunkEnd = chunk.firstValueIndex + chunk.valueCount;
            if (chunkEnd <= begin || chunk.firstValueIndex >= last) continue;

            byte[] stored = source.read(chunk.filePosition, chunk.storedByteLength);
            V10Cursor c = new V10Cursor(stored);
            c.magic("H10V");
            require(c.byteValue() == chunk.codec && c.byteValue() == chunk.transform,
                    "vector chunk codec/transform mismatch");
            c.zero(2);
            require(c.word() == chunk.uncompressedByteLength && c.word() == chunk.valueCount,
                    "vector chunk size mismatch");
            byte[] data = V10Zstd.decompress(c, chunk.uncompressedByteLength);
            c.done();

            V10Cursor values = new V10Cursor(data);
            long previousBits = 0;
            for (int k = 0; k < chunk.valueCount; k++) {
                long bits;
                if (chunk.transform == V10.TRANSFORM_BYTE_SHUFFLE) {
                    bits = 0;
                    for (int lane = 0; lane < 4; lane++) {
                        bits |= (long) (data[lane * chunk.valueCount + k] & 0xFF) << (8 * lane);
                    }
                } else {
                    bits = values.word();
                    if (chunk.transform == V10.TRANSFORM_XOR32 && k > 0) {
                        bits ^= previousBits;
                    }
                }
                previousBits = bits;
                long index = chunk.firstValueIndex + k;
                if (index >= begin && index < last) {
                    result[(int) (index - begin)] = V10.bitsToFloat(bits & 0xFFFFFFFFL);
                }
            }
        }
        return result;
    }
}
