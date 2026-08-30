package javastraw.reader.v10;

/**
 * One independently addressable, Zstandard-compressed run of f32 values from a
 * normalization or expected-value vector (Section J.1).
 */
public class V10Chunk {

    public final long firstValueIndex;
    public final int valueCount;
    public final int transform;
    public final int codec;
    public final long filePosition;
    public final long storedByteLength;
    public final int uncompressedByteLength;

    V10Chunk(long firstValueIndex, int valueCount, int transform, int codec,
             long filePosition, long storedByteLength, int uncompressedByteLength) {
        this.firstValueIndex = firstValueIndex;
        this.valueCount = valueCount;
        this.transform = transform;
        this.codec = codec;
        this.filePosition = filePosition;
        this.storedByteLength = storedByteLength;
        this.uncompressedByteLength = uncompressedByteLength;
    }
}
