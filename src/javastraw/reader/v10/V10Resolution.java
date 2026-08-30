package javastraw.reader.v10;

/**
 * A logical-resolution record from the variable header (Section C.3).
 */
public class V10Resolution {

    public final int binSize;
    public final int storageMode;
    public final int aggregation;
    /**
     * Index into the same unit's resolution list, or {@link V10#NO_SOURCE_RESOLUTION}.
     */
    public final long sourceResolutionIndex;

    V10Resolution(int binSize, int storageMode, int aggregation, long sourceResolutionIndex) {
        this.binSize = binSize;
        this.storageMode = storageMode;
        this.aggregation = aggregation;
        this.sourceResolutionIndex = sourceResolutionIndex;
    }

    public boolean isDerived() {
        return storageMode == V10.DERIVED;
    }
}
