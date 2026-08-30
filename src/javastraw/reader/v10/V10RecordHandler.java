package javastraw.reader.v10;

/**
 * Receives decoded cells of a V10 logical block.
 * <p>
 * Counts are delivered as exact unsigned integers so that derived-resolution
 * aggregation never rounds; scores are delivered as their exact f32 value.
 */
public interface V10RecordHandler {

    /**
     * @param binColumn global bin column (chromosome 1 of the canonical pair)
     * @param binRow    global bin row (chromosome 2 of the canonical pair)
     * @param count     exact contact count, meaningful when {@code isScore} is false
     * @param score     exact stored score, meaningful when {@code isScore} is true
     * @param isScore   true for SCORE_FLOAT32 matrices, false for COUNT_UINT
     */
    void record(int binColumn, int binRow, long count, float score, boolean isScore);
}
