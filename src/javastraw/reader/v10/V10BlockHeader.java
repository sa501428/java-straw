package javastraw.reader.v10;

/**
 * The fixed 40-byte header of a V10 logical block (Section I.1).
 * <p>
 * A stored block carries explicit genomic bin offsets, width and height, which
 * readers must use when reconstructing cells: for the rotated cis grid the
 * block's region is rotated in genomic coordinates, so its extent cannot be
 * inferred from the block number alone.
 */
public class V10BlockHeader {

    public final int representation;
    public final int valueMode;
    public final int valueType;
    public final int flags;
    public final long binColumnOffset;
    public final long binRowOffset;
    public final long blockWidth;
    public final long blockHeight;
    public final long occupiedCellCount;
    public final long positionStreamBytes;
    public final long valueStreamBytes;

    V10BlockHeader(int representation, int valueMode, int valueType, int flags,
                   long binColumnOffset, long binRowOffset, long blockWidth, long blockHeight,
                   long occupiedCellCount, long positionStreamBytes, long valueStreamBytes) {
        this.representation = representation;
        this.valueMode = valueMode;
        this.valueType = valueType;
        this.flags = flags;
        this.binColumnOffset = binColumnOffset;
        this.binRowOffset = binRowOffset;
        this.blockWidth = blockWidth;
        this.blockHeight = blockHeight;
        this.occupiedCellCount = occupiedCellCount;
        this.positionStreamBytes = positionStreamBytes;
        this.valueStreamBytes = valueStreamBytes;
    }

    public long lastBinColumn() {
        return binColumnOffset + blockWidth - 1;
    }

    public long lastBinRow() {
        return binRowOffset + blockHeight - 1;
    }
}
