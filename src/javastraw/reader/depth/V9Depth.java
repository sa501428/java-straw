package javastraw.reader.depth;

public abstract class V9Depth {
    protected final int blockBinCount;
    protected double BASE;

    V9Depth(int blockBinCount) {
        this.blockBinCount = blockBinCount;
    }

    public static V9Depth setDepthMethod(int depthBase, int blockBinCount) {
        if (depthBase > 1) {
            return new LogDepth(depthBase, blockBinCount);
        } else if (depthBase < 0) {
            return new ConstantDepth(-depthBase, blockBinCount);
        }

        // Default
        return new LogDepth(2, blockBinCount);
    }

    public int getDepth(int val1, int val2) {
        return logBase(Math.abs(val1 - val2) / Math.sqrt(2) / blockBinCount);
    }

    public int getDepth(long val1, long val2) {
        return logBase(Math.abs(val1 - val2) / Math.sqrt(2) / blockBinCount);
    }

    protected abstract int logBase(double value);
}
