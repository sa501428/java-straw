package javastraw.reader;

public class ConstantDepth extends V9Depth {

    public ConstantDepth(int base, int blockBinCount) {
        super(blockBinCount);
        this.BASE = Math.abs(base);
    }

    @Override
    protected int logBase(double v) {
        return (int) Math.abs(v / BASE);
    }
}
