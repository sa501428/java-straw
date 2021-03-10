package javastraw.reader.depth;

public class LogDepth extends V9Depth {

    public LogDepth(int base, int blockBinCount) {
        super(blockBinCount);
        this.BASE = Math.log(base);
    }

    @Override
    protected int logBase(double v) {
        return (int) (Math.log(1 + v) / BASE);
    }
}
