package javastraw.expected;

public class Zscore {

    private final double mean;
    private final double stdDev;

    public Zscore(double mean, double stdDev) {
        this.mean = mean;
        this.stdDev = stdDev;
    }

    public float getZscore(double value) {
        return (float) ((value - mean) / stdDev);
    }
}
