package javastraw.expected;

public class Zscore {

    private final double mean;
    private final double stdDev;

    public Zscore(double mean, double stdDev) {
        this.mean = mean;
        this.stdDev = stdDev;
    }

    public double getZscore(double value) {
        return (value - mean) / stdDev;
    }

    public double getValForZscore(double z) {
        return mean + (z * stdDev);
    }
}
