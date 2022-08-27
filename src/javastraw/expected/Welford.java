package javastraw.expected;

public class Welford {
    private long counts = 0;
    private double mu = 0;
    private double aggSquaredDiffs = 0;
    private double min = Double.MAX_VALUE;
    private double max = 0;

    public void addValue(double x) {
        counts++;
        double nextMu = mu + ((x - mu) / counts);
        aggSquaredDiffs += (x - mu) * (x - nextMu);
        mu = nextMu;
        if (x < min) min = x;
        if (x > max) max = x;
    }

    public double getMean() {
        return mu;
    }

    public double getStdDev() {
        if (counts > 2) {
            return Math.sqrt(aggSquaredDiffs / (counts - 1));
        }
        return 0;
    }

    public double getRange() {
        return max - min;
    }

    public Zscore getZscore() {
        return new Zscore(getMean(), getStdDev());
    }

    public long getCounts() {
        return counts;
    }

    public void addZeroIfBelow(int length) {
        int numToAdd = (int) (length - counts);
        for (int z = 0; z < numToAdd; z++) {
            addValue(0);
        }
    }

    public String getSummary() {
        return " mu " + getMean() + " sigma " + getStdDev() + " min " + min + " max " + max;
    }

    public double getMin() {
        return min;
    }

    public double getMax() {
        return max;
    }
}



