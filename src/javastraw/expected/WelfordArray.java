package javastraw.expected;

public class WelfordArray {
    private final long[] counts;
    private final double[] mu;
    private final double[] aggSquaredDiffs;

    public WelfordArray(int n) {
        counts = new long[n];
        mu = new double[n];
        aggSquaredDiffs = new double[n];
    }

    public void addValue(int i, double x) {
        counts[i]++;
        double nextMu = mu[i] + ((x - mu[i]) / counts[i]);
        aggSquaredDiffs[i] += (x - mu[i]) * (x - nextMu);
        mu[i] = nextMu;
    }

    public double[] getMean() {
        return mu;
    }

    public double[] getStdDev() {
        double[] std = new double[counts.length];
        for (int i = 0; i < counts.length; i++) {
            if (counts[i] > 2) {
                std[i] = Math.sqrt(aggSquaredDiffs[i] / (counts[i] - 1));
            }
        }
        return std;
    }

    public ZScoreArray getZscores() {
        return new ZScoreArray(getMean(), getStdDev());
    }
}



