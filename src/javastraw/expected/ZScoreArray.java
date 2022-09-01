package javastraw.expected;

import java.util.Arrays;

public class ZScoreArray {

    private final double[] mean;
    private final double[] stdDev;

    public ZScoreArray(double[] mean, double[] stdDev) {
        this.mean = mean;
        this.stdDev = stdDev;
    }

    public double getZscore(int index, double value) {
        return (value - mean[index]) / stdDev[index];
    }

    public double getValForZscore(int index, double z) {
        return mean[index] + (z * stdDev[index]);
    }

    public void print() {
        System.out.println(Arrays.toString(mean));
        System.out.println(Arrays.toString(stdDev));
    }
}
