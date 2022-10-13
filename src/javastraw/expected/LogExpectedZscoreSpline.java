package javastraw.expected;

import javastraw.reader.basics.Chromosome;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.NormalizationType;
import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;

import java.util.Iterator;

public class LogExpectedZscoreSpline extends ExpectedModel {

    private final PolynomialSplineFunction mean, std;
    private final double nearDiagonalSignal;
    private final float maxX;

    public LogExpectedZscoreSpline(MatrixZoomData zd, NormalizationType norm, Chromosome chrom, int res) {
        int maxBin = (int) (chrom.getLength() / res + 1);
        int[] maxDist = new int[1];
        PolynomialSplineFunction[] functions = fitDataToFunction(zd, norm, maxBin, maxDist);
        mean = functions[0];
        std = functions[1];
        maxX = logp1i(Math.min(maxBin, maxDist[0]));
        nearDiagonalSignal = Math.expm1(mean.value(0.5));
    }

    private PolynomialSplineFunction[] fitDataToFunction(MatrixZoomData zd, NormalizationType norm, int maxBin, int[] maxDist) {

        WelfordArray values = new WelfordArray(logp1i(maxBin) + 2);
        WelfordArray distances = new WelfordArray(logp1i(maxBin) + 2);
        populateWithCounts(zd, norm, values, distances, maxBin, maxDist);

        double[] avgExpected = values.getMean();
        double[] avgExpectedStds = values.getStdDev();
        double[] avgDistances = distances.getMean();
        maxDist[0] = (int) Math.expm1(avgDistances[avgDistances.length - 1]);

        return new PolynomialSplineFunction[]{new SplineInterpolator().interpolate(avgDistances, avgExpected),
                new SplineInterpolator().interpolate(avgDistances, avgExpectedStds)};
    }

    private void populateWithCounts(MatrixZoomData zd, NormalizationType norm, WelfordArray values,
                                    WelfordArray distances, int maxBin, int[] maxDist) {
        Iterator<ContactRecord> records = getIterator(zd, norm);
        while (records.hasNext()) {
            ContactRecord record = records.next();
            int dist = getDist(record);
            if (dist == 0) {
                values.addValue(0, logp1(record.getCounts()));
                distances.addValue(0, logp1(dist));
            } else if (dist < maxBin) {
                maxDist[0] = Math.max(maxDist[0], dist);
                values.addValue(logp1i(dist) + 1, logp1(record.getCounts()));
                distances.addValue(logp1i(dist) + 1, logp1(dist));
            }
        }
    }

    @Override
    public double getExpectedFromUncompressedBin(int dist0) {
        double dist = Math.min(Math.max(0, logp1(dist0)), maxX);
        return Math.expm1(mean.value(dist));
    }

    public double getExpectedValForZFromUncompressedBin(int dist0, float z) {
        double dist = Math.min(Math.max(0, logp1(dist0)), maxX);
        return Math.expm1(mean.value(dist) + z * std.value(dist));
    }

    public double getZscoreForObservedUncompressedBin(int dist0, float observed) {
        double dist = Math.min(Math.max(0, logp1(dist0)), maxX);
        return (logp1(observed) - mean.value(dist)) / std.value(dist);
    }

    @Override
    public double getNearDiagonalSignal() {
        return nearDiagonalSignal;
    }
}
