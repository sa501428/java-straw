package javastraw.expected;

import javastraw.feature2D.Feature2D;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.NormalizationType;

import java.util.Iterator;

public class ExpectedUtils {
    public static double[] calculateExpected(MatrixZoomData zd, NormalizationType norm, int maxBinDist,
                                             boolean useLog) {
        double[] expected = new double[maxBinDist];
        long[] counts = new long[maxBinDist];

        Iterator<ContactRecord> iterator = zd.getNormalizedIterator(norm);
        updateExpectedAndCounts(maxBinDist, useLog, expected, counts, iterator);

        normalizeByCounts(expected, counts);

        if (useLog) {
            expm1(expected);
        }

        return expected;
    }

    private static void normalizeByCounts(double[] vector, long[] counts) {
        for (int z = 0; z < vector.length; z++) {
            if (counts[z] > 0) {
                vector[z] /= counts[z];
            }
        }
    }

    private static void expm1(double[] vector) {
        for (int z = 0; z < vector.length; z++) {
            vector[z] = Math.expm1(vector[z]);
        }
    }

    private static void updateExpectedAndCounts(int maxBinDist, boolean useLog, double[] expected, long[] counts, Iterator<ContactRecord> iterator) {
        while (iterator.hasNext()) {
            ContactRecord record = iterator.next();
            int dist = getDist(record);
            if (dist < maxBinDist) {
                if (useLog) {
                    expected[dist] += Math.log(1 + record.getCounts());
                } else {
                    expected[dist] += record.getCounts();
                }
                counts[dist]++;
            }
        }
    }

    public static int getDist(Feature2D loop, int resolution) {
        int binXStart = (int) (loop.getMidPt1() / resolution);
        int binYStart = (int) (loop.getMidPt2() / resolution);
        return Math.abs(binXStart - binYStart);
    }

    public static int getDist(ContactRecord record) {
        return Math.abs(record.getBinX() - record.getBinY());
    }

    public static Iterator<ContactRecord> getIterator(MatrixZoomData zd, NormalizationType norm) {
        if (norm.getLabel().equalsIgnoreCase("none")) {
            return zd.getDirectIterator();
        } else {
            return zd.getNormalizedIterator(norm);
        }
    }
}
