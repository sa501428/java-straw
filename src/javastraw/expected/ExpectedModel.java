package javastraw.expected;

import javastraw.reader.block.ContactRecord;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.NormalizationType;

import java.util.Iterator;

public abstract class ExpectedModel {

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

    public static float getP(double obs, double expected, double superDiagonal) {
        // P = (O - E)/(SD - E)
        return (float) ((obs - expected) / (superDiagonal - expected));
    }

    public static double logp1(double x) {
        return Math.log(1 + x);
    }

    public abstract double getExpectedFromUncompressedBin(int dist);

    public abstract double getNearDiagonalSignal();

    public float getPercentContact(ContactRecord cr) {
        return getPercentContact(ExpectedUtils.getDist(cr), cr.getCounts());
    }

    public float getPercentContact(int dist, float counts) {
        double baseline = getExpectedFromUncompressedBin(dist);
        return getP(counts, baseline, getNearDiagonalSignal());
    }

    public double[] expm1(double[] input) {
        double[] vec = new double[input.length];
        for (int k = 0; k < vec.length; k++) {
            vec[k] = Math.expm1(input[k]);
        }
        return vec;
    }

    public int logp1i(int x) {
        return (int) Math.log(1 + x);
    }

    public boolean isReasonablePercentContact(ContactRecord cr, float minPC, float maxPC) {
        float percentContact = getPercentContact(cr);
        return percentContact > minPC && percentContact < maxPC;
    }

    public boolean isReasonableEnrichment(ContactRecord cr, float enrichmentVSExpected) {
        return cr.getCounts() / getExpectedFromUncompressedBin(ExpectedUtils.getDist(cr)) > enrichmentVSExpected;
    }

    public void print(int maxBin) {
        for (int k = 0; k < maxBin; k++) {
            System.out.println(k + "  " + getExpectedFromUncompressedBin(k));
        }
    }
}
