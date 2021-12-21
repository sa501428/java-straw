package javastraw.reader.pearsons;

import javastraw.matrices.BasicMatrix;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.expected.ExpectedValueFunction;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import java.util.Arrays;
import java.util.BitSet;
import java.util.Iterator;

public class PearsonsManager {

    /**
     * Computes eigenvector from Pearson's.
     *
     * @return Eigenvector
     */
    public static double[] computeEigenvector(BasicMatrix pearsons, int which) {

        int dim = pearsons.getRowDimension();
        BitSet bitSet = new BitSet(dim);
        for (int i = 0; i < dim; i++) {
            float tmp = pearsons.getEntry(i, i);
            if (tmp > .999) { // checking if diagonal is 1
                bitSet.set(i);
            }
        }

        int[] newPosToOrig = getMapNewPosToOriginal(dim, bitSet);

        RealMatrix subMatrix = getSubsetOfMatrix(newPosToOrig, bitSet.cardinality(), pearsons);
        RealVector rv = (new EigenDecomposition(subMatrix)).getEigenvector(which);
        double[] ev = rv.toArray();

        int size = pearsons.getColumnDimension();
        double[] eigenvector = new double[size];
        Arrays.fill(eigenvector, Double.NaN);
        for (int i = 0; i < newPosToOrig.length; i++) {
            int oldI = newPosToOrig[i];
            eigenvector[oldI] = ev[i];
        }
        return eigenvector;
    }

    private static int[] getMapNewPosToOriginal(int dim, BitSet bitSet) {
        int[] newPosToOrig = new int[bitSet.cardinality()];
        int count = 0;
        for (int i = 0; i < dim; i++) {
            if (bitSet.get(i)) {
                newPosToOrig[count++] = i;
            }
        }
        return newPosToOrig;
    }

    private static RealMatrix getSubsetOfMatrix(int[] newPosToOrig, int subsetN, BasicMatrix pearsons) {
        double[][] data = new double[subsetN][subsetN];

        for (int newI = 0; newI < newPosToOrig.length; newI++) {
            int oldI = newPosToOrig[newI];
            for (int newJ = newI; newJ < newPosToOrig.length; newJ++) {
                int oldJ = newPosToOrig[newJ];
                float tmp = pearsons.getEntry(oldI, oldJ);
                data[newI][newJ] = tmp;
                data[newJ][newI] = tmp;
            }
        }
        return new Array2DRowRealMatrix(data);
    }

    public static BasicMatrix computePearsons(ExpectedValueFunction df, Iterator<ContactRecord> iterator,
                                              Chromosome chr1, int binSize) {

        int dim = (int) (chr1.getLength() / binSize) + 1;

        // Compute O/E column vectors
        double[][] oeMatrix = new double[dim][dim];
        BitSet bitSet = new BitSet(dim);
        populateOEMatrixAndBitset(oeMatrix, bitSet, df, iterator, chr1.getIndex());
        return PearsonsUtils.computeParallelizedPearsons(oeMatrix, dim, bitSet);
    }

    private static void populateOEMatrixAndBitset(double[][] oeMatrix, BitSet bitSet, ExpectedValueFunction df,
                                                  Iterator<ContactRecord> iterator, int chr1Index) {
        while (iterator.hasNext()) {
            ContactRecord record = iterator.next();
            float counts = record.getCounts();
            if (Float.isNaN(counts)) continue;

            int i = record.getBinX();
            int j = record.getBinY();
            int dist = Math.abs(i - j);

            double expected = df.getExpectedValue(chr1Index, dist);
            double oeValue = counts / expected;

            oeMatrix[i][j] = oeValue;
            oeMatrix[j][i] = oeValue;

            bitSet.set(i);
            bitSet.set(j);
        }
    }
}
