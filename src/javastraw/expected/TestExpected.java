package javastraw.expected;

import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationHandler;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.HiCFileTools;
import javastraw.tools.MatrixTools;

public class TestExpected {

    public static void main(String[] args) {

        Dataset ds = HiCFileTools.extractDatasetForCLT(
                "/Users/muhammad/Desktop/hicfiles/subsample/chr10_subsample0.25.hic",
                //"",
                false, false, false);
        Chromosome chrom = ds.getChromosomeHandler().getChromosomeFromName("chr10");

        NormalizationType[] norms = new NormalizationType[]{NormalizationHandler.SCALE,
                NormalizationHandler.VC
        };

        for (int res : new int[]{50000, 10000, 5000, 2000, 1000}) {
            MatrixZoomData zd = ds.getMatrix(chrom, chrom).getZoomData(new HiCZoom(res));

            for (NormalizationType norm : norms) {

                int maxBin = (int) (chrom.getLength() / res);
                LogExpectedSpline polynomial2 = new LogExpectedSpline(zd, norm, chrom, res);
                LogExpectedZscoreSpline polyZ = new LogExpectedZscoreSpline(zd, norm, chrom, res);

                System.out.println("Got all points");

                double[][] data = new double[5][maxBin];
                double[] evec = ds.getExpectedValues(new HiCZoom(res), norm, false).getExpectedValuesWithNormalization(chrom.getIndex()).getValues().get(0);
                for (int x = 0; x < maxBin; x++) {
                    data[0][x] = x;
                    data[1][x] = polynomial2.getExpectedFromUncompressedBin(x);
                    data[2][x] = evec[x];
                    data[3][x] = polyZ.getExpectedFromUncompressedBin(x);
                    data[4][x] = polyZ.getExpectedValForZFromUncompressedBin(x, 1);
                }
                MatrixTools.saveMatrixTextNumpy("multi_interp_" + norm.getLabel() + "_" + res + ".npy", data);
            }
        }
    }
}
