package javastraw.tools;

import javastraw.reader.Dataset;
import javastraw.reader.Matrix;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.block.Block;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.datastructures.ListOfDoubleArrays;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.norm.NormalizationVector;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationType;

import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

public class IterationUtils {
    public static void processDataset(String path, int resolution, int id1) {
        Dataset ds1 = HiCFileTools.extractDatasetForCLT(path, true, false, false);
        Chromosome[] array = ds1.getChromosomeHandler().getChromosomeArrayWithoutAllByAll();
        List<NormalizationType> norms = ds1.getNormalizationTypes();
        HiCZoom zoom = ds1.getZoomForBPResolution(resolution);

        for (Chromosome chrom : array) {
            for (NormalizationType norm : norms) {
                if (norm.getLabel().toLowerCase().contains("vc")) continue;
                NormalizationVector nv1 = ds1.getNormalizationVector(chrom.getIndex(), zoom, norm);
                if (nv1 == null) continue;
                System.out.println(chrom.getName() + "_" + norm.getLabel());
                double[] vals;
                if (id1 == 0) {
                    vals = getRowSums(ds1, chrom, norm, zoom, nv1.getData()); //norm
                } else {
                    vals = getRowSumsOld(ds1, chrom, norm, zoom, nv1.getData()); //norm
                }

                System.out.println(Arrays.toString(vals));
            }
        }
    }

    // 4 seconds
    public static double[] getRowSums(Dataset ds1, Chromosome chrom, NormalizationType norm,
                                      HiCZoom zoom, ListOfDoubleArrays nv1) {
        Matrix matrix = ds1.getMatrix(chrom, chrom);
        if (matrix == null) return null;
        MatrixZoomData zd = matrix.getZoomData(zoom);
        if (zd == null) return null;

        double[] vals = new double[3];
        Iterator<ContactRecord> it = zd.getNormalizedIterator(norm);
        while (it.hasNext()) {
            updateCounters(it.next(), vals);
        }
        return vals;
    }

    public static void updateCounters(ContactRecord cr, double[] vals) {
        vals[0]++;
        if (cr.getCounts() > 0) {
            vals[1]++;
            vals[2] += cr.getCounts();
            if (cr.getBinX() != cr.getBinY()) {
                vals[2] += cr.getCounts();
            }
        }
    }

    public static double[] getRowSumsOld(Dataset ds1, Chromosome chrom, NormalizationType norm,
                                         HiCZoom zoom, ListOfDoubleArrays nv1) {
        Matrix matrix = ds1.getMatrix(chrom, chrom);
        if (matrix == null) return null;
        MatrixZoomData zd = matrix.getZoomData(zoom);
        if (zd == null) return null;

        long maxBin = nv1.getLength();
        double[] vals = new double[3];

        List<Block> blocks = zd.getNormalizedBlocksOverlapping(0, 0, maxBin, maxBin,
                norm, false);
        for (Block block : blocks) {
            for (ContactRecord cr : block.getContactRecords()) {
                updateCounters(cr, vals);
            }
        }
        blocks.clear();

        return vals;
    }
}
