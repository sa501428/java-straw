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

public class Tester {

    public static void main(String[] args) {
        String path = "/Desktop/chr10.hic";
        processDataset(path, 50);
    }

    private static void processDataset(String path, int resolution) {
        Dataset ds1 = HiCFileTools.extractDatasetForCLT(path, true, false);
        Chromosome[] array = ds1.getChromosomeHandler().getChromosomeArrayWithoutAllByAll();
        List<NormalizationType> norms = ds1.getNormalizationTypes();
        HiCZoom zoom = ds1.getZoomForBPResolution(resolution);
        //NormalizationType none = ds1.getNormalizationHandler().getNormTypeFromString("NONE");

        for (Chromosome chrom : array) {
            for (NormalizationType norm : norms) {
                if (norm.getLabel().toLowerCase().contains("vc")) continue;
                NormalizationVector nv1 = ds1.getNormalizationVector(chrom.getIndex(), zoom, norm);
                if (nv1 == null) continue;
                System.out.println(chrom.getName() + "_" + norm.getLabel());
                long[] vals = getRowSums(ds1, chrom, norm, zoom, nv1.getData()); //norm
                System.out.println(Arrays.toString(vals));
            }
        }
    }


    public static long[] getRowSums(Dataset ds1, Chromosome chrom, NormalizationType norm,
                                    HiCZoom zoom, ListOfDoubleArrays nv1) {
        Matrix matrix = ds1.getMatrix(chrom, chrom);
        if (matrix == null) return null;
        System.out.println("M");
        MatrixZoomData zd = matrix.getZoomData(zoom);
        if (zd == null) return null;
        System.out.println("Z");

        long maxBin = nv1.getLength();
        System.out.println("Max " + maxBin);
        long[] vals = new long[2];

        Iterator<ContactRecord> it = zd.getNormalizedIterator(norm);
        while (it.hasNext()) {
            ContactRecord cr = it.next();
            vals[0]++;
            if (cr.getCounts() > 0) {
                vals[1] += cr.getCounts();
                if (cr.getBinX() != cr.getBinY()) {
                    vals[1] += cr.getCounts();
                }
            }
        }

        return vals;
    }

    public static long[] getRowSumsOld(Dataset ds1, Chromosome chrom, NormalizationType norm,
                                       HiCZoom zoom, ListOfDoubleArrays nv1) {
        Matrix matrix = ds1.getMatrix(chrom, chrom);
        if (matrix == null) return null;
        System.out.println("M");
        MatrixZoomData zd = matrix.getZoomData(zoom);
        if (zd == null) return null;
        System.out.println("Z");

        long maxBin = nv1.getLength();
        System.out.println("Max " + maxBin);
        long[] vals = new long[2];

        List<Block> blocks = zd.getNormalizedBlocksOverlapping(0, 0, maxBin, maxBin,
                norm, false, false);
        System.out.println("Num boxes " + blocks.size());
        for (Block block : blocks) {
            System.out.println("B# " + block.getNumber());
            vals[0] += block.getContactRecords().size();
            for (ContactRecord cr : block.getContactRecords()) {
                if (cr.getCounts() > 0) {
                    vals[1] += cr.getCounts();
                    if (cr.getBinX() != cr.getBinY()) {
                        vals[1] += cr.getCounts();
                    }
                }
            }
        }
        blocks.clear();

        return vals;
    }
}
