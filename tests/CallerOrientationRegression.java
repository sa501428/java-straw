import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.block.Block;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.mzd.Matrix;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.NormalizationHandler;
import javastraw.tools.HiCFileTools;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

public final class CallerOrientationRegression {
    public static void main(String[] args) {
        if (args.length != 1) throw new IllegalArgumentException("usage: test FILE.hic");
        Dataset dataset = HiCFileTools.extractDatasetForCLT(args[0], false, true, false);
        Chromosome[] chromosomes = dataset.getChromosomeHandler().getChromosomeArrayWithoutAllByAll();
        if (chromosomes.length < 2) throw new AssertionError("fixture needs two chromosomes");
        Chromosome first = chromosomes[0];
        Chromosome second = chromosomes[1];

        Matrix forward = dataset.getMatrix(first, second);
        Matrix reverse = dataset.getMatrix(second, first);
        if (forward == null || reverse == null) throw new AssertionError("missing inter matrix");
        MatrixZoomData a = forward.getFirstZoomData();
        MatrixZoomData b = reverse.getFirstZoomData();
        if (a.getChr1Idx() != first.getIndex() || a.getChr2Idx() != second.getIndex() ||
            b.getChr1Idx() != second.getIndex() || b.getChr2Idx() != first.getIndex()) {
            throw new AssertionError("matrix axes do not follow caller order");
        }

        long aX = first.getLength() / a.getBinSize() + 1;
        long aY = second.getLength() / a.getBinSize() + 1;
        Map<String, Float> expected = collect(a.getNormalizedBlocksOverlapping(
            0, 0, aX, aY, NormalizationHandler.NONE, false), false);
        Map<String, Float> actual = collect(b.getNormalizedBlocksOverlapping(
            0, 0, aY, aX, NormalizationHandler.NONE, false), true);
        if (expected.isEmpty() || !expected.equals(actual)) {
            throw new AssertionError("reverse query did not transpose every contact");
        }
    }

    private static Map<String, Float> collect(List<Block> blocks, boolean transpose) {
        Map<String, Float> records = new HashMap<>();
        for (Block block : blocks) for (ContactRecord record : block.getContactRecords()) {
            int x = transpose ? record.getBinY() : record.getBinX();
            int y = transpose ? record.getBinX() : record.getBinY();
            records.put(x + ":" + y, record.getCounts());
        }
        return records;
    }
}
