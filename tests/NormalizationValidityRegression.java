import javastraw.reader.*;
import javastraw.reader.block.*;
import javastraw.reader.mzd.*;
import javastraw.reader.norm.*;
import javastraw.reader.type.*;
import javastraw.reader.datastructures.ListOfDoubleArrays;
import java.util.*;

public class NormalizationValidityRegression {
    public static void main(String[] args) throws Exception {
        DatasetReaderV10 reader = new DatasetReaderV10(args[0], true) {
            @Override public NormalizationVector getNormalizationVector(int chr, HiCZoom zoom, NormalizationType type) {
                ListOfDoubleArrays values = new ListOfDoubleArrays(8, 2);
                values.set(0, 0);
                values.set(1, Double.POSITIVE_INFINITY);
                return new NormalizationVector(type, chr, zoom.getUnit(), zoom.getBinSize(), values);
            }
        };
        Dataset ds = reader.read();
        MatrixZoomData zd = ds.getMatrix(ds.getChromosomeHandler().getChromosomeFromIndex(0), ds.getChromosomeHandler().getChromosomeFromIndex(0)).getZoomData(new HiCZoom(10));
        int cells = 0;
        for (Block block : zd.getNormalizedBlocksOverlapping(0, 0, 7, 7, NormalizationHandler.VC, false)) {
            for (ContactRecord record : block.getContactRecords()) {
                if (record.getCounts() != 1.25f) throw new AssertionError("invalid normalized record");
                cells++;
            }
        }
        if (cells != 1) throw new AssertionError("zero and infinite vectors were not excluded");
        cells = 0;
        for (Iterator<ContactRecord> it = zd.getNormalizedIterator(NormalizationHandler.VC); it.hasNext();) {
            if (it.next().getCounts() != 1.25f) throw new AssertionError("invalid iterator record");
            cells++;
        }
        if (cells != 1) throw new AssertionError("iterator did not exclude invalid vectors");
        System.out.println("Normalization validity: blocks and iterators skip zero and infinite denominators");
    }
}
