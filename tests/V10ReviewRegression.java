import javastraw.reader.*;
import javastraw.reader.mzd.*;
import javastraw.reader.block.*;
import javastraw.reader.type.*;
import java.util.*;

/** Cross-block score ordering and bounded normalization at fine resolutions. */
public class V10ReviewRegression {
    public static void main(String[] args) throws Exception {
        DatasetReaderV10 reader = new DatasetReaderV10(args[0] + "/score-order-v10.hic", true);
        Dataset ds = reader.read();
        MatrixZoomData zoom = ds.getMatrix(ds.getChromosomeHandler().getChromosomeFromIndex(0),
                ds.getChromosomeHandler().getChromosomeFromIndex(1)).getZoomData(new HiCZoom(10000000));
        int cells = 0;
        for (Block block : zoom.getNormalizedBlocksOverlapping(0, 0, 0, 0, NormalizationHandler.NONE, false))
            for (ContactRecord record : block.getContactRecords()) {
                if (record.getCounts() != 1f) throw new AssertionError("score order in regional query");
                cells++;
            }
        if (cells != 1) throw new AssertionError("score cell missing");
        Iterator<ContactRecord> iterator = zoom.getDirectIterator();
        if (!iterator.hasNext() || iterator.next().getCounts() != 1f || iterator.hasNext())
            throw new AssertionError("score order in iterator");

        reader = new DatasetReaderV10(args[0] + "/large-vector-v10.hic", true);
        ds = reader.read();
        javastraw.reader.basics.Chromosome chr = ds.getChromosomeHandler().getChromosomeFromIndex(0);
        zoom = ds.getMatrix(chr, chr).getZoomData(new HiCZoom(1));
        cells = 0;
        for (Block block : zoom.getNormalizedBlocksOverlapping(0, 0, 0, 0, NormalizationHandler.VC, false))
            for (ContactRecord record : block.getContactRecords()) {
                if (record.getCounts() != 1.25f) throw new AssertionError("chunked normalization");
                cells++;
            }
        if (cells != 1) throw new AssertionError("normalized cell missing");
        iterator = zoom.getNormalizedIterator(NormalizationHandler.VC);
        if (!iterator.hasNext() || iterator.next().getCounts() != 1.25f || iterator.hasNext())
            throw new AssertionError("chunked normalized iterator");
        if (reader.readNormalizationVector(NormalizationHandler.VC_SQRT, 0, HiCZoom.HiCUnit.BP, 1) != null ||
            reader.getNormalizationVectorForQuery(0, new HiCZoom(1), NormalizationHandler.VC_SQRT) != null)
            throw new AssertionError("inferred VC_SQRT capability");
        System.out.println("V10 review: deterministic scores, bounded vectors, explicit capabilities passed");
    }
}
