import javastraw.reader.DatasetReaderV10;
import javastraw.reader.Dataset;
import javastraw.reader.v10.V10AllOverview;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationHandler;
import javastraw.reader.block.ContactRecord;
import java.math.BigInteger;
import java.util.Iterator;

public class ScaledAllRegression {
    public static void main(String[] args) throws Exception {
        DatasetReaderV10 reader = new DatasetReaderV10(args[0], true);
        Dataset ds = reader.read();
        V10AllOverview all = reader.getHeader().allOverview;
        assert all != null && all.chromosomeIndex == 2 && all.finestResolution == 1;
        assert all.sourceResolution == 10 && all.groupFactor == 1 && all.binCount == 18;
        assert all.chromosomeForCoordinate(0) == 0;
        assert all.chromosomeForCoordinate(10) == 1;
        assert all.positionForCoordinate(1, 11) == 10;
        assert ds.getChromosomeHandler().getChromosomeArrayWithoutAllByAll().length == 2;
        assert ds.getChromosomeHandler().getChromosomeArrayWithoutAllByAll()[0].getName().equals("chrA");
        MatrixZoomData zd = reader.readMatrix("2_2", -1).getZoomData(new HiCZoom(HiCZoom.HiCUnit.BP, 1));
        Iterator<ContactRecord> records = zd.getNormalizedIterator(NormalizationHandler.NONE);
        BigInteger sum = BigInteger.ZERO;
        int count = 0;
        while (records.hasNext()) {
            ContactRecord record = records.next();
            if (record.getBinX() == 1 && record.getBinY() == 2) assert record.getExactCount().equals(BigInteger.ONE.shiftLeft(53).add(BigInteger.valueOf(17)));
            sum = sum.add(record instanceof javastraw.reader.block.V10ContactRecord ? record.getExactCount() : BigInteger.valueOf((long) record.getCounts()));
            count++;
        }
        assert count == 3 && sum.equals(BigInteger.ONE.shiftLeft(53).add(BigInteger.valueOf(49)));
        assert reader.readMatrix("0_2", -1) == null;
        System.out.println("Scaled All: metadata, exact counts, appended index, layout and omitted trans passed");
    }
}
