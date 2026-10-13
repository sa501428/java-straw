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
        if (args.length > 1) {
            DatasetReaderV10 kbReader = new DatasetReaderV10(args[1], true); kbReader.read();
            V10AllOverview kb = kbReader.getHeader().allOverview;
            assert kb.version == 2 && kb.syntheticLength == 170 && kb.storageStride == 17;
            assert kb.displayResolution == 17 && kb.effectiveResolution == 17000 && kb.binCount == 10;
            assert kb.syntheticBoundary(0) == 95 && kb.syntheticBoundary(1) == 170;
            assert kb.chromosomeForCoordinate(100) == 1 && kb.positionForCoordinate(1, 100) == 5000;
            Iterator<ContactRecord> raw = kbReader.readMatrix("2_2", -1)
                .getZoomData(new HiCZoom(HiCZoom.HiCUnit.BP, 1)).getNormalizedIterator(NormalizationHandler.NONE);
            sum = BigInteger.ZERO; count = 0;
            while (raw.hasNext()) {
                ContactRecord record = raw.next();
                assert record.getBinX() % 17 == 0 && record.getBinY() % 17 == 0;
                sum = sum.add(record.getExactCount()); ++count;
            }
            assert count == 2 && sum.equals(BigInteger.ONE.shiftLeft(53).add(BigInteger.valueOf(49)));
        }
        System.out.println("Scaled All: metadata, exact counts, appended index, layout and omitted trans passed");
    }
}
