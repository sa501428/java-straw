package javastraw.reader.expected;

import javastraw.reader.Dataset;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.block.Block;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.HiCFileTools;

import java.util.List;

public class MedianExpectedVector {
    public static ExpectedValueFunction generateVector(Dataset dataset, HiCZoom zoom, NormalizationType type) {

        ChromosomeHandler handler = dataset.getChromosomeHandler();
        ExpectedValueCalculation ev = new ExpectedValueCalculation(handler, zoom.getBinSize(), type);

        // Loop through chromosomes
        for (Chromosome chr : handler.getChromosomeArrayWithoutAllByAll()) {

            MatrixZoomData zd = HiCFileTools.getMatrixZoomData(dataset, chr, chr, zoom);
            if (zd == null) continue;
            int lengthChr1 = (int) (chr.getLength() / zoom.getBinSize() + 1);

            try {
                List<Block> blocks = HiCFileTools.getAllRegionBlocks(zd, 0, lengthChr1, 0, lengthChr1,
                        type, false);
                for (Block block : blocks) {
                    for (ContactRecord cr : block.getContactRecords()) {
                        ev.addDistance(chr.getIndex(), cr);
                    }
                }
            } catch (Exception e) {
                e.printStackTrace();
            }

            zd.clearCache();
        }
        return ev.getExpectedValueFunction();

    }
}
