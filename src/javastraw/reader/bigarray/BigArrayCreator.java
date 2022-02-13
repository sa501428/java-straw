package javastraw.reader.bigarray;

import javastraw.reader.Dataset;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.iterators.GenomeWideIterator;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.HiCZoom;

import java.util.Iterator;

public class BigArrayCreator {
    public static BigArray createFromZD(MatrixZoomData zd) {
        return populateBigArrayFromSingleIterator(zd.getDirectIterator(), 10000000);
    }

    public static BigArray createForWholeGenome(Dataset dataset, ChromosomeHandler handler,
                                                HiCZoom zoom, boolean includeIntra) {
        int limit = 10000000;
        if (includeIntra) {
            limit = 50000000;
        }
        return populateBigArrayFromSingleIterator(new GenomeWideIterator(dataset, handler, zoom, includeIntra), limit);
    }

    public static BigArray populateBigArrayFromSingleIterator(Iterator<ContactRecord> iterator, int limit) {
        BigArray allRecords = new BigArray(limit);
        int[] x = new int[limit];
        int[] y = new int[limit];
        float[] c = new float[limit];
        int counter = 0;
        while (iterator.hasNext()) {
            ContactRecord cr = iterator.next();
            x[counter] = cr.getBinX();
            y[counter] = cr.getBinY();
            c[counter] = cr.getCounts();
            counter++;
            if (counter > limit) {
                allRecords.addSubList(x, y, c);
                x = new int[limit];
                y = new int[limit];
                c = new float[limit];
                counter = 0;
            }
        }
        if (counter > 0) {
            allRecords.addSubList(x, y, c, counter);
        }
        return allRecords;
    }
}
