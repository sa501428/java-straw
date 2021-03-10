package javastraw.reader.iterators;

import javastraw.reader.Dataset;
import javastraw.reader.DatasetReader;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.block.Block;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.HiCZoom;
import org.broad.igv.util.collections.LRUCache;

import java.util.Iterator;
import java.util.List;

public class IteratorContainer {

    private final static int USE_ZD = 1, USE_LIST = 2, USE_GW = 3;
    private final LRUCache<String, Block> blockCache;
    private final DatasetReader reader;
    private final MatrixZoomData zd;
    private final boolean useCache;
    private final int type;
    private final List<ContactRecord> readList;
    private final Dataset dataset;
    private final ChromosomeHandler handler;
    private final HiCZoom zoom;
    private final boolean includeIntra;
    private final long matrixSize;
    private long numberOfContactRecords = -1;

    public IteratorContainer(DatasetReader reader, MatrixZoomData zd, LRUCache<String, Block> blockCache, boolean useCache) {
        type = USE_ZD;
        this.reader = reader;
        this.zd = zd;
        this.matrixSize = zd.getChr1().getLength() / zd.getZoom().getBinSize() + 1;
        this.blockCache = blockCache;
        this.readList = null;
        this.useCache = useCache;

        this.dataset = null;
        this.handler = null;
        this.zoom = null;
        this.includeIntra = false;
    }

    public IteratorContainer(List<ContactRecord> readList, long matrixSize) {
        type = USE_LIST;
        this.reader = null;
        this.zd = null;
        this.blockCache = null;
        this.readList = readList;
        this.matrixSize = matrixSize;
        this.useCache = false;

        this.dataset = null;
        this.handler = null;
        this.zoom = null;
        this.includeIntra = false;
    }

    public IteratorContainer(Dataset dataset, ChromosomeHandler handler,
                             HiCZoom zoom, boolean includeIntra) {
        type = USE_GW;
        this.reader = null;
        this.zd = null;
        this.blockCache = null;
        this.readList = null;
        this.useCache = false;

        this.dataset = dataset;
        this.handler = handler;
        this.zoom = zoom;
        this.includeIntra = includeIntra;

        long totalSize = 0;
        for (Chromosome c1 : handler.getChromosomeArrayWithoutAllByAll()) {
            totalSize += (c1.getLength() / zoom.getBinSize()) + 1;
        }
        this.matrixSize = totalSize;

    }

    public Iterator<ContactRecord> getNewContactRecordIterator() {
        if (type == USE_ZD) {
            if (reader != null && zd != null && blockCache != null) {
                return new ContactRecordIterator(reader, zd, blockCache, useCache);
            }
        }
        if (type == USE_LIST) {
            if (readList != null) {
                return readList.iterator();
            }
        }
        if (type == USE_GW) {
            if (dataset != null && handler != null && zoom != null) {
                return new GenomeWideIterator(dataset, handler, zoom, includeIntra);
            }
        }

        System.err.println("Null Contact Record Iterator");
        return null;
    }

    public long getNumberOfContactRecords() {
        if (numberOfContactRecords > 0) return numberOfContactRecords;

        numberOfContactRecords = 0;
        Iterator<ContactRecord> iterator = getNewContactRecordIterator();
        while (iterator.hasNext()) {
            iterator.next();
            numberOfContactRecords++;
        }

        return numberOfContactRecords;
    }

    public long getMatrixSize() {
        return matrixSize;
    }
}
