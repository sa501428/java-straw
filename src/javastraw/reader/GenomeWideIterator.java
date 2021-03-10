package javastraw.reader;

import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ContactRecord;
import javastraw.type.HiCZoom;

import java.util.Iterator;

public class GenomeWideIterator implements Iterator<ContactRecord> {

    private final Chromosome[] chromosomes;
    private final boolean includeIntra;
    private final HiCZoom zoom;
    private final Dataset dataset;
    private final boolean iterationIsDone = false;
    private Iterator<ContactRecord> currentIterator = null;

    private int recentAddX = 0;
    private int recentAddY = 0;
    private int c1i = 0, c2i = 0;

    public GenomeWideIterator(Dataset dataset, ChromosomeHandler handler,
                              HiCZoom zoom, boolean includeIntra) {
        this.chromosomes = handler.getChromosomeArrayWithoutAllByAll();
        this.includeIntra = includeIntra;
        this.zoom = zoom;
        this.dataset = dataset;
        getNextIterator();
    }

    @Override
    public boolean hasNext() {
        if (currentIterator.hasNext()) {
            return true;
        } else {
            recentAddY += chromosomes[c2i].getLength() / zoom.getBinSize() + 1;
            c2i++;
        }
        return getNextIterator();
    }

    private boolean getNextIterator() {
        while (c1i < chromosomes.length) {
            Chromosome c1 = chromosomes[c1i];
            while (c2i < chromosomes.length) {
                Chromosome c2 = chromosomes[c2i];

                if (c1.getIndex() < c2.getIndex() || (c1.equals(c2) && includeIntra)) {
                    MatrixZoomData zd = HiCFileTools.getMatrixZoomData(dataset, c1, c2, zoom);
                    if (zd != null) {
                        currentIterator = zd.getIteratorContainer().getNewContactRecordIterator();
                        if (currentIterator != null && currentIterator.hasNext()) {
                            return true;
                        }
                    }
                }
                recentAddY += c2.getLength() / zoom.getBinSize() + 1;
                c2i++;
            }
            recentAddX += c1.getLength() / zoom.getBinSize() + 1;
            recentAddY = 0;
            c1i++;
            c2i = 0;
        }
        return false;
    }

    @Override
    public ContactRecord next() {
        ContactRecord cr = currentIterator.next();
        int binX = cr.getBinX() + recentAddX;
        int binY = cr.getBinY() + recentAddY;
        return new ContactRecord(binX, binY, cr.getCounts());
    }
}
