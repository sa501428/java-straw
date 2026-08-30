package javastraw.reader.mzd;

import javastraw.reader.basics.Chromosome;
import javastraw.reader.block.Block;
import javastraw.reader.block.BlockModifier;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.type.NormalizationType;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/** A caller-oriented view of a canonical inter-chromosomal matrix zoom. */
final class OrientedMatrixZoomData extends MatrixZoomData {
    private final MatrixZoomData delegate;

    OrientedMatrixZoomData(MatrixZoomData delegate) {
        super(delegate);
        if (delegate.getChr1Idx() == delegate.getChr2Idx()) {
            throw new IllegalArgumentException("An intra-chromosomal matrix cannot be reoriented");
        }
        this.delegate = delegate;
    }

    @Override public Chromosome getChr1() { return delegate.getChr2(); }
    @Override public Chromosome getChr2() { return delegate.getChr1(); }
    @Override public int getChr1Idx() { return delegate.getChr2Idx(); }
    @Override public int getChr2Idx() { return delegate.getChr1Idx(); }
    @Override public long getMatrixSize() { return getChr1().getLength() / getBinSize() + 1; }
    @Override public String getKey() {
        return triKey(getChr1().getName(), getChr2().getName(), getZoom().getKey());
    }
    @Override public String getDescription() {
        return getChr1().getName() + " - " + getChr2().getName() + " - " + getZoom();
    }

    @Override
    public List<Block> getNormalizedBlocksOverlapping(long binX1, long binY1, long binX2, long binY2,
                                                      NormalizationType norm, boolean fillUnderDiagonal) {
        return transpose(delegate.getNormalizedBlocksOverlapping(
                binY1, binX1, binY2, binX2, norm, fillUnderDiagonal));
    }

    @Override
    public List<Block> getNormalizedBlocksOverlapping(long binX1, long binY1, long binX2, long binY2,
                                                      NormalizationType norm, boolean fillUnderDiagonal,
                                                      BlockModifier modifier) {
        List<Block> oriented = getNormalizedBlocksOverlapping(
                binX1, binY1, binX2, binY2, norm, fillUnderDiagonal);
        List<Block> modified = new ArrayList<>(oriented.size());
        for (Block block : oriented) {
            modified.add(modifier.modify(block, getKey(), getBinSize(), getChr1(), getChr2()));
        }
        return modified;
    }

    @Override public Iterator<ContactRecord> getDirectIterator() {
        return transpose(delegate.getDirectIterator());
    }
    @Override public Iterator<ContactRecord> getNormalizedIterator(NormalizationType normType) {
        return transpose(delegate.getNormalizedIterator(normType));
    }
    @Override public double getAverageCount() { return delegate.getAverageCount(); }
    @Override public Integer getBlockSize(int blockNum) { return delegate.getBlockSize(blockNum); }
    @Override void clearCache() { delegate.clearCache(); }

    private List<Block> transpose(List<Block> blocks) {
        List<Block> result = new ArrayList<>(blocks.size());
        for (Block block : blocks) {
            List<ContactRecord> records = new ArrayList<>(block.getContactRecords().size());
            for (ContactRecord record : block.getContactRecords()) records.add(record.transpose());
            result.add(new Block(block.getNumber(), records, getKey()));
        }
        return result;
    }

    private Iterator<ContactRecord> transpose(final Iterator<ContactRecord> source) {
        return new Iterator<ContactRecord>() {
            @Override public boolean hasNext() { return source.hasNext(); }
            @Override public ContactRecord next() { return source.next().transpose(); }
            @Override public void remove() { source.remove(); }
        };
    }
}
