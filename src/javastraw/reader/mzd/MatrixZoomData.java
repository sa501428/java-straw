/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2011-2020 Broad Institute, Aiden Lab, Rice University, Baylor College of Medicine
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 *  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 *  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 *  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 *  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 *  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 *  THE SOFTWARE.
 */


package javastraw.reader.mzd;

import javastraw.matrices.BasicMatrix;
import javastraw.reader.DatasetReader;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.block.Block;
import javastraw.reader.block.BlockModifier;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.block.IdentityModifier;
import javastraw.reader.depth.LogDepth;
import javastraw.reader.depth.V9Depth;
import javastraw.reader.expected.ExpectedValueFunction;
import javastraw.reader.iterators.IteratorContainer;
import javastraw.reader.iterators.ListOfListGenerator;
import javastraw.reader.iterators.ZDIteratorContainer;
import javastraw.reader.pearsons.PearsonsManager;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationType;

import java.util.*;

public class MatrixZoomData {

    protected final Chromosome chr1;  // Chromosome on the X axis
    protected final Chromosome chr2;  // Chromosome on the Y axis
    protected final boolean isIntra;
    protected final HiCZoom zoom;    // Unit and bin size
    // Observed values are organized into sub-matrices ("blocks")
    protected final int blockBinCount;   // block size in bins
    protected final int blockColumnCount;     // number of block columns
    protected final long correctedBinCount;
    // Cache the last 20 blocks loaded
    protected final BlockCache blockCache;
    protected final V9Depth v9Depth;
    protected final Map<NormalizationType, BasicMatrix> pearsonsMap;
    protected DatasetReader reader;
    protected final Map<String, double[]> eigenvectorMap;
    protected final BlockModifier identity = new IdentityModifier();
    protected double averageCount = -1;
    protected IteratorContainer iteratorContainer = null;
    public static boolean useIteratorDontPutAllInRAM = false;
    public static boolean shouldCheckRAMUsage = false;

    public MatrixZoomData(Chromosome chr1, Chromosome chr2, HiCZoom zoom, int blockBinCount, int blockColumnCount,
                          int[] chr1Sites, int[] chr2Sites, DatasetReader reader) {
        this.chr1 = chr1;
        this.chr2 = chr2;
        this.zoom = zoom;
        this.isIntra = chr1.getIndex() == chr2.getIndex();
        this.reader = reader;
        this.blockBinCount = blockBinCount;
        blockCache = new BlockCache();
        if (reader.getVersion() > 8) {
            v9Depth = V9Depth.setDepthMethod(reader.getDepthBase(), blockBinCount);
        } else {
            v9Depth = new LogDepth(2, blockBinCount);
        }
        this.blockColumnCount = blockColumnCount;
        if (!(this instanceof DynamicMatrixZoomData)) {
            if (reader.getVersion() < 8 && chr1.getLength() < chr2.getLength()) {
                boolean isFrag = zoom.getUnit() == HiCZoom.HiCUnit.FRAG;
                long len1 = chr1.getLength();
                long len2 = chr2.getLength();
                if (chr1Sites != null && chr2Sites != null && isFrag) {
                    len1 = chr1Sites.length + 1;
                    len2 = chr2Sites.length + 1;
                }
                long nBinsX = Math.max(len1, len2) / zoom.getBinSize() + 1;
                correctedBinCount = nBinsX / blockColumnCount + 1;
            } else {
                correctedBinCount = blockBinCount;
            }
        } else {
            correctedBinCount = blockBinCount;
        }
        pearsonsMap = new HashMap<>();
        eigenvectorMap = new HashMap<>();
    }

    protected MatrixZoomData(MatrixZoomData zd0) {
        this.chr1 = zd0.chr1;
        this.chr2 = zd0.chr2;
        this.isIntra = zd0.isIntra;
        this.zoom = zd0.zoom;
        this.blockBinCount = zd0.blockBinCount;
        this.blockColumnCount = zd0.blockColumnCount;
        this.correctedBinCount = zd0.correctedBinCount;
        this.blockCache = zd0.blockCache;
        this.v9Depth = zd0.v9Depth;
        this.averageCount = zd0.averageCount;
        this.reader = zd0.reader;
        this.iteratorContainer = zd0.iteratorContainer;
        this.pearsonsMap = zd0.pearsonsMap;
        this.eigenvectorMap = zd0.eigenvectorMap;
    }

    public void setUseCache(boolean useCache) {
        blockCache.setUseCache(useCache);
    }

    public Chromosome getChr1() {
        return chr1;
    }

    public Chromosome getChr2() {
        return chr2;
    }

    public int getBinSize() {
        return zoom.getBinSize();
    }

    public int getChr1Idx() {
        return chr1.getIndex();
    }

    public int getChr2Idx() {
        return chr2.getIndex();
    }

    public long getMatrixSize() {
        return chr1.getLength() / zoom.getBinSize() + 1;
    }

    public long getCorrectedBinCount() {
        return correctedBinCount;
    }

    public int getBlockBinCount() {
        return blockBinCount;
    }

    public HiCZoom getZoom() {
        return zoom;
    }

    public int getBlockColumnCount() {
        return blockColumnCount;
    }

    public static String triKey(String s1, String s2, String s3) {
        return s1 + "_" + s2 + "_" + s3;
    }

    public String getKey() {
        return triKey(chr1.getName(), chr2.getName(), zoom.getKey());
    }

    public String getKey(int chr1, int chr2) {
        return triKey("" + chr1, "" + chr2, zoom.getKey());
    }

    public String getBlockKey(int blockNumber, NormalizationType no) {
        return triKey(getKey(), "" + blockNumber, "" + no);
    }

    public String getNormLessBlockKey(Block block) {
        return triKey(getKey(), "" + block.getNumber(), block.getUniqueRegionID());
    }

    public List<Block> getNormalizedBlocksOverlapping(long binX1, long binY1, long binX2, long binY2,
                                                      final NormalizationType no, boolean isImportant,
                                                      boolean fillUnderDiagonal) {
        return getNormalizedBlocksOverlapping(binX1, binY1, binX2, binY2,
                no, isImportant, fillUnderDiagonal, identity);
    }

    public List<Block> getNormalizedBlocksOverlapping(long binX1, long binY1, long binX2, long binY2,
                                                      final NormalizationType no, boolean isImportant,
                                                      boolean fillUnderDiagonal, BlockModifier modifier) {
        final List<Block> blockList = Collections.synchronizedList(new ArrayList<>());
        if (reader.getVersion() > 8 && isIntra) {
            return V9IntraBlockReader.addNormalizedBlocksToListV9(blockList, (int) binX1, (int) binY1,
                    (int) binX2, (int) binY2, no, modifier, blockBinCount, v9Depth,
                    blockColumnCount, blockCache, getKey(), chr1, chr2, zoom, reader);
        } else {
            return LegacyVersionBlockReader.addNormalizedBlocksToList(blockList, (int) binX1, (int) binY1, (int) binX2,
                    (int) binY2, no, fillUnderDiagonal, modifier, blockBinCount, blockColumnCount, blockCache,
                    getKey(), chr1, chr2, zoom, reader);
        }
    }

    /**
     * Utility for printing description of this matrix.
     */
    public String getDescription() {
        return chr1.getName() + " - " + chr2.getName() + " - " + getZoom();
    }

    public void printFullDescription() {
        System.out.println("Chromosomes: " + chr1.getName() + " - " + chr2.getName());
        System.out.println("unit: " + zoom.getUnit());
        System.out.println("binSize (bp): " + zoom.getBinSize());
        System.out.println("blockBinCount (bins): " + blockBinCount);
        System.out.println("blockColumnCount (columns): " + blockColumnCount);
        System.out.println("Block size (bp): " + blockBinCount * zoom.getBinSize());
        System.out.println();
    }

    /**
     * For a specified region, select the block numbers corresponding to it
     */
    protected List<Integer> getBlockNumbersForRegionFromGenomePosition(long[] regionIndices) {
        int resolution = zoom.getBinSize();
        long[] regionBinIndices = new long[4];
        for (int i = 0; i < regionBinIndices.length; i++) {
            regionBinIndices[i] = regionIndices[i] / resolution;
        }
        return getBlockNumbersForRegionFromBinPosition(regionBinIndices);
    }

    private List<Integer> getBlockNumbersForRegionFromBinPosition(long[] regionBinIndices) {
        if (reader.getVersion() > 8 && isIntra) {
            return V9IntraBlockReader.getBlockNumbersForRegionFromBinPosition(regionBinIndices,
                    blockBinCount, blockColumnCount, isIntra, v9Depth);
        } else {
            return LegacyVersionBlockReader.getBlockNumbersForRegionFromBinPosition(regionBinIndices,
                    blockBinCount, blockColumnCount, isIntra);
        }
    }

    /**
     * Returns the average count
     *
     * @return Average count
     */
    public double getAverageCount() {
        return averageCount;
    }

    /**
     * Sets the average count
     *
     * @param averageCount Average count to set
     */
    public void setAverageCount(double averageCount) {
        this.averageCount = averageCount;
    }

    public void clearCache() {
        blockCache.clear();
        pearsonsMap.clear();
        eigenvectorMap.clear();
    }

    public long getNumberOfContactRecords() {
        return getIteratorContainer().getNumberOfContactRecords();
    }

    protected Iterator<ContactRecord> getNewContactRecordIterator() {
        return getIteratorContainer().getNewContactRecordIterator();
        //return new ContactRecordIterator(reader, this, blockCache);
    }

    public IteratorContainer getIteratorContainer() {
        if (iteratorContainer == null) {
            iteratorContainer = ListOfListGenerator.createFromZD(reader, this, blockCache,
                    useIteratorDontPutAllInRAM, shouldCheckRAMUsage);
        }
        return iteratorContainer;
    }

    public IteratorContainer getFromFileIteratorContainer() {
        return new ZDIteratorContainer(reader, this, blockCache);
    }

    /*
    public BigListOfContactRecords getContactRecords(){
        new ContactRecordIterator(reader, getKey(), blockCache,
                getChr1Idx(), getChr2Idx(), getZoom());
    }
    */


    public BasicMatrix getPearsons(ExpectedValueFunction df) {
        if (chr1.getIndex() != chr2.getIndex()) {
            throw new RuntimeException("Cannot compute pearsons for non-diagonal matrices");
        }

        BasicMatrix pearsons = pearsonsMap.get(df.getNormalizationType());
        if (pearsons != null) {
            return pearsons;
        }

        pearsons = PearsonsManager.computePearsons(df, getNewContactRecordIterator(), chr1, zoom.getBinSize());
        pearsonsMap.put(df.getNormalizationType(), pearsons);
        return pearsonsMap.get(df.getNormalizationType());
    }

    protected BasicMatrix getPearsons(NormalizationType type) {
        return pearsonsMap.get(type);
    }

    public float getPearsonValue(int binX, int binY, NormalizationType type) {
        BasicMatrix pearsons = pearsonsMap.get(type);
        if (pearsons != null) {
            return pearsons.getEntry(binX, binY);
        } else {
            return 0;
        }
    }

    public double[] getEigenvector(ExpectedValueFunction df, int which) {
        if (chr1.getIndex() != chr2.getIndex()) {
            throw new RuntimeException("Cannot compute eigenvector for non-diagonal matrices");
        }

        String eigKey = getEigenvectorKey(df.getNormalizationType(), which);

        double[] eig = eigenvectorMap.get(eigKey);
        if (eig != null) {
            return eig;
        }

        BasicMatrix pearsons = getPearsons(df);
        if (pearsons == null) {
            return null;
        }

        eig = PearsonsManager.computeEigenvector(pearsons, which);
        eigenvectorMap.put(eigKey, eig);
        return eigenvectorMap.get(eigKey);
    }

    protected String getEigenvectorKey(NormalizationType normalizationType, int which) {
        return normalizationType.getLabel() + "_" + which;
    }

    protected double[] getEigenvector(NormalizationType normalizationType, int which) {
        return eigenvectorMap.get(getEigenvectorKey(normalizationType, which));
    }
}
