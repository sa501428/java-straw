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
import javastraw.tools.ParallelizationTools;
import org.broad.igv.util.collections.LRUCache;

import java.io.IOException;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.atomic.AtomicInteger;

public class MatrixZoomData {

    final Chromosome chr1;  // Chromosome on the X axis
    final Chromosome chr2;  // Chromosome on the Y axis
    private final boolean isIntra;
    final HiCZoom zoom;    // Unit and bin size
    // Observed values are organized into sub-matrices ("blocks")
    protected final int blockBinCount;   // block size in bins
    protected final int blockColumnCount;     // number of block columns
    private final long correctedBinCount;
    // Cache the last 20 blocks loaded
    protected final LRUCache<String, Block> blockCache = new LRUCache<>(500);
    private final V9Depth v9Depth;
    private double averageCount = -1;
    protected DatasetReader reader;
    private IteratorContainer iteratorContainer = null;
    private boolean useCache = true;
    private final Map<NormalizationType, BasicMatrix> pearsonsMap;
    private final Map<String, double[]> eigenvectorMap;
    public static boolean useIteratorDontPutAllInRAM = false;
    public static boolean shouldCheckRAMUsage = false;
    private final BlockModifier identity = new IdentityModifier();

    public MatrixZoomData(Chromosome chr1, Chromosome chr2, HiCZoom zoom, int blockBinCount, int blockColumnCount,
                          int[] chr1Sites, int[] chr2Sites, DatasetReader reader) {
        this.chr1 = chr1;
        this.chr2 = chr2;
        this.zoom = zoom;
        this.isIntra = chr1.getIndex() == chr2.getIndex();
        this.reader = reader;
        this.blockBinCount = blockBinCount;
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

    public void setUseCache(boolean useCache) {
        this.useCache = useCache;
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

    public String getKey() {
        return chr1.getName() + "_" + chr2.getName() + "_" + zoom.getKey();
    }

    public String getKey(int chr1, int chr2) {
        return chr1 + "_" + chr2 + "_" + zoom.getKey();
    }

    public String getBlockKey(int blockNumber, NormalizationType no) {
        return getKey() + "_" + blockNumber + "_" + no;
    }

    public String getNormLessBlockKey(Block block) {
        return getKey() + "_" + block.getNumber() + "_" + block.getUniqueRegionID();
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
            return addNormalizedBlocksToListV9(blockList, (int) binX1, (int) binY1, (int) binX2, (int) binY2, no, modifier);
        } else {
            return addNormalizedBlocksToList(blockList, (int) binX1, (int) binY1, (int) binX2, (int) binY2, no, fillUnderDiagonal, modifier);
        }
    }

    public int getBlockNumberVersion9FromPADAndDepth(int positionAlongDiagonal, int depth) {
        return depth * blockColumnCount + positionAlongDiagonal;
    }

    protected void populateBlocksToLoadV9(int positionAlongDiagonal, int depth, NormalizationType no, List<Block> blockList, Set<Integer> blocksToLoad) {
        int blockNumber = getBlockNumberVersion9FromPADAndDepth(positionAlongDiagonal, depth);
        String key = getBlockKey(blockNumber, no);
        Block b;
        if (useCache && blockCache.containsKey(key)) {
            b = blockCache.get(key);
            blockList.add(b);
        } else {
            blocksToLoad.add(blockNumber);
        }
    }

    protected List<Block> addNormalizedBlocksToListV9(final List<Block> blockList, int binX1, int binY1, int binX2, int binY2,
                                                      final NormalizationType norm, BlockModifier modifier) {

        Set<Integer> blocksToLoad = new HashSet<>();

        // PAD = positionAlongDiagonal (~projected)
        // Depth is axis perpendicular to diagonal; nearer means closer to diagonal
        int translatedLowerPAD = (binX1 + binY1) / 2 / blockBinCount;
        int translatedHigherPAD = (binX2 + binY2) / 2 / blockBinCount + 1;
        int translatedNearerDepth = v9Depth.getDepth(binX1, binY2);
        int translatedFurtherDepth = v9Depth.getDepth(binX2, binY1);

        // because code above assume above diagonal; but we could be below diagonal
        int nearerDepth = Math.min(translatedNearerDepth, translatedFurtherDepth);
        if ((binX1 > binY2 && binX2 < binY1) || (binX2 > binY1 && binX1 < binY2)) {
            nearerDepth = 0;
        }
        int furtherDepth = Math.max(translatedNearerDepth, translatedFurtherDepth) + 1; // +1; integer divide rounds down


        for (int depth = nearerDepth; depth <= furtherDepth; depth++) {
            for (int pad = translatedLowerPAD; pad <= translatedHigherPAD; pad++) {
                populateBlocksToLoadV9(pad, depth, norm, blockList, blocksToLoad);
            }
        }

        actuallyLoadGivenBlocks(blockList, blocksToLoad, norm, modifier);

        return new ArrayList<>(new HashSet<>(blockList));
    }

    protected void populateBlocksToLoad(int r, int c, NormalizationType no, List<Block> blockList, Set<Integer> blocksToLoad) {
        int blockNumber = r * getBlockColumnCount() + c;
        String key = getBlockKey(blockNumber, no);
        Block b;
        if (useCache && blockCache.containsKey(key)) {
            b = blockCache.get(key);
            blockList.add(b);
        } else {
            blocksToLoad.add(blockNumber);
        }
    }

    /**
     * Return the blocks of normalized, observed values overlapping the rectangular region specified.
     *
     * @param binY1 leftmost position in "bins"
     * @param binX2 rightmost position in "bins"
     * @param binY2 bottom position in "bins"
     * @param norm  normalization type
     * @return List of overlapping blocks, normalized
     */
    protected List<Block> addNormalizedBlocksToList(final List<Block> blockList, int binX1, int binY1, int binX2, int binY2,
                                                    final NormalizationType norm, boolean getBelowDiagonal, BlockModifier modifier) {

        Set<Integer> blocksToLoad = new HashSet<>();

        // have to do this regardless (just in case)
        int col1 = binX1 / blockBinCount;
        int row1 = binY1 / blockBinCount;
        int col2 = binX2 / blockBinCount;
        int row2 = binY2 / blockBinCount;

        for (int r = row1; r <= row2; r++) {
            for (int c = col1; c <= col2; c++) {
                populateBlocksToLoad(r, c, norm, blockList, blocksToLoad);
            }
        }

        if (getBelowDiagonal && binY1 < binX2) {
            for (int r = row1; r <= row2; r++) {
                for (int c = col1; c <= col2; c++) {
                    populateBlocksToLoad(c, r, norm, blockList, blocksToLoad);
                }
            }
        }

        actuallyLoadGivenBlocks(blockList, blocksToLoad, norm, modifier);

        return new ArrayList<>(new HashSet<>(blockList));
    }

    protected void actuallyLoadGivenBlocks(final List<Block> blockList, Set<Integer> blocksToLoad,
                                           final NormalizationType no, BlockModifier modifier) {
        final AtomicInteger errorCounter = new AtomicInteger();

        ExecutorService service = Executors.newFixedThreadPool(200);
        for (final int blockNumber : blocksToLoad) {
            String key = getBlockKey(blockNumber, no);
            readBlockUpdateListAndCache(blockNumber, reader, no, blockList, key, errorCounter, service, modifier);
        }
        ParallelizationTools.shutDownServiceAndWait(service, errorCounter);
    }

    protected void readBlockUpdateListAndCache(int blockNumber, DatasetReader reader, NormalizationType no,
                                               List<Block> blockList, String key, final AtomicInteger errorCounter,
                                               ExecutorService service, BlockModifier modifier) {
        Runnable loader = () -> {
            try {
                Block b = reader.readNormalizedBlock(blockNumber, MatrixZoomData.this, no);
                if (b == null) {
                    b = new Block(blockNumber, key);
                }
                b = modifier.modify(b);
                if (useCache) {
                    blockCache.put(key, b);
                }
                blockList.add(b);
            } catch (IOException e) {
                errorCounter.incrementAndGet();
            }
        };
        service.submit(loader);
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

    // todo V9 needs a diff method
    protected List<Integer> getBlockNumbersForRegionFromBinPosition(long[] regionIndices) {

        // cast should be fine - this is for V8
        int col1 = (int) (regionIndices[0] / blockBinCount);
        int col2 = (int) ((regionIndices[1] + 1) / blockBinCount);
        int row1 = (int) (regionIndices[2] / blockBinCount);
        int row2 = (int) ((regionIndices[3] + 1) / blockBinCount);

        // first check the upper triangular matrix
        Set<Integer> blocksSet = new HashSet<>();
        for (int r = row1; r <= row2; r++) {
            for (int c = col1; c <= col2; c++) {
                int blockNumber = r * getBlockColumnCount() + c;
                blocksSet.add(blockNumber);
            }
        }
        // check region part that overlaps with lower left triangle
        // but only if intrachromosomal
        if (chr1.getIndex() == chr2.getIndex()) {
            for (int r = col1; r <= col2; r++) {
                for (int c = row1; c <= row2; c++) {
                    int blockNumber = r * getBlockColumnCount() + c;
                    blocksSet.add(blockNumber);
                }
            }
        }

        List<Integer> blocksToIterateOver = new ArrayList<>(blocksSet);
        Collections.sort(blocksToIterateOver);
        return blocksToIterateOver;
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
            iteratorContainer = ListOfListGenerator.createFromZD(reader, this, blockCache, useCache,
                    useIteratorDontPutAllInRAM, shouldCheckRAMUsage);
        }
        return iteratorContainer;
    }

    public IteratorContainer getFromFileIteratorContainer() {
        return new ZDIteratorContainer(reader, this, blockCache, useCache);
    }


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
}
