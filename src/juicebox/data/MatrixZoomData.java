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


package juicebox.data;

import htsjdk.tribble.util.LittleEndianOutputStream;
import juicebox.HiCGlobals;
import juicebox.data.basics.Chromosome;
import juicebox.windowui.HiCZoom;
import juicebox.windowui.MatrixType;
import juicebox.windowui.NormalizationType;
import org.broad.igv.util.collections.LRUCache;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;



/**
 * @author jrobinso
 * @since Aug 10, 2010
 */
public class MatrixZoomData {

    final Chromosome chr1;  // Chromosome on the X axis
    final Chromosome chr2;  // Chromosome on the Y axis
    private final boolean isIntra;
    final HiCZoom zoom;    // Unit and bin size
    // Observed values are organized into sub-matrices ("blocks")
    protected final int blockBinCount;   // block size in bins
    protected final int blockColumnCount;     // number of block columns
    // Cache the last 20 blocks loaded
    protected final LRUCache<String, Block> blockCache = new LRUCache<>(500);
    DatasetReader reader;
    private double averageCount = -1;
    private List<List<ContactRecord>> localCacheOfRecords = null;
    private long numberOfContactRecords = 0;

    /**
     * Constructor, sets the grid axes.  Called when read from file.
     *
     * @param chr1             Chromosome 1
     * @param chr2             Chromosome 2
     * @param zoom             Zoom (bin size and BP or FRAG)
     * @param blockBinCount    Number of bins divided by number of columns (around 1000)
     * @param blockColumnCount Number of bins divided by 1000 (BLOCK_SIZE)
     * @param chr1Sites        Used for looking up fragment
     * @param chr2Sites        Used for looking up fragment
     * @param reader           Pointer to file reader
     */
    public MatrixZoomData(Chromosome chr1, Chromosome chr2, HiCZoom zoom, int blockBinCount, int blockColumnCount,
                          int[] chr1Sites, int[] chr2Sites, DatasetReader reader) {
        this.chr1 = chr1;
        this.chr2 = chr2;
        this.zoom = zoom;
        this.isIntra = chr1.getIndex() == chr2.getIndex();
        this.reader = reader;
        this.blockBinCount = blockBinCount;
        this.blockColumnCount = blockColumnCount;

        long correctedBinCount = blockBinCount;
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
            }
        }
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

    public HiCZoom getZoom() {
        return zoom;
    }

    private int getBlockColumnCount() {
        return blockColumnCount;
    }

    public String getKey() {
        return chr1.getName() + "_" + chr2.getName() + "_" + zoom.getKey();
    }

    // i think this is how it should be? todo sxxgrc please confirm use case
    private String getKey(int chr1, int chr2) {
        return chr1 + "_" + chr2 + "_" + zoom.getKey();
    }

    public String getBlockKey(int blockNumber, NormalizationType no) {
        return getKey() + "_" + blockNumber + "_" + no;
    }

    public String getNormLessBlockKey(Block block) {
        return getKey() + "_" + block.getNumber() + "_" + block.getUniqueRegionID();
    }

    private String getBlockKey(int blockNumber, NormalizationType no, int chr1, int chr2) {
        return getKey(chr1, chr2) + "_" + blockNumber + "_" + no;
    }
    
    private static int log2(double v) {
        return (int) (Math.log(v) / Math.log(2));
    }
    
    /**
     * Return the blocks of normalized, observed values overlapping the rectangular region specified.
     * The units are "bins"
     *
     * @param binY1       leftmost position in "bins"
     * @param binX2       rightmost position in "bins"
     * @param binY2       bottom position in "bins"
     * @param no          normalization type
     * @param isImportant used for debugging
     * @return List of overlapping blocks, normalized
     */
    public List<Block> getNormalizedBlocksOverlapping(long binX1, long binY1, long binX2, long binY2, final NormalizationType no,
                                                      boolean isImportant, boolean fillUnderDiagonal) {
        
        final List<Block> blockList = Collections.synchronizedList(new ArrayList<>());
        if (reader.getVersion() > 8 && isIntra) {
            return addNormalizedBlocksToListV9(blockList, (int) binX1, (int) binY1, (int) binX2, (int) binY2, no);
        } else {
            return addNormalizedBlocksToList(blockList, (int) binX1, (int) binY1, (int) binX2, (int) binY2, no, fillUnderDiagonal);
        }
    }
    
    public int getBlockNumberVersion9FromPADAndDepth(int positionAlongDiagonal, int depth) {
        return depth * blockColumnCount + positionAlongDiagonal;
    }
    
    private void populateBlocksToLoadV9(int positionAlongDiagonal, int depth, NormalizationType no, List<Block> blockList, Set<Integer> blocksToLoad) {
        int blockNumber = getBlockNumberVersion9FromPADAndDepth(positionAlongDiagonal, depth);
        String key = getBlockKey(blockNumber, no);
        Block b;
        if (HiCGlobals.useCache && blockCache.containsKey(key)) {
            b = blockCache.get(key);
            blockList.add(b);
        } else {
            blocksToLoad.add(blockNumber);
        }
    }
    
    // for reference
    public int getBlockNumberVersion9(int binI, int binJ) {
        //int numberOfBlocksOnDiagonal = numberOfBinsInThisIntraMatrixAtResolution / blockSizeInBinCount + 1;
        // assuming number of blocks on diagonal is blockClolumnSize
        int depth = log2(1 + Math.abs(binI - binJ) / Math.sqrt(2) / blockBinCount);
        int positionAlongDiagonal = ((binI + binJ) / 2 / blockBinCount);
        return getBlockNumberVersion9FromPADAndDepth(positionAlongDiagonal, depth);
    }

    public int[] getBlockBoundsFromNumberVersion9Up(int blockNumber) {
        int positionAlongDiagonal = blockNumber % blockColumnCount;
        int depth = blockNumber / blockColumnCount;
        int avgPosition1 = positionAlongDiagonal * blockBinCount;
        int avgPosition2 = (positionAlongDiagonal + 1) * blockBinCount;
        double difference1 = (Math.pow(2, depth) - 1) * blockBinCount * Math.sqrt(2);
        double difference2 = (Math.pow(2, depth + 1) - 1) * blockBinCount * Math.sqrt(2);
        int c1 = avgPosition1 + (int) difference1 / 2 - 1;
        int c2 = avgPosition2 + (int) difference2 / 2 + 1;
        int r1 = avgPosition1 - (int) difference2 / 2 - 1;
        int r2 = avgPosition2 - (int) difference1 / 2 + 1;
        return new int[]{c1, c2, r1, r2};
    }
    
    public int[] getBlockBoundsFromNumberVersion8Below(int blockNumber) {
        int c = (blockNumber % blockColumnCount);
        int r = blockNumber / blockColumnCount;
        int c1 = c * blockBinCount;
        int c2 = c1 + blockBinCount - 1;
        int r1 = r * blockBinCount;
        int r2 = r1 + blockBinCount - 1;
        return new int[]{c1, c2, r1, r2};
    }
    
    private List<Block> addNormalizedBlocksToListV9(final List<Block> blockList, int binX1, int binY1, int binX2, int binY2,
                                                    final NormalizationType norm) {
        
        Set<Integer> blocksToLoad = new HashSet<>();
        
        // PAD = positionAlongDiagonal (~projected)
        // Depth is axis perpendicular to diagonal; nearer means closer to diagonal
        int translatedLowerPAD = (binX1 + binY1) / 2 / blockBinCount;
        int translatedHigherPAD = (binX2 + binY2) / 2 / blockBinCount + 1;
        int translatedNearerDepth = log2(1 + Math.abs(binX1 - binY2) / Math.sqrt(2) / blockBinCount);
        int translatedFurtherDepth = log2(1 + Math.abs(binX2 - binY1) / Math.sqrt(2) / blockBinCount);
        
        // because code above assume above diagonal; but we could be below diagonal
        int nearerDepth = Math.min(translatedNearerDepth, translatedFurtherDepth);
        int furtherDepth = Math.max(translatedNearerDepth, translatedFurtherDepth) + 1; // +1; integer divide rounds down
        
        
        for (int depth = nearerDepth; depth <= furtherDepth; depth++) {
            for (int pad = translatedLowerPAD; pad <= translatedHigherPAD; pad++) {
                populateBlocksToLoadV9(pad, depth, norm, blockList, blocksToLoad);
            }
        }
        
        actuallyLoadGivenBlocks(blockList, blocksToLoad, norm);
        
        return new ArrayList<>(new HashSet<>(blockList));
    }
    
    private void populateBlocksToLoad(int r, int c, NormalizationType no, List<Block> blockList, Set<Integer> blocksToLoad) {
        int blockNumber = r * getBlockColumnCount() + c;
        String key = getBlockKey(blockNumber, no);
        Block b;
        if (HiCGlobals.useCache && blockCache.containsKey(key)) {
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
    private List<Block> addNormalizedBlocksToList(final List<Block> blockList, int binX1, int binY1, int binX2, int binY2,
                                                  final NormalizationType norm, boolean getBelowDiagonal) {

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

        actuallyLoadGivenBlocks(blockList, blocksToLoad, norm);

        return new ArrayList<>(new HashSet<>(blockList));
    }

    private List<Block> addNormalizedBlocksToList(final List<Block> blockList, int binX1, int binY1, int binX2, int binY2,
                                                  final NormalizationType no, int chr1, int chr2) {

        Set<Integer> blocksToLoad = new HashSet<>();
    
        // for V8 - these will always be ints
        // have to do this regardless (just in case)
        int col1 = binX1 / blockBinCount;
        int row1 = binY1 / blockBinCount;
        int col2 = binX2 / blockBinCount;
        int row2 = binY2 / blockBinCount;

        for (int r = row1; r <= row2; r++) {
            for (int c = col1; c <= col2; c++) {
                populateBlocksToLoad(r, c, no, blockList, blocksToLoad);
            }
        }

        actuallyLoadGivenBlocks(blockList, blocksToLoad, no, chr1, chr2);
        System.out.println("I am block size: " + blockList.size());
        System.out.println("I am first block: " + blockList.get(0).getNumber());
        return new ArrayList<>(new HashSet<>(blockList));
    }

    private void actuallyLoadGivenBlocks(final List<Block> blockList, Set<Integer> blocksToLoad,
                                         final NormalizationType no) {
        final AtomicInteger errorCounter = new AtomicInteger();

        ExecutorService service = Executors.newFixedThreadPool(200);

        final int binSize = getBinSize();
        final int chr1Index = chr1.getIndex();
        final int chr2Index = chr2.getIndex();

        for (final int blockNumber : blocksToLoad) {
            Runnable loader = new Runnable() {
                @Override
                public void run() {
                    try {
                        String key = getBlockKey(blockNumber, no);
                        Block b = reader.readNormalizedBlock(blockNumber, MatrixZoomData.this, no);
                        if (b == null) {
                            b = new Block(blockNumber, key);   // An empty block
                        }
                        //Run out of memory if do it here
                        if (HiCGlobals.useCache) {
                            blockCache.put(key, b);
                        }
                        blockList.add(b);
                    } catch (IOException e) {
                        errorCounter.incrementAndGet();
                    }
                }
            };

            service.submit(loader);
        }

        // done submitting all jobs
        service.shutdown();

        // wait for all to finish
        try {
            service.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
        } catch (InterruptedException e) {
            System.err.println("Error loading mzd data " + e.getLocalizedMessage());
            if (HiCGlobals.printVerboseComments) {
                e.printStackTrace();
            }
        }

        // error printing
        if (errorCounter.get() > 0) {
            System.err.println(errorCounter.get() + " errors while reading blocks");
        }
    }

    private void actuallyLoadGivenBlocks(final List<Block> blockList, Set<Integer> blocksToLoad,
                                         final NormalizationType no, final int chr1Id, final int chr2Id) {
        final AtomicInteger errorCounter = new AtomicInteger();

        ExecutorService service = Executors.newFixedThreadPool(200);

        final int binSize = getBinSize();

        for (final int blockNumber : blocksToLoad) {
            Runnable loader = new Runnable() {
                @Override
                public void run() {
                    try {
                        String key = getBlockKey(blockNumber, no, chr1Id, chr2Id);
                        Block b = reader.readNormalizedBlock(blockNumber, MatrixZoomData.this, no);
                        if (b == null) {
                            b = new Block(blockNumber, key);   // An empty block
                        }
                        //Run out of memory if do it here
                        if (HiCGlobals.useCache) {
                            blockCache.put(key, b);
                        }
                        blockList.add(b);
                    } catch (IOException e) {
                        errorCounter.incrementAndGet();
                    }
                }
            };

            service.submit(loader);
        }

        // done submitting all jobs
        service.shutdown();

        // wait for all to finish
        try {
            service.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
        } catch (InterruptedException e) {
            System.err.println("Error loading mzd data " + e.getLocalizedMessage());
            if (HiCGlobals.printVerboseComments) {
                e.printStackTrace();
            }
        }

        // error printing
        if (errorCounter.get() > 0) {
            System.err.println(errorCounter.get() + " errors while reading blocks");
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
     *
     * @param regionIndices
     * @return
     */
    List<Integer> getBlockNumbersForRegionFromGenomePosition(long[] regionIndices) {
        int resolution = zoom.getBinSize();
        long[] regionBinIndices = new long[4];
        for (int i = 0; i < regionBinIndices.length; i++) {
            regionBinIndices[i] = regionIndices[i] / resolution;
        }
        return getBlockNumbersForRegionFromBinPosition(regionBinIndices);
    }
    
    // todo V9 needs a diff method
    private List<Integer> getBlockNumbersForRegionFromBinPosition(long[] regionIndices) {
        
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
    
    
    public void dump(PrintWriter printWriter, LittleEndianOutputStream les, NormalizationType norm, MatrixType matrixType,
                     boolean useRegionIndices, long[] regionIndices, ExpectedValueFunction df, boolean dense) throws IOException {
        
        // determine which output will be used
        if (printWriter == null && les == null) {
            printWriter = new PrintWriter(System.out);
        }
        boolean usePrintWriter = printWriter != null && les == null;
        boolean isIntraChromosomal = chr1.getIndex() == chr2.getIndex();
        
        // Get the block index keys, and sort
        List<Integer> blocksToIterateOver;
        if (useRegionIndices) {
            blocksToIterateOver = getBlockNumbersForRegionFromGenomePosition(regionIndices);
        } else {
            blocksToIterateOver = reader.getBlockNumbers(this);
            Collections.sort(blocksToIterateOver);
        }

        if (!dense) {
            for (Integer blockNumber : blocksToIterateOver) {
                Block b = reader.readNormalizedBlock(blockNumber, MatrixZoomData.this, norm);
                if (b != null) {
                    for (ContactRecord rec : b.getContactRecords()) {
                        float counts = rec.getCounts();
                        int x = rec.getBinX();
                        int y = rec.getBinY();
                        int xActual = x * zoom.getBinSize();
                        int yActual = y * zoom.getBinSize();
                        float oeVal = 0f;
                        if (matrixType == MatrixType.OE) {
                            double expected = 0;
                            if (chr1 == chr2) {
                                if (df != null) {
                                    int dist = Math.abs(x - y);
                                    expected = df.getExpectedValue(chr1.getIndex(), dist);
                                }
                            } else {
                                expected = (averageCount > 0 ? averageCount : 1);
                            }

                            double observed = rec.getCounts(); // Observed is already normalized
                            oeVal = (float) (observed / expected);
                        }
                        if (!useRegionIndices || // i.e. use full matrix
                                // or check regions that overlap with upper left
                                (xActual >= regionIndices[0] && xActual <= regionIndices[1] &&
                                        yActual >= regionIndices[2] && yActual <= regionIndices[3]) ||
                                // or check regions that overlap with lower left
                                (isIntraChromosomal && yActual >= regionIndices[0] && yActual <= regionIndices[1] &&
                                        xActual >= regionIndices[2] && xActual <= regionIndices[3])) {
                            // but leave in upper right triangle coordinates
                            if (usePrintWriter) {
                                if (matrixType == MatrixType.OBSERVED) {
                                    printWriter.println(xActual + "\t" + yActual + "\t" + counts);
                                } else if (matrixType == MatrixType.OE) {
                                    printWriter.println(xActual + "\t" + yActual + "\t" + oeVal);
                                }
                            } else {
                                // TODO I suspect this is wrong - should be writing xActual - but this is for binary dumping and we never use it
                                if (matrixType == MatrixType.OBSERVED) {
                                    les.writeInt(x);
                                    les.writeInt(y);
                                    les.writeFloat(counts);
                                } else if (matrixType == MatrixType.OE) {
                                    les.writeInt(x);
                                    les.writeInt(y);
                                    les.writeFloat(oeVal);
                                }
                            }
                        }
                    }
                }
            }

            if (usePrintWriter) {
                printWriter.close();
            }
            else {
                les.close();
            }
        }
        else {
            int maxX = 0;
            int maxY = 0;
            for (Integer blockNumber : blocksToIterateOver) {
                Block b = reader.readNormalizedBlock(blockNumber, MatrixZoomData.this, norm);
                if (b != null) {
                    for (ContactRecord rec : b.getContactRecords()) {
                        int x = rec.getBinX();
                        int y = rec.getBinY();
                        if (maxX < x) maxX = x;
                        if (maxY < y) maxY = y;
                    }
                }
            }
            if (isIntraChromosomal) {
                if (maxX < maxY) {
                    maxX = maxY;
                } else {
                    maxY = maxX;
                }
            }

            maxX++;
            maxY++;
            float[][] matrix = new float[maxX][maxY];  // auto initialized to 0

            for (Integer blockNumber : blocksToIterateOver) {
                Block b = reader.readNormalizedBlock(blockNumber, MatrixZoomData.this, norm);
                if (b != null) {
                    for (ContactRecord rec : b.getContactRecords()) {
                        float counts = rec.getCounts();
                        int x = rec.getBinX();
                        int y = rec.getBinY();

                        int xActual = x * zoom.getBinSize();
                        int yActual = y * zoom.getBinSize();
                        float oeVal = 0f;
                        if (matrixType == MatrixType.OE) {
                            int dist = Math.abs(x - y);
                            double expected = 0;
                            try {
                                expected = df.getExpectedValue(chr1.getIndex(), dist);
                            } catch (Exception e) {
                                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                            }
                            double observed = rec.getCounts(); // Observed is already normalized
                            oeVal = (float) (observed / expected);
                        }
                        if (!useRegionIndices || // i.e. use full matrix
                                // or check regions that overlap with upper left
                                (xActual >= regionIndices[0] && xActual <= regionIndices[1] &&
                                        yActual >= regionIndices[2] && yActual <= regionIndices[3]) ||
                                // or check regions that overlap with lower left
                                (isIntraChromosomal && yActual >= regionIndices[0] && yActual <= regionIndices[1] &&
                                        xActual >= regionIndices[2] && xActual <= regionIndices[3])) {

                            if (matrixType == MatrixType.OBSERVED) {
                                matrix[x][y] = counts;
                                if (isIntraChromosomal) {
                                    matrix[y][x] = counts;
                                }
                                // printWriter.println(xActual + "\t" + yActual + "\t" + counts);
                            } else if (matrixType == MatrixType.OE) {
                                matrix[x][y] = oeVal;
                                if (isIntraChromosomal) {
                                    matrix[y][x] = oeVal;
                                }
                                // printWriter.println(xActual + "\t" + yActual + "\t" + oeVal);
                            }
                        }
                    }
                }
            }
            if (usePrintWriter) {
                for (int i = 0; i < maxX; i++) {
                    for (int j = 0; j < maxY; j++) {
                        printWriter.print(matrix[i][j] + "\t");
                    }
                    printWriter.println();
                }
            } else {
                for (int i = 0; i < maxX; i++) {
                    for (int j = 0; j < maxY; j++) {
                        les.writeFloat(matrix[i][j]);

                    }
    
                }
            }
    
            if (usePrintWriter) {
                printWriter.close();
            } else {
                les.close();
            }
        }
    }
    
    public void dump1DTrackFromCrossHairAsWig(PrintWriter printWriter, long binStartPosition,
                                              boolean isIntraChromosomal, long[] regionBinIndices,
                                              NormalizationType norm, MatrixType matrixType) {
        
        if (!MatrixType.isObservedOrControl(matrixType)) {
            System.out.println("This feature is only available for Observed or Control views");
            return;
        }
        
        int binCounter = 0;
        
        // Get the block index keys, and sort
        List<Integer> blocksToIterateOver = getBlockNumbersForRegionFromBinPosition(regionBinIndices);
        Collections.sort(blocksToIterateOver);

        for (Integer blockNumber : blocksToIterateOver) {
            Block b = null;
            try {
                b = reader.readNormalizedBlock(blockNumber, MatrixZoomData.this, norm);
            } catch (Exception e) {
                System.err.println("Skipping block " + blockNumber);
            }
            if (b != null) {
                for (ContactRecord rec : b.getContactRecords()) {
                    float counts = rec.getCounts();
                    int x = rec.getBinX();
                    int y = rec.getBinY();

                    if (    //check regions that overlap with upper left
                            (x >= regionBinIndices[0] && x <= regionBinIndices[1] &&
                                    y >= regionBinIndices[2] && y <= regionBinIndices[3]) ||
                                    // or check regions that overlap with lower left
                                    (isIntraChromosomal && x >= regionBinIndices[2] && x <= regionBinIndices[3] &&
                                            y >= regionBinIndices[0] && y <= regionBinIndices[1])) {
                        // but leave in upper right triangle coordinates

                        if (x == binStartPosition) {
                            while (binCounter < y) {
                                printWriter.println("0");
                                binCounter++;
                            }
                        } else if (y == binStartPosition) {
                            while (binCounter < x) {
                                printWriter.println("0");
                                binCounter++;
                            }
                        } else {
                            System.err.println("Something went wrong while generating 1D track");
                            System.err.println("Improper input was likely provided");
                        }

                        printWriter.println(counts);
                        binCounter++;

                    }
                }
            }
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

    public List<List<ContactRecord>> getContactRecordList() {
        if (localCacheOfRecords == null || localCacheOfRecords.size() < 1) {
            numberOfContactRecords = 0;
            localCacheOfRecords = new ArrayList<>(3);
            List<ContactRecord> currentList = new ArrayList<>(1000000);
            int localCounter = 0;
            int maxAllowed = 9 * (Integer.MAX_VALUE / 10);
            Iterator<ContactRecord> iterator = new ContactRecordIterator(reader, this, blockCache);
            while (iterator.hasNext()) {
                ContactRecord cr = iterator.next();
                currentList.add(cr);
                numberOfContactRecords++;
                localCounter++;
                if (localCounter > maxAllowed) {
                    // make new ArrayList
                    localCacheOfRecords.add(currentList);
                    localCounter = 0;
                    currentList = new ArrayList<>();
                }
            }
            if (localCounter > 0) { // add last list
                localCacheOfRecords.add(currentList);
            }
        }
        return localCacheOfRecords;
    }

    public void clearCache() {
        blockCache.clear();
    }

    public long getNumberOfContactRecords() {
        if (numberOfContactRecords == 0) {
            getContactRecordList();
        }
        return numberOfContactRecords;
    }
}
