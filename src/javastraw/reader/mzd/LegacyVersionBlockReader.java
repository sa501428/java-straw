package javastraw.reader.mzd;

import javastraw.reader.DatasetReader;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.block.Block;
import javastraw.reader.block.BlockModifier;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationType;
import org.broad.igv.util.collections.LRUCache;

import java.util.*;

public class LegacyVersionBlockReader {

    /**
     * Return the blocks of normalized, observed values overlapping the rectangular region specified.
     *
     * @param binY1         leftmost position in "bins"
     * @param binX2         rightmost position in "bins"
     * @param binY2         bottom position in "bins"
     * @param norm          normalization type
     * @param blockBinCount
     * @return List of overlapping blocks, normalized
     */
    public static List<Block> addNormalizedBlocksToList(final List<Block> blockList, int binX1, int binY1,
                                                        int binX2, int binY2, final NormalizationType norm,
                                                        boolean getBelowDiagonal, BlockModifier modifier,
                                                        int blockBinCount, int blockColumnCount, boolean useCache,
                                                        LRUCache<String, Block> blockCache, String zdKey,
                                                        Chromosome chrom1, Chromosome chrom2, HiCZoom zoom,
                                                        DatasetReader reader) {

        Set<Integer> blocksToLoad = new HashSet<>();

        // have to do this regardless (just in case)
        int col1 = binX1 / blockBinCount;
        int row1 = binY1 / blockBinCount;
        int col2 = binX2 / blockBinCount;
        int row2 = binY2 / blockBinCount;

        for (int r = row1; r <= row2; r++) {
            for (int c = col1; c <= col2; c++) {
                populateBlocksToLoad(r, c, norm, blockList, blocksToLoad, blockColumnCount,
                        useCache, blockCache, zdKey);
            }
        }

        if (getBelowDiagonal && binY1 < binX2) {
            for (int r = row1; r <= row2; r++) {
                for (int c = col1; c <= col2; c++) {
                    populateBlocksToLoad(c, r, norm, blockList, blocksToLoad, blockColumnCount,
                            useCache, blockCache, zdKey);
                }
            }
        }

        BlockLoader.actuallyLoadGivenBlocks(blockList, blocksToLoad, norm, modifier, zdKey,
                chrom1, chrom2, zoom, useCache, blockCache, reader);

        return new ArrayList<>(new HashSet<>(blockList));
    }

    protected static void populateBlocksToLoad(int r, int c, NormalizationType no, List<Block> blockList,
                                               Set<Integer> blocksToLoad, int blockColumnCount,
                                               boolean useCache, LRUCache<String, Block> blockCache,
                                               String zdKey) {
        int blockNumber = r * blockColumnCount + c;
        String key = BlockLoader.getBlockKey(zdKey, blockNumber, no);
        Block b;
        if (useCache && blockCache.containsKey(key)) {
            b = blockCache.get(key);
            blockList.add(b);
        } else {
            blocksToLoad.add(blockNumber);
        }
    }

    public static List<Integer> getBlockNumbersForRegionFromBinPosition(long[] regionIndices, int blockBinCount,
                                                                        int blockColumnCount, boolean isIntra) {

        // cast should be fine - this is for V8
        int row1 = (int) (regionIndices[0] / blockBinCount);
        int row2 = (int) ((regionIndices[1] + 1) / blockBinCount);
        int col1 = (int) (regionIndices[2] / blockBinCount);
        int col2 = (int) ((regionIndices[3] + 1) / blockBinCount);

        // first check the upper triangular matrix
        Set<Integer> blocksSet = new HashSet<>();
        for (int r = row1; r <= row2; r++) {
            for (int c = col1; c <= col2; c++) {
                int blockNumber = r * blockColumnCount + c;
                blocksSet.add(blockNumber);
            }
        }
        // check region part that overlaps with lower left triangle
        // but only if intrachromosomal
        if (isIntra) {
            for (int r = col1; r <= col2; r++) {
                for (int c = row1; c <= row2; c++) {
                    int blockNumber = r * blockColumnCount + c;
                    blocksSet.add(blockNumber);
                }
            }
        }

        List<Integer> blocksToIterateOver = new ArrayList<>(blocksSet);
        Collections.sort(blocksToIterateOver);
        return blocksToIterateOver;
    }
}
