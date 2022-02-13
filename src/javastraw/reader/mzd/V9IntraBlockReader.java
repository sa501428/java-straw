package javastraw.reader.mzd;

import javastraw.reader.DatasetReader;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.block.Block;
import javastraw.reader.block.BlockModifier;
import javastraw.reader.depth.V9Depth;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationType;

import java.util.*;

public class V9IntraBlockReader {
    public static List<Block> addNormalizedBlocksToListV9(final List<Block> blockList, int binX1, int binY1, int binX2, int binY2,
                                                          final NormalizationType norm, BlockModifier modifier,
                                                          int blockBinCount, V9Depth v9Depth,
                                                          int blockColumnCount, BlockCache blockCache, String zdKey,
                                                          Chromosome chrom1, Chromosome chrom2, HiCZoom zoom,
                                                          DatasetReader reader) {

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
                populateBlocksToLoadV9(pad, depth, norm, blockList, blocksToLoad, blockColumnCount, blockCache, zdKey);
            }
        }

        BlockLoader.actuallyLoadGivenBlocks(blockList, blocksToLoad, norm, modifier, zdKey,
                chrom1, chrom2, zoom, blockCache, reader);

        return blockList;
    }


    protected static void populateBlocksToLoadV9(int positionAlongDiagonal, int depth, NormalizationType no,
                                                 List<Block> blockList, Set<Integer> blocksToLoad,
                                                 int blockColumnCount, BlockCache blockCache, String zdKey) {
        int blockNumber = getBlockNumberVersion9FromPADAndDepth(positionAlongDiagonal, depth, blockColumnCount);
        String key = BlockLoader.getBlockKey(zdKey, blockNumber, no);
        Block b;
        if (blockCache.containsKey(key)) {
            b = blockCache.get(key);
            blockList.add(b);
        } else {
            blocksToLoad.add(blockNumber);
        }
    }

    public static int getBlockNumberVersion9FromPADAndDepth(int positionAlongDiagonal, int depth, int blockColumnCount) {
        return depth * blockColumnCount + positionAlongDiagonal;
    }

    public static List<Integer> getBlockNumbersForRegionFromBinPosition(long[] regionBinIndices, int blockBinCount,
                                                                        int blockColumnCount, boolean isIntra, V9Depth v9Depth) {
        long binX1 = regionBinIndices[0];
        long binX2 = regionBinIndices[1];
        long binY1 = regionBinIndices[2];
        long binY2 = regionBinIndices[3];
        Set<Integer> blocksSet = new HashSet<>();

        int translatedLowerPAD = (int) ((binX1 + binY1) / 2 / blockBinCount);
        int translatedHigherPAD = (int) ((binX2 + binY2) / 2 / blockBinCount + 1);
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
                blocksSet.add(getBlockNumberVersion9FromPADAndDepth(pad, depth, blockColumnCount));
            }
        }

        List<Integer> blocksToIterateOver = new ArrayList<>(blocksSet);
        Collections.sort(blocksToIterateOver);
        return blocksToIterateOver;
    }
}
