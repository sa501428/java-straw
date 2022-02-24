package javastraw.reader.mzd;

import javastraw.reader.DatasetReader;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.block.Block;
import javastraw.reader.block.BlockModifier;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.ParallelizationTools;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.concurrent.atomic.AtomicInteger;

public class BlockLoader {
    public static void actuallyLoadGivenBlocks(final List<Block> globalBlockList, Set<Integer> blocksToLoad,
                                               final NormalizationType no, BlockModifier modifier,
                                               final String zdKey, Chromosome chrom1, Chromosome chrom2, HiCZoom zoom,
                                               BlockCache globalBlockCache, DatasetReader reader) {
        final AtomicInteger errorCounter = new AtomicInteger();
        final Object listLock = new Object();
        final Object cacheLock = new Object();
        final Set<Block> globalBlockSet = new HashSet<>();

        int numJobs = Math.min(blocksToLoad.size(), 100);
        List<Integer> blockIds = new ArrayList<>(blocksToLoad);
        AtomicInteger index = new AtomicInteger(0);

        ParallelizationTools.launchParallelizedCode(1, () -> { //numJobs
            List<Block> blockList = new ArrayList<>();
            BlockCache blockCache = new BlockCache();

            int i = index.getAndIncrement();
            while (i < (blockIds).size()) {
                int blockNumber = blockIds.get(i);
                String key = getBlockKey(zdKey, blockNumber, no);
                try {
                    getBlockFromReader(blockList, no, modifier, zdKey, chrom1, chrom2, zoom, blockCache,
                            reader, blockNumber, key);
                } catch (IOException e) {
                    errorCounter.incrementAndGet();
                }
                i = index.getAndIncrement();
            }

            synchronized (listLock) {
                globalBlockSet.addAll(blockList);
            }
            blockList.clear();
            synchronized (cacheLock) {
                globalBlockCache.addAll(blockCache);
            }
            blockCache.clear();
        });

        globalBlockList.addAll(globalBlockSet);
    }

    private static void getBlockFromReader(List<Block> blockList, NormalizationType no, BlockModifier modifier,
                                           String zdKey, Chromosome chrom1, Chromosome chrom2, HiCZoom zoom,
                                           BlockCache blockCache, DatasetReader reader,
                                           int blockNumber, String key) throws IOException {
        Block b = reader.readNormalizedBlock(blockNumber, zdKey, no,
                chrom1.getIndex(), chrom2.getIndex(), zoom);
        if (b == null) {
            b = new Block(blockNumber, key);
        }
        b = modifier.modify(b, key, zoom.getBinSize(), chrom1, chrom2);
        blockCache.put(key, b);
        blockList.add(b);
    }

    public static String getBlockKey(String zdKey, int blockNumber, NormalizationType no) {
        return MatrixZoomData.triKey(zdKey, "" + blockNumber, "" + no);
    }
}
