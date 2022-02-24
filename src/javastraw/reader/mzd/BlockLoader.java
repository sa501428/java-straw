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
import java.util.List;
import java.util.Set;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.atomic.AtomicInteger;

public class BlockLoader {
    public static void actuallyLoadGivenBlocks(final List<Block> globalBlockList, Set<Integer> blocksToLoad,
                                               final NormalizationType no, BlockModifier modifier,
                                               final String zdKey, Chromosome chrom1, Chromosome chrom2, HiCZoom zoom,
                                               BlockCache globalBlockCache, DatasetReader reader) {
        final AtomicInteger errorCounter = new AtomicInteger();

        int numJobs = Math.min(blocksToLoad.size(), 200);
        List<Integer> blockIds = new ArrayList<>(blocksToLoad);
        AtomicInteger index = new AtomicInteger(0);

        ParallelizationTools.launchParallelizedCode(numJobs, () -> {
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

            synchronized (globalBlockList) {
                globalBlockList.addAll(blockList);
            }
            blockList.clear();
            synchronized (globalBlockCache) {
                globalBlockCache.addAll(blockCache);
            }
            blockCache.clear();
        });
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

    private static void readBlockUpdateListAndCache(int blockNumber, DatasetReader reader, NormalizationType no,
                                                    List<Block> blockList, String key, final AtomicInteger errorCounter,
                                                    ExecutorService service, BlockModifier modifier,
                                                    Chromosome chrom1, Chromosome chrom2, HiCZoom zoom,
                                                    String zdKey, BlockCache blockCache) {
        Runnable loader = () -> {
            try {
                getBlockFromReader(blockList, no, modifier, zdKey, chrom1, chrom2, zoom,
                        blockCache, reader, blockNumber, key);
            } catch (IOException e) {
                errorCounter.incrementAndGet();
            }
        };
        service.submit(loader);
    }

    public static String getBlockKey(String zdKey, int blockNumber, NormalizationType no) {
        return MatrixZoomData.triKey(zdKey, "" + blockNumber, "" + no);
    }
}
