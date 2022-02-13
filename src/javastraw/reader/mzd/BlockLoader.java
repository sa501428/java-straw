package javastraw.reader.mzd;

import javastraw.reader.DatasetReader;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.block.Block;
import javastraw.reader.block.BlockModifier;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.ParallelizationTools;

import java.io.IOException;
import java.util.List;
import java.util.Set;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.atomic.AtomicInteger;

public class BlockLoader {
    public static void actuallyLoadGivenBlocks(final List<Block> blockList, Set<Integer> blocksToLoad,
                                               final NormalizationType no, BlockModifier modifier,
                                               final String zdKey,
                                               Chromosome chrom1, Chromosome chrom2, HiCZoom zoom,
                                               BlockCache blockCache, DatasetReader reader) {
        final AtomicInteger errorCounter = new AtomicInteger();

        ExecutorService service = Executors.newFixedThreadPool(200);
        for (final int blockNumber : blocksToLoad) {
            String key = getBlockKey(zdKey, blockNumber, no);
            readBlockUpdateListAndCache(blockNumber, reader, no, blockList, key, errorCounter, service, modifier,
                    chrom1, chrom2, zoom, zdKey, blockCache);
        }
        ParallelizationTools.shutDownServiceAndWait(service, errorCounter);
    }

    public static String getBlockKey(String zdKey, int blockNumber, NormalizationType no) {
        return MatrixZoomData.triKey(zdKey, "" + blockNumber, "" + no);
    }

    private static void readBlockUpdateListAndCache(int blockNumber, DatasetReader reader, NormalizationType no,
                                                    List<Block> blockList, String key, final AtomicInteger errorCounter,
                                                    ExecutorService service, BlockModifier modifier,
                                                    Chromosome chrom1, Chromosome chrom2, HiCZoom zoom,
                                                    String zdKey, BlockCache blockCache) {
        Runnable loader = () -> {
            try {
                Block b = reader.readNormalizedBlock(blockNumber, zdKey, no,
                        chrom1.getIndex(), chrom2.getIndex(), zoom);
                if (b == null) {
                    b = new Block(blockNumber, key);
                }
                b = modifier.modify(b, key, zoom.getBinSize(), chrom1, chrom2);
                blockCache.put(key, b);
                blockList.add(b);
            } catch (IOException e) {
                errorCounter.incrementAndGet();
            }
        };
        service.submit(loader);
    }
}
