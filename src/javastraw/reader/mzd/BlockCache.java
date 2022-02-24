package javastraw.reader.mzd;

import javastraw.reader.block.Block;
import org.broad.igv.util.collections.LRUCache;

public class BlockCache {
    LRUCache<String, Block> cache = new LRUCache<>(500);
    private boolean useCache = true;

    public void clear() {
        cache.clear();
    }

    public void put(String key, Block b) {
        if (useCache) {
            cache.put(key, b);
        }
    }

    public boolean containsKey(String key) {
        return useCache && cache.containsKey(key);
    }

    public Block get(String key) {
        return cache.get(key);
    }

    public void setUseCache(boolean useCache) {
        this.useCache = useCache;
        if (!useCache) {
            cache.clear();
        }
    }

    public void addAll(BlockCache other) {
        for (String key : other.cache.keySet()) {
            cache.put(key, other.cache.get(key));
        }
    }
}
