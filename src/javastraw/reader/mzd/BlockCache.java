package javastraw.reader.mzd;

import javastraw.reader.block.Block;
import org.broad.igv.util.collections.LRUCache;

import java.util.HashSet;
import java.util.Set;

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
        if (useCache) {
            return cache.get(key);
        }
        System.err.println("Invalid situation - cache is not set");
        System.exit(9);
        return null;
    }

    public void setUseCache(boolean useCache) {
        this.useCache = useCache;
        if (!useCache) {
            cache.clear();
        }
    }

    public void addAll(BlockCache other) {
        if (useCache) {
            Set<String> keys = new HashSet<>(other.cache.keySet());
            for (String key : keys) {
                cache.put(key, other.cache.get(key));
            }
        }
    }

    public boolean getUseCache() {
        return useCache;
    }
}
