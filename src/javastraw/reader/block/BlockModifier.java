package javastraw.reader.block;

import javastraw.reader.basics.Chromosome;

public interface BlockModifier {
    Block modify(Block b, String key, int binSize, Chromosome chr1, Chromosome chr2);
}
