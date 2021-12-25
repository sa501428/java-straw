package javastraw.reader.block;

import javastraw.reader.basics.Chromosome;

public class IdentityModifier implements BlockModifier {
    @Override
    public Block modify(Block b, String key, int binSize, Chromosome chr1, Chromosome chr2) {
        return b;
    }
}
