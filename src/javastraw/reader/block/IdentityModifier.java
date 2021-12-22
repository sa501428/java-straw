package javastraw.reader.block;

public class IdentityModifier implements BlockModifier {
    @Override
    public Block modify(Block b) {
        return b;
    }
}
