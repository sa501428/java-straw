package javastraw.reader.basics;

import java.util.Arrays;
import java.util.Objects;

public class ChromosomeArrayPair {
    public final Chromosome[] a, b;

    public ChromosomeArrayPair(Chromosome[] a, Chromosome[] b) {
        this.a = a;
        this.b = b;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o instanceof ChromosomeArrayPair) {
            ChromosomeArrayPair other = (ChromosomeArrayPair) o;

            if (a.length != other.a.length || b.length != other.b.length) {
                return false;
            }

            for (int i = 0; i < a.length; i++) {
                if (notEqual(a[i], other.a[i])) {
                    return false;
                }
            }

            for (int i = 0; i < b.length; i++) {
                if (notEqual(b[i], other.b[i])) {
                    return false;
                }
            }
        }
        return true;
    }

    private boolean notEqual(Chromosome c, Chromosome d) {
        boolean areEqual = c.getName().equalsIgnoreCase(d.getName()) &&
                c.getIndex() == d.getIndex() &&
                c.getLength() == d.getLength();
        return !areEqual;
    }

    @Override
    public int hashCode() {
        return Objects.hash(Arrays.hashCode(a), Arrays.hashCode(b));
    }

}
