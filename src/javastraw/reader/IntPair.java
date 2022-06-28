package javastraw.reader;

import java.util.Objects;

public class IntPair {
    public final int a, b;

    public IntPair(int a, int b) {
        this.a = a;
        this.b = b;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o instanceof IntPair) {
            IntPair other = (IntPair) o;
            return a == other.a && b == other.b;
        }
        return false;
    }

    @Override
    public int hashCode() {
        return Objects.hash(a, b);
    }

}
