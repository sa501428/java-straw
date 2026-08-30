package javastraw.reader.v10;

import javastraw.reader.datastructures.ListOfDoubleArrays;
import javastraw.reader.expected.ExpectedValueFunction;
import javastraw.reader.expected.ExpectedValueFunctionImpl;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationType;

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

/**
 * Expected-value vector backed by a V10 EVI0 or NEVI entry.
 * <p>
 * Values are fetched from the chunked, Zstandard-compressed store on demand, so
 * a viewport query never has to materialize a whole chromosome-scale vector.
 * Per-chromosome scale factors are applied exactly as Section J.7 requires:
 * {@code effectiveExpected(d) = expectedVector[d] / scaleFactor}.
 */
public class V10ExpectedValueFunction implements ExpectedValueFunction {

    private static final int WINDOW = 500000;

    private final NormalizationType type;
    private final HiCZoom.HiCUnit unit;
    private final int binSize;
    private final V10VectorEntry entry;
    private final V10Source source;

    private double[] window;
    private long windowStart = -1;
    private long windowEnd = -1;

    public V10ExpectedValueFunction(NormalizationType type, HiCZoom.HiCUnit unit, int binSize,
                                    V10VectorEntry entry, V10Source source) {
        this.type = type;
        this.unit = unit;
        this.binSize = binSize;
        this.entry = entry;
        this.source = source;
    }

    public V10VectorEntry getEntry() {
        return entry;
    }

    /**
     * Chromosome scale factor; 1.0 when the file records none for this chromosome.
     */
    public double getScaleFactor(int chrIdx) {
        return entry.getScaleFactor(chrIdx);
    }

    @Override
    public double getExpectedValue(int chrIdx, long distance) {
        double scale = entry.getScaleFactor(chrIdx);
        long length = entry.vectorValueCount;
        if (length <= 0) {
            System.err.println("Expected values array is empty");
            return -1;
        }
        long index = Math.min(Math.max(distance, 0), length - 1);
        return valueAt(index) / scale;
    }

    private synchronized double valueAt(long index) {
        if (window == null || index < windowStart || index >= windowEnd) {
            long start = Math.max(0, index - WINDOW / 2);
            long end = Math.min(entry.vectorValueCount, start + WINDOW);
            try {
                window = entry.read(source, start, end);
            } catch (IOException e) {
                throw new V10FormatException("could not read expected vector", e);
            }
            windowStart = start;
            windowEnd = start + window.length;
        }
        return window[(int) (index - windowStart)];
    }

    @Override
    public long getLength() {
        return entry.vectorValueCount;
    }

    @Override
    public NormalizationType getNormalizationType() {
        return type;
    }

    @Override
    public HiCZoom.HiCUnit getUnit() {
        return unit;
    }

    @Override
    public int getBinSize() {
        return binSize;
    }

    @Override
    public ListOfDoubleArrays getExpectedValuesNoNormalization() {
        return readAll();
    }

    @Override
    public ListOfDoubleArrays getExpectedValuesWithNormalization(int chrIdx) {
        ListOfDoubleArrays values = readAll();
        double scale = entry.getScaleFactor(chrIdx);
        if (scale != 1.0) {
            values.multiplyEverythingBy(1.0 / scale);
        }
        return values;
    }

    private ListOfDoubleArrays readAll() {
        long n = entry.vectorValueCount;
        ListOfDoubleArrays values = new ListOfDoubleArrays(n);
        long position = 0;
        try {
            while (position < n) {
                long end = Math.min(n, position + WINDOW);
                double[] part = entry.read(source, position, end);
                for (int i = 0; i < part.length; i++) {
                    values.set(position + i, part[i]);
                }
                position = end;
            }
        } catch (IOException e) {
            throw new V10FormatException("could not read expected vector", e);
        }
        return values;
    }

    @Override
    public ExpectedValueFunction getCorrectedVersion(int window) {
        ListOfDoubleArrays smoothed = readAll();
        smoothed.doRollingMedian(window / binSize);
        Map<Integer, Double> normFactors = new HashMap<>();
        for (Map.Entry<Integer, Float> factor : entry.scaleFactors.entrySet()) {
            normFactors.put(factor.getKey(), (double) factor.getValue());
        }
        return new ExpectedValueFunctionImpl(type, unit, binSize, smoothed, normFactors);
    }
}
