package javastraw.reader.v10;

import javastraw.reader.basics.Chromosome;
import java.util.ArrayList;
import java.util.List;

/** Optional scaled All convention, represented entirely by existing attributes. */
public final class V10AllOverview {
    public final int chromosomeIndex, finestResolution, binCount, sourceResolution;
    public final long groupFactor, effectiveResolution;
    public final int[] chromosomeIndices;
    public final long[] sourceOffsets, chromosomeBins;

    private V10AllOverview(V10Header h) {
        List<Chromosome> real = new ArrayList<>();
        int all = -1;
        for (Chromosome c : h.chromosomes) {
            if (c.getName().equalsIgnoreCase("all")) all = c.getIndex(); else real.add(c);
        }
        chromosomeIndex = all;
        finestResolution = Math.toIntExact(number(h, "finestResolution"));
        sourceResolution = Math.toIntExact(number(h, "sourceResolution"));
        groupFactor = number(h, "groupFactor");
        binCount = Math.toIntExact(number(h, "binCount"));
        effectiveResolution = Math.multiplyExact(sourceResolution, groupFactor);
        V10.require(all >= 0 && !h.resolutions[0].isEmpty() &&
                finestResolution == h.resolutions[0].get(0).binSize &&
                h.chromosomes.get(all).getLength() == Math.multiplyExact((long) binCount, finestResolution),
                "inconsistent scaled All chromosome");
        String[] bins = h.attributes.getOrDefault("hicAll.chromosomeBins", "").split(",");
        V10.require(bins.length == real.size(), "inconsistent scaled All chromosome layout");
        chromosomeIndices = new int[real.size()];
        sourceOffsets = new long[real.size()];
        chromosomeBins = new long[real.size()];
        long total = 0;
        for (int i = 0; i < real.size(); i++) {
            Chromosome c = real.get(i);
            long n = Long.parseLong(bins[i]);
            V10.require(n > 0 && n == c.getLength() / sourceResolution + (c.getLength() % sourceResolution != 0 ? 1 : 0),
                    "inconsistent scaled All chromosome bins");
            chromosomeIndices[i] = c.getIndex(); sourceOffsets[i] = total; chromosomeBins[i] = n;
            total = Math.addExact(total, n);
        }
        V10.require(total / groupFactor + (total % groupFactor != 0 ? 1 : 0) == binCount,
                "inconsistent scaled All bin count");
    }
    private static long number(V10Header h, String key) {
        String text = h.attributes.get("hicAll." + key);
        V10.require(text != null && text.matches("[0-9]+"), "invalid scaled All " + key);
        long n = Long.parseLong(text); V10.require(n > 0, "invalid scaled All " + key); return n;
    }
    public static V10AllOverview from(V10Header h) {
        return "1".equals(h.attributes.get("hicAll.version")) ? new V10AllOverview(h) : null;
    }
    public double syntheticBoundary(int ordinal) {
        return (sourceOffsets[ordinal] + chromosomeBins[ordinal]) * (double) finestResolution / groupFactor;
    }
    public int chromosomeForCoordinate(double coordinate) {
        double bin = coordinate * groupFactor / finestResolution;
        for (int i = 0; i < chromosomeIndices.length; i++)
            if (bin < sourceOffsets[i] + chromosomeBins[i]) return chromosomeIndices[i];
        return chromosomeIndices[chromosomeIndices.length - 1];
    }
    public long positionForCoordinate(int chromosome, double coordinate) {
        for (int i = 0; i < chromosomeIndices.length; i++) if (chromosomeIndices[i] == chromosome)
            return Math.max(0L, (long) ((coordinate * groupFactor / finestResolution - sourceOffsets[i]) * sourceResolution));
        throw new IllegalArgumentException("Not a real chromosome in the overview");
    }
}
