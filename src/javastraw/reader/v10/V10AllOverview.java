package javastraw.reader.v10;

import javastraw.reader.basics.Chromosome;
import java.util.ArrayList;
import java.util.List;

/** Optional scaled All convention, represented entirely by existing attributes. */
public final class V10AllOverview {
    public final int version, chromosomeIndex, finestResolution, binCount, sourceResolution, displayResolution;
    public final long coordinateScale, storageStride, syntheticLength;
    public final long[] chromosomeLengths;
    public final long groupFactor, effectiveResolution;
    public final int[] chromosomeIndices;
    public final long[] sourceOffsets, chromosomeBins;

    private V10AllOverview(V10Header h) {
        List<Chromosome> real = new ArrayList<>();
        int all = -1;
        for (Chromosome c : h.chromosomes) {
            if (c.getName().equalsIgnoreCase("all")) all = c.getIndex(); else real.add(c);
        }
        version = Integer.parseInt(h.attributes.get("hicAll.version"));
        chromosomeIndex = all;
        syntheticLength = all < 0 ? 0 : h.chromosomes.get(all).getLength();
        finestResolution = Math.toIntExact(number(h, "finestResolution"));
        sourceResolution = Math.toIntExact(number(h, "sourceResolution"));
        groupFactor = number(h, "groupFactor");
        binCount = Math.toIntExact(number(h, "binCount"));
        coordinateScale = version == 2 ? number(h, "coordinateScale") : 0;
        storageStride = version == 2 ? number(h, "storageStride") : 1;
        displayResolution = Math.toIntExact(Math.multiplyExact(finestResolution, storageStride));
        effectiveResolution = version == 2 ? Math.multiplyExact(displayResolution, coordinateScale)
                : Math.multiplyExact(sourceResolution, groupFactor);
        V10.require(all >= 0 && !h.resolutions[0].isEmpty() &&
                finestResolution == h.resolutions[0].get(0).binSize &&
                (version == 2 || syntheticLength == Math.multiplyExact((long) binCount, finestResolution)),
                "inconsistent scaled All chromosome");
        String[] bins = h.attributes.getOrDefault("hicAll.chromosomeBins", "").split(",");
        V10.require(bins.length == real.size(), "inconsistent scaled All chromosome layout");
        chromosomeIndices = new int[real.size()];
        sourceOffsets = new long[real.size()];
        chromosomeBins = new long[real.size()];
        chromosomeLengths = new long[real.size()];
        long total = 0;
        for (int i = 0; i < real.size(); i++) {
            Chromosome c = real.get(i);
            long n = Long.parseLong(bins[i]);
            V10.require(n > 0 && n == c.getLength() / sourceResolution + (c.getLength() % sourceResolution != 0 ? 1 : 0),
                    "inconsistent scaled All chromosome bins");
            chromosomeIndices[i] = c.getIndex(); sourceOffsets[i] = total; chromosomeBins[i] = n;
            chromosomeLengths[i] = c.getLength();
            total = Math.addExact(total, version == 2 ? c.getLength() : n);
        }
        V10.require(version == 2
                ? syntheticLength == Math.max(1, total / coordinateScale) &&
                  binCount == syntheticLength / displayResolution + (syntheticLength % displayResolution != 0 ? 1 : 0)
                : total / groupFactor + (total % groupFactor != 0 ? 1 : 0) == binCount,
                "inconsistent scaled All bin count");
    }
    private static long number(V10Header h, String key) {
        String text = h.attributes.get("hicAll." + key);
        V10.require(text != null && text.matches("[0-9]+"), "invalid scaled All " + key);
        long n = Long.parseLong(text); V10.require(n > 0, "invalid scaled All " + key); return n;
    }
    public static V10AllOverview from(V10Header h) {
        String version = h.attributes.get("hicAll.version");
        return "1".equals(version) || "2".equals(version) ? new V10AllOverview(h) : null;
    }
    public double syntheticBoundary(int ordinal) {
        if (version == 2) return Math.min(syntheticLength, (sourceOffsets[ordinal] + chromosomeLengths[ordinal]) / (double) coordinateScale);
        return (sourceOffsets[ordinal] + chromosomeBins[ordinal]) * (double) finestResolution / groupFactor;
    }
    public int chromosomeForCoordinate(double coordinate) {
        double bin = version == 2 ? coordinate * coordinateScale : coordinate * groupFactor / finestResolution;
        for (int i = 0; i < chromosomeIndices.length; i++)
            if (bin < sourceOffsets[i] + (version == 2 ? chromosomeLengths[i] : chromosomeBins[i])) return chromosomeIndices[i];
        return chromosomeIndices[chromosomeIndices.length - 1];
    }
    public long positionForCoordinate(int chromosome, double coordinate) {
        for (int i = 0; i < chromosomeIndices.length; i++) if (chromosomeIndices[i] == chromosome)
            if (version == 2) return Math.max(0L, (long) (coordinate * coordinateScale - sourceOffsets[i]));
            else return Math.max(0L, (long) ((coordinate * groupFactor / finestResolution - sourceOffsets[i]) * sourceResolution));
        throw new IllegalArgumentException("Not a real chromosome in the overview");
    }
}
