package javastraw.reader.v10;

import javastraw.reader.basics.Chromosome;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import static javastraw.reader.v10.V10.require;

/**
 * The V10 fixed prefix plus variable header (Sections C.1-C.4).
 */
public class V10Header {

    public long headerByteLength;
    public V10Locator footer;
    public V10Locator normVectorIndex;
    public V10Locator expectedValueIndex;
    public V10Locator normExpectedValueIndex;

    public String genomeId;
    /**
     * Attribute pairs in stored order, duplicates preserved (Section C.2).
     */
    public final List<String[]> attributeList = new ArrayList<>();
    public final Map<String, String> attributes = new LinkedHashMap<>();

    public final List<Chromosome> chromosomes = new ArrayList<>();
    /**
     * Index 0 holds BP resolutions, index 1 holds FRAG resolutions.
     */
    @SuppressWarnings("unchecked")
    public final List<V10Resolution>[] resolutions = new List[]{
            new ArrayList<V10Resolution>(), new ArrayList<V10Resolution>()};
    /**
     * Number of restriction cut sites per chromosome; a chromosome with n sites
     * has n + 1 fragments.
     */
    public int[] fragmentSiteCounts;
    public final List<String> normalizations = new ArrayList<>();

    /**
     * Number of bins for a chromosome at a given unit and resolution index.
     */
    public long bins(int chrIndex, int unit, int resolutionIndex) {
        long length = unit == V10.UNIT_FRAG
                ? (long) fragmentSiteCounts[chrIndex] + 1
                : chromosomes.get(chrIndex).getLength();
        long binSize = resolutions[unit].get(resolutionIndex).binSize;
        return length / binSize + (length % binSize != 0 ? 1 : 0);
    }

    public int resolutionIndex(int unit, int binSize) {
        List<V10Resolution> list = resolutions[unit];
        for (int i = 0; i < list.size(); i++) {
            if (list.get(i).binSize == binSize) return i;
        }
        return -1;
    }

    public int totalResolutionCount() {
        return resolutions[V10.UNIT_BP].size() + resolutions[V10.UNIT_FRAG].size();
    }

    public int chromosomeIndex(String name) {
        for (int i = 0; i < chromosomes.size(); i++) {
            if (chromosomes.get(i).getName().equals(name)) return i;
        }
        return -1;
    }

    public int[] binSizes(int unit) {
        List<V10Resolution> list = resolutions[unit];
        int[] out = new int[list.size()];
        for (int i = 0; i < out.length; i++) {
            out[i] = list.get(i).binSize;
        }
        return out;
    }

    /**
     * Parses the complete header. {@code bytes} must be exactly
     * {@code headerByteLength} bytes starting at file position zero.
     */
    public static V10Header parse(byte[] bytes) {
        V10Cursor c = new V10Cursor(bytes);
        V10Header h = new V10Header();
        c.magic("HIC\0");
        require(c.word() == V10.VERSION, "unsupported file version");
        h.headerByteLength = c.wide();
        require(h.headerByteLength == bytes.length && h.headerByteLength >= 88, "header length mismatch");
        h.footer = V10Locator.read(c);
        require(h.footer.isPresent(), "missing footer");
        h.normVectorIndex = V10Locator.read(c);
        h.expectedValueIndex = V10Locator.read(c);
        h.normExpectedValueIndex = V10Locator.read(c);
        c.zero(8); // fileFlags + reserved

        h.genomeId = c.str();

        long nAttributes = c.word();
        require(nAttributes <= c.left() / 2, "attribute count out of bounds");
        for (long i = 0; i < nAttributes; i++) {
            String key = c.str();
            String value = c.str();
            h.attributeList.add(new String[]{key, value});
            h.attributes.put(key, value);
        }

        long nChromosomes = c.word();
        require(nChromosomes > 0 && nChromosomes <= c.left() / 10 && nChromosomes <= Integer.MAX_VALUE,
                "chromosome count out of bounds");
        Set<String> names = new HashSet<>();
        for (int i = 0; i < (int) nChromosomes; i++) {
            String name = c.str();
            long length = c.wide();
            require(name.length() > 0 && names.add(name) && length > 0, "invalid chromosome");
            h.chromosomes.add(new Chromosome(i, name, length));
        }

        for (int unit = 0; unit < 2; unit++) {
            List<V10Resolution> list = h.resolutions[unit];
            long n = c.word();
            require(n <= c.left() / 12, "resolution count out of bounds");
            for (int i = 0; i < (int) n; i++) {
                long binSize = c.word();
                int storageMode = c.byteValue();
                int aggregation = c.byteValue();
                c.zero(2);
                long source = c.word();
                require(binSize > 0 && binSize <= Integer.MAX_VALUE
                        && (list.isEmpty() || binSize > list.get(list.size() - 1).binSize), "invalid resolution order");
                require(storageMode <= 1 && aggregation <= 1, "unknown resolution enumeration");
                if (storageMode == V10.MATERIALIZED) {
                    require(source == 0xFFFFFFFFL, "materialized source must be absent");
                } else {
                    require(aggregation == V10.AGGREGATION_SUM && source < i
                                    && list.get((int) source).storageMode == V10.MATERIALIZED
                                    && binSize % list.get((int) source).binSize == 0,
                            "invalid derived source");
                }
                list.add(new V10Resolution((int) binSize, storageMode, aggregation, source));
            }
        }
        validateMandatoryDerivationPolicy(h.resolutions[V10.UNIT_BP]);

        h.fragmentSiteCounts = new int[h.chromosomes.size()];
        if (!h.resolutions[V10.UNIT_FRAG].isEmpty()) {
            for (int i = 0; i < h.chromosomes.size(); i++) {
                long nSites = c.word();
                require(nSites <= c.left() / 8 && nSites <= Integer.MAX_VALUE - 1,
                        "fragment-site count out of bounds");
                long previous = 0;
                for (long j = 0; j < nSites; j++) {
                    long site = c.wide();
                    require(site > previous && site < h.chromosomes.get(i).getLength(), "invalid fragment site");
                    previous = site;
                }
                h.fragmentSiteCounts[i] = (int) nSites;
            }
        }

        long nNorms = c.word();
        require(nNorms <= c.left() / 2, "normalization count out of bounds");
        names.clear();
        for (long i = 0; i < nNorms; i++) {
            String name = c.str();
            require(name.length() > 0 && !name.equals("NONE") && names.add(name), "invalid normalization name");
            h.normalizations.add(name);
        }
        c.done();

        for (int unit = 0; unit < 2; unit++) {
            for (int ri = 0; ri < h.resolutions[unit].size(); ri++) {
                for (int chr = 0; chr < h.chromosomes.size(); chr++) {
                    require(h.bins(chr, unit, ri) <= 0xFFFFFFFFL, "chromosome bin count exceeds uint32");
                }
            }
        }
        return h;
    }

    /**
     * The five mandatory BP derivations and the materialized 500 kb rule
     * (Section C.3).
     */
    private static void validateMandatoryDerivationPolicy(List<V10Resolution> bp) {
        for (int i = 0; i < bp.size(); i++) {
            V10Resolution r = bp.get(i);
            int requiredSource = requiredSourceBinSize(r.binSize);
            if (requiredSource > 0) {
                int sourceIndex = -1;
                for (int j = 0; j < bp.size(); j++) {
                    if (bp.get(j).binSize == requiredSource) {
                        sourceIndex = j;
                        break;
                    }
                }
                require(sourceIndex >= 0 && r.isDerived() && r.sourceResolutionIndex == sourceIndex
                                && !bp.get(sourceIndex).isDerived(),
                        "mandatory BP derivation policy is not satisfied");
            }
            require(r.binSize != 500000 || !r.isDerived(), "500 kb must be materialized");
        }
    }

    private static int requiredSourceBinSize(int binSize) {
        if (binSize == 20 || binSize == 50) return 10;
        if (binSize == 200 || binSize == 500) return 100;
        if (binSize == 2000) return 1000;
        return 0;
    }
}
