/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2011-2020 Broad Institute, Aiden Lab, Rice University, Baylor College of Medicine
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 *  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 *  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 *  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 *  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 *  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 *  THE SOFTWARE.
 */

package javastraw.reader;

import com.google.common.primitives.Ints;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.expected.ExpectedValueFunction;
import javastraw.reader.expected.ExpectedValueFunctionImpl;
import javastraw.reader.norm.NormalizationVector;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationHandler;
import javastraw.reader.type.NormalizationType;
import org.broad.igv.util.Pair;
import org.broad.igv.util.collections.LRUCache;

import java.io.IOException;
import java.util.*;

public class Dataset {

    public static final String V9_DEPTH_BASE = "v9-depth-base";
    public static final int DEFAULT_V9_DEPTH_BASE = 2;
    public static final String HIC_FILE_SCALING = "hicFileScalingFactor";
    public static final String STATISTICS = "statistics";
    public static final String GRAPHS = "graphs";
    public static final String SOFTWARE = "software";
    public static final String NVI_INDEX = "nviIndex";
    public static final String NVI_LENGTH = "nviLength";

    private final Map<String, Matrix> matrices = new HashMap<>(625);
    private final DatasetReader reader;
    private final LRUCache<String, double[]> eigenvectorCache;
    private final LRUCache<String, NormalizationVector> normalizationVectorCache;
    private final String restrictionEnzyme = null;
    protected String genomeId;
    protected List<HiCZoom> bpZooms, dynamicZooms, fragZooms;
    private Map<String, ExpectedValueFunction> correctedExpectedValueFunctionMap = null;
    private Map<String, ExpectedValueFunction> expectedValueFunctionMap = null;
    private int v9DepthBase = 0;
    private List<Integer> bpZoomResolutions;
    private Map<String, String> attributes;
    private Map<String, Integer> fragmentCounts;
    protected NormalizationHandler normalizationHandler = new NormalizationHandler();
    private List<NormalizationType> normalizationTypes;
    private ChromosomeHandler chromosomeHandler;

    public Dataset(DatasetReader reader) {
        this.reader = reader;
        eigenvectorCache = new LRUCache<>(25);
        normalizationVectorCache = new LRUCache<>(25);
        normalizationTypes = new ArrayList<>();
    }

    public void clearCache() {
        for (Matrix matrix : matrices.values()) {
            try {
                matrix.clearCache();
            } catch (Exception e) {
            }
        }
        eigenvectorCache.clear();
        normalizationVectorCache.clear();
        normalizationTypes.clear();
    }

    public Matrix getMatrix(Chromosome chr1, Chromosome chr2) {

        // order is arbitrary, convention is lower # chr first
        if (chr1 == null || chr2 == null) return null;

        //System.out.println("from dataset");
        String key = Matrix.generateKey(chr1, chr2);
        Matrix m = matrices.get(key);

        if (m == null && reader != null) {
            try {
                m = reader.readMatrix(key);
                matrices.put(key, m);
            } catch (Exception e) {
                System.err.println("Error fetching matrix for: " + chr1.getName() + "-" + chr2.getName());
                e.printStackTrace();
            }
        }

        return m;
    }

    public void addDynamicResolution(int newRes) {

        int highRes = -1;
        for (int potentialRes : bpZoomResolutions) {
            if (potentialRes < newRes && potentialRes > highRes && newRes % potentialRes == 0) {
                highRes = potentialRes;
            }
        }
        if (highRes < 0) {
            System.err.println("No suitable higher resolution found");
            return;
        }

        for (Matrix matrix : matrices.values()) {
            matrix.createDynamicResolutionMZD(new Pair<>(newRes, highRes), true);
        }
        dynamicZooms.add(new HiCZoom(HiCZoom.HiCUnit.BP, newRes));
    }



    public void setAttributes(Map<String, String> map) {
        this.attributes = map;
        try {
            v9DepthBase = Integer.parseInt(attributes.get(V9_DEPTH_BASE));
        } catch (Exception e) {
            v9DepthBase = 0;
        }
    }

    public List<NormalizationType> getNormalizationTypes() {
        return normalizationTypes;
    }

    public Map<String, NormalizationType> getNormalizationTypesMap() {
        Map<String, NormalizationType> normsForDatasetMap = new HashMap<>();
        for (NormalizationType nType : normalizationTypes) {
            normsForDatasetMap.put(nType.getLabel(), nType);
        }
        return normsForDatasetMap;
    }

    public void setNormalizationTypes(List<NormalizationType> normalizationTypes) {
        this.normalizationTypes = normalizationTypes;
    }

    public void addNormalizationType(NormalizationType type) {
        if (!normalizationTypes.contains(type)) normalizationTypes.add(type);
    }

    public int getNumberZooms(HiCZoom.HiCUnit unit) {
        return unit == HiCZoom.HiCUnit.BP ? bpZooms.size() + dynamicZooms.size() : fragZooms.size();
    }

    // todo deprecate
    public HiCZoom getZoom(HiCZoom.HiCUnit unit, int index) {
        return unit == HiCZoom.HiCUnit.BP ? bpZooms.get(index) : fragZooms.get(index);
    }

    public HiCZoom getZoomForBPResolution(Integer resolution) {
        for (HiCZoom zoom : bpZooms) {
            if (zoom.getBinSize() == resolution) {
                return zoom;
            }
        }
        for (HiCZoom zoom : dynamicZooms) {
            if (zoom.getBinSize() == resolution) {
                return zoom;
            }
        }
        return null;
    }

    public ExpectedValueFunction getExpectedValues(HiCZoom zoom, NormalizationType type, boolean getCorrectedVersion) {
        if (expectedValueFunctionMap == null || zoom == null || type == null) return null;
        String key = ExpectedValueFunctionImpl.getKey(zoom, type);
        if (getCorrectedVersion) {
            return getCorrectedVersionOfExpectedVector(zoom, type);
        }
        return expectedValueFunctionMap.get(key);
    }

    private ExpectedValueFunction getCorrectedVersionOfExpectedVector(HiCZoom zoom, NormalizationType type) {
        if (correctedExpectedValueFunctionMap == null || expectedValueFunctionMap == null
                || zoom == null || type == null) {
            return null;
        }
        String key = ExpectedValueFunctionImpl.getKey(zoom, type);
        if (!correctedExpectedValueFunctionMap.containsKey(key)) {
            ExpectedValueFunction evf = expectedValueFunctionMap.get(key);
            if (evf == null) return null;
            correctedExpectedValueFunctionMap.put(key, evf.getCorrectedVersion());
        }
        return correctedExpectedValueFunctionMap.get(key);
    }

    public ExpectedValueFunction getExpectedValuesOrExit(HiCZoom zoom, NormalizationType type, Chromosome chromosome,
                                                         boolean isIntra, boolean getCorrectedVersion) {
        ExpectedValueFunction df = getExpectedValues(zoom, type, getCorrectedVersion);
        if (isIntra && df == null) {
            System.err.println("O/E data not available at " + chromosome.getName() + " " + zoom + " " + type);
            System.exit(14);
        }
        return df;
    }

    public Map<String, ExpectedValueFunction> getExpectedValueFunctionMap() {
        return expectedValueFunctionMap;
    }

    public void setExpectedValueFunctionMap(Map<String, ExpectedValueFunction> df) {
        this.expectedValueFunctionMap = df;
        correctedExpectedValueFunctionMap = new HashMap<>();
    }

    public ChromosomeHandler getChromosomeHandler() {
        return chromosomeHandler;
    }

    public void setChromosomeHandler(ChromosomeHandler chromosomeHandler) {
        this.chromosomeHandler = chromosomeHandler;
    }

    public int getVersion() {
        return reader.getVersion();
    }

    public String getGenomeId() {
        return genomeId;
    }

    public void setGenomeId(String genomeId) {
        if (genomeId.equals("GRCm38"))
            genomeId = "mm10";
        this.genomeId = genomeId;
    }

    public String getSoftware() {
        if (attributes != null) return attributes.get(SOFTWARE);
        else return null;
    }

    public String getHiCFileScalingFactor() {
        if (attributes != null) return attributes.get(HIC_FILE_SCALING);
        else return null;
    }

    public String getStatistics() {
        String stats = null;
        if (attributes != null) stats = attributes.get(STATISTICS);
        if (stats == null) {
            try {
                attributes.put(STATISTICS, reader.readStats());
            } catch (IOException error) {
                return null;
            }
        }
        return attributes.get(STATISTICS);
    }

    public String getGraphs() {
        if (attributes == null) return null;
        return attributes.get(GRAPHS);
    }

    public List<HiCZoom> getBpZooms() {
        List<HiCZoom> zooms = new ArrayList<>(bpZooms);
        zooms.addAll(dynamicZooms);
        zooms.sort(Collections.reverseOrder());
        return zooms;
    }

    public void setBpZooms(int[] bpBinSizes) {

        bpZoomResolutions = Ints.asList(bpBinSizes);

        bpZooms = new ArrayList<>(bpBinSizes.length);
        for (int bpBinSize : bpZoomResolutions) {
            bpZooms.add(new HiCZoom(HiCZoom.HiCUnit.BP, bpBinSize));
        }
        dynamicZooms = new ArrayList<>();
    }

    public List<HiCZoom> getFragZooms() {
        return fragZooms;
    }

    public void setFragZooms(int[] fragBinSizes) {

        // Don't show fragments in restricted mode
//        if (MainWindow.isRestricted()) return;

        this.fragZooms = new ArrayList<>(fragBinSizes.length);
        for (int fragBinSize : fragBinSizes) {
            fragZooms.add(new HiCZoom(HiCZoom.HiCUnit.FRAG, fragBinSize));
        }
    }

    public boolean hasFrags() {
        return fragZooms != null && fragZooms.size() > 0;
    }

    public Map<String, Integer> getFragmentCounts() {
        return fragmentCounts;
    }

    public void setFragmentCounts(Map<String, Integer> map) {
        fragmentCounts = map;
    }
    
    /**
     * Return the "next" zoom level, relative to the current one, in the direction indicated
     *
     * @param zoom               - current zoom level
     * @param useIncreasingOrder -- direction, true == increasing resolution, false decreasing
     * @return Next zoom level
     */

    public HiCZoom getNextZoom(HiCZoom zoom, boolean useIncreasingOrder) {
        final HiCZoom.HiCUnit currentUnit = zoom.getUnit();
        List<HiCZoom> zoomList = currentUnit == HiCZoom.HiCUnit.BP ? getBpZooms() : fragZooms;

        // TODO MSS - is there a reason not to just rewrite this using indexOf? cleaner?
        if (useIncreasingOrder) {
            for (int i = 0; i < zoomList.size() - 1; i++) {
                if (zoom.equals(zoomList.get(i))) return zoomList.get(i + 1);
            }
            return zoomList.get(zoomList.size() - 1);

        } else {
            // Decreasing
            for (int i = zoomList.size() - 1; i > 0; i--) {
                if (zoom.equals(zoomList.get(i))) {
                    return zoomList.get(i - 1);
                }
            }
            return zoomList.get(0);
        }
    }

    public NormalizationVector getNormalizationVector(int chrIdx, HiCZoom zoom, NormalizationType type) {

        String key = NormalizationVector.getKey(type, chrIdx, zoom.getUnit().toString(), zoom.getBinSize());

        if (type.equals(NormalizationHandler.NONE)) {
            return null;
        }  else if (!normalizationVectorCache.containsKey(key)) {
            try {
                NormalizationVector nv = reader.readNormalizationVector(type, chrIdx, zoom.getUnit(), zoom.getBinSize());
                normalizationVectorCache.put(key, nv);
            } catch (IOException e) {
                normalizationVectorCache.put(key, null);
            }
        }

        return normalizationVectorCache.get(key);
    }

    public NormalizationVector getPartNormalizationVector(int chrIdx, HiCZoom zoom, NormalizationType type, int bound1, int bound2) {
        String key = NormalizationVector.getKey(type, chrIdx, zoom.getUnit().toString(), zoom.getBinSize());
        NormalizationVector nv;

        if (type.equals(NormalizationHandler.NONE)) {
            return null;
        } else {
            try {
                nv = reader.readNormalizationVectorPart(type, chrIdx, zoom.getUnit(), zoom.getBinSize(), bound1, bound2);
            } catch (IOException e) {
                return null;
            }
        }

        return nv;
    }

    public List<HiCZoom> getAllPossibleResolutions() {
        List<HiCZoom> resolutions = new ArrayList<>();
        resolutions.addAll(bpZooms);
        resolutions.addAll(dynamicZooms);
        resolutions.addAll(fragZooms);
        return resolutions;
    }

    public NormalizationHandler getNormalizationHandler() {
        return normalizationHandler;
    }

    public int getDepthBase() {
        return v9DepthBase;
    }

    public String getPath() {
        return reader.getPath();
    }
}
