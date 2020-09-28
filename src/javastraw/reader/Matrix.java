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

import javastraw.reader.basics.Chromosome;
import javastraw.type.HiCZoom;
import org.broad.igv.util.Pair;

import java.util.*;

/**
 * @author jrobinso
 * @since Aug 12, 2010
 */
public class Matrix {

    private final int chr1;
    private final int chr2;
    private final static Set<Pair<Integer, Integer>> dynamicZoomResolutions = new HashSet<>();
    protected List<MatrixZoomData> bpZoomData = new ArrayList<>();
    protected List<MatrixZoomData> fragZoomData = new ArrayList<>();
    protected List<MatrixZoomData> dynamicBPZoomData = new ArrayList<>();
    private final Comparator<MatrixZoomData> comparator = new Comparator<MatrixZoomData>() {
        @Override
        public int compare(MatrixZoomData o1, MatrixZoomData o2) {
            return o2.getBinSize() - o1.getBinSize();
        }
    };

    /**
     * Constructor for creating a matrix from precomputed data.
     *
     * @param chr1
     * @param chr2
     * @param zoomDataList
     */
    public Matrix(int chr1, int chr2, List<MatrixZoomData> zoomDataList) {
        this.chr1 = chr1;
        this.chr2 = chr2;
        initZoomDataMap(zoomDataList);
    }

    private void initZoomDataMap(List<MatrixZoomData> zoomDataList) {
        for (MatrixZoomData zd : zoomDataList) {
            if (zd.getZoom().getUnit() == HiCZoom.HiCUnit.BP) {
                bpZoomData.add(zd);
            } else {
                fragZoomData.add(zd);
            }
    
            // Zooms should be sorted, but in case they are not...
    
            bpZoomData.sort(comparator);
            fragZoomData.sort(comparator);
        }

        for (Pair<Integer, Integer> resPair : dynamicZoomResolutions) {
            try {
                createDynamicResolutionMZD(resPair, false);
            } catch (Exception e) {
                System.err.println("Dynamic resolution could not be made");
            }
        }
        dynamicBPZoomData.sort(comparator);

    }

    public static String generateKey(int chr1, int chr2) {
        if (chr2 < chr1) return "" + chr2 + "_" + chr1;
        return "" + chr1 + "_" + chr2;
    }

    public static Matrix createAssemblyChromosomeMatrix(ChromosomeHandler handler,
                                                        final Map<String, Matrix> matrices, DatasetReader reader) {
        Map<HiCZoom, MatrixZoomData> assemblyZDs = new HashMap<>();

        Matrix matrix = null;
        int numAttempts = 0;
        while (matrix == null && numAttempts < 3) {
            try {
                matrix = reader.readMatrix("1_1");
            } catch (Exception ignored) {
                numAttempts++;
            }
        }

        long length = handler.getChromosomeFromName("pseudoassembly").getLength(); // TODO: scaling; also maybe chromosome ends need to shift to start with new bin at every zoom?
        for (MatrixZoomData zd : matrix.bpZoomData) {
            // todo @dudcha is this done for resolutions where conversion will be lossy?
            assemblyZDs.put(zd.getZoom(), new MatrixZoomData(handler.getChromosomeFromName("pseudoassembly"), handler.getChromosomeFromName("pseudoassembly"), zd.getZoom(),
                    (int) (length / zd.getBinSize()), (int) (length / zd.getBinSize()), null, null, reader));
        }


        //TODO: assumption that we are doing this before modifying the handler

//        for (Chromosome i : handler.getChromosomeArrayWithoutAllByAll()) {
//            for (Chromosome j : handler.getChromosomeArrayWithoutAllByAll()) {
//
//                //System.out.println("from mtrx");
//                String key = Matrix.generateKey(i, j);
//                try {
//                    Matrix m = matrices.get(key);
//                    if (m == null) {
//                        // TODO sometimes this fails once or twice, but later succeeds -
//                        // TODO high priority, needs to be fixed??????
//                        int numAttempts = 0;
//                        while (m == null && numAttempts < 3) {
//                            try {
//                                m = reader.readMatrix(key);
//                            } catch (Exception ignored) {
//                                numAttempts++;
//                            }
//                        }
//
//                        for(MatrixZoomData tempMatrixZoomData : m.bpZoomData){
//                            tempMatrixZoomData.
//                        }
//
//
//                        // modify m for each zoom
////                        matrices.put(key, m); //perhaps move it to the end
//                    }
//                    for (MatrixZoomData zd : m.bpZoomData) {
//                        updateCustomZoomDataRegions(newChr1, newChr2, handler, key, zd, assemblyZDs, reader);
//                    }
////                    for (MatrixZoomData zd : m.fragZoomData) {
////                        updateCustomZoomDataRegions(newChr1, newChr2, handler, key, zd, customZDs, reader);
////                    }
//                } catch (Exception ee) {
//                    System.err.println("Everything failed in creatingAssemblyChromosomeMatrix " + key);
//                    ee.printStackTrace();
//                }
//            }
//        }

        Matrix m = new Matrix(handler.size(), handler.size(), new ArrayList<>(assemblyZDs.values()));
        matrices.put(generateKey(handler.size(), handler.size()), m);
        return m;
    }

    public void createDynamicResolutionMZD(Pair<Integer, Integer> resPair, boolean addToSet) {
        int newRes = resPair.getFirst();
        int highRes = resPair.getSecond();

        MatrixZoomData highMZD = getZoomData(new HiCZoom(HiCZoom.HiCUnit.BP, highRes));
        MatrixZoomData newMZD = new DynamicMatrixZoomData(new HiCZoom(HiCZoom.HiCUnit.BP, newRes), highMZD);
        if (addToSet) {
            dynamicZoomResolutions.add(resPair);
        }
        dynamicBPZoomData.add(newMZD);
    }



    public static String generateKey(Chromosome chr1, Chromosome chr2) {
        int t1 = Math.min(chr1.getIndex(), chr2.getIndex());
        int t2 = Math.max(chr1.getIndex(), chr2.getIndex());
        return generateKey(t1, t2);
    }

    public String getKey() {
        return generateKey(chr1, chr2);
    }

    public MatrixZoomData getFirstZoomData() {
        if (bpZoomData != null && bpZoomData.size() > 0) {
            return getFirstZoomData(HiCZoom.HiCUnit.BP);
        } else {
            return getFirstZoomData(HiCZoom.HiCUnit.FRAG);
        }
    }

    public MatrixZoomData getFirstZoomData(HiCZoom.HiCUnit unit) {
        if (unit == HiCZoom.HiCUnit.BP) {
            return bpZoomData != null && bpZoomData.size() > 0 ? bpZoomData.get(0) : null;
        } else {
            return fragZoomData != null && fragZoomData.size() > 0 ? fragZoomData.get(0) : null;
        }
    }

    public MatrixZoomData getFirstPearsonZoomData(HiCZoom.HiCUnit unit) {
        if (unit == HiCZoom.HiCUnit.BP) {
            return bpZoomData != null ? bpZoomData.get(2) : null;
        } else {
            return fragZoomData != null ? fragZoomData.get(2) : null;
        }

    }

    public MatrixZoomData getZoomData(HiCZoom zoom) {
        int targetZoom = zoom.getBinSize();
        List<MatrixZoomData> zdList = (zoom.getUnit() == HiCZoom.HiCUnit.BP) ? bpZoomData : fragZoomData;
        //linear search for bin size, the lists are not large
        for (MatrixZoomData zd : zdList) {
            if (zd.getBinSize() == targetZoom) {
                return zd;
            }
        }

        for (MatrixZoomData zd : dynamicBPZoomData) {
            if (zd.getBinSize() == targetZoom) {
                return zd;
            }
        }

        // special exception for all by all
        if (chr1 == 0 && chr2 == 0) {

            MatrixZoomData closestValue = zdList.get(0);
            int distance = Math.abs(closestValue.getBinSize() - targetZoom);
            for (MatrixZoomData zd : zdList) {
                int cdistance = Math.abs(zd.getBinSize() - targetZoom);
                if (cdistance < distance) {
                    closestValue = zd;
                    distance = cdistance;
                }
            }

            return closestValue;
        }

        return null;
    }

    public int getNumberOfZooms(HiCZoom.HiCUnit unit) {
        return (unit == HiCZoom.HiCUnit.BP) ? bpZoomData.size() : fragZoomData.size();
    }

    public boolean isNotIntra() {
        return chr1 != chr2;
    }
}