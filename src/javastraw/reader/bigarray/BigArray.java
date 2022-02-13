package javastraw.reader.bigarray;

import javastraw.reader.datastructures.ListOfDoubleArrays;
import javastraw.reader.datastructures.ListOfFloatArrays;
import javastraw.reader.iterators.IteratorContainer;
import javastraw.tools.ParallelizationTools;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;

public class BigArray {

    private final int limit; // 10 million
    private final List<int[]> binXs = new ArrayList<>();
    private final List<int[]> binYs = new ArrayList<>();
    private final List<float[]> binVals = new ArrayList<>();
    private long numOfContactRecords = 0;

    public BigArray(int limit) {
        this.limit = limit;
    }

    public static void matrixVectorMult(ListOfFloatArrays vector,
                                        ListOfDoubleArrays sumVector, int x, int y, float c) {
        double counts = c;
        if (x == y) {
            counts *= .5;
        }
        sumVector.addTo(x, counts * vector.get(y));
        sumVector.addTo(y, counts * vector.get(x));
    }

    public void addSubList(int[] x, int[] y, float[] c) {
        binXs.add(x);
        binYs.add(y);
        binVals.add(c);
        numOfContactRecords += x.length;
    }

    public void addSubList(int[] x, int[] y, float[] c, int counter) {
        int[] x2 = new int[counter];
        int[] y2 = new int[counter];
        float[] c2 = new float[counter];
        System.arraycopy(x, 0, x2, 0, counter);
        System.arraycopy(y, 0, y2, 0, counter);
        System.arraycopy(c, 0, c2, 0, counter);
        addSubList(x2, y2, c2);
    }

    public void addAllSubLists(BigArray other) {
        binXs.addAll(other.binXs);
        binYs.addAll(other.binYs);
        binVals.addAll(other.binVals);
        for (float[] records : other.binVals) {
            numOfContactRecords += records.length;
        }
    }

    public long getTotalSize() {
        return numOfContactRecords;
    }

    public int getNumLists() {
        return binXs.size();
    }

    public void clear() {
        binXs.clear();
        binYs.clear();
        binVals.clear();
        numOfContactRecords = 0;
    }

    public ListOfFloatArrays sparseMultiplyAcrossLists(ListOfFloatArrays vector, long vectorLength) {
        final ListOfDoubleArrays totalSumVector = new ListOfDoubleArrays(vectorLength);

        AtomicInteger index = new AtomicInteger(0);
        ParallelizationTools.launchParallelizedCode(IteratorContainer.numCPUMatrixThreads, () -> {
            int sIndx = index.getAndIncrement();
            ListOfDoubleArrays sumVector = new ListOfDoubleArrays(vectorLength);
            while (sIndx < binXs.size()) {
                int[] subBinXs = binXs.get(sIndx);
                int[] subBinYs = binYs.get(sIndx);
                float[] subBinVals = binVals.get(sIndx);

                for (int z = 0; z < subBinXs.length; z++) {
                    matrixVectorMult(vector, sumVector,
                            subBinXs[z], subBinYs[z], subBinVals[z]);
                }
                sIndx = index.getAndIncrement();
            }

            synchronized (totalSumVector) {
                totalSumVector.addValuesFrom(sumVector);
            }
        });

        return totalSumVector.convertToFloats();
    }
}
