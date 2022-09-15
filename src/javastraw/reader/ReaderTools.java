package javastraw.reader;

import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.tribble.util.LittleEndianInputStream;
import javastraw.StrawGlobals;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.block.*;
import javastraw.reader.datastructures.ListOfDoubleArrays;
import javastraw.reader.expected.ExpectedValueFunction;
import javastraw.reader.expected.ExpectedValueFunctionImpl;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.norm.NormFactorMapReader;
import javastraw.reader.norm.NormalizationVector;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationType;
import org.broad.igv.util.CompressionUtils;
import org.broad.igv.util.stream.IGVSeekableStreamFactory;

import java.io.BufferedInputStream;
import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.SequenceInputStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;

public class ReaderTools {

    private final static int dynamicResolutionLimit = 50;
    private static final IGVSeekableStreamFactory streamFactory = IGVSeekableStreamFactory.getInstance();
    private static final int maxLengthEntryName = 100;
    private static final int MAX_BYTE_READ_SIZE = Integer.MAX_VALUE - 10;

    public static SeekableStream getValidStream(String path) throws IOException {
        SeekableStream stream;
        do {
            stream = streamFactory.getStreamFor(path);
        } while (stream == null);
        return stream;
    }

    public static SeekableStream getValidStream(String path, long position) throws IOException {
        SeekableStream stream = getValidStream(path);
        stream.seek(position);
        return stream;
    }


    public static LittleEndianInputStream createStreamFromSeveralBuffers(LargeIndexEntry idx, String path) throws IOException {
        List<byte[]> buffer = seekAndFullyReadLargeCompressedBytes(idx, path);
        List<ByteArrayInputStream> disList = new ArrayList<>();
        for (int i = 0; i < buffer.size(); i++) {
            disList.add(new ByteArrayInputStream(buffer.get(i)));
        }
        return new LittleEndianInputStream(new SequenceInputStream(Collections.enumeration(disList)));
    }

    public static byte[] seekAndFullyReadCompressedBytes(IndexEntry idx, String path) throws IOException {
        byte[] compressedBytes = new byte[idx.size];
        SeekableStream stream = ReaderTools.getValidStream(path, idx.position);
        stream.readFully(compressedBytes);
        stream.close();
        return compressedBytes;
    }

    public static List<byte[]> seekAndFullyReadLargeCompressedBytes(LargeIndexEntry idx, String path) throws IOException {
        List<byte[]> compressedBytes = new ArrayList<>();
        long counter = idx.size;
        while (counter > MAX_BYTE_READ_SIZE) {
            compressedBytes.add(new byte[MAX_BYTE_READ_SIZE]);
            counter = counter - MAX_BYTE_READ_SIZE;
        }
        compressedBytes.add(new byte[(int) counter]);

        SeekableStream stream = ReaderTools.getValidStream(path, idx.position);
        for (int i = 0; i < compressedBytes.size(); i++) {
            stream.readFully(compressedBytes.get(i));
        }
        stream.close();
        return compressedBytes;
    }

    public static MatrixZoomData readMatrixZoomData(Chromosome chr1, Chromosome chr2, int[] chr1Sites, int[] chr2Sites,
                                                    long filePointer, String path, boolean useCache,
                                                    DatasetReader reader, int specificResolution,
                                                    boolean allowDynamicBlockIndex, long[] storeFilePosition) throws IOException {
        SeekableStream stream = ReaderTools.getValidStream(path, filePointer);
        LittleEndianInputStream dis = new LittleEndianInputStream(new BufferedInputStream(stream, StrawGlobals.bufferSize));

        String hicUnitStr = dis.readString();
        HiCZoom.HiCUnit unit = HiCZoom.valueOfUnit(hicUnitStr);
        dis.readInt();                // Old "zoom" index -- not used

        // Stats.  Not used yet, but we need to read them anyway
        double sumCounts = dis.readFloat();
        float occupiedCellCount = dis.readFloat();
        float stdDev = dis.readFloat();
        float percent95 = dis.readFloat();

        int binSize = dis.readInt();
        HiCZoom zoom = new HiCZoom(unit, binSize);
        // TODO: Default binSize value for "ALL" is 6197...
        //  (actually (genomeLength/1000)/500; depending on bug fix, could be 6191 for hg19);
        //  We need to make sure our maps hold a valid binSize value as default.

        int blockBinCount = dis.readInt();
        int blockColumnCount = dis.readInt();

        int nBlocks = dis.readInt();

        BlockIndices blockIndices;
        long currentFilePointer = filePointer + (9 * 4) + hicUnitStr.getBytes().length + 1; // i think 1 byte for 0 terminated string?
        if (specificResolution > 0) {
            if (binSize != specificResolution) {
                int maxPossibleBlockNumber = blockColumnCount * blockColumnCount - 1;
                blockIndices = new DynamicBlockIndices(ReaderTools.getValidStream(path), nBlocks, maxPossibleBlockNumber, currentFilePointer);
            } else {
                blockIndices = new BlockIndices(nBlocks);
                blockIndices.populateBlocks(dis);
            }
        } else {
            if (allowDynamicBlockIndex && binSize < dynamicResolutionLimit) {
                int maxPossibleBlockNumber = blockColumnCount * blockColumnCount - 1;
                blockIndices = new DynamicBlockIndices(ReaderTools.getValidStream(path), nBlocks, maxPossibleBlockNumber, currentFilePointer);
            } else {
                blockIndices = new BlockIndices(nBlocks);
                blockIndices.populateBlocks(dis);
            }
        }
        currentFilePointer += (nBlocks * 16L);

        MatrixZoomData zd = new MatrixZoomData(chr1, chr2, zoom, blockBinCount, blockColumnCount, chr1Sites, chr2Sites,
                reader, blockIndices, useCache, sumCounts);

        stream.close();
        storeFilePosition[0] = currentFilePointer;
        return zd;
    }

    public static long readExpectedVectorInFooter(long currentPosition,
                                                  Map<String, ExpectedValueFunction> expectedValuesMap,
                                                  NormalizationType norm, int version, String path, DatasetReader reader) throws IOException {
        SeekableStream stream = ReaderTools.getValidStream(path, currentPosition);
        LittleEndianInputStream dis = new LittleEndianInputStream(new BufferedInputStream(stream, 50));
        String unitString = dis.readString();
        currentPosition += (unitString.length() + 1);
        HiCZoom.HiCUnit unit = HiCZoom.valueOfUnit(unitString);
        int binSize = dis.readInt();
        currentPosition += 4;

        long[] nValues = new long[1];
        currentPosition += readVectorLength(dis, nValues, version);
        /* todo time
        if (binSize >= 500) {
            currentPosition = ReaderTools.readWholeNormalizationVector(currentPosition, dis, expectedValuesMap, unit, binSize,
                    nValues[0], norm, version);
        } else {
            currentPosition = ReaderTools.setUpPartialVectorStreaming(currentPosition, expectedValuesMap, unit, binSize,
                    nValues[0], norm, version, path, reader);
        }
         */
        currentPosition = ReaderTools.setUpPartialVectorStreaming(currentPosition, expectedValuesMap, unit, binSize,
                nValues[0], norm, version, path, reader);
        stream.close();
        return currentPosition;
    }

    public static long readVectorLength(LittleEndianInputStream dis, long[] nValues, int version) throws IOException {
        if (version > 8) {
            nValues[0] = dis.readLong();
            return 8;
        } else {
            nValues[0] = dis.readInt();
            return 4;
        }
    }

    public static long setUpPartialVectorStreaming(long currentPosition, Map<String, ExpectedValueFunction> expectedValuesMap,
                                                   HiCZoom.HiCUnit unit, int binSize, long nValues,
                                                   NormalizationType norm, int version, String path, DatasetReader reader) throws IOException {
        long skipPosition = currentPosition;
        long expectedVectorIndexPosition = currentPosition;
        if (version > 8) {
            skipPosition += (nValues * 4);
        } else {
            skipPosition += (nValues * 8);
        }

        SeekableStream stream = ReaderTools.getValidStream(path, skipPosition);
        int nNormalizationFactors = ReaderTools.readIntFromBytes(stream);
        currentPosition = skipPosition + 4;

        if (nNormalizationFactors > 0) {
            NormFactorMapReader hmReader = new NormFactorMapReader(nNormalizationFactors, version, currentPosition, path);
            currentPosition += hmReader.getOffset();
            ExpectedValueFunction df = new ExpectedValueFunctionImpl(norm, unit, binSize, nValues,
                    expectedVectorIndexPosition, hmReader.getNormFactors(), reader);
            String key = ExpectedValueFunction.getKey(unit, binSize, norm);
            expectedValuesMap.put(key, df);
        }

        stream.close();
        return currentPosition;
    }

    /*
    public long readWholeNormalizationVector(long currentPosition, LittleEndianInputStream dis,
                                             Map<String, ExpectedValueFunction> expectedValuesMap,
                                             HiCZoom.HiCUnit unit, int binSize, long nValues, NormalizationType norm,
                                             int version) throws IOException {
        ListOfDoubleArrays values = new ListOfDoubleArrays(nValues);
        if (version > 8) {
            currentPosition += readVectorOfFloats(dis, nValues, values);
        } else {
            currentPosition += readVectorOfDoubles(dis, nValues, values);
        }

        int nNormalizationFactors = dis.readInt();
        currentPosition += 4;

        NormFactorMapReader hmReader = new NormFactorMapReader(nNormalizationFactors, version,
                getValidStream(path, currentPosition));
        currentPosition += hmReader.getOffset();

        String key = ExpectedValueFunction.getKey(unit, binSize, norm);
        ExpectedValueFunction df = new ExpectedValueFunctionImpl(norm, unit, binSize, values, hmReader.getNormFactors());
        expectedValuesMap.put(key, df);
        return currentPosition;
    }
    */

    public static long readVectorOfFloats(LittleEndianInputStream dis, long nValues,
                                          ListOfDoubleArrays values) throws IOException {
        for (long j = 0; j < nValues; j++) {
            values.set(j, dis.readFloat());
        }
        return 4 * nValues;
    }

    public static long readVectorOfDoubles(LittleEndianInputStream dis, long nValues,
                                           ListOfDoubleArrays values) throws IOException {
        for (long j = 0; j < nValues; j++) {
            values.set(j, dis.readDouble());
        }
        return 8 * nValues;
    }

    public static NormalizationVector createNormalizationVector(NormalizationType type, int chrIdx,
                                                                HiCZoom.HiCUnit unit, int binSize,
                                                                boolean useVCForVCSQRT, LittleEndianInputStream dis,
                                                                long nValues, int version) throws IOException {
        ListOfDoubleArrays values = new ListOfDoubleArrays(nValues);
        boolean allNaN = true;
        for (long i = 0; i < nValues; i++) {
            double val = version > 8 ? (double) dis.readFloat() : dis.readDouble();
            if (!useVCForVCSQRT) {
                values.set(i, val);
            } else {
                values.set(i, Math.sqrt(val));
            }
            if (!Double.isNaN(val)) {
                allNaN = false;
            }
        }
        if (allNaN) return null;
        else return new NormalizationVector(type, chrIdx, unit, binSize, values);
    }

    public static void populateContactRecordsColShort(LittleEndianInputStream dis, List<ContactRecord> records, int binXOffset, boolean useShort, int binY) throws IOException {
        int colCount = dis.readShort();
        for (int j = 0; j < colCount; j++) {
            int binX = binXOffset + dis.readShort();
            float counts = useShort ? dis.readShort() : dis.readFloat();
            records.add(new ContactRecord(binX, binY, counts));
        }
    }

    public static void populateContactRecordsColInt(LittleEndianInputStream dis, List<ContactRecord> records, int binXOffset, boolean useShort, int binY) throws IOException {
        int colCount = dis.readInt();
        for (int j = 0; j < colCount; j++) {
            int binX = binXOffset + dis.readInt();
            float counts = useShort ? dis.readShort() : dis.readFloat();
            records.add(new ContactRecord(binX, binY, counts));
        }
    }

    public static byte[] decompress(byte[] compressedBytes) {
        CompressionUtils compressionUtils = new CompressionUtils();
        return compressionUtils.decompress(compressedBytes);
    }

    public static int[] readSites(long position, int nSites, String path) throws IOException {
        IndexEntry idx = new IndexEntry(position, 4 + nSites * 4);
        byte[] buffer = ReaderTools.seekAndFullyReadCompressedBytes(idx, path);
        LittleEndianInputStream les = new LittleEndianInputStream(new ByteArrayInputStream(buffer));
        int[] sites = new int[nSites];
        for (int s = 0; s < nSites; s++) {
            sites[s] = les.readInt();
        }
        return sites;
    }

    public static int readIntFromBytes(SeekableStream stream) throws IOException {
        byte[] buffer = new byte[4];
        int actualBytes = stream.read(buffer);
        if (actualBytes == 4) {
            LittleEndianInputStream dis = new LittleEndianInputStream(new ByteArrayInputStream(buffer));
            return dis.readInt();
        } else {
            System.err.println("Actually read " + actualBytes + " bytes instead of " + 4);
            System.exit(110);
            return 0;
        }
    }
}
