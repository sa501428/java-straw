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


import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.tribble.util.LittleEndianInputStream;
import javastraw.StrawGlobals;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.block.*;
import javastraw.reader.datastructures.ListOfDoubleArrays;
import javastraw.reader.expected.ExpectedValueFunction;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.norm.NormalizationVector;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationHandler;
import javastraw.reader.type.NormalizationType;
import org.broad.igv.Globals;
import org.broad.igv.exceptions.HttpResponseException;
import org.broad.igv.util.Pair;
import org.broad.igv.util.ParsingUtils;

import java.io.*;
import java.util.*;


public class DatasetReaderV2 extends AbstractDatasetReader {

    /**
     * Cache of chromosome name -> array of restriction sites
     */
    private final Map<String, int[]> fragmentSitesCache = new HashMap<>();
    private final Map<String, IndexEntry> masterIndex = Collections.synchronizedMap(new HashMap<>());
    private Map<String, LargeIndexEntry> normVectorIndex;
    private final Dataset dataset;
    private int version = -1;
    private Map<String, FragIndexEntry> fragmentSitesIndex;
    private final Map<String, BlockIndex> blockIndexMap = Collections.synchronizedMap(new HashMap<>());
    private long masterIndexPos;
    private long normVectorFilePosition;
    private boolean activeStatus = true;
    public static double[] globalTimeDiffThings = new double[5];
    private final boolean useCache;
    private long nviHeaderPosition;

    public DatasetReaderV2(String path, boolean useCache) {
        super(path);
        this.useCache = useCache;
        dataset = new Dataset(this);
    }

    @Override
    public Dataset read() throws IOException {

        try {
            SeekableStream stream = ReaderTools.getValidStream(path);
            long position = 0L;

            // Read the header
            LittleEndianInputStream dis = new LittleEndianInputStream(new BufferedInputStream(stream, StrawGlobals.bufferSize));

            String magicString = dis.readString();
            position += magicString.length() + 1;
            if (!magicString.equals("HIC")) {
                throw new IOException("Magic string is not HIC, this does not appear to be a hic file.");
            }

            version = dis.readInt();
            position += 4;

            System.out.println("HiC file version: " + version);
            masterIndexPos = dis.readLong();

            position += 8;

            // will set genomeId below
            String genomeId = dis.readString();
            position += genomeId.length() + 1;

            if (version > 8) {
                // read NVI todo
                nviHeaderPosition = position;
                long nvi = dis.readLong();
                long nviSize = dis.readLong();
                System.err.println(nvi + " " + nviSize);
                position += 16;
            }

            Map<String, String> attributes = new HashMap<>();
            // Attributes  (key-value pairs)
            if (version > 4) {
                int nAttributes = dis.readInt();
                position += 4;

                for (int i = 0; i < nAttributes; i++) {
                    String key = dis.readString();
                    position += key.length() + 1;

                    String value = dis.readString();
                    position += value.length() + 1;
                    attributes.put(key, value);
                }
            }

            dataset.setAttributes(attributes);

            // Read chromosome dictionary
            int nchrs = dis.readInt();
            position += 4;

            List<Chromosome> chromosomes = new ArrayList<>(nchrs);
            for (int i = 0; i < nchrs; i++) {
                String name = dis.readString();
                position += name.length() + 1;

                long size;
                if (version > 8) {
                    size = dis.readLong();
                    position += 8;
                } else {
                    size = dis.readInt();
                    position += 4;
                }

                chromosomes.add(new Chromosome(i, name, size));
            }

            ChromosomeHandler chromosomeHandler = new ChromosomeHandler(chromosomes, genomeId, false);

            dataset.setChromosomeHandler(chromosomeHandler);
            // guess genomeID from chromosomes
            String genomeId1 = chromosomeHandler.getGenomeID();
            // if cannot find matching genomeID, set based on file
            dataset.setGenomeId(genomeId1);

            int nBpResolutions = dis.readInt();
            position += 4;

            int[] bpBinSizes = new int[nBpResolutions];
            for (int i = 0; i < nBpResolutions; i++) {
                bpBinSizes[i] = dis.readInt();
                position += 4;
            }
            dataset.setBpZooms(bpBinSizes);

            int nFragResolutions = dis.readInt();
            position += 4;

            int[] fragBinSizes = new int[nFragResolutions];
            for (int i = 0; i < nFragResolutions; i++) {
                fragBinSizes[i] = dis.readInt();
                position += 4;
            }
            dataset.setFragZooms(fragBinSizes);

            // Now we need to skip  through stream reading # fragments, stream on buffer is not needed so null it to
            // prevent accidental use
            dis = null;
            if (nFragResolutions > 0) {
                stream.seek(position);
                fragmentSitesIndex = new HashMap<>();
                Map<String, Integer> map = new HashMap<>();
                String firstChrName = null;
                for (int i = 0; i < nchrs; i++) {
                    String chr = chromosomes.get(i).getName();
                    if (!chr.equals(Globals.CHR_ALL)) {
                        firstChrName = chr;
                    }
                    byte[] buffer = new byte[4];
                    stream.readFully(buffer);
                    int nSites = (new LittleEndianInputStream(new ByteArrayInputStream(buffer))).readInt();
                    position += 4;

                    FragIndexEntry entry = new FragIndexEntry(position, nSites);
                    fragmentSitesIndex.put(chr, entry);
                    map.put(chr, nSites);

                    stream.skip(nSites * 4L);
                    position += nSites * 4L;
                }
                dataset.setFragmentCounts(map);
            }

            readFooter(masterIndexPos);

            stream.close();
        } catch (IOException e) {
            System.err.println("Error reading dataset : " + e.getLocalizedMessage());
            e.printStackTrace();
        }

        return dataset;
    }

    public String readStats() throws IOException {
        String statsFileName = path.substring(0, path.lastIndexOf('.')) + "_stats.html";
        String stats;
        BufferedReader reader = null;
        try {
            StringBuilder builder = new StringBuilder();
            reader = ParsingUtils.openBufferedReader(statsFileName);
            String nextLine;
            int count = 0; // if there is an big text file that happens to be named the same, don't read it forever
            while ((nextLine = reader.readLine()) != null && count < 1000) {
                builder.append(nextLine);
                builder.append("\n");
                count++;
            }
            stats = builder.toString();
        } finally {
            if (reader != null) {
                reader.close();
            }
        }

        return stats;
    }

    @Override
    public NormalizationVector getNormalizationVector(int chr1Idx, HiCZoom zoom, NormalizationType normalizationType) {
        return dataset.getNormalizationVector(chr1Idx, zoom, normalizationType);
    }

    @Override
    public int getDepthBase() {
        return dataset.getDepthBase();
    }

    private String checkGraphs(String graphs) {
        boolean reset = false;
        if (graphs == null) reset = true;
        else {
            Scanner scanner = new Scanner(graphs);
            try {
                while (!scanner.next().equals("[")) ;

                for (int idx = 0; idx < 2000; idx++) {
                    scanner.nextLong();
                }

                while (!scanner.next().equals("[")) ;
                for (int idx = 0; idx < 201; idx++) {
                    scanner.nextInt();
                    scanner.nextInt();
                    scanner.nextInt();
                }

                while (!scanner.next().equals("[")) ;
                for (int idx = 0; idx < 100; idx++) {
                    scanner.nextInt();
                    scanner.nextInt();
                    scanner.nextInt();
                    scanner.nextInt();
                }
            } catch (NoSuchElementException exception) {
                reset = true;
            }
        }

/*        if (reset) {
            try {
                graphs = readGraphs(null);
            } catch (IOException e) {
                graphs = null;
            }
        }*/
        return graphs;

    }

    private String readGraphs(String graphFileName) throws IOException {
        String graphs;
        try (BufferedReader reader = ParsingUtils.openBufferedReader(graphFileName)) {
            if (reader == null) return null;
            StringBuilder builder = new StringBuilder();
            String nextLine;
            while ((nextLine = reader.readLine()) != null) {
                builder.append(nextLine);
                builder.append("\n");
            }
            graphs = builder.toString();
        } catch (IOException e) {
            System.err.println("Error while reading graphs file: " + e);
            graphs = null;
        }
        return graphs;
    }


    @Override
    public boolean isActive() {
        return activeStatus;
    }

    @Override
    public void setActive(boolean status) {
        activeStatus = status;
    }

    @Override
    public int getVersion() {
        return version;
    }

    public long getNviHeaderPosition() {
        return nviHeaderPosition;
    }

    private void readFooter(long position) throws IOException {

        long currentPosition;
        if (version > 8) {
            currentPosition = determineNormVectorFilePosition(8, position);
        } else {
            currentPosition = determineNormVectorFilePosition(4, position);
        }

        currentPosition = populateMasterIndex(currentPosition);
        currentPosition = readExpectedValuesMapForNone(currentPosition);

        // Normalized expected values (v6 and greater only)
        if (version >= 6) {
            currentPosition = normVectorFilePosition;
            SeekableStream stream = ReaderTools.getValidStream(path);
            stream.seek(currentPosition);
            LittleEndianInputStream dis = new LittleEndianInputStream(new BufferedInputStream(stream, StrawGlobals.bufferSize));


            int nNormExpectedValueVectors;
            try {
                nNormExpectedValueVectors = dis.readInt();
                currentPosition += 4;
            } catch (EOFException | HttpResponseException e) {
                if (StrawGlobals.printVerboseComments) {
                    System.out.println("No normalization vectors");
                }
                return;
            }

            for (int i = 0; i < nNormExpectedValueVectors; i++) {
                stream.seek(currentPosition);
                dis = new LittleEndianInputStream(new BufferedInputStream(stream, StrawGlobals.bufferSize));

                String typeString = dis.readString();
                NormalizationType norm = dataset.getNormalizationHandler().getNormTypeFromString(typeString);
                currentPosition += (typeString.length() + 1);
                currentPosition = ReaderTools.readExpectedVectorInFooter(currentPosition, dataset.getExpectedValueFunctionMap(), norm,
                        version, path, this);
            }

            // Normalization vectors (indexed)
            if (StrawGlobals.printVerboseComments) {
                System.out.println("NVI " + currentPosition);
            }

            stream.seek(currentPosition);
            dis = new LittleEndianInputStream(new BufferedInputStream(stream, StrawGlobals.bufferSize));

            int nNormVectors = dis.readInt();
            normVectorIndex = new HashMap<>(nNormVectors * 2);
            for (int i = 0; i < nNormVectors; i++) {

                NormalizationType type = dataset.getNormalizationHandler().getNormTypeFromString(dis.readString());
                int chrIdx = dis.readInt();
                String unit = dis.readString();
                int resolution = dis.readInt();
                long filePosition = dis.readLong();
                long sizeInBytes = version > 8 ? dis.readLong() : dis.readInt();

                String key = NormalizationVector.getKey(type, chrIdx, unit, resolution);

                dataset.addNormalizationType(type);

                normVectorIndex.put(key, new LargeIndexEntry(filePosition, sizeInBytes));
            }
            stream.close();
        }
    }

    private long readExpectedValuesMapForNone(long currentPosition) throws IOException {
        Map<String, ExpectedValueFunction> expectedValuesMap = new LinkedHashMap<>();
        SeekableStream stream = ReaderTools.getValidStream(path);
        stream.seek(currentPosition);
        LittleEndianInputStream dis = new LittleEndianInputStream(new BufferedInputStream(stream, StrawGlobals.bufferSize));
        int nExpectedValues = dis.readInt();
        currentPosition += 4;
        for (int i = 0; i < nExpectedValues; i++) {
            NormalizationType norm = NormalizationHandler.NONE;
            currentPosition = ReaderTools.readExpectedVectorInFooter(currentPosition, expectedValuesMap, norm,
                    version, path, this);
        }
        dataset.setExpectedValueFunctionMap(expectedValuesMap);
        stream.close();
        return currentPosition;
    }

    private long populateMasterIndex(long currentPosition) throws IOException {
        SeekableStream stream = ReaderTools.getValidStream(path);
        stream.seek(currentPosition);
        LittleEndianInputStream dis = new LittleEndianInputStream(new BufferedInputStream(stream, StrawGlobals.bufferSize));

        int nEntries = dis.readInt();
        currentPosition += 4;

        for (int i = 0; i < nEntries; i++) {
            String key = dis.readString();
            currentPosition += (key.length() + 1);
            long filePosition = dis.readLong();
            int sizeInBytes = dis.readInt();
            currentPosition += 12;
            masterIndex.put(key, new IndexEntry(filePosition, sizeInBytes));
        }
        stream.close();
        return currentPosition;
    }

    private long determineNormVectorFilePosition(int numBytesInVar, long position) throws IOException {
        SeekableStream stream = ReaderTools.getValidStream(path);
        stream.seek(position);
        byte[] buffer = new byte[numBytesInVar];
        int actualBytes = stream.read(buffer);
        if (numBytesInVar == actualBytes) {
            LittleEndianInputStream dis = new LittleEndianInputStream(new ByteArrayInputStream(buffer));
            long nBytes;
            if (numBytesInVar == 4) {
                nBytes = dis.readInt();
            } else {
                nBytes = dis.readLong();
            }
            normVectorFilePosition = masterIndexPos + nBytes + numBytesInVar;
        } else {
            System.err.println("Actually read " + actualBytes + " bytes instead of " + numBytesInVar);
            System.exit(10);
        }
        stream.close();
        return position + actualBytes;
    }

    @Override
    public Matrix readMatrix(String key) throws IOException {
        IndexEntry idx = masterIndex.get(key);
        if (idx == null) {
            return null;
        }

        byte[] buffer = ReaderTools.seekAndFullyReadCompressedBytes(idx, path);
        LittleEndianInputStream dis = new LittleEndianInputStream(new ByteArrayInputStream(buffer));

        int c1 = dis.readInt();
        int c2 = dis.readInt();

        // TODO weird bug
        // interesting bug with local files; difficult to reliably repeat, but just occurs on loading a region
        // indices that are read (c1, c2) seem to be excessively large / wrong
        // maybe some int overflow is occurring?
        // uncomment next 2 lines to help debug
        // System.err.println("read in mtrx indcs "+c1+ "  " +c2+"  key  "+key+"    idx "+idx.position+
        //         " sz  "+idx.size+ " "+stream.getSource()+" "+stream.position()+" "+stream );
        if (c1 < 0 || c1 > dataset.getChromosomeHandler().getChromosomeArray().length ||
                c2 < 0 || c2 > dataset.getChromosomeHandler().getChromosomeArray().length) {
            System.err.println("WEIRD BUG HAPPENED AGAIN!!");
            return null;
        }

        Chromosome chr1 = dataset.getChromosomeHandler().getChromosomeFromIndex(c1);
        Chromosome chr2 = dataset.getChromosomeHandler().getChromosomeFromIndex(c2);

        // # of resolution levels (bp and frags)
        int nResolutions = dis.readInt();
        long currentFilePosition = idx.position + 12;
        List<MatrixZoomData> zdList = new ArrayList<>();
        int[] chr1Sites = retrieveFragmentSitesFromCache(chr1);
        int[] chr2Sites = retrieveFragmentSitesFromCache(chr2);

        for (int i = 0; i < nResolutions; i++) {
            try {
                Pair<MatrixZoomData, Long> result = ReaderTools.readMatrixZoomData(chr1, chr2, chr1Sites, chr2Sites,
                        currentFilePosition, path, useCache, blockIndexMap, this);
                zdList.add(result.getFirst());
                currentFilePosition = result.getSecond();
            } catch (Exception ee) {
                System.err.println("Weird error happened with trying to read MZD at currentFilePosition: " + currentFilePosition);
                ee.printStackTrace();
            }
        }

        return new Matrix(c1, c2, zdList);
    }

    int getFragCount(Chromosome chromosome) {
        FragIndexEntry entry = null;
        if (fragmentSitesIndex != null)
            entry = fragmentSitesIndex.get(chromosome.getName());

        if (entry != null) {
            return entry.nSites;
        } else return -1;
    }

    private synchronized int[] retrieveFragmentSitesFromCache(Chromosome chromosome) throws IOException {
        int[] chrSites = fragmentSitesCache.get(chromosome.getName());
        if (chrSites == null && fragmentSitesIndex != null) {
            FragIndexEntry entry = fragmentSitesIndex.get(chromosome.getName());
            if (entry != null && entry.nSites > 0) {
                chrSites = ReaderTools.readSites(entry.position, entry.nSites, path);
            }
            fragmentSitesCache.put(chromosome.getName(), chrSites);
        }
        return chrSites;
    }

    @Override
    public List<Integer> getBlockNumbers(MatrixZoomData zd) {
        BlockIndex blockIndex = blockIndexMap.get(zd.getKey());
        return blockIndex == null ? null : blockIndex.getBlockNumbers();
    }

    public Map<String, LargeIndexEntry> getNormVectorIndex() {
        return normVectorIndex;
    }

    public long getNormFilePosition() {
        return version <= 5 ? (new File(this.path)).length() : normVectorFilePosition;
    }

    static class FragIndexEntry {
        final long position;
        final int nSites;

        FragIndexEntry(long position, int nSites) {
            this.position = position;
            this.nSites = nSites;
        }
    }

    @Override
    public NormalizationVector readNormalizationVector(NormalizationType type, int chrIdx, HiCZoom.HiCUnit unit, int binSize) throws IOException {
        String key = NormalizationVector.getKey(type, chrIdx, unit.toString(), binSize);
        if (normVectorIndex == null) return null;
        LargeIndexEntry idx = normVectorIndex.get(key);
        boolean useVCForVCSQRT = false;
        if (idx == null && type.equals(NormalizationHandler.VC_SQRT)) {
            key = NormalizationVector.getKey(NormalizationHandler.VC, chrIdx, unit.toString(), binSize);
            idx = normVectorIndex.get(key);
            useVCForVCSQRT = true;
        }
        if (idx == null) return null;


        LittleEndianInputStream dis = ReaderTools.createStreamFromSeveralBuffers(idx, path);

        long nValues;
        if (version > 8) {
            nValues = dis.readLong();
        } else {
            nValues = dis.readInt();
        }
        return ReaderTools.createNormalizationVector(type, chrIdx, unit, binSize, useVCForVCSQRT, dis, nValues, version);
    }

    @Override
    public NormalizationVector readNormalizationVectorPart(NormalizationType type, int chrIdx, HiCZoom.HiCUnit unit, int binSize, int bound1, int bound2) throws IOException {
        String key = NormalizationVector.getKey(type, chrIdx, unit.toString(), binSize);
        if (normVectorIndex == null) return null;
        LargeIndexEntry idx = normVectorIndex.get(key);
        boolean useVCForVCSQRT = false;
        if (idx == null && type.equals(NormalizationHandler.VC_SQRT)) {
            key = NormalizationVector.getKey(NormalizationHandler.VC, chrIdx, unit.toString(), binSize);
            idx = normVectorIndex.get(key);
            useVCForVCSQRT = true;
        }
        if (idx == null) return null;

        long partPosition = version > 8 ? idx.position + 8 + 4L * bound1 : idx.position + 4 + 8L * bound1;
        long partSize = version > 8 ? (bound2 - bound1 + 1) * 4L : (bound2 - bound1 + 1) * 8L;

        LittleEndianInputStream dis = ReaderTools.createStreamFromSeveralBuffers(new LargeIndexEntry(partPosition, partSize), path);
        long nValues = bound2 - bound1 + 1;
        return ReaderTools.createNormalizationVector(type, chrIdx, unit, binSize, useVCForVCSQRT, dis, nValues, version);
    }

    @Override
    public ListOfDoubleArrays readExpectedVectorPart(long position, long nVals) throws IOException {
        long size = version > 8 ? nVals * 4 : nVals * 8;
        LargeIndexEntry idx = new LargeIndexEntry(position, size);
        LittleEndianInputStream dis = ReaderTools.createStreamFromSeveralBuffers(idx, path);
        ListOfDoubleArrays values = new ListOfDoubleArrays(nVals);
        for (int i = 0; i < nVals; i++) {
            double val = version > 8 ? dis.readFloat() : dis.readDouble();
            values.set(i, val);
        }
        return values;
    }

    @Override
    public Block readNormalizedBlock(int blockNumber, MatrixZoomData zd, NormalizationType no) throws IOException {

        if (no == null) {
            throw new IOException("Norm " + no + " is null");
        } else if (no.equals(NormalizationHandler.NONE)) {
            return readBlock(blockNumber, zd);
        } else {
            long[] timeDiffThings = new long[4];
            timeDiffThings[0] = System.currentTimeMillis();

            NormalizationVector nv1 = dataset.getNormalizationVector(zd.getChr1Idx(), zd.getZoom(), no);
            NormalizationVector nv2 = dataset.getNormalizationVector(zd.getChr2Idx(), zd.getZoom(), no);

            if (nv1 == null || nv2 == null) {
                if (StrawGlobals.printVerboseComments) { // todo should this print an error always instead?
                    System.err.println("Norm " + no + " missing for: " + zd.getDescription());
                    System.err.println(nv1 + " - " + nv2);
                }
                return null;
            }
            ListOfDoubleArrays nv1Data = nv1.getData();
            ListOfDoubleArrays nv2Data = nv2.getData();
            timeDiffThings[1] = System.currentTimeMillis();
            Block rawBlock = readBlock(blockNumber, zd);
            timeDiffThings[2] = System.currentTimeMillis();
            if (rawBlock == null) return null;

            Collection<ContactRecord> records = rawBlock.getContactRecords();
            List<ContactRecord> normRecords = new ArrayList<>(records.size());
            for (ContactRecord rec : records) {
                int x = rec.getBinX();
                int y = rec.getBinY();
                double denominator = nv1Data.get(x) * nv2Data.get(y);
                float counts = (float) (rec.getCounts() / denominator);
                if (!Float.isNaN(counts)) {
                    normRecords.add(new ContactRecord(x, y, counts));
                }
            }
            timeDiffThings[3] = System.currentTimeMillis();

            return new Block(blockNumber, normRecords, zd.getBlockKey(blockNumber, no));
        }
    }

    private Block readBlock(int blockNumber, MatrixZoomData zd) throws IOException {

        long[] timeDiffThings = new long[6];
        timeDiffThings[0] = System.currentTimeMillis();

        Block b = null;
        BlockIndex blockIndex = blockIndexMap.get(zd.getKey());
        if (blockIndex != null) {

            IndexEntry idx = blockIndex.getBlock(blockNumber);
            if (idx != null) {

                //System.out.println(" blockIndexPosition:" + idx.position);
                timeDiffThings[1] = System.currentTimeMillis();
                byte[] compressedBytes = ReaderTools.seekAndFullyReadCompressedBytes(idx, path);
                timeDiffThings[2] = System.currentTimeMillis();
                byte[] buffer;

                try {
                    buffer = ReaderTools.decompress(compressedBytes);
                    timeDiffThings[3] = System.currentTimeMillis();

                } catch (Exception e) {
                    throw new RuntimeException("Block read error: " + e.getMessage());
                }

                LittleEndianInputStream dis = new LittleEndianInputStream(new ByteArrayInputStream(buffer));
                int nRecords = dis.readInt();
                List<ContactRecord> records = new ArrayList<>(nRecords);
                timeDiffThings[4] = System.currentTimeMillis();

                if (version < 7) {
                    for (int i = 0; i < nRecords; i++) {
                        int binX = dis.readInt();
                        int binY = dis.readInt();
                        float counts = dis.readFloat();
                        records.add(new ContactRecord(binX, binY, counts));
                    }
                } else {

                    int binXOffset = dis.readInt();
                    int binYOffset = dis.readInt();

                    boolean useShort = dis.readByte() == 0;
                    boolean useShortBinX = true, useShortBinY = true;
                    if (version > 8) {
                        useShortBinX = dis.readByte() == 0;
                        useShortBinY = dis.readByte() == 0;
                    }

                    byte type = dis.readByte();

                    switch (type) {
                        case 1:
                            if (useShortBinX && useShortBinY) {
                                // List-of-rows representation
                                int rowCount = dis.readShort();
                                for (int i = 0; i < rowCount; i++) {
                                    int binY = binYOffset + dis.readShort();
                                    ReaderTools.populateContactRecordsColShort(dis, records, binXOffset, useShort, binY);
                                }
                            } else if (useShortBinX) { // && !useShortBinY
                                // List-of-rows representation
                                int rowCount = dis.readInt();
                                for (int i = 0; i < rowCount; i++) {
                                    int binY = binYOffset + dis.readInt();
                                    ReaderTools.populateContactRecordsColShort(dis, records, binXOffset, useShort, binY);
                                }
                            } else if (useShortBinY) { // && !useShortBinX
                                // List-of-rows representation
                                int rowCount = dis.readShort();
                                for (int i = 0; i < rowCount; i++) {
                                    int binY = binYOffset + dis.readShort();
                                    ReaderTools.populateContactRecordsColInt(dis, records, binXOffset, useShort, binY);
                                }
                            } else {
                                // List-of-rows representation
                                int rowCount = dis.readInt();
                                for (int i = 0; i < rowCount; i++) {
                                    int binY = binYOffset + dis.readInt();
                                    ReaderTools.populateContactRecordsColInt(dis, records, binXOffset, useShort, binY);
                                }
                            }
                            break;
                        case 2:

                            int nPts = dis.readInt();
                            int w = dis.readShort();

                            for (int i = 0; i < nPts; i++) {
                                //int idx = (p.y - binOffset2) * w + (p.x - binOffset1);
                                int row = i / w;
                                int col = i - row * w;
                                int bin1 = binXOffset + col;
                                int bin2 = binYOffset + row;

                                if (useShort) {
                                    short counts = dis.readShort();
                                    if (counts != Short.MIN_VALUE) {
                                        records.add(new ContactRecord(bin1, bin2, counts));
                                    }
                                } else {
                                    float counts = dis.readFloat();
                                    if (!Float.isNaN(counts)) {
                                        records.add(new ContactRecord(bin1, bin2, counts));
                                    }
                                }
                            }

                            break;
                        default:
                            throw new RuntimeException("Unknown block type: " + type);
                    }
                }
                b = new Block(blockNumber, records, zd.getBlockKey(blockNumber, NormalizationHandler.NONE));
                timeDiffThings[5] = System.currentTimeMillis();
                for (int ii = 0; ii < timeDiffThings.length - 1; ii++) {
                    globalTimeDiffThings[ii] += (timeDiffThings[ii + 1] - timeDiffThings[ii]) / 1000.0;
                }
            }
        }

        // If no block exists, mark with an "empty block" to prevent further attempts
        if (b == null) {
            b = new Block(blockNumber, zd.getBlockKey(blockNumber, NormalizationHandler.NONE));
        }
        return b;
    }
}
