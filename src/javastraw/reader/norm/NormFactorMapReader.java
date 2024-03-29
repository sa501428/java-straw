package javastraw.reader.norm;

import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.tribble.util.LittleEndianInputStream;
import javastraw.reader.ReaderTools;

import java.io.BufferedInputStream;
import java.io.IOException;
import java.util.LinkedHashMap;
import java.util.Map;

public class NormFactorMapReader {
    private final int version, nFactors;
    private final Map<Integer, Double> normFactors = new LinkedHashMap<>();

    public NormFactorMapReader(int nFactors, int version, long position, String path)
            throws IOException {
        this.version = version;
        this.nFactors = nFactors;

        SeekableStream stream = ReaderTools.getValidStream(path, position);
        LittleEndianInputStream dis = new LittleEndianInputStream(new BufferedInputStream(stream, getOffset()));

        for (int j = 0; j < nFactors; j++) {
            int chrIdx = dis.readInt();
            if (version > 8) {
                normFactors.put(chrIdx, (double) dis.readFloat());
            } else {
                normFactors.put(chrIdx, dis.readDouble());
            }
        }
    }

    public Map<Integer, Double> getNormFactors() {
        return normFactors;
    }

    public int getOffset() {
        if (version > 8) {
            return 8 * nFactors;
        } else {
            return 12 * nFactors;
        }
    }
}
