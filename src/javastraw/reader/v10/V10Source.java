package javastraw.reader.v10;

import htsjdk.samtools.seekablestream.SeekableStream;
import javastraw.reader.ReaderTools;

import java.io.IOException;

import static javastraw.reader.v10.V10.require;

/**
 * Random-access byte source for a local or remote .hic file. Every read is
 * validated against the file length (when known) and against the per-record
 * allocation ceiling before any buffer is allocated (Section N).
 */
public class V10Source {

    private final String path;
    private final long fileLength;

    public V10Source(String path) throws IOException {
        this.path = path;
        long length = 0;
        SeekableStream stream = null;
        try {
            stream = ReaderTools.getValidStream(path);
            length = stream.length();
        } catch (Exception e) {
            // Some remote sources do not advertise a length; interval checks are then
            // deferred to the read itself.
            length = 0;
        } finally {
            if (stream != null) {
                try {
                    stream.close();
                } catch (IOException ignored) {
                }
            }
        }
        this.fileLength = Math.max(length, 0);
    }

    public String getPath() {
        return path;
    }

    public long getFileLength() {
        return fileLength;
    }

    /**
     * Validates that [position, position + length) lies within the file.
     */
    public void interval(long position, long length) {
        require(length > 0 && position >= 0, "file interval out of bounds");
        require(V10.add(position, length) >= 0, "file interval out of bounds");
        if (fileLength > 0) {
            require(position <= fileLength && length <= fileLength - position, "file interval out of bounds");
        }
    }

    public byte[] read(long position, long length) throws IOException {
        require(length > 0 && length <= V10.ALLOCATION_LIMIT, "record exceeds allocation limit");
        interval(position, length);
        byte[] buffer = new byte[(int) length];
        SeekableStream stream = ReaderTools.getValidStream(path, position);
        try {
            stream.readFully(buffer);
        } finally {
            stream.close();
        }
        return buffer;
    }
}
