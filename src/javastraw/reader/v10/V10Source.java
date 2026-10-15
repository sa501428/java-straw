package javastraw.reader.v10;

import htsjdk.samtools.seekablestream.SeekableStream;
import javastraw.reader.ReaderTools;

import java.io.IOException;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.ArrayDeque;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.FutureTask;

import static javastraw.reader.v10.V10.require;

/**
 * Random-access byte source for a local or remote .hic file. Every read is
 * validated against the file length (when known) and against the per-record
 * allocation ceiling before any buffer is allocated (Section N).
 */
public class V10Source implements AutoCloseable {

    private final String path;
    private final long fileLength;
    // Opening an HTTP SeekableStream probes Content-Length. Reuse idle stream
    // objects so subsequent ranges do not pay for another metadata request.
    private final ArrayDeque<SeekableStream> remoteStreams = new ArrayDeque<>();
    private final boolean remote;
    private final Map<V10Chunk, float[]> vectorChunks = new LinkedHashMap<V10Chunk, float[]>(16, .75f, true);
    private final ConcurrentHashMap<V10Chunk, FutureTask<float[]>> pendingVectorChunks = new ConcurrentHashMap<>();
    private long vectorBytes;
    private static final long VECTOR_CACHE_BYTES = 32L * 1024 * 1024;

    /** Immutable decoded chunks shared by normalization and expected queries. */
    public float[] readVectorChunk(final V10Chunk chunk) throws IOException {
        synchronized (vectorChunks) {
            float[] cached = vectorChunks.get(chunk);
            if (cached != null) return cached;
        }
        FutureTask<float[]> task = new FutureTask<>(() -> {
            // Recheck after registering: an earlier load may have just finished.
            synchronized (vectorChunks) {
                float[] cached = vectorChunks.get(chunk);
                if (cached != null) return cached;
            }
            float[] values = V10VectorEntry.decodeChunk(this, chunk);
            long bytes = (long) values.length * Float.BYTES;
            if (bytes <= VECTOR_CACHE_BYTES) synchronized (vectorChunks) {
                while (!vectorChunks.isEmpty() && (vectorChunks.size() >= 64 || vectorBytes + bytes > VECTOR_CACHE_BYTES)) {
                    V10Chunk oldest = vectorChunks.keySet().iterator().next();
                    vectorBytes -= (long) vectorChunks.remove(oldest).length * Float.BYTES;
                }
                vectorChunks.put(chunk, values);
                vectorBytes += bytes;
            }
            return values;
        });
        FutureTask<float[]> existing = pendingVectorChunks.putIfAbsent(chunk, task);
        if (existing != null) task = existing;
        else task.run();
        try {
            return task.get();
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            throw new IOException("Interrupted reading vector chunk", e);
        } catch (ExecutionException e) {
            Throwable cause = e.getCause();
            if (cause instanceof IOException) throw (IOException) cause;
            if (cause instanceof RuntimeException) throw (RuntimeException) cause;
            throw new IOException("Could not read vector chunk", cause);
        } finally {
            if (task.isDone()) pendingVectorChunks.remove(chunk, task);
        }
    }

    public V10Source(String path) throws IOException {
        this.path = path;
        this.remote = path.startsWith("http://") || path.startsWith("https://");
        long length = 0;
        SeekableStream stream = null;
        try {
            stream = ReaderTools.getValidStream(path);
            length = stream.length();
            if (remote) {
                remoteStreams.add(stream);
                stream = null;
            }
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
        SeekableStream stream = null;
        if (remote) synchronized (remoteStreams) { stream = remoteStreams.poll(); }
        if (stream == null) stream = ReaderTools.getValidStream(path);
        boolean success = false;
        try {
            stream.seek(position);
            stream.readFully(buffer);
            success = true;
        } finally {
            boolean retained = false;
            if (remote && success) synchronized (remoteStreams) {
                if (remoteStreams.size() < 4) {
                    remoteStreams.add(stream);
                    retained = true;
                }
            }
            if (!retained) stream.close();
        }
        return buffer;
    }

    @Override public void close() throws IOException {
        synchronized (remoteStreams) {
            while (!remoteStreams.isEmpty()) remoteStreams.remove().close();
        }
        synchronized (vectorChunks) {
            vectorChunks.clear();
            vectorBytes = 0;
        }
    }
}
