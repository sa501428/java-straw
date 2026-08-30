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

import htsjdk.samtools.seekablestream.SeekableHTTPStream;
import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.tribble.util.LittleEndianInputStream;

import java.io.BufferedInputStream;
import java.io.FileInputStream;
import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URL;

public class DatasetReaderFactory {

    /**
     * Returns the reader matching the file's declared version: {@link DatasetReaderV10}
     * for version 10 files, {@link DatasetReaderV2} for versions 6 through 9.
     * Prefer this over {@link #getReaderForFile} unless a V2 reader is required.
     */
    public static DatasetReader getReader(String file, boolean useCache, boolean useDynamicBlockIndex)
            throws IOException {
        String magicString = getMagicString(file);
        if (magicString == null || !magicString.equals("HIC")) {
            System.err.println("This version is deprecated and is no longer supported.");
            return null;
        }
        if (getVersion(file) == javastraw.reader.v10.V10.VERSION) {
            return new DatasetReaderV10(file, useCache);
        }
        return new DatasetReaderV2(file, useCache, useDynamicBlockIndex);
    }

    /**
     * The version declared in the four bytes following the magic string.
     */
    public static int getVersion(String path) throws IOException {
        SeekableStream stream = ReaderTools.getValidStream(path);
        try {
            byte[] buffer = new byte[8];
            stream.readFully(buffer);
            return ((buffer[4] & 0xFF)) | ((buffer[5] & 0xFF) << 8)
                    | ((buffer[6] & 0xFF) << 16) | ((buffer[7] & 0xFF) << 24);
        } finally {
            stream.close();
        }
    }

    public static DatasetReaderV2 getReaderForFile(String file, boolean useCache, boolean useDynamicBlockIndex) throws IOException {
        String magicString = getMagicString(file);

        if (magicString != null) {
            if (magicString.equals("HIC")) {
                return new DatasetReaderV2(file, useCache, useDynamicBlockIndex);
            } else {
                System.err.println("This version is deprecated and is no longer supported.");
                //reader = new DatasetReaderV1(file);
                // file not actually read, usually canceled the read of password-protected file
                //if (reader.getVersion() == -1)
            }
        }
        return null;
    }

    public static String getMagicString(String path) throws IOException {

        SeekableStream stream = null;
        LittleEndianInputStream dis = null;

        try {
            stream = new SeekableHTTPStream(new URL(path)); // IGVSeekableStreamFactory.getStreamFor(path);
            dis = new LittleEndianInputStream(new BufferedInputStream(stream));
        } catch (MalformedURLException e) {
            try {
                dis = new LittleEndianInputStream(new FileInputStream(path));
            } catch (Exception e2) {
                System.out.println("File could not be found\n(" + path + ")");
                System.out.println(e2.getLocalizedMessage());
            }
        } finally {
            if (stream != null) stream.close();

        }
        if (dis != null) {
            return dis.readString();
        }
        return null;
    }

}
