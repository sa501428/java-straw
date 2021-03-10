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
import org.broad.igv.ui.util.MessageUtils;

import java.io.BufferedInputStream;
import java.io.FileInputStream;
import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URL;

/**
 * @author jrobinso
 *         Date: 12/22/12
 *         Time: 1:06 PM
 */
public class DatasetReaderFactory {

    public static DatasetReaderV2 getReaderForFile(String file, boolean useCache) throws IOException {
        String magicString = getMagicString(file);

        if (magicString != null) {
            if (magicString.equals("HIC")) {
                return new DatasetReaderV2(file, useCache);
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
                MessageUtils.showErrorMessage("File could not be found\n(" + path + ")", e2);
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
