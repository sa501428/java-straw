package javastraw.igv;


import org.broad.igv.util.HttpUtils;
import org.broad.igv.util.ResourceLocator;

import java.io.*;
import java.net.URL;
import java.net.URLDecoder;
import java.nio.charset.StandardCharsets;
import java.util.zip.GZIPInputStream;

/**
 * @author jrobinso
 */
public class ParsingUtils {


    /**
     * Open a BufferedReader on the path, which might be
     * a local file or URL, and might be gzipped or not.
     *
     * @param pathOrUrl
     * @return
     * @throws IOException
     */
    public static BufferedReader openBufferedReader(String pathOrUrl) throws IOException {
        BufferedReader reader = new BufferedReader(new InputStreamReader(openInputStream(pathOrUrl)));
        return reader;
    }

    public static InputStream openInputStream(String path) throws IOException {
        return openInputStreamGZ(new ResourceLocator(path));
    }

    /**
     * Open an InputStream on the resource.  Wrap it in a GZIPInputStream if necessary.
     *
     * @param locator
     * @return
     * @throws IOException
     */
    public static InputStream openInputStreamGZ(ResourceLocator locator) throws IOException {

        InputStream inputStream = null;

        if (isDataURL(locator.getPath())) {
            String dataURL = locator.getPath();
            int commaIndex = dataURL.indexOf(',');
            if (commaIndex < 0) {
                throw new Error("dataURL missing commas");
            }
            // TODO -- check optional media type - only text/plain supported
            String contents = URLDecoder.decode(dataURL.substring(commaIndex + 1), String.valueOf(StandardCharsets.UTF_8));
            return new ByteArrayInputStream(contents.getBytes(StandardCharsets.UTF_8));
        } else if (HttpUtils.isRemoteURL(locator.getPath())) {
            URL url = HttpUtils.createURL(locator.getPath());
            inputStream = HttpUtils.getInstance().openConnectionStream(url);
        } else {
            String path = locator.getPath();
            if (path.startsWith("file://")) {
                path = path.substring(7);
            }
            File file = new File(path);
            inputStream = new FileInputStream(file);
        }

        if (locator.getPath().endsWith("gz")) {
            return new GZIPInputStream(inputStream);
        } else {
            return inputStream;
        }
    }

    public static boolean isDataURL(String url) {
        return url != null && url.startsWith("data:") && !((new File(url)).exists());
    }
}
