import com.sun.net.httpserver.HttpServer;
import javastraw.reader.*;
import javastraw.reader.type.*;
import javastraw.reader.block.*;
import java.net.InetSocketAddress;
import java.nio.file.*;
import java.util.concurrent.atomic.AtomicInteger;

public class HttpRangeRegression {
    public static void main(String[] args) throws Exception {
        byte[] bytes = Files.readAllBytes(Paths.get(args[0]));
        AtomicInteger ranges = new AtomicInteger();
        HttpServer server = HttpServer.create(new InetSocketAddress("127.0.0.1", 0), 0);
        server.createContext("/cis.hic", exchange -> {
            exchange.getResponseHeaders().set("Content-Length", Integer.toString(bytes.length));
            if (exchange.getRequestMethod().equals("HEAD")) {
                exchange.sendResponseHeaders(200, -1);
                exchange.close();
                return;
            }
            String range = exchange.getRequestHeaders().getFirst("Range");
            int begin = 0, end = bytes.length - 1;
            if (range != null) {
                String[] interval = range.substring(6).split("-", -1);
                begin = Integer.parseInt(interval[0]);
                if (!interval[1].isEmpty()) end = Math.min(end, Integer.parseInt(interval[1]));
                ranges.incrementAndGet();
                exchange.getResponseHeaders().set("Content-Range", "bytes " + begin + "-" + end + "/" + bytes.length);
            }
            exchange.sendResponseHeaders(range == null ? 200 : 206, end - begin + 1);
            exchange.getResponseBody().write(bytes, begin, end - begin + 1);
            exchange.close();
        });
        server.start();
        try {
            String url = "http://127.0.0.1:" + server.getAddress().getPort() + "/cis.hic";
            DatasetReader reader = DatasetReaderFactory.getReader(url, true, false);
            if (!(reader instanceof DatasetReaderV10)) throw new AssertionError("HTTP factory chose wrong reader");
            Dataset dataset = reader.read();
            int count = 0;
            for (Block block : dataset.getMatrix(dataset.getChromosomeHandler().getChromosomeFromIndex(0), dataset.getChromosomeHandler().getChromosomeFromIndex(0)).getZoomData(new HiCZoom(20)).getNormalizedBlocksOverlapping(0, 0, 3, 3, NormalizationHandler.NONE, false)) {
                count += block.getContactRecords().size();
            }
            if (count != 3 || ranges.get() == 0) throw new AssertionError("HTTP derived query failed");
            System.out.println("Java HTTP ranges: factory, header, vectors, blocks, and derived query passed");
        } finally { server.stop(0); }
    }
}
