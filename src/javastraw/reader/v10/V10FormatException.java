package javastraw.reader.v10;

/**
 * Raised whenever a V10 byte stream violates the normative requirements of
 * hic-format/HiCFormatV10.md. A malformed optional structure is never silently
 * treated as absent; corruption is reported instead of returning a plausible
 * partial result (Section N).
 */
public class V10FormatException extends RuntimeException {

    public V10FormatException(String message) {
        super("V10: " + message);
    }

    public V10FormatException(String message, Throwable cause) {
        super("V10: " + message, cause);
    }
}
