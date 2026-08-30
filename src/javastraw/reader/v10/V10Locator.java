package javastraw.reader.v10;

import static javastraw.reader.v10.V10.require;

/**
 * A (position, length) pair addressing a record in the file. A zero position
 * and zero length mean the optional structure is absent; a pair with only one
 * zero value is invalid (Section A.3).
 */
public class V10Locator {

    public final long position;
    public final long length;

    public V10Locator(long position, long length) {
        require((position == 0) == (length == 0), "incomplete locator");
        this.position = position;
        this.length = length;
    }

    public static V10Locator read(V10Cursor cursor) {
        return new V10Locator(cursor.wide(), cursor.wide());
    }

    public boolean isPresent() {
        return length > 0;
    }
}
