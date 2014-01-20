package de.helmholtz_muenchen.ibis.utils.datatypes.file;

import org.knime.core.data.DataCell;
import org.knime.core.data.DataType;
import org.knime.core.node.NodeLogger;

public class FileCellFactory {
    /**
     * Minimum size for blobs in bytes. That is, if a given string is at least
     * as large as this value, it will be represented by a blob cell
     */
    public static final int DEF_MIN_BLOB_SIZE_IN_BYTES = 8 * 1024;

    private static final int MIN_BLOB_SIZE_IN_BYTES;

    static {
        int size = DEF_MIN_BLOB_SIZE_IN_BYTES;
        String envVar = "org.knime.Fileminblobsize";
        String property = System.getProperty(envVar);
        if (property != null) {
            String s = property.trim();
            int multiplier = 1;
            if (s.endsWith("m") || s.endsWith("M")) {
                s = s.substring(0, s.length() - 1);
                multiplier = 1024 * 1024;
            } else if (s.endsWith("k") || s.endsWith("K")) {
                s = s.substring(0, s.length() - 1);
                multiplier = 1024;
            }
            try {
                int newSize = Integer.parseInt(s);
                if (newSize < 0) {
                    throw new NumberFormatException("Size < 0" + newSize);
                }
                size = newSize * multiplier;
                NodeLogger.getLogger(FileCellFactory.class).debug(
                        "Setting min blob size for File cells to " + size
                                + " bytes");
            } catch (NumberFormatException e) {
                NodeLogger.getLogger(FileCellFactory.class).warn(
                        "Unable to parse property " + envVar
                                + ", using default", e);
            }
        }
        MIN_BLOB_SIZE_IN_BYTES = size;
    }

    /** Type for CML cells. */
    public static final DataType TYPE = FileCell.TYPE;

    /** Don't instantiate this class. */
    private FileCellFactory() {

    }

    /**
     * Factory method to create {@link DataCell} representing a File path.
     * The returned cell is either of type {@link FileCell} (for small strings)
     * or {@link FileBlobCell} (otherwise, default threshold is
     * {@value #DEF_MIN_BLOB_SIZE_IN_BYTES} bytes or larger).
     *
     * @param string String representing the File path.
     * @return DataCell representing the File path
     * @throws NullPointerException if argument is null
     */
    public static DataCell create(final String string) {
        if (string == null) {
            throw new NullPointerException("File must not be null");
        }
        if (string.length() >= MIN_BLOB_SIZE_IN_BYTES) {
            return new FileBlobCell(string);
        } else {
            return new FileCell(string);
        }
    }
}
