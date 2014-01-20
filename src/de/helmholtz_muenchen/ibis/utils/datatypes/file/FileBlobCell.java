package de.helmholtz_muenchen.ibis.utils.datatypes.file;

import java.io.IOException;

import org.knime.core.data.DataCell;
import org.knime.core.data.DataCellDataInput;
import org.knime.core.data.DataCellDataOutput;
import org.knime.core.data.DataCellSerializer;
import org.knime.core.data.DataValue;
import org.knime.core.data.StringValue;
import org.knime.core.data.container.BlobDataCell;

public class FileBlobCell extends BlobDataCell implements StringValue, FileValue {
    /**
     * Returns the preferred value class of this cell implementation. This
     * method is called per reflection to determine which is the preferred
     * renderer, comparator, etc.
     *
     * @return CMLValue.class;
     */
    public static final Class<? extends DataValue> getPreferredValueClass() {
        return FileValue.class;
    }

    private static final FileBlobSerializer SERIALIZER = new FileBlobSerializer();

    /**
     * Returns the factory to read/write DataCells of this class from/to a
     * DataInput/DataOutput. This method is called via reflection.
     *
     * @return a serializer for reading/writing cells of this kind
     * @see DataCell
     */
    public static final FileBlobSerializer getCellSerializer() {
        return SERIALIZER;
    }

    private final String m_FileString;

    /**
     * Creates a new CML Cell based on the given String value. This constructor
     * is used from the {@link CMLCellFactory#create(String) create} method of
     * the accompanying {@link CMLCellFactory}.
     *
     * @param str the String value to store
     * @throws NullPointerException if the given String value is
     *             <code>null</code>
     */
    FileBlobCell(final String str) {
        if (str == null) {
            throw new NullPointerException("CML value must not be null.");
        }
        m_FileString = str;
    }

    /** {@inheritDoc} */
    public String getStringValue() {
        return m_FileString;
    }

    /**
     * {@inheritDoc}
     */
    public String getCMLValue() {
        return m_FileString;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public String toString() {
        return getStringValue();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected boolean equalsDataCell(final DataCell dc) {
        return m_FileString.equals(((FileBlobCell)dc).m_FileString);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int hashCode() {
        return m_FileString.hashCode();
    }

    /** Factory for (de-)serializing a CMLBlobCell. */
    private static class FileBlobSerializer implements
            DataCellSerializer<FileBlobCell> {
        /**
         * {@inheritDoc}
         */
        public void serialize(final FileBlobCell cell,
                final DataCellDataOutput output) throws IOException {
            output.writeUTF(cell.getStringValue());
        }

        /**
         * {@inheritDoc}
         */
        public FileBlobCell deserialize(final DataCellDataInput input)
                throws IOException {
            String s = input.readUTF();
            return new FileBlobCell(s);
        }
    }

	@Override
	public String getFilePath() {
		// TODO Auto-generated method stub
		return m_FileString;
	}
}
