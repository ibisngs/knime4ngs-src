package de.helmholtz_muenchen.ibis.utils.datatypes.file;

import java.io.IOException;

import org.knime.core.data.DataCell;
import org.knime.core.data.DataCellDataInput;
import org.knime.core.data.DataCellDataOutput;
import org.knime.core.data.DataCellSerializer;
import org.knime.core.data.DataType;
import org.knime.core.data.StringValue;

@SuppressWarnings("serial")
public class FileCell extends DataCell implements FileValue, StringValue{

	public static final DataType TYPE = DataType.getType(FileCell.class);
	
	private final String m_file;
	
	FileCell(String f) {
		this.m_file = f;
	}

	public static class FileSerializer implements DataCellSerializer<FileCell> {
		
		 public void serialize(final FileCell cell,
	                final DataCellDataOutput output) throws IOException {
	            output.writeUTF(cell.getStringValue());
	        }
		
		public FileCell deserialize(final DataCellDataInput input)
                throws IOException {
            String s = input.readUTF();
            return new FileCell(s);
        }
	}
	
	@Override
	public String toString() {
		return m_file;
	}

	@Override
	protected boolean equalsDataCell(DataCell dc) {
		FileCell fc = (FileCell) dc;
		return this.m_file.equals(fc.m_file);
	}

	@Override
	public int hashCode() {
		return this.m_file.hashCode();
	}

	@Override
	public String getStringValue() {
		return m_file;
	}

	@Override
	public String getFilePath() {
		return m_file;
	}
}
