package de.helmholtz_muenchen.ibis.utils.datatypes.file;

import java.io.IOException;

import org.knime.core.data.DataCellDataInput;
import org.knime.core.data.DataCellDataOutput;
import org.knime.core.data.DataCellSerializer;
import org.knime.core.data.DataType;

@SuppressWarnings("serial")
public class FastQCell extends FileCell {

	public static final DataType TYPE = DataType.getType(FastQCell.class);
	
	FastQCell(String str) {
		super(str);
	}
	
	public static class FastQSerializer implements DataCellSerializer<FastQCell> {
		
		public void serialize(final FastQCell cell,
                final DataCellDataOutput output) throws IOException {
            output.writeUTF(cell.getStringValue());
        }
		
		public FastQCell deserialize(final DataCellDataInput input)
               throws IOException {
           String s = input.readUTF();
           return new FastQCell(s);
       }
	}
	
	public static class FastQCellFactory extends FileCellFactory {
		
		@Override
		public DataType getDataType() {
			return FastQCell.TYPE;
		}
	}
}
