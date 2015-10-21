package de.helmholtz_muenchen.ibis.utils.datatypes.file;

import java.io.IOException;

import org.knime.core.data.DataCellDataInput;
import org.knime.core.data.DataCellDataOutput;
import org.knime.core.data.DataCellSerializer;
import org.knime.core.data.DataType;

@SuppressWarnings("serial")
public class SAMCell extends FileCell {

	public static final DataType TYPE = DataType.getType(SAMCell.class);
	
	SAMCell(String str) {
		super(str);
	}
	
	public static class SAMSerializer implements DataCellSerializer<SAMCell> {
		
		public void serialize(final SAMCell cell,
                final DataCellDataOutput output) throws IOException {
            output.writeUTF(cell.getStringValue());
        }
		
		public SAMCell deserialize(final DataCellDataInput input)
               throws IOException {
           String s = input.readUTF();
           return new SAMCell(s);
       }
	}
	
	public static class SAMCellFactory extends FileCellFactory {
		
		@Override
		public DataType getDataType() {
			return SAMCell.TYPE;
		}
	}
}
