package de.helmholtz_muenchen.ibis.utils.datatypes.file;

import java.io.IOException;

import org.knime.core.data.DataCellDataInput;
import org.knime.core.data.DataCellDataOutput;
import org.knime.core.data.DataCellSerializer;
import org.knime.core.data.DataType;

@SuppressWarnings("serial")
public class BAMCell extends FileCell {

	public static final DataType TYPE = DataType.getType(BAMCell.class);
	
	BAMCell(String str) {
		super(str);
	}
	
	public static class BAMSerializer implements DataCellSerializer<BAMCell> {
		
		public void serialize(final BAMCell cell,
                final DataCellDataOutput output) throws IOException {
            output.writeUTF(cell.getStringValue());
        }
		
		public BAMCell deserialize(final DataCellDataInput input)
               throws IOException {
           String s = input.readUTF();
           return new BAMCell(s);
       }
	}
	
	public static class BAMCellFactory extends FileCellFactory {
		
		@Override
		public DataType getDataType() {
			return BAMCell.TYPE;
		}
	}
}
