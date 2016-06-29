package de.helmholtz_muenchen.ibis.utils.datatypes.file;

import java.io.IOException;

import org.knime.core.data.DataCellDataInput;
import org.knime.core.data.DataCellDataOutput;
import org.knime.core.data.DataCellSerializer;
import org.knime.core.data.DataType;

@SuppressWarnings("serial")
public class FastACell extends FileCell {

	public static final DataType TYPE = DataType.getType(FastACell.class);
	
	FastACell(String str) {
		super(str);
	}
	
	public static class FastASerializer implements DataCellSerializer<FastACell> {
		
		public void serialize(final FastACell cell,
                final DataCellDataOutput output) throws IOException {
            output.writeUTF(cell.getStringValue());
        }
		
		public FastACell deserialize(final DataCellDataInput input)
               throws IOException {
           String s = input.readUTF();
           return new FastACell(s);
       }
	}
	
	public static class FastACellFactory extends FileCellFactory {
		
		@Override
		public DataType getDataType() {
			return FastACell.TYPE;
		}
	}
}