package de.helmholtz_muenchen.ibis.utils.datatypes.file;

import java.io.IOException;

import org.knime.core.data.DataCellDataInput;
import org.knime.core.data.DataCellDataOutput;
import org.knime.core.data.DataCellSerializer;
import org.knime.core.data.DataType;

@SuppressWarnings("serial")
public class GVCFCell extends FileCell {

	public static final DataType TYPE = DataType.getType(GVCFCell.class);
	
	GVCFCell(String str) {
		super(str);
	}
	
	public static class GVCFSerializer implements DataCellSerializer<GVCFCell> {
		
		public void serialize(final GVCFCell cell,
                final DataCellDataOutput output) throws IOException {
            output.writeUTF(cell.getStringValue());
        }
		
		public GVCFCell deserialize(final DataCellDataInput input)
               throws IOException {
           String s = input.readUTF();
           return new GVCFCell(s);
       }
	}
	
	public static class GVCFCellFactory extends FileCellFactory {
		
		@Override
		public DataType getDataType() {
			return GVCFCell.TYPE;
		}
	}
}
