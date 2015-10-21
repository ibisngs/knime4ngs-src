package de.helmholtz_muenchen.ibis.utils.datatypes.file;

import java.io.IOException;

import org.knime.core.data.DataCellDataInput;
import org.knime.core.data.DataCellDataOutput;
import org.knime.core.data.DataCellSerializer;
import org.knime.core.data.DataType;

@SuppressWarnings("serial")
public class VCFCell extends FileCell{
	
	public static final DataType TYPE = DataType.getType(VCFCell.class);
	
	VCFCell(String f) {
		super(f);
	}
	
	public static class VCFSerializer implements DataCellSerializer<VCFCell> {
		
		public void serialize(final VCFCell cell,
                final DataCellDataOutput output) throws IOException {
            output.writeUTF(cell.getStringValue());
        }
		
		public VCFCell deserialize(final DataCellDataInput input)
               throws IOException {
           String s = input.readUTF();
           return new VCFCell(s);
       }
	}
	
	public static class VCFCellFactory extends FileCellFactory {
		
		@Override
		public DataType getDataType() {
			return VCFCell.TYPE;
		}
	}
}
