/**
 *  Copyright (C) 2016 the Knime4NGS contributors.
 *  Website: http://ibisngs.github.io/knime4ngs
 *  
 *  This file is part of the KNIME4NGS KNIME extension.
 *  
 *  The KNIME4NGS extension is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
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
