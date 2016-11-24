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

import org.knime.core.data.DataCell;
import org.knime.core.data.DataCellFactory.FromSimpleString;
import org.knime.core.data.DataType;

public class FileCellFactory implements FromSimpleString {

	@Override
	public DataType getDataType() {
		return FileCell.TYPE;
	}

	@Override
	public DataCell createCell(String input) {
		return create(input);
	}
	
	public static FileCell create(final String s) {
        if (s == null) {
            throw new NullPointerException("File path must not be null!");
        }
        return new FileCell(s);
    }
}
