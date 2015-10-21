package de.helmholtz_muenchen.ibis.utils.datatypes.file;

import org.knime.core.data.DataCell;
import org.knime.core.data.DataType;
import org.knime.core.data.DataCellFactory.FromSimpleString;

public class FileCellFactory implements FromSimpleString {

	@Override
	public DataType getDataType() {
		return FileCell.TYPE;
	}

	@Override
	public DataCell createCell(String input) {
		return create(input);
	}
	
	public static DataCell create(final String s) {
        if (s == null) {
            throw new NullPointerException("File path must not be null!");
        }
        return new FileCell(s);
    }
}
