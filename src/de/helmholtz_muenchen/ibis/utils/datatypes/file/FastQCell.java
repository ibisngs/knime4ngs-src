package de.helmholtz_muenchen.ibis.utils.datatypes.file;

import org.knime.core.data.DataType;

@SuppressWarnings("serial")
public class FastQCell extends FileCell {

	public static final DataType TYPE = DataType.getType(FastQCell.class);
	
	FastQCell(String str) {
		super(str);
	}
}
