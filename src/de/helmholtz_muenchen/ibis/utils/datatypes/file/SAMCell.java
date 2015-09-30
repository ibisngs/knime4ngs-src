package de.helmholtz_muenchen.ibis.utils.datatypes.file;

import org.knime.core.data.DataType;

@SuppressWarnings("serial")
public class SAMCell extends FileCell {

	public static final DataType TYPE = DataType.getType(BAMCell.class);
	
	SAMCell(String str) {
		super(str);
	}
}
