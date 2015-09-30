package de.helmholtz_muenchen.ibis.utils.datatypes.file;

import org.knime.core.data.DataType;

@SuppressWarnings("serial")
public class BAMCell extends FileCell {

	public static final DataType TYPE = DataType.getType(BAMCell.class);
	
	BAMCell(String str) {
		super(str);
	}
}
