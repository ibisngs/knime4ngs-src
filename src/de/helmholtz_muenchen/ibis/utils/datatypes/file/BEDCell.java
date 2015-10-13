package de.helmholtz_muenchen.ibis.utils.datatypes.file;

import org.knime.core.data.DataType;

@SuppressWarnings("serial")
public class BEDCell extends FileCell {

	public static final DataType TYPE = DataType.getType(BEDCell.class);
	
	BEDCell(String str) {
		super(str);
	}

}
