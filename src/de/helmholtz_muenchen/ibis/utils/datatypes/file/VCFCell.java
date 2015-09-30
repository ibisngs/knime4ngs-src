package de.helmholtz_muenchen.ibis.utils.datatypes.file;

import org.knime.core.data.DataType;

@SuppressWarnings("serial")
public class VCFCell extends FileCell {

	public static final DataType TYPE = DataType.getType(VCFCell.class);
	
	VCFCell(String str) {
		super(str);
	}

}
